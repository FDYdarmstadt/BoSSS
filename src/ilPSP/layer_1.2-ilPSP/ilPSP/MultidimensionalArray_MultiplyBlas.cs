using System;
using System.Numerics;
using ilPSP.Utils;
using SimdVector = System.Numerics.Vector<double>;

namespace ilPSP {
    public partial class MultidimensionalArray {
        // DIAGNOSTIC (temporary): BOSSS_MULTIPLY_FASTPATH=0 disables all fast paths,
        // =gemm keeps only the direct (non-packing) DGEMM/DGEMV paths.
        static readonly string s_FastPathMode = (System.Environment.GetEnvironmentVariable("BOSSS_MULTIPLY_FASTPATH") ?? "all").ToLowerInvariant();
        static readonly bool s_FastPathOff = s_FastPathMode == "0";
        static readonly bool s_FastPathGemmOnly = s_FastPathMode == "gemm";
        // =check: run reference kernel into a copy and compare (slow!)
        static readonly bool s_FastPathCheck = s_FastPathMode == "check";
        [ThreadStatic] static string s_LastPath;
        static long s_CheckCount, s_CheckFail;
        static int s_ExitHooked;
        static readonly System.Collections.Concurrent.ConcurrentDictionary<string, long[]> s_PathStats = new System.Collections.Concurrent.ConcurrentDictionary<string, long[]>();

        static unsafe void FastPathCompare(double[] r, double[] t, string path, int DT,
            int* lenRun, int* ct, int* ca, int* cb, int* lenSum, int* sa, int* sb, double alpha, double beta) {
            double maxDiff = 0, maxRef = 0; int iMax = -1;
            for(int i = 0; i < r.Length; i++) {
                double ri = r[i], ti = t[i];
                if(double.IsNaN(ri) && double.IsNaN(ti)) continue;
                maxRef = Math.Max(maxRef, Math.Abs(ri));
                if(ri == ti) continue;
                double d = Math.Abs(ri - ti);
                if(!(d <= maxDiff)) { maxDiff = double.IsNaN(d) ? double.PositiveInfinity : d; iMax = i; }
            }
            if(System.Threading.Interlocked.Exchange(ref s_ExitHooked, 1) == 0)
                AppDomain.CurrentDomain.ProcessExit += (o, e) => {
                    var sb2 = new System.Text.StringBuilder($"FASTPATH-SUMMARY compared={s_CheckCount} mismatches={s_CheckFail}:");
                    foreach(var kv in s_PathStats) sb2.Append($" {kv.Key}={kv.Value[0]} (maxRel={BitConverter.Int64BitsToDouble(kv.Value[1]):g3})");
                    Console.Error.WriteLine(sb2.ToString());
                    try { System.IO.File.AppendAllText("FASTPATH-SUMMARY.txt", sb2.ToString() + System.Environment.NewLine); } catch(Exception) { }
                };
            {
                var st = s_PathStats.GetOrAdd(path ?? "?", _ => new long[2]);
                System.Threading.Interlocked.Increment(ref st[0]);
                double rel = maxDiff / Math.Max(maxRef, 1e-300);
                lock(st) { if(rel > BitConverter.Int64BitsToDouble(st[1])) st[1] = BitConverter.DoubleToInt64Bits(rel); }
            }
            long n = System.Threading.Interlocked.Increment(ref s_CheckCount);
            if(maxDiff > 1e-11 * Math.Max(maxRef, 1e-300)) {
                long f = System.Threading.Interlocked.Increment(ref s_CheckFail);
                if(f <= 200) {
                    var sw = new System.Text.StringBuilder();
                    sw.Append($"FASTPATH-MISMATCH #{f}/{n} path={path} maxDiff={maxDiff:g4} maxRef={maxRef:g4} at {iMax} ref={(iMax>=0?r[iMax]:0):g6} got={(iMax>=0?t[iMax]:0):g6} alpha={alpha} beta={beta} len=[");
                    for(int d = 0; d < DT; d++) sw.Append(lenRun[d] + " ");
                    sw.Append("] ct=["); for(int d = 0; d < DT; d++) sw.Append(ct[d] + " ");
                    sw.Append("] ca=["); for(int d = 0; d < DT; d++) sw.Append(ca[d] + " ");
                    sw.Append("] cb=["); for(int d = 0; d < DT; d++) sw.Append(cb[d] + " ");
                    sw.Append($"] K={lenSum[0]} sa={sa[0]} sb={sb[0]}");
                    Console.Error.WriteLine(sw.ToString());
                }
            }
            if((n & ((1 << 18) - 1)) == 0)
                Console.Error.WriteLine($"FASTPATH-CHECK {n} calls compared, {s_CheckFail} mismatches");
        }

        // Selection order: direct GEMM, direct GEMV, packed GEMM, SIMD, fallback.
        static unsafe bool TryMultiplyBlas(int runDim, int sumDim,
            double* t, double* a, double* b, int* lengths,
            int* ct, int* ca, int* cb, int* sums, int* sa, int* sb,
            double alpha, double beta, bool directAllowed) {
            // In-place Multiply has established traversal semantics when a
            // shallow view aliases the destination. Preserve that reference
            // behavior instead of silently changing it to snapshot semantics.
            if(!directAllowed || s_FastPathOff)
                return false;
            if(sumDim != 1)
                { bool ok = TryMultiplySimd(runDim, sumDim, t, a, b, lengths, ct, ca, cb, sums, sa, sb, alpha, beta); if(ok) s_LastPath = "Simd"; return ok; }
            if(sums[0] <= 0 || alpha == 0.0)
                return false;
            if(TryDirectGemmRank4(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta)) { s_LastPath = "TryDirectGemmRank4"; return true; }
            if(directAllowed && TryDirectGemm(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta)) { s_LastPath = "TryDirectGemm"; return true; }
            if(directAllowed && TryDirectGemv(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta)) { s_LastPath = "TryDirectGemv"; return true; }
            if(s_FastPathGemmOnly) return false;
            if(TryPackedGemm(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta)) { s_LastPath = "TryPackedGemm"; return true; }
            { bool ok = TryMultiplySimd(runDim, sumDim, t, a, b, lengths, ct, ca, cb, sums, sa, sb, alpha, beta); if(ok) s_LastPath = "Simd"; return ok; }
        }

        // Handles contractions such as T[i,k,n,a] = sum_b A[i,k,a,b] B[k,n,b].
        // The two operand-exclusive axes form the GEMM dimensions; all other
        // running axes are independent batches.  The strict stride checks keep
        // this path limited to the row/column-major layouts BLAS can consume
        // without packing.
        static unsafe bool TryDirectGemmRank4(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta) {
            if(runDim != 4 || k < 8) return false;
            int p = -1, q = -1;
            for(int d = 0; d < runDim; d++) {
                if(ca[d] == 0 && cb[d] != 0) { if(p >= 0) return false; p = d; }
                if(ca[d] != 0 && cb[d] == 0) { if(q >= 0) return false; q = d; }
            }
            if(p < 0 || q < 0 || p == q) return false;
            int m = lengths[p], n = lengths[q];
            if(m < 4 || n < 4 || (long)m * n * k < 4096) return false;
            if(ct[q] != 1 || ct[p] < n || sa[0] != 1 || sb[0] != 1) return false;
            if(ca[q] < k || cb[p] < k) return false;

            int* batch = stackalloc int[2];
            int batchCount = 1;
            int batchDimCount = 0;
            for(int d = 0; d < runDim; d++) {
                if(d == p || d == q) continue;
                if(ca[d] == 0 || cb[d] == 0 || ct[d] == 0) return false;
                batch[batchDimCount++] = d;
                batchCount *= lengths[d];
            }
            int ta = 't', tb = 'n', M = n, N = m, K = k;
            int lda = ca[q], ldb = cb[p], ldc = ct[p];
            for(int z = 0; z < batchCount; z++) {
                int rem = z;
                int ot = 0, oa = 0, ob = 0;
                for(int j = 0; j < batchDimCount; j++) {
                    int d = batch[j], idx = rem % lengths[d]; rem /= lengths[d];
                    ot += idx * ct[d]; oa += idx * ca[d]; ob += idx * cb[d];
                }
                double* tp = t + ot, ap = a + oa, bp = b + ob;
                double nativeBeta = beta;
                if(beta == 0.0) {
                    for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) tp[i * ldc + j] *= beta;
                    nativeBeta = 1.0;
                }
                BLAS.m_seq_BLAS.DGEMM(ref ta, ref tb, ref M, ref N, ref K,
                    ref alpha, ap, ref lda, bp, ref ldb, ref nativeBeta, tp, ref ldc);
            }
            return true;
        }

        // Fast path for the common transformed layout: the transformed index is
        // the outer run index and each mapped slice is an ordinary matrix
        // product.  The generated transformed kernels currently bypass the
        // selector above, so handle this case before dispatching to them.
        static unsafe bool TryMultiplyWTrafoGemm(int runDim, int sumDim,
            double* t, double* a, double* b, int* lengths,
            int* ct, int* ca, int* cb, int* sums, int* sa, int* sb,
            double alpha, double beta, int* trfT, int* trfA, int* trfB,
            int trfT0Sw, int trfASw, int trfBSw,
            int preA, int cycleA, int postA, int preB, int cycleB, int postB,
            bool directAllowed) {
            if(!directAllowed || s_FastPathOff || runDim < 3 || sumDim != 1 || lengths[0] < 1 || sums[0] < 8)
                return false;
            int innerRun = runDim - 1;
            if(innerRun != 2) return false;
            // Validate every mapped slice before modifying the destination.
            for(int i = 0; i < lengths[0]; i++) {
                int ti = (trfT[i] + 0) * trfT0Sw;
                int ai = (trfA[i * cycleA + preA] + postA) * trfASw;
                int bi = (trfB[i * cycleB + preB] + postB) * trfBSw;
                if(ti < 0 || ai < 0 || bi < 0) return false;
                double* ts = t + ct[0] * (i * (1 - trfT0Sw) + ti);
                double* ap = a + ca[0] * (i * (1 - trfASw) + ai);
                double* bp = b + cb[0] * (i * (1 - trfBSw) + bi);
                if(!CanDirectGemm(innerRun, ap, bp, lengths + 1, ct + 1, ca + 1, cb + 1, sums[0], sa, sb))
                    return false;
            }
            for(int i = 0; i < lengths[0]; i++) {
                int ti = trfT[i] * trfT0Sw;
                int ai = (trfA[i * cycleA + preA] + postA) * trfASw;
                int bi = (trfB[i * cycleB + preB] + postB) * trfBSw;
                double* ts = t + ct[0] * (i * (1 - trfT0Sw) + ti);
                double* ap = a + ca[0] * (i * (1 - trfASw) + ai);
                double* bp = b + cb[0] * (i * (1 - trfBSw) + bi);
                TryDirectGemm(innerRun, ts, ap, bp, lengths + 1, ct + 1, ca + 1, cb + 1, sums[0], sa, sb, alpha, beta);
            }
            return true;
        }

        static unsafe bool CanDirectGemm(int runDim, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb) {
            if(runDim != 2) return false;
            int m = lengths[0], n = lengths[1];
            if(m < 8 || n < 8 || k < 8 || (long)m * n * k < 4096) return false;
            if(ca[1] != 0 || cb[0] != 0 || ct[1] != 1 || ct[0] < n) return false;
            if(!((sa[0] == 1 && ca[0] >= k) || (ca[0] == 1 && sa[0] >= m))) return false;
            if(!((cb[1] == 1 && sb[0] >= n) || (sb[0] == 1 && cb[1] >= k))) return false;
            return true;
        }

        static unsafe bool TryDirectGemm(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta) {
            if(runDim != 2 && runDim != 3) return false;
            int row = runDim - 2, col = row + 1;
            int m = lengths[row], n = lengths[col];
            if(m < 8 || n < 8 || k < 8 || (long)m * n * k < 4096) return false;
            if(ca[col] != 0 || cb[row] != 0 || ct[col] != 1 || ct[row] < n) return false;
            int transB, ldb;
            if(sa[0] == 1 && ca[row] >= k) { transB = 'n'; ldb = ca[row]; }
            else if(ca[row] == 1 && sa[0] >= m) { transB = 't'; ldb = sa[0]; }
            else return false;
            int transA, lda;
            if(cb[col] == 1 && sb[0] >= n) { transA = 'n'; lda = sb[0]; }
            else if(sb[0] == 1 && cb[col] >= k) { transA = 't'; lda = cb[col]; }
            else return false;
            int batches = row == 0 ? 1 : lengths[0], ldc = ct[row];
            for(int q = 0; q < batches; q++) {
                double nativeBeta = beta;
                if(beta == 0.0) {
                    for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) t[i * ldc + j] *= beta;
                    nativeBeta = 1.0;
                }
                BLAS.m_seq_BLAS.DGEMM(ref transA, ref transB, ref n, ref m, ref k,
                    ref alpha, b, ref lda, a, ref ldb, ref nativeBeta, t, ref ldc);
                if(row != 0) { t += ct[0]; a += ca[0]; b += cb[0]; }
            }
            return true;
        }

        // Pack arbitrary leading dimensions. The caller excludes aliases so the
        // established in-place traversal semantics remain unchanged.
        static unsafe bool TryPackedGemm(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta) {
            if(runDim != 2 && runDim != 3) return false;
            int row = runDim - 2, col = row + 1, m = lengths[row], n = lengths[col];
            if(k < 1 || m < 4 || n < 4 || (long)m * n * k < 4096 || ca[col] != 0 || cb[row] != 0) return false;
            int batches = row == 0 ? 1 : lengths[0];
            double[] aa = new double[m * k], bb = new double[k * n], cc = new double[m * n];
            fixed(double* pa = aa, pb = bb, pc = cc) {
                for(int q = 0; q < batches; q++) {
                    for(int i = 0; i < m; i++) for(int l = 0; l < k; l++) aa[i * k + l] = a[i * ca[row] + l * sa[0]];
                    for(int l = 0; l < k; l++) for(int j = 0; j < n; j++) bb[l * n + j] = b[l * sb[0] + j * cb[col]];
                    for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) cc[i * n + j] = t[i * ct[row] + j * ct[col]];
                    int ta = 'n', tb = 'n', M = n, N = m, K = k, lda = n, ldb = k, ldc = n;
                    double nativeBeta = beta;
                    if(beta == 0.0) { for(int i = 0; i < m * n; i++) cc[i] *= beta; nativeBeta = 1.0; }
                    BLAS.m_seq_BLAS.DGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha,
                        pb, ref lda, pa, ref ldb, ref nativeBeta, pc, ref ldc);
                    for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) t[i * ct[row] + j * ct[col]] = cc[i * n + j];
                    if(row != 0) { t += ct[0]; a += ca[0]; b += cb[0]; }
                }
            }
            return true;
        }

        static unsafe bool TryDirectGemv(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta) {
            if(runDim != 1) return false;
            int m = lengths[0];
            if(m < 8 || k < 8 || (long)m * k < 512) return false;
            bool left = ca[0] != 0 && cb[0] == 0, right = cb[0] != 0 && ca[0] == 0;
            if(!left && !right) return false;
            if(ct[0] != 1) return false;
            // BoSSS is row-major while BLAS is column-major: a row-major
            // m-by-k matrix is presented as its k-by-m transpose.
            int trans = 't', M = k, N = m, lda = left ? ca[0] : cb[0], inc = 1;
            double nativeBeta = beta;
            if(beta == 0.0) {
                for(int i = 0; i < m; i++) t[i] *= beta;
                nativeBeta = 1.0;
            }
            if(left && sa[0] == 1 && sb[0] == 1) {
                BLAS.m_seq_BLAS.DGEMV(ref trans, ref M, ref N, ref alpha, a, ref lda, b, ref inc, ref nativeBeta, t, ref inc); return true;
            }
            if(right && sb[0] == 1 && sa[0] == 1) {
                BLAS.m_seq_BLAS.DGEMV(ref trans, ref M, ref N, ref alpha, b, ref lda, a, ref inc, ref nativeBeta, t, ref inc); return true;
            }
            return false;
        }

        // Vectorize across adjacent output entries. This targets short (2/3/4)
        // reductions where vectorizing the reduction itself has too few terms.
        static unsafe bool TryMultiplySimd(int runDim, int sumDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int* sums, int* sa, int* sb, double alpha, double beta) {
            if(sumDim != 1 || runDim < 1 || runDim > 4 || sums[0] > 4) return false;
            int lanes = SimdVector.Count, last = runDim - 1, count = lengths[last];
            if(count < lanes || ct[last] != 1 || ca[last] != 1 || (cb[last] != 0 && cb[last] != 1)) return false;
            int outer = 1; for(int d = 0; d < last; d++) outer *= lengths[d];
            int k = sums[0];
            double[] avStorage = new double[lanes];
            double[] bvStorage = new double[lanes];
            double[] result = new double[lanes];
            for(int q = 0; q < outer; q++) {
                // Flattening the outer dimensions into one stride is only valid
                // for contiguous arrays.  Multiply also accepts padded views;
                // derive each operand's actual outer offset instead.
                int rem = q, ot = 0, oa = 0, ob = 0;
                for(int d = 0; d < last; d++) {
                    int idx = rem % lengths[d];
                    rem /= lengths[d];
                    ot += idx * ct[d];
                    oa += idx * ca[d];
                    ob += idx * cb[d];
                }
                double* tq = t + ot;
                double* aq = a + oa;
                double* bq = b + ob;
                int vectorEnd = count - count % lanes;
                for(int i = 0; i < vectorEnd; i += lanes) {
                    SimdVector acc = SimdVector.Zero;
                    for(int l = 0; l < k; l++) {
                        for(int lane = 0; lane < lanes; lane++) {
                            avStorage[lane] = aq[i + lane + l * sa[0]];
                            bvStorage[lane] = cb[last] == 0 ? bq[l * sb[0]] : bq[i + lane + l * sb[0]];
                        }
                        SimdVector av = new SimdVector(avStorage);
                        SimdVector bv = new SimdVector(bvStorage);
                        acc += av * bv;
                    }
                    for(int lane = 0; lane < lanes; lane++) result[lane] = tq[i + lane];
                    acc = acc * new SimdVector(alpha) + new SimdVector(result) * new SimdVector(beta);
                    acc.CopyTo(result);
                    for(int lane = 0; lane < lanes; lane++) tq[i + lane] = result[lane];
                }
                // SIMD covers only complete vectors.  The scalar tail is part
                // of the same selected kernel; returning success without it
                // silently left valid output entries unchanged.
                for(int i = vectorEnd; i < count; i++) {
                    double sum = 0.0;
                    for(int l = 0; l < k; l++) {
                        double av = aq[i + l * sa[0]];
                        double bv = cb[last] == 0 ? bq[l * sb[0]] : bq[i + l * sb[0]];
                        sum += av * bv;
                    }
                    tq[i] = alpha * sum + beta * tq[i];
                }
            }
            return true;
        }
    }
}
