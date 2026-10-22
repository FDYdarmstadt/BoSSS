using System;
using System.Numerics;
using ilPSP.Utils;
using SimdVector = System.Numerics.Vector<double>;

namespace ilPSP {
    public partial class MultidimensionalArray {
        // ==========================================================================================
        // DIAGNOSTIC SWITCH (inactive by default; kept on purpose, owner decision 2026-10-01)
        // ------------------------------------------------------------------------------------------
        // To activate, define the compilation symbol MULTIPLY_FASTPATH_DIAGNOSTICS for the ilPSP
        // project (e.g. <DefineConstants>$(DefineConstants);MULTIPLY_FASTPATH_DIAGNOSTICS</DefineConstants>
        // in ilPSP.csproj). Then, the environment variable BOSSS_MULTIPLY_FASTPATH selects:
        //   - "0":     all fast paths below are disabled, i.e. the generated kernels are used
        //              (behaviour before the tensor-product optimization);
        //   - "gemm":  only the direct (non-packing) DGEMM/DGEMV paths are used;
        //   - "nogett": the blocked SIMD path (MultidimensionalArray_MultiplyGett.cs) is disabled;
        //   - "census": all fast paths active; per (program, extents, strides, path) call counts and
        //              times are recorded and written at exit (see CensusRecord below);
        //   - "check": every fast-path call is additionally computed by the generated reference
        //              kernel into a copy of the destination storage and compared (slow!);
        //              mismatches are printed as FASTPATH-MISMATCH lines, and a summary per path
        //              is written at process exit to stderr and to FASTPATH-SUMMARY.txt;
        //   - unset:   all fast paths active (same as without the compilation symbol).
        // Used for the analysis in .agent/kummer/tensorproductoptimization/FINDINGS.md (2026-10-01).
        // ==========================================================================================
#if MULTIPLY_FASTPATH_DIAGNOSTICS
        static readonly string s_FastPathMode = (System.Environment.GetEnvironmentVariable("BOSSS_MULTIPLY_FASTPATH") ?? "all").ToLowerInvariant();
        static readonly bool s_FastPathOff = s_FastPathMode == "0";
        static readonly bool s_FastPathGemmOnly = s_FastPathMode == "gemm";
        // =nogett: blocked SIMD path (MultidimensionalArray_MultiplyGett.cs) off, all other fast paths on
        static readonly bool s_GettOff = s_FastPathMode == "nogett" || s_FastPathGemmOnly;
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

        // ------------------------------------------------------------------------------------------
        // census mode (BOSSS_MULTIPLY_FASTPATH=census; PLAN-4 Phase C): per (program, extents,
        // strides, trafo, alpha/beta class, aliasing, selected path) record the number of calls and
        // the accumulated time of the core Multiply (after program compilation). Callers are sampled
        // (stack trace on the 1st and every 4096th call of an entry). Fast paths stay active.
        // Written at process exit to the file given by BOSSS_MULTIPLY_CENSUS_FILE
        // (default MULTIPLY-CENSUS.tsv), one line per entry.
        // ------------------------------------------------------------------------------------------
        static readonly bool s_Census = s_FastPathMode == "census";

        sealed class CensusEntry {
            public long Count;
            public long Ticks;
            public string Desc;
            public readonly System.Collections.Concurrent.ConcurrentDictionary<string, int> Callers = new System.Collections.Concurrent.ConcurrentDictionary<string, int>();
        }

        static readonly System.Collections.Concurrent.ConcurrentDictionary<ulong, CensusEntry> s_CensusEntries = new System.Collections.Concurrent.ConcurrentDictionary<ulong, CensusEntry>();
        static int s_CensusExitHooked;

        static ulong CensusMix(ulong h, long v) {
            h ^= (ulong)v;
            h *= 0x100000001b3UL;
            return h;
        }

        static unsafe string CensusList(int* v, int n) {
            var sb = new System.Text.StringBuilder();
            for(int i = 0; i < n; i++) { if(i > 0) sb.Append(' '); sb.Append(v[i]); }
            return sb.ToString();
        }

        static string CensusCaller() {
            var st = new System.Diagnostics.StackTrace(2, false);
            var sb = new System.Text.StringBuilder();
            int found = 0;
            for(int i = 0; i < st.FrameCount && found < 3; i++) {
                var m = st.GetFrame(i).GetMethod();
                if(m == null || m.DeclaringType == typeof(MultidimensionalArray))
                    continue;
                if(found > 0) sb.Append(" < ");
                sb.Append((m.DeclaringType?.FullName ?? "?") + "." + m.Name);
                found++;
            }
            return sb.ToString();
        }

        static unsafe void CensusRecord(long t0, ref MultiplyProgram mp, int* lenRun, int* ct, int* ca, int* cb,
            int* lenSum, int* sa, int* sb, double alpha, double beta, bool noAlias) {
            long dt = System.Diagnostics.Stopwatch.GetTimestamp() - t0;
            int DT = mp.DT, NS = mp.NoOfSumCycles;
            int ab = (alpha == 1.0 ? 0 : alpha == 0.0 ? 1 : 2) + 3 * (beta == 0.0 ? 0 : beta == 1.0 ? 1 : 2);
            string path = s_LastPath ?? "?";

            // the transformed loop has been swapped to position 0 by the driver; report original order
            int k = mp.iTrafoIdx;
            int* l = stackalloc int[4 * DT];
            int* t_ = l + DT, a_ = t_ + DT, b_ = a_ + DT;
            for(int d = 0; d < DT; d++) {
                int s = (k > 0 && d == 0) ? k : (k > 0 && d == k) ? 0 : d;
                l[d] = lenRun[s]; t_[d] = ct[s]; a_[d] = ca[s]; b_[d] = cb[s];
            }

            ulong h = 0xcbf29ce484222325UL;
            h = CensusMix(h, mp.DT); h = CensusMix(h, mp.DA); h = CensusMix(h, mp.DB); h = CensusMix(h, NS);
            h = CensusMix(h, mp.iTrafoIdx); h = CensusMix(h, mp.TrfA0Sw); h = CensusMix(h, mp.TrfB0Sw);
            h = CensusMix(h, ab); h = CensusMix(h, noAlias ? 1 : 0); h = CensusMix(h, path.GetHashCode());
            for(int d = 0; d < DT; d++) { h = CensusMix(h, l[d]); h = CensusMix(h, t_[d]); h = CensusMix(h, a_[d]); h = CensusMix(h, b_[d]); }
            for(int d = 0; d < NS; d++) { h = CensusMix(h, lenSum[d]); h = CensusMix(h, sa[d]); h = CensusMix(h, sb[d]); }

            if(!s_CensusEntries.TryGetValue(h, out var e)) {
                var mpLoc = mp;
                string desc = string.Join("\t", mpLoc.Describe(), path, mp.iTrafoIdx >= 0 ? "1" : "0", ab.ToString(), noAlias ? "1" : "0",
                    CensusList(l, DT), CensusList(t_, DT), CensusList(a_, DT), CensusList(b_, DT),
                    CensusList(lenSum, NS), CensusList(sa, NS), CensusList(sb, NS));
                e = s_CensusEntries.GetOrAdd(h, _ => new CensusEntry() { Desc = desc });
                if(System.Threading.Interlocked.Exchange(ref s_CensusExitHooked, 1) == 0)
                    AppDomain.CurrentDomain.ProcessExit += (o, ea) => CensusWrite();
            }
            long cnt = System.Threading.Interlocked.Increment(ref e.Count);
            System.Threading.Interlocked.Add(ref e.Ticks, dt);
            if((cnt & 4095) == 1)
                e.Callers.AddOrUpdate(CensusCaller(), 1, (_, c) => c + 1);
        }

        static void CensusWrite() {
            string file = System.Environment.GetEnvironmentVariable("BOSSS_MULTIPLY_CENSUS_FILE");
            if(string.IsNullOrWhiteSpace(file)) file = "MULTIPLY-CENSUS.tsv";
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("seconds\tcalls\tprogram\tpath\ttrafo\talphabeta\tnoalias\tlenRun\tcycT\tcycA\tcycB\tlenSum\tsumA\tsumB\tcallers");
            double f = 1.0 / System.Diagnostics.Stopwatch.Frequency;
            long totCalls = 0; double totSec = 0;
            var all = new System.Collections.Generic.List<CensusEntry>(s_CensusEntries.Values);
            all.Sort((x, y) => y.Ticks.CompareTo(x.Ticks));
            foreach(var e in all) {
                totCalls += e.Count; totSec += e.Ticks * f;
                var callers = new System.Collections.Generic.List<string>();
                foreach(var kv in e.Callers) callers.Add(kv.Value + "x " + kv.Key);
                sb.Append((e.Ticks * f).ToString("g6", System.Globalization.CultureInfo.InvariantCulture)).Append('\t')
                  .Append(e.Count).Append('\t').Append(e.Desc).Append('\t').Append(string.Join(" | ", callers)).AppendLine();
            }
            Console.Error.WriteLine($"MULTIPLY-CENSUS: {all.Count} entries, {totCalls} calls, {totSec:g4} s -> {file}");
            try { System.IO.File.WriteAllText(file, sb.ToString()); } catch(Exception ex) { Console.Error.WriteLine("MULTIPLY-CENSUS: write failed: " + ex.Message); }
        }

#endif

        /// <summary>
        /// Records which fast path was selected (only in diagnostic builds) and returns true.
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        static bool FastPathSelected(string path) {
#if MULTIPLY_FASTPATH_DIAGNOSTICS
            s_LastPath = path;
#endif
            return true;
        }

        // Selection order: direct GEMM (rank 4, rank 2/3), direct GEMV, packed GEMM, fallback.
        // The short-reduction SIMD path (TryMultiplySimd) is currently disabled, see there.
        //
        // BLAS vs. blocked SIMD path (MultidimensionalArray_MultiplyGett.cs): if the latter is
        // applicable, DGEMM is only used for GEMMs with m*n*k above the following limits. Measured on pcmit36
        // (i7-6700, MKL, 1 thread; PLAN-4 Phase E, benchmarks/gett-f/results-crossover.csv): the SIMD path
        // runs at ~20-25 GFlop/s; direct DGEMM is 0-17% faster up to 20^3, 1.4x at 24^3, 1.6-2.4x for
        // 32^3...128^3. Below the limit, the SIMD path is preferred: comparable speed, and its results are
        // bitwise identical to the generated kernels. Packed DGEMM additionally copies all operands
        // (4.7x slower than the SIMD path for a census shape 20x49/6), hence the higher limit.
        // Keep the limits in one place, for re-tuning per machine/library (PLAN-5).
        internal const long GETT_PREFERRED_BELOW_MNK = 12000;
        internal const long GETT_PREFERRED_BELOW_MNK_PACKED = 64 * 64 * 64;
        static unsafe bool TryMultiplyBlas(int runDim, int sumDim,
            double* t, double* a, double* b, int* lengths,
            int* ct, int* ca, int* cb, int* sums, int* sa, int* sb,
            double alpha, double beta, bool directAllowed) {
            // In-place Multiply has established traversal semantics when a
            // shallow view aliases the destination. Preserve that reference
            // behavior instead of silently changing it to snapshot semantics.
            if(!directAllowed)
                return false;
#if MULTIPLY_FASTPATH_DIAGNOSTICS
            if(s_FastPathOff)
                return false;
#endif
            if(sumDim < 1 || sumDim > 2 || alpha == 0.0)
                return false;
            for(int s = 0; s < sumDim; s++)
                if(sums[s] <= 0) return false;
            // small GEMMs: MKL JIT is faster than both the blocked SIMD path and DGEMM (PLAN-5, RESULTS-5.md)
            if(sumDim == 1 && s_MklJitGemm && TryJitGemm(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta)) return FastPathSelected("TryJitGemm");
            long gettBelow = 0, gettBelowPacked = 0;
            {
                // every BLAS path needs m*n*k >= 4096 per matrix; only then, the SIMD path plan is evaluated
                long work = sums[0] * (sumDim > 1 ? (long)sums[1] : 1);
                for(int d = 0; d < runDim; d++) work *= lengths[d];
                if(work < 4096)
                    return false;
                if(GettWouldRun(runDim, sumDim, lengths, ct, ca, cb, sums, sa, sb, 1)) {
                    gettBelow = GETT_PREFERRED_BELOW_MNK;
                    gettBelowPacked = GETT_PREFERRED_BELOW_MNK_PACKED;
                }
            }
            if(sumDim != 1)
                goto general;
            if(TryDirectGemmRank4(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta, gettBelow)) return FastPathSelected("TryDirectGemmRank4");
            if(TryDirectGemm(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta, gettBelow)) return FastPathSelected("TryDirectGemm");
            if(TryDirectGemv(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta)) return FastPathSelected("TryDirectGemv");
#if MULTIPLY_FASTPATH_DIAGNOSTICS
            if(s_FastPathGemmOnly)
                return false;
#endif
            if(s_BlisGett && TryBlisGett(runDim, sumDim, t, a, b, lengths, ct, ca, cb, sums, sa, sb, alpha, beta, Math.Max(s_BlisGettMinMnk, gettBelow))) return FastPathSelected("TryBlisGett");
            if(TryPackedGemm(runDim, t, a, b, lengths, ct, ca, cb, sums[0], sa, sb, alpha, beta, gettBelowPacked)) return FastPathSelected("TryPackedGemm");
            goto packedGeneral;
        general:
#if MULTIPLY_FASTPATH_DIAGNOSTICS
            if(s_FastPathGemmOnly)
                return false;
#endif
            if(s_BlisGett && TryBlisGett(runDim, sumDim, t, a, b, lengths, ct, ca, cb, sums, sa, sb, alpha, beta, Math.Max(s_BlisGettMinMnk, gettBelow))) return FastPathSelected("TryBlisGett");
        packedGeneral:
            if(TryPackedGemmGeneral(runDim, sumDim, t, a, b, lengths, ct, ca, cb, sums, sa, sb, alpha, beta, gettBelowPacked)) return FastPathSelected("TryPackedGemmGeneral");
            return false;
        }

        /// <summary>
        /// Packed DGEMM for general contractions (PLAN-4 Phase E): the running loops are classified as
        /// M (varies in A only), N (B only) or batch (A and B); M, N and the summation loops (K) may each
        /// consist of several loops with arbitrary strides (e.g. T[i,j,k] = sum_{p,q} A[p,j,q] B[i,p,k,q]:
        /// M = j, N = (i,k), K = (p,q)). A, B and the T block are packed into contiguous buffers (offset
        /// tables, TempBuffer), one DGEMM per batch entry, T copied back.
        /// </summary>
        static unsafe bool TryPackedGemmGeneral(int runDim, int sumDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int* sums, int* sa, int* sb, double alpha, double beta, long gettBelow) {
            const int MAXD = 8;
            if(!s_PackedGemmGeneral || runDim > MAXD || sumDim < 1 || sumDim > 2) return false;
            int* dM = stackalloc int[MAXD]; int* dN = stackalloc int[MAXD]; int* dB = stackalloc int[MAXD];
            int nM = 0, nN = 0, nB = 0;
            long m = 1, n = 1, nb = 1, k = 1;
            for(int d = 0; d < runDim; d++) {
                if(lengths[d] <= 0) return false;
                if(lengths[d] == 1) continue;
                if(ca[d] != 0 && cb[d] == 0) { dM[nM++] = d; m *= lengths[d]; }
                else if(ca[d] == 0 && cb[d] != 0) { dN[nN++] = d; n *= lengths[d]; }
                else if(ca[d] != 0 && cb[d] != 0) { dB[nB++] = d; nb *= lengths[d]; }
                else return false; // T-only (broadcast) loop
            }
            for(int s = 0; s < sumDim; s++) k *= sums[s];
            if(m < 8 || n < 8 || k < 8) return false;
            long mnk = m * n * k;
            if(mnk < GEMM_PACKED_GENERAL_MIN_MNK || mnk < gettBelow) return false;
            if(m * k + k * n + m * n > int.MaxValue / 2 || m + n + 2 * k > (1 << 24)) return false;
            if(!GettTNonOverlapping(runDim, lengths, ct)) return false; // T entries pairwise distinct

            int M = (int)m, N = (int)n, K = (int)k;
            // offset tables: A over M and K, B over K and N, T over M and N
            int[] offs = new int[2 * M + 2 * N + 2 * K];
            fixed(int* pOffs = offs) {
                int* oAM = pOffs, oTM = oAM + M, oBN = oTM + M, oTN = oBN + N, oAK = oTN + N, oBK = oAK + K;
                FillOffsets(nM, dM, lengths, ca, ct, oAM, oTM);
                FillOffsets(nN, dN, lengths, cb, ct, oBN, oTN);
                for(int l0 = 0, l = 0; l0 < sums[0]; l0++)
                    for(int l1 = 0; l1 < (sumDim > 1 ? sums[1] : 1); l1++, l++) {
                        oAK[l] = l0 * sa[0] + (sumDim > 1 ? l1 * sa[1] : 0);
                        oBK[l] = l0 * sb[0] + (sumDim > 1 ? l1 * sb[1] : 0);
                    }

                double[] scratch = TempBuffer.GetTempBuffer(out int iScratch, M * K + K * N + M * N);
                try {
                    fixed(double* pScratch = scratch) {
                        double* pa = pScratch, pb = pa + M * K, pc = pb + K * N;
                        int* idx = stackalloc int[MAXD];
                        for(int j = 0; j < nB; j++) idx[j] = 0;
                        for(long q = 0; q < nb; q++) {
                            int ot = 0, oa = 0, ob = 0;
                            for(int j = 0; j < nB; j++) { int d = dB[j]; ot += idx[j] * ct[d]; oa += idx[j] * ca[d]; ob += idx[j] * cb[d]; }
                            double* tq = t + ot, aq = a + oa, bq = b + ob;
                            for(int i = 0; i < M; i++) { double* ai = aq + oAM[i]; double* pai = pa + i * K; for(int l = 0; l < K; l++) pai[l] = ai[oAK[l]]; }
                            for(int l = 0; l < K; l++) { double* bl = bq + oBK[l]; double* pbl = pb + l * N; for(int j = 0; j < N; j++) pbl[j] = bl[oBN[j]]; }
                            for(int i = 0; i < M; i++) { double* ti = tq + oTM[i]; double* pci = pc + i * N; for(int j = 0; j < N; j++) pci[j] = ti[oTN[j]]; }
                            int ta = 'n', tb = 'n', lda = N, ldb = K, ldc = N;
                            double nativeBeta = beta;
                            if(beta == 0.0) { for(int i = 0; i < M * N; i++) pc[i] *= beta; nativeBeta = 1.0; } // keeps NaN/Inf propagation
                            BLAS.SeqDGEMM(ref ta, ref tb, ref N, ref M, ref K, ref alpha,
                                pb, ref lda, pa, ref ldb, ref nativeBeta, pc, ref ldc);
                            for(int i = 0; i < M; i++) { double* ti = tq + oTM[i]; double* pci = pc + i * N; for(int j = 0; j < N; j++) ti[oTN[j]] = pci[j]; }
                            // next batch entry (odometer, innermost = last)
                            for(int j = nB - 1; j >= 0; j--) { if(++idx[j] < lengths[dB[j]]) break; idx[j] = 0; }
                        }
                    }
                } finally {
                    TempBuffer.FreeTempBuffer(iScratch);
                }
            }
            return true;
        }

        /// <summary>
        /// minimal GEMM size (m*n*k per batch entry) for <see cref="TryPackedGemmGeneral"/>; below, the
        /// packing costs dominate (it is additionally subject to GETT_PREFERRED_BELOW_MNK_PACKED)
        /// </summary>
        internal const long GEMM_PACKED_GENERAL_MIN_MNK = 64 * 64 * 64;

        /// <summary>
        /// Switch for <see cref="TryPackedGemmGeneral"/>; ON (owner decision 2026-10-03): 1.7x-2.2x faster
        /// than the blocked SIMD path for the large multi-index contractions of the volume quadratures
        /// (e.g. `ijk,pjq,ipkq`, 3D p4); it changes the summation order (results differ at round-off level).
        /// </summary>
        internal static bool s_PackedGemmGeneral = true;

        /// <summary> offsets of all multi-indices over the loops dims[0..nd-1] (last loop fastest), for two operands </summary>
        static unsafe void FillOffsets(int nd, int* dims, int* lengths, int* c1, int* c2, int* o1, int* o2) {
            int cnt = 1;
            o1[0] = 0; o2[0] = 0;
            for(int jj = nd - 1; jj >= 0; jj--) {
                // prepend loop dims[jj] as the slower index of the ones already processed
                int d = dims[jj], L = lengths[d];
                for(int i = L - 1; i >= 0; i--)
                    for(int r = 0; r < cnt; r++) {
                        o1[i * cnt + r] = o1[r] + i * c1[d];
                        o2[i * cnt + r] = o2[r] + i * c2[d];
                    }
                cnt *= L;
            }
        }

        // Handles contractions such as T[i,k,n,a] = sum_b A[i,k,a,b] B[k,n,b].
        // The two operand-exclusive axes form the GEMM dimensions; all other
        // running axes are independent batches.  The strict stride checks keep
        // this path limited to the row/column-major layouts BLAS can consume
        // without packing.
        static unsafe bool TryDirectGemmRank4(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta, long gettBelow) {
            if(runDim != 4 || k < 8) return false;
            int p = -1, q = -1;
            for(int d = 0; d < runDim; d++) {
                if(ca[d] == 0 && cb[d] != 0) { if(p >= 0) return false; p = d; }
                if(ca[d] != 0 && cb[d] == 0) { if(q >= 0) return false; q = d; }
            }
            if(p < 0 || q < 0 || p == q) return false;
            int m = lengths[p], n = lengths[q];
            if(m < 4 || n < 4 || (long)m * n * k < 4096 || (long)m * n * k < gettBelow) return false;
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
                BLAS.SeqDGEMM(ref ta, ref tb, ref M, ref N, ref K,
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
            if(!directAllowed || runDim < 3 || sumDim != 1 || lengths[0] < 1 || sums[0] < 1)
                return false;
#if MULTIPLY_FASTPATH_DIAGNOSTICS
            if(s_FastPathOff)
                return false;
#endif
            int innerRun = runDim - 1;
            if(innerRun != 2) return false;
            // small slices: MKL JIT kernel (PLAN-5), all slices have the same shape and strides
            if(s_MklJitGemm && JitGemmKernel(lengths + 1, ct + 1, ca + 1, cb + 1, sums[0], sa, sb, alpha, beta, out var jitKern)) {
                for(int i = 0; i < lengths[0]; i++) {
                    if(trfT[i] * trfT0Sw < 0 || (trfA[i * cycleA + preA] + postA) * trfASw < 0 || (trfB[i * cycleB + preB] + postB) * trfBSw < 0)
                        return false; // validate every mapped slice before modifying the destination
                }
                if(beta == 0.0) {
                    // all slices first: scaling directly before the kernel of the same slice stalls the kernel's
                    // vector loads (store-to-load forwarding; 13x slower on AMD EPYC 9654, RESULTS-5.md)
                    for(int i = 0; i < lengths[0]; i++) {
                        double* ts = t + ct[0] * (i * (1 - trfT0Sw) + trfT[i] * trfT0Sw);
                        JitGemmScale(ts, lengths[1], lengths[2], ct[1], beta);
                    }
                }
                for(int i = 0; i < lengths[0]; i++) {
                    int ti = trfT[i] * trfT0Sw;
                    int ai = (trfA[i * cycleA + preA] + postA) * trfASw;
                    int bi = (trfB[i * cycleB + preB] + postB) * trfBSw;
                    double* ts = t + ct[0] * (i * (1 - trfT0Sw) + ti);
                    double* ap = a + ca[0] * (i * (1 - trfASw) + ai);
                    double* bp = b + cb[0] * (i * (1 - trfBSw) + bi);
                    jitKern.Invoke(bp, ap, ts);
                }
                return FastPathSelected("WTrafoJit");
            }
            if(sums[0] < 8) return false;
            // the blocked SIMD path is preferred for small GEMMs, see GETT_PREFERRED_BELOW_MNK
            if((long)lengths[1] * lengths[2] * sums[0] < GETT_PREFERRED_BELOW_MNK
                && GettWouldRun(innerRun, sumDim, lengths + 1, ct + 1, ca + 1, cb + 1, sums, sa, sb, lengths[0])
                && GettTNonOverlapping(runDim, lengths, ct))
                return false;
            // Validate every mapped slice before modifying the destination.
            for(int i = 0; i < lengths[0]; i++) {
                int ti = trfT[i] * trfT0Sw;
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
                TryDirectGemm(innerRun, ts, ap, bp, lengths + 1, ct + 1, ca + 1, cb + 1, sums[0], sa, sb, alpha, beta, 0);
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
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta, long gettBelow) {
            if(runDim != 2 && runDim != 3) return false;
            int row = runDim - 2, col = row + 1;
            int m = lengths[row], n = lengths[col];
            if(m < 8 || n < 8 || k < 8 || (long)m * n * k < 4096 || (long)m * n * k < gettBelow) return false;
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
            if(s_MklBatchDirect && row != 0
                && DirectGemmMklBatch(batches, transA, transB, n, m, k, alpha, b, lda, cb[0], a, ldb, ca[0], beta, t, ldc, ct[0], m, n))
                return true;
            for(int q = 0; q < batches; q++) {
                double nativeBeta = beta;
                if(beta == 0.0) {
                    for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) t[i * ldc + j] *= beta;
                    nativeBeta = 1.0;
                }
                BLAS.SeqDGEMM(ref transA, ref transB, ref n, ref m, ref k,
                    ref alpha, b, ref lda, a, ref ldb, ref nativeBeta, t, ref ldc);
                if(row != 0) { t += ct[0]; a += ca[0]; b += cb[0]; }
            }
            return true;
        }

        // Pack arbitrary leading dimensions. The caller excludes aliases so the
        // established in-place traversal semantics remain unchanged.
        static unsafe bool TryPackedGemm(int runDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta, long gettBelow) {
            if(runDim != 2 && runDim != 3) return false;
            int row = runDim - 2, col = row + 1, m = lengths[row], n = lengths[col];
            if(k < 1 || m < 4 || n < 4 || (long)m * n * k < 4096 || (long)m * n * k < gettBelow || ca[col] != 0 || cb[row] != 0) return false;
            int batches = row == 0 ? 1 : lengths[0];
            // one temporary buffer for the packed A (m x k), B (k x n) and T (m x n) blocks
            double[] scratch = TempBuffer.GetTempBuffer(out int iScratch, m * k + k * n + m * n);
            try {
                fixed(double* pScratch = scratch) {
                    double* pa = pScratch, pb = pa + m * k, pc = pb + k * n;
                    for(int q = 0; q < batches; q++) {
                        for(int i = 0; i < m; i++) for(int l = 0; l < k; l++) pa[i * k + l] = a[i * ca[row] + l * sa[0]];
                        for(int l = 0; l < k; l++) for(int j = 0; j < n; j++) pb[l * n + j] = b[l * sb[0] + j * cb[col]];
                        for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) pc[i * n + j] = t[i * ct[row] + j * ct[col]];
                        int ta = 'n', tb = 'n', M = n, N = m, K = k, lda = n, ldb = k, ldc = n;
                        double nativeBeta = beta;
                        if(beta == 0.0) { for(int i = 0; i < m * n; i++) pc[i] *= beta; nativeBeta = 1.0; }
                        BLAS.SeqDGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha,
                            pb, ref lda, pa, ref ldb, ref nativeBeta, pc, ref ldc);
                        for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) t[i * ct[row] + j * ct[col]] = pc[i * n + j];
                        if(row != 0) { t += ct[0]; a += ca[0]; b += cb[0]; }
                    }
                }
            } finally {
                TempBuffer.FreeTempBuffer(iScratch);
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
            if(sa[0] != 1 || sb[0] != 1) return false;
            // BoSSS is row-major while BLAS is column-major: a row-major
            // m-by-k matrix is presented as its k-by-m transpose.
            int trans = 't', M = k, N = m, lda = left ? ca[0] : cb[0], inc = 1;
            // BLAS requires lda >= max(1, M); otherwise, DGEMV rejects the call (parameter error)
            // and leaves T unchanged, although this path would report success.
            if(lda < k) return false;
            // all eligibility checks done; only now, the destination may be modified
            double nativeBeta = beta;
            if(beta == 0.0) {
                for(int i = 0; i < m; i++) t[i] *= beta;
                nativeBeta = 1.0;
            }
            double* mtx = left ? a : b;
            double* vec = left ? b : a;
            BLAS.SeqDGEMV(ref trans, ref M, ref N, ref alpha, mtx, ref lda, vec, ref inc, ref nativeBeta, t, ref inc);
            return true;
        }

        // Vectorize across adjacent output entries. This targets short (2/3/4)
        // reductions where vectorizing the reduction itself has too few terms.
        //
        // CURRENTLY NOT USED (disabled 2026-10-01, PLAN-3 A3): benchmarked against the generated
        // kernels for all accepted shapes (T[j,n] = sum_d A[j,d,n]*B[j,d] resp. *B[j,d,n],
        // J=64, D=2..4, N=4..64), this implementation was 3x to 7x SLOWER, because of the
        // element-wise copies through managed arrays (see
        // .agent/kummer/tensorproductoptimization/benchmarks/simd-a3/results.csv).
        // A prototype with direct vector loads was 1.1x-2.4x faster for N >= 16, but slower
        // for N = 7, 10 (scalar tail). A redesign is part of PLAN-4 (Phase F).
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
