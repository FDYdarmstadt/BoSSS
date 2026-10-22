using System;
using System.Linq;
using System.Reflection;
using ilPSP.Utils;
using NUnit.Framework;

namespace ilPSP {

    /// <summary>
    /// Tests for the optional native BLAS extensions (<see cref="MKLBatch"/>, <see cref="MKLJit"/>, <see cref="BLIS"/>, <see cref="LIBXSMM"/>)
    /// and the Multiply paths which use them, against naive reference loops.
    /// Unavailable capabilities are reported and skipped; the availability reason must then be non-empty.
    /// </summary>
    [TestFixture]
    public unsafe class BLASExtensions_Tests {

        static double[] Rnd(Random r, int n) => Enumerable.Range(0, n).Select(_ => r.NextDouble() - 0.5).ToArray();

        static void SkipIfUnavailable(bool avail, string reason, string name) {
            if(!avail) {
                Assert.IsFalse(string.IsNullOrEmpty(reason));
                Assert.Ignore(name + " unavailable: " + reason);
            }
        }

        /// <summary>
        /// column-major reference: C = alpha*op(A)*op(B) + beta*C
        /// </summary>
        static void RefGemm(bool ta, bool tb, int m, int n, int k, double alpha, double[] A, int oA, int lda, double[] B, int oB, int ldb, double beta, double[] C, int oC, int ldc) {
            for(int i = 0; i < m; i++)
                for(int j = 0; j < n; j++) {
                    double s = 0;
                    for(int l = 0; l < k; l++)
                        s += (ta ? A[oA + l + i * lda] : A[oA + i + l * lda]) * (tb ? B[oB + j + l * ldb] : B[oB + l + j * ldb]);
                    C[oC + i + j * ldc] = alpha * s + beta * C[oC + i + j * ldc];
                }
        }

        static void AssertClose(double[] expected, double[] actual, string msg) {
            Assert.AreEqual(expected.Length, actual.Length);
            for(int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i], 1e-12 * (1 + Math.Abs(expected[i])), msg + " at " + i);
        }

        [Test]
        public static void Availability() {
            Console.WriteLine($"MKL batch: {(MKLBatch.Available ? "available" : MKLBatch.UnavailableReason)}");
            Console.WriteLine($"MKL JIT:   {(MKLJit.Available ? "available" : MKLJit.UnavailableReason)}");
            Console.WriteLine($"BLIS:      {(BLIS.Available ? "available" : BLIS.UnavailableReason)}");
            Console.WriteLine($"LIBXSMM:   {(LIBXSMM.Available ? "available, " + LIBXSMM.LibraryPath : LIBXSMM.UnavailableReason)}");
            foreach(var (a, r) in new[] { (MKLBatch.Available, MKLBatch.UnavailableReason), (MKLJit.Available, MKLJit.UnavailableReason), (BLIS.Available, BLIS.UnavailableReason), (LIBXSMM.Available, LIBXSMM.UnavailableReason) })
                Assert.IsTrue(a == (r == null));
        }

        [Test]
        public static void MklBatch([Values(false, true)] bool ta, [Values(false, true)] bool tb) {
            SkipIfUnavailable(MKLBatch.Available, MKLBatch.UnavailableReason, "MKL batch");
            var rnd = new Random(11);
            int m = 7, n = 5, k = 9, batch = 6;
            int lda = (ta ? k : m) + 1, ldb = (tb ? n : k) + 2, ldc = m + 3;
            int sA = lda * (ta ? m : k) + 1, sB = ldb * (tb ? k : n), sC = ldc * n + 2;
            double[] A = Rnd(rnd, sA * batch), B = Rnd(rnd, sB * batch), C = Rnd(rnd, sC * batch);
            double alpha = 1.25, beta = -0.5;

            // strided
            double[] Cref = C.CloneAs(), C1 = C.CloneAs();
            for(int q = 0; q < batch; q++)
                RefGemm(ta, tb, m, n, k, alpha, A, q * sA, lda, B, q * sB, ldb, beta, Cref, q * sC, ldc);
            fixed(double* pA = A, pB = B, pC = C1) {
                MKLBatch.DgemmBatchStrided(MKLBatch.ColMajor, ta ? MKLBatch.Trans : MKLBatch.NoTrans, tb ? MKLBatch.Trans : MKLBatch.NoTrans,
                    m, n, k, alpha, pA, lda, sA, pB, ldb, sB, beta, pC, ldc, sC, batch);
            }
            AssertClose(Cref, C1, "dgemm_batch_strided");

            // grouped: one group, entries in reversed order (pointer arrays)
            double[] C2 = C.CloneAs();
            fixed(double* pA = A, pB = B, pC = C2) {
                double** Ap = stackalloc double*[batch], Bp = stackalloc double*[batch], Cp = stackalloc double*[batch];
                for(int q = 0; q < batch; q++) { int r = batch - 1 - q; Ap[q] = pA + r * sA; Bp[q] = pB + r * sB; Cp[q] = pC + r * sC; }
                int _ta = ta ? MKLBatch.Trans : MKLBatch.NoTrans, _tb = tb ? MKLBatch.Trans : MKLBatch.NoTrans, _m = m, _n = n, _k = k, _lda = lda, _ldb = ldb, _ldc = ldc, gs = batch;
                double _al = alpha, _be = beta;
                MKLBatch.DgemmBatch(MKLBatch.ColMajor, &_ta, &_tb, &_m, &_n, &_k, &_al, Ap, &_lda, Bp, &_ldb, &_be, Cp, &_ldc, 1, &gs);
            }
            AssertClose(Cref, C2, "dgemm_batch");

            // gemv strided: y = alpha*op(A)*x + beta*y, A is m x k with lda
            int lx = ta ? m : k, ly = ta ? k : m, incx = 2, incy = 3, sx = lx * incx + 1, sy = ly * incy;
            double[] X = Rnd(rnd, sx * batch), Y = Rnd(rnd, sy * batch), Yref = Y.CloneAs(), Y1 = Y.CloneAs();
            int ldA2 = m + 1, sA2 = ldA2 * k;
            double[] A2 = Rnd(rnd, sA2 * batch);
            for(int q = 0; q < batch; q++)
                for(int i = 0; i < ly; i++) {
                    double s = 0;
                    for(int l = 0; l < lx; l++)
                        s += (ta ? A2[q * sA2 + l + i * ldA2] : A2[q * sA2 + i + l * ldA2]) * X[q * sx + l * incx];
                    Yref[q * sy + i * incy] = alpha * s + beta * Y[q * sy + i * incy];
                }
            fixed(double* pA = A2, px = X, py = Y1) {
                MKLBatch.DgemvBatchStrided(MKLBatch.ColMajor, ta ? MKLBatch.Trans : MKLBatch.NoTrans, m, k, alpha, pA, ldA2, sA2, px, incx, sx, beta, py, incy, sy, batch);
            }
            AssertClose(Yref, Y1, "dgemv_batch_strided");
        }

        [Test]
        public static void MklJit([Values(false, true)] bool ta, [Values(false, true)] bool tb, [Values(0.0, 1.0, -0.75)] double beta) {
            SkipIfUnavailable(MKLJit.Available, MKLJit.UnavailableReason, "MKL JIT");
            var rnd = new Random(12);
            int m = 12, n = 10, k = 7;
            int lda = (ta ? k : m) + 1, ldb = (tb ? n : k), ldc = m + 2;
            double[] A = Rnd(rnd, lda * (ta ? m : k)), B = Rnd(rnd, ldb * (tb ? k : n)), C = Rnd(rnd, ldc * n);
            double alpha = 0.5;
            double[] Cref = C.CloneAs(), C1 = C.CloneAs();
            RefGemm(ta, tb, m, n, k, alpha, A, 0, lda, B, 0, ldb, beta, Cref, 0, ldc);
            bool ok = MKLJit.TryGetKernel(MKLBatch.ColMajor, ta ? MKLBatch.Trans : MKLBatch.NoTrans, tb ? MKLBatch.Trans : MKLBatch.NoTrans, m, n, k, alpha, lda, ldb, beta, ldc, out var kern);
            if(!ok)
                Assert.Ignore("MKL declined JIT for this shape");
            int nk = MKLJit.NumKernels;
            Assert.IsTrue(MKLJit.TryGetKernel(MKLBatch.ColMajor, ta ? MKLBatch.Trans : MKLBatch.NoTrans, tb ? MKLBatch.Trans : MKLBatch.NoTrans, m, n, k, alpha, lda, ldb, beta, ldc, out var kern2));
            Assert.AreEqual(nk, MKLJit.NumKernels, "cache hit must not create a kernel");
            fixed(double* pA = A, pB = B, pC = C1) {
                kern2.Invoke(pA, pB, pC);
            }
            AssertClose(Cref, C1, "MKL JIT");

            // concurrent use of one kernel
            ilPSP.Environment.ParallelFor(0, 16, delegate (int i) {
                double[] Ci = C.CloneAs();
                fixed(double* pA = A, pB = B, pC = Ci) {
                    kern.Invoke(pA, pB, pC);
                }
                for(int j = 0; j < Ci.Length; j++)
                    if(Math.Abs(Ci[j] - Cref[j]) > 1e-12 * (1 + Math.Abs(Cref[j])))
                        throw new Exception("MKL JIT concurrent mismatch");
            });
        }

        [Test]
        public static void Blis([Values(false, true)] bool ta, [Values(false, true)] bool tb) {
            SkipIfUnavailable(BLIS.Available, BLIS.UnavailableReason, "BLIS");
            var rnd = new Random(13);
            int m = 11, n = 9, k = 13;
            // general strides: A stored row-major with padding, B with non-unit strides in both directions, C row-major
            int rsA = 3 * k + 1, csA = 3, rsB = 2, csB = 2 * k + 5, rsC = n + 4, csC = 1;
            // opA(i,l) = A(i,l) if !ta, else A^T: the stored matrix has (rows, cols) = ta ? (k, m) : (m, k)
            int rowsA = ta ? k : m, colsA = ta ? m : k, rowsB = tb ? n : k, colsB = tb ? k : n;
            double[] A = Rnd(rnd, rowsA * rsA + colsA * csA), B = Rnd(rnd, rowsB * rsB + colsB * csB), C = Rnd(rnd, m * rsC + n);
            double alpha = -1.5, beta = 0.25;
            double[] Cref = C.CloneAs();
            for(int i = 0; i < m; i++)
                for(int j = 0; j < n; j++) {
                    double s = 0;
                    for(int l = 0; l < k; l++)
                        s += (ta ? A[l * rsA + i * csA] : A[i * rsA + l * csA]) * (tb ? B[j * rsB + l * csB] : B[l * rsB + j * csB]);
                    Cref[i * rsC + j * csC] = alpha * s + beta * C[i * rsC + j * csC];
                }
            double[] C1 = C.CloneAs();
            fixed(double* pA = A, pB = B, pC = C1) {
                BLIS.Dgemm(ta, tb, m, n, k, alpha, pA, rsA, csA, pB, rsB, csB, beta, pC, rsC, csC);
            }
            AssertClose(Cref, C1, "BLIS dgemm");

            // gemv: y = alpha*op(A)*x + beta*y, stored A is (rowsA x colsA)
            int lx = ta ? rowsA : colsA, ly = ta ? colsA : rowsA, incx = 3, incy = 2;
            double[] x = Rnd(rnd, lx * incx), y = Rnd(rnd, ly * incy), yref = y.CloneAs();
            for(int i = 0; i < ly; i++) {
                double s = 0;
                for(int l = 0; l < lx; l++) s += (ta ? A[l * rsA + i * csA] : A[i * rsA + l * csA]) * x[l * incx];
                yref[i * incy] = alpha * s + beta * y[i * incy];
            }
            fixed(double* pA = A, px = x, py = y) {
                BLIS.Dgemv(ta, rowsA, colsA, alpha, pA, rsA, csA, px, incx, beta, py, incy);
            }
            AssertClose(yref, y, "BLIS dgemv");
        }

        [Test]
        public static void Libxsmm([Values(true, false)] bool beta0) {
            SkipIfUnavailable(LIBXSMM.Available, LIBXSMM.UnavailableReason, "LIBXSMM");
            var rnd = new Random(14);
            int m = 9, n = 6, k = 11, lda = m + 1, ldb = k + 2, ldc = m + 3;
            double[] A = Rnd(rnd, lda * k), B = Rnd(rnd, ldb * n), C = Rnd(rnd, ldc * n), Cref = C.CloneAs(), C1 = C.CloneAs();
            RefGemm(false, false, m, n, k, 1.0, A, 0, lda, B, 0, ldb, beta0 ? 0.0 : 1.0, Cref, 0, ldc);
            if(!LIBXSMM.TryGetKernel(m, n, k, lda, ldb, ldc, beta0, out var kern))
                Assert.Ignore("LIBXSMM declined the shape");
            fixed(double* pA = A, pB = B, pC = C1) {
                kern.Invoke(pA, pB, pC);
            }
            AssertClose(Cref, C1, "LIBXSMM");
        }

        delegate bool BlisGettDel(int runDim, int sumDim, double* t, double* a, double* b,
            int* lengths, int* ct, int* ca, int* cb, int* sums, int* sa, int* sb, double alpha, double beta, long minMnk);

        static readonly BlisGettDel s_BlisGett = (BlisGettDel)typeof(MultidimensionalArray).GetMethod("TryBlisGett", BindingFlags.NonPublic | BindingFlags.Static).CreateDelegate(typeof(BlisGettDel));

        /// <summary>
        /// <see cref="MultidimensionalArray"/>`.TryBlisGett` on random layouts of T[i,j,k] = alpha*sum_{p,q} A[..] B[..] + beta*T,
        /// with operand index orders permuted, padding, M/N/batch classification and one or two summation loops.
        /// </summary>
        [Test]
        public static void BlisGettRandomLayouts() {
            SkipIfUnavailable(BLIS.Available, BLIS.UnavailableReason, "BLIS");
            var rnd = new Random(15);
            int tested = 0;
            for(int trial = 0; trial < 200; trial++) {
                int runDim = 3, sumDim = rnd.Next(1, 3);
                int[] len = { rnd.Next(1, 9), rnd.Next(4, 12), rnd.Next(4, 12) }, sums = { rnd.Next(4, 12), rnd.Next(1, 4) };
                // role of each running loop: 0 = batch (A and B), 1 = M (A only), 2 = N (B only); ensure at least one M and one N
                int[] role = new int[3];
                int perm = rnd.Next(3);
                role[perm] = 1; role[(perm + 1) % 3] = 2; role[(perm + 2) % 3] = rnd.Next(3);

                // strides: random index order and padding per operand
                int[] ct = new int[3], ca = new int[3], cb = new int[3], sa = new int[2], sb = new int[2];
                int Strides(int[] cRun, int[] cSum, bool[] runMember) {
                    var dims = Enumerable.Range(0, 3).Where(d => runMember[d]).Select(d => d).Concat(Enumerable.Range(0, sumDim).Select(s => -1 - s)).OrderBy(_ => rnd.Next()).ToList();
                    int stride = 1;
                    foreach(int d in dims) {
                        int l = d >= 0 ? len[d] : sums[-1 - d];
                        if(d >= 0) cRun[d] = stride; else cSum[-1 - d] = stride;
                        stride = stride * l + (rnd.NextDouble() < 0.3 ? rnd.Next(1, 3) : 0);
                    }
                    return stride + 1;
                }
                int sizeT = Strides(ct, new int[2], new[] { true, true, true });
                int sizeA = Strides(ca, sa, role.Select(r => r != 2).ToArray());
                int sizeB = Strides(cb, sb, role.Select(r => r != 1).ToArray());

                double[] T = Rnd(rnd, sizeT), A = Rnd(rnd, sizeA), B = Rnd(rnd, sizeB);
                double alpha = rnd.NextDouble() + 0.5, beta = new[] { 0.0, 1.0, -0.5 }[rnd.Next(3)];
                if(beta == 0.0) T[0] = double.NaN; // beta = 0 must keep NaN*0 = NaN semantics of the generated kernels
                double[] Tref = T.CloneAs();
                for(int i = 0; i < len[0]; i++)
                    for(int j = 0; j < len[1]; j++)
                        for(int k = 0; k < len[2]; k++) {
                            double s = 0;
                            for(int p = 0; p < sums[0]; p++)
                                for(int q = 0; q < (sumDim > 1 ? sums[1] : 1); q++)
                                    s += A[i * ca[0] + j * ca[1] + k * ca[2] + p * sa[0] + q * sa[1]] * B[i * cb[0] + j * cb[1] + k * cb[2] + p * sb[0] + q * sb[1]];
                            int o = i * ct[0] + j * ct[1] + k * ct[2];
                            Tref[o] = Tref[o] * beta + alpha * s;
                        }

                bool done;
                fixed(double* pT = T, pA = A, pB = B)
                fixed(int* pLen = len, pct = ct, pca = ca, pcb = cb, pSums = sums, psa = sa, psb = sb) {
                    done = s_BlisGett(runDim, sumDim, pT, pA, pB, pLen, pct, pca, pcb, pSums, psa, psb, alpha, beta, 1);
                }
                if(!done)
                    continue;
                tested++;
                for(int i = 0; i < T.Length; i++) {
                    if(double.IsNaN(Tref[i])) {
                        Assert.IsTrue(double.IsNaN(T[i]), $"trial {trial}: NaN expected at {i}");
                        continue;
                    }
                    Assert.AreEqual(Tref[i], T[i], 1e-12 * (1 + Math.Abs(Tref[i])), $"trial {trial} at {i}");
                }
            }
            Console.WriteLine($"TryBlisGett applied in {tested} of 200 random layouts");
            Assert.Greater(tested, 100);
        }

        /// <summary>
        /// Multiply through the public API with the MKL JIT path (on by default) against the generated kernels (JIT switched off, which
        /// selects the bitwise-identical SIMD path or the generated kernels for these sizes): equal up to round-off, including
        /// transformed programs, beta = 0 with NaN in T, and repeated calls (cached kernels).
        /// </summary>
        [Test]
        public static void MultiplyJitPrograms() {
            SkipIfUnavailable(MKLJit.Available, MKLJit.UnavailableReason, "MKL JIT");
            var fJit = typeof(MultidimensionalArray).GetField("s_MklJitGemm", BindingFlags.NonPublic | BindingFlags.Static);
            object oldJit = fJit.GetValue(null);
            var rnd = new Random(17);
            var cases = new (string T, string A, string B, int[] lenT, int lenP)[] {
                ("ij", "ip", "pj", new[] { 13, 10 }, 16),
                ("ij", "ip", "pj", new[] { 5, 5 }, 5),
                ("ij", "ip", "jp", new[] { 20, 49 }, 6),
                ("ij", "pi", "pj", new[] { 7, 9 }, 11),
                ("ijk", "ijp", "ipk", new[] { 4, 15, 15 }, 15),
                ("ijk", "ijp", "pk", new[] { 4, 6, 6 }, 6),
                ("ijk", "T(i)pj", "ipk", new[] { 20, 21, 21 }, 9),
                ("ijk", "T(i)pk", "ipj", new[] { 7, 10, 6 }, 5),
            };
            int jitBefore = MKLJit.NumKernels;
            try {
                foreach(var c in cases) {
                    foreach(double beta in new[] { 0.0, 1.0, 0.25 }) {
                        var idx = new System.Collections.Generic.Dictionary<char, int>();
                        for(int d = 0; d < c.T.Length; d++) idx[c.T[d]] = c.lenT[d];
                        idx['p'] = c.lenP;
                        int[] Shape(string s) => s.Replace("T(", "").Replace(")", "").Select(ch => idx[ch]).ToArray();
                        var T = MultidimensionalArray.Create(Shape(c.T));
                        var A = MultidimensionalArray.Create(Shape(c.A));
                        var B = MultidimensionalArray.Create(Shape(c.B));
                        for(int i = 0; i < T.Storage.Length; i++) T.Storage[i] = rnd.NextDouble() - .5;
                        for(int i = 0; i < A.Storage.Length; i++) A.Storage[i] = rnd.NextDouble() - .5;
                        for(int i = 0; i < B.Storage.Length; i++) B.Storage[i] = rnd.NextDouble() - .5;
                        if(beta == 0.0) T.Storage[1] = double.NaN; // NaN*0 = NaN as in the generated kernels
                        bool trafo = c.A.Contains("T(") || c.B.Contains("T(");
                        int[] trf = trafo ? Enumerable.Range(0, c.lenT[0]).Select(i => (i * 5 + 3) % c.lenT[0]).ToArray() : null;
                        var mp = MultidimensionalArray.MultiplyProgram.Compile(c.T, c.A, c.B, trafo);

                        fJit.SetValue(null, false);
                        var Tr = T.CloneAs();
                        Tr.Multiply(0.75, A, B, beta, ref mp, trf);
                        fJit.SetValue(null, true);
                        for(int rep = 0; rep < 2; rep++) {
                            var Tj = T.CloneAs();
                            Tj.Multiply(0.75, A, B, beta, ref mp, trf);
                            for(int i = 0; i < Tr.Storage.Length; i++) {
                                if(double.IsNaN(Tr.Storage[i])) {
                                    Assert.IsTrue(double.IsNaN(Tj.Storage[i]), $"{c.T},{c.A},{c.B} beta={beta}: NaN expected at {i}");
                                    continue;
                                }
                                Assert.AreEqual(Tr.Storage[i], Tj.Storage[i], 1e-13 * (1 + Math.Abs(Tr.Storage[i])), $"{c.T},{c.A},{c.B} beta={beta} at {i}");
                            }
                        }
                    }
                }
            } finally {
                fJit.SetValue(null, oldJit);
            }
            Console.WriteLine($"MKL JIT kernels: {jitBefore} -> {MKLJit.NumKernels}");
            Assert.Greater(MKLJit.NumKernels, jitBefore, "the JIT path was not used");
        }

        /// <summary>
        /// Multiply through the public API with the BLIS path and the MKL batch path enabled, against the reference with all fast paths
        /// of this phase disabled
        /// </summary>
        [Test]
        public static void MultiplyWithNativePaths([Values("ijk,pjq,ipkq", "ijk,ijp,ipk", "ij,ip,pj", "ijk,ipk,ijp")] string prog) {
            if(!BLIS.Available && !MKLBatch.Available)
                Assert.Ignore("neither BLIS nor MKL batch available");
            var fBlis = typeof(MultidimensionalArray).GetField("s_BlisGett", BindingFlags.NonPublic | BindingFlags.Static);
            var fMin = typeof(MultidimensionalArray).GetField("s_BlisGettMinMnk", BindingFlags.NonPublic | BindingFlags.Static);
            var fBatch = typeof(MultidimensionalArray).GetField("s_MklBatchDirect", BindingFlags.NonPublic | BindingFlags.Static);
            object oldBlis = fBlis.GetValue(null), oldMin = fMin.GetValue(null), oldBatch = fBatch.GetValue(null);

            var rnd = new Random(16);
            // extents by index letter
            int L(char c) => c switch { 'i' => 7, 'j' => 24, 'k' => 20, 'p' => 30, 'q' => 3, _ => throw new ArgumentException() };
            string[] parts = prog.Split(',');
            MultidimensionalArray Make(string idx) {
                var M = MultidimensionalArray.Create(idx.Select(L).ToArray());
                M.ApplyAll(_ => rnd.NextDouble() - 0.5);
                return M;
            }
            var A = Make(parts[1]); var B = Make(parts[2]);
            var T0 = Make(parts[0]);
            try {
                foreach(double beta in new[] { 0.0, 1.0, 0.5 }) {
                    fBlis.SetValue(null, false); fBatch.SetValue(null, false);
                    var Tref = T0.CloneAs();
                    Tref.Multiply(1.5, A, B, beta, parts[0], parts[1], parts[2]);

                    fBlis.SetValue(null, true); fMin.SetValue(null, 1L); fBatch.SetValue(null, true);
                    var T = T0.CloneAs();
                    T.Multiply(1.5, A, B, beta, parts[0], parts[1], parts[2]);

                    double[] t1 = T.ResizeShallow(T.Length).To1DArray(), r1 = Tref.ResizeShallow(Tref.Length).To1DArray();
                    double err = t1.Zip(r1, (u, v) => Math.Abs(u - v)).Max();
                    double nrm = r1.Max(Math.Abs);
                    Assert.LessOrEqual(err, 1e-12 * (1 + nrm), $"{prog} beta={beta}");
                }
            } finally {
                fBlis.SetValue(null, oldBlis); fMin.SetValue(null, oldMin); fBatch.SetValue(null, oldBatch);
            }
        }
    }
}
