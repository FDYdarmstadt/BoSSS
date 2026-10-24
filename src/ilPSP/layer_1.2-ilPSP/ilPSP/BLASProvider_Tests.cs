using System;
using System.Linq;
using ilPSP.Utils;
using NUnit.Framework;

namespace ilPSP {

    /// <summary>
    /// Tests for the MKL/OpenBLAS provider selection in <see cref="BLAS"/>: all wrapped BLAS operations for each available provider,
    /// the OpenBLAS thread-count policy at parallel-section boundaries, concurrent sequential tensor calls and the calibration.
    /// If OpenBLAS is unavailable, the tests check that MKL remains active and that the reason is reported.
    /// </summary>
    [TestFixture]
    public unsafe class BLASProvider_Tests {

        static BLASProvider[] AvailableProviders() {
            return BLAS.OpenBLASAvailable ? new[] { BLASProvider.MKL, BLASProvider.OpenBLAS } : new[] { BLASProvider.MKL };
        }

        static void WithProvider(BLASProvider p, Action a) {
            var old = BLAS.Provider;
            Assert.IsTrue(BLAS.TryActivateProvider(p, out string reason), reason);
            try {
                Assert.AreEqual(p, BLAS.Provider);
                Assert.AreEqual(p, BLAS.F77_BLAS.Provider);
                a();
            } finally {
                Assert.IsTrue(BLAS.TryActivateProvider(old, out reason), reason);
            }
        }

        static double[] Rnd(Random r, int n) => Enumerable.Range(0, n).Select(_ => r.NextDouble() - 0.5).ToArray();

        /// <summary>
        /// reports the provider state; fails if OpenBLAS is unavailable but MKL is not active
        /// </summary>
        [Test]
        public static void ProviderState() {
            Console.WriteLine($"BLAS provider: {BLAS.Provider} (requested: {BLAS.RequestedProvider}); {BLAS.ProviderSelectionReason}");
            Console.WriteLine($"MKL: {BLAS.MKLLibraryPath}");
            if(BLAS.OpenBLASAvailable) {
                Console.WriteLine($"OpenBLAS: {BLAS.OpenBLASLoadInfo}; {BLAS.OpenBLASConfig}; threads now {BLAS.OpenBLASNumThreads}");
            } else {
                Console.WriteLine($"OpenBLAS unavailable: {BLAS.OpenBLASUnavailableReason}");
                Assert.IsNotNull(BLAS.OpenBLASUnavailableReason);
                Assert.AreEqual(BLASProvider.MKL, BLAS.Provider);
                Assert.IsFalse(BLAS.TryActivateProvider(BLASProvider.OpenBLAS, out string reason));
                Assert.IsNotNull(reason);
                Assert.AreEqual(BLASProvider.MKL, BLAS.Provider);
            }
        }

        /// <summary>
        /// all BLAS operations which are wrapped by <see cref="UnsafeDBLAS"/>, against managed reference results
        /// </summary>
        [Test]
        public static void AllOperations([Values(BLASProvider.MKL, BLASProvider.OpenBLAS)] BLASProvider p) {
            if(!AvailableProviders().Contains(p)) {
                Assert.Ignore("OpenBLAS unavailable: " + BLAS.OpenBLASUnavailableReason);
            }
            WithProvider(p, delegate () {
                var rnd = new Random(1234);
                const double tol = 1e-12;
                int n = 37;

                // ddot, dnrm2
                double[] x = Rnd(rnd, n), y = Rnd(rnd, n);
                double dotRef = 0, nrmRef = 0;
                for(int i = 0; i < n; i++) { dotRef += x[i] * y[i]; nrmRef += x[i] * x[i]; }
                Assert.AreEqual(dotRef, BLAS.ddot(n, x, 1, y, 1), tol);
                Assert.AreEqual(Math.Sqrt(nrmRef), BLAS.dnrm2(n, x, 1), tol);

                // daxpy, dscal, dswap
                double[] y2 = y.CloneAs();
                BLAS.daxpy(n, 0.75, x, 1, y2, 1);
                for(int i = 0; i < n; i++) Assert.AreEqual(y[i] + 0.75 * x[i], y2[i], tol);
                BLAS.dscal(n, -2.0, y2, 1);
                for(int i = 0; i < n; i++) Assert.AreEqual(-2.0 * (y[i] + 0.75 * x[i]), y2[i], tol);
                double[] xs = x.CloneAs(), ys = y.CloneAs();
                BLAS.dswap(n, xs, 1, ys, 1);
                Assert.IsTrue(xs.SequenceEqual(y) && ys.SequenceEqual(x));

                // dgemm, all transpositions, alpha/beta (column-major)
                foreach(char ta in new[] { 'n', 't' }) {
                    foreach(char tb in new[] { 'n', 't' }) {
                        int M = 13, N = 7, K = 11;
                        double[] A = Rnd(rnd, M * K), B = Rnd(rnd, K * N), C = Rnd(rnd, M * N);
                        int lda = ta == 'n' ? M : K, ldb = tb == 'n' ? K : N;
                        double Aij(int i, int l) => ta == 'n' ? A[i + l * lda] : A[l + i * lda];
                        double Bij(int l, int j) => tb == 'n' ? B[l + j * ldb] : B[j + l * ldb];
                        double[] Cref = C.CloneAs();
                        for(int i = 0; i < M; i++)
                            for(int j = 0; j < N; j++) {
                                double s = 0;
                                for(int l = 0; l < K; l++) s += Aij(i, l) * Bij(l, j);
                                Cref[i + j * M] = 1.5 * s + 0.25 * C[i + j * M];
                            }
                        fixed(double* pA = A, pB = B, pC = C) {
                            BLAS.dgemm(ta, tb, M, N, K, 1.5, pA, lda, pB, ldb, 0.25, pC, M);
                        }
                        for(int i = 0; i < C.Length; i++) Assert.AreEqual(Cref[i], C[i], tol, $"dgemm {ta}{tb}");

                        // sequential tensor entry
                        double[] C3 = Rnd(new Random(5), M * N), C3ref = C3.CloneAs();
                        for(int i = 0; i < M; i++)
                            for(int j = 0; j < N; j++) {
                                double s = 0;
                                for(int l = 0; l < K; l++) s += Aij(i, l) * Bij(l, j);
                                C3ref[i + j * M] = -1.0 * s + 1.0 * C3[i + j * M];
                            }
                        int _ta = ta, _tb = tb, _M = M, _N = N, _K = K, _lda = lda, _ldb = ldb, ldc = M;
                        double alpha = -1.0, beta = 1.0;
                        fixed(double* pA = A, pB = B, pC = C3) {
                            BLAS.SeqDGEMM(ref _ta, ref _tb, ref _M, ref _N, ref _K, ref alpha, pA, ref _lda, pB, ref _ldb, ref beta, pC, ref ldc);
                        }
                        for(int i = 0; i < C3.Length; i++) Assert.AreEqual(C3ref[i], C3[i], tol, $"SeqDGEMM {ta}{tb}");
                    }
                }

                // dgemv, sgemv
                foreach(char t in new[] { 'n', 't' }) {
                    int M = 9, N = 14;
                    double[] A = Rnd(rnd, M * N), xv = Rnd(rnd, t == 'n' ? N : M), yv = Rnd(rnd, t == 'n' ? M : N);
                    double[] yref = yv.CloneAs();
                    for(int i = 0; i < yv.Length; i++) {
                        double s = 0;
                        for(int l = 0; l < xv.Length; l++) s += (t == 'n' ? A[i + l * M] : A[l + i * M]) * xv[l];
                        yref[i] = 2.0 * s - 1.0 * yv[i];
                    }
                    double[] y4 = yv.CloneAs();
                    fixed(double* pA = A, px = xv, py = yv, py4 = y4) {
                        BLAS.dgemv(t, M, N, 2.0, pA, M, px, 1, -1.0, py, 1);
                        int _t = t, _M = M, _N = N, lda = M, inc = 1;
                        double alpha = 2.0, beta = -1.0;
                        BLAS.SeqDGEMV(ref _t, ref _M, ref _N, ref alpha, pA, ref lda, px, ref inc, ref beta, py4, ref inc);
                    }
                    for(int i = 0; i < yv.Length; i++) Assert.AreEqual(yref[i], yv[i], tol, $"dgemv {t}");
                    for(int i = 0; i < yv.Length; i++) Assert.AreEqual(yref[i], y4[i], tol, $"SeqDGEMV {t}");

                    float[] Af = A.Select(v => (float)v).ToArray(), xf = xv.Select(v => (float)v).ToArray(), yf = yv.Select(_ => 0.0f).ToArray();
                    fixed(float* pA = Af, px = xf, py = yf) {
                        BLAS.sgemv(t, M, N, 1.0f, pA, M, px, 1, 0.0f, py, 1);
                    }
                    for(int i = 0; i < yf.Length; i++) {
                        double s = 0;
                        for(int l = 0; l < xf.Length; l++) s += (t == 'n' ? Af[i + l * M] : Af[l + i * M]) * xf[l];
                        Assert.AreEqual(s, yf[i], 1e-5, $"sgemv {t}");
                    }
                }

                // sgemm, via the raw interface
                {
                    int M = 5, N = 6, K = 4, ta = 'n', tb = 'n';
                    float alpha = 1.0f, beta = 0.0f;
                    float[] A = Rnd(rnd, M * K).Select(v => (float)v).ToArray(), B = Rnd(rnd, K * N).Select(v => (float)v).ToArray(), C = new float[M * N];
                    fixed(float* pA = A, pB = B, pC = C) {
                        BLAS.F77_BLAS.SGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha, pA, ref M, pB, ref K, ref beta, pC, ref M);
                    }
                    for(int i = 0; i < M; i++)
                        for(int j = 0; j < N; j++) {
                            double s = 0;
                            for(int l = 0; l < K; l++) s += A[i + l * M] * B[l + j * K];
                            Assert.AreEqual(s, C[i + j * M], 1e-5, "sgemm");
                        }
                }
            });
        }

        /// <summary>
        /// OpenBLAS: thread count is 1 within a parallel section, the budget outside (if threading is enabled);
        /// concurrent sequential DGEMMs within the section give the reference result.
        /// </summary>
        [Test]
        public static void OpenBLASThreadPolicy() {
            if(!BLAS.OpenBLASAvailable)
                Assert.Ignore("OpenBLAS unavailable: " + BLAS.OpenBLASUnavailableReason);

            WithProvider(BLASProvider.OpenBLAS, delegate () {
                int outside = BLAS.OpenBLASNumThreads;
                Console.WriteLine($"NumThreads {ilPSP.Environment.NumThreads}, MaxNumOpenMPthreads {ilPSP.Environment.MaxNumOpenMPthreads}, OpenMP enabled {ilPSP.Environment.OpenMPenabled}, OpenBLAS threads outside {outside}");
                if(ilPSP.Environment.OpenMPenabled)
                    Assert.AreEqual(Math.Max(1, Math.Min(ilPSP.Environment.MaxNumOpenMPthreads, ilPSP.Environment.NumThreads)), outside);
                else
                    Assert.AreEqual(1, outside);

                int M = 24, N = 24, K = 24;
                var rnd = new Random(99);
                double[] A = Rnd(rnd, M * K), B = Rnd(rnd, K * N);
                double[] Cref = new double[M * N];
                fixed(double* pA = A, pB = B, pC = Cref) {
                    BLAS.m_seq_BLAS.DGEMM(ref Ntr, ref Ntr, ref M, ref N, ref K, ref One, pA, ref M, pB, ref K, ref Zero, pC, ref M);
                }

                int nItems = 64;
                int[] countsInside = new int[nItems];
                bool[] ok = new bool[nItems];
                long changesBefore = BLAS.OpenBLASThreadCountChanges;
                ilPSP.Environment.ParallelFor(0, nItems, delegate (int i) {
                    countsInside[i] = BLAS.OpenBLASNumThreads;
                    double[] C = new double[M * N];
                    int _M = M, _N = N, _K = K, lda = M, ldb = K, ldc = M;
                    double alpha = 1.0, beta = 0.0;
                    for(int rep = 0; rep < 50; rep++) {
                        fixed(double* pA = A, pB = B, pC = C) {
                            BLAS.SeqDGEMM(ref Ntr2, ref Ntr2, ref _M, ref _N, ref _K, ref alpha, pA, ref lda, pB, ref ldb, ref beta, pC, ref ldc);
                        }
                    }
                    ok[i] = C.SequenceEqual(Cref) || C.Zip(Cref, (u, v) => Math.Abs(u - v)).Max() < 1e-13;
                });
                long changes = BLAS.OpenBLASThreadCountChanges - changesBefore;
                Console.WriteLine($"OpenBLAS thread count changes during one parallel section: {changes}");

                Assert.IsTrue(ok.All(b => b), "SeqDGEMM result mismatch in parallel section");
                if(ilPSP.Environment.NumThreads > 1)
                    Assert.IsTrue(countsInside.All(c => c == 1), "OpenBLAS thread count inside parallel section: " + countsInside.Distinct().ToConcatString("", ",", ""));
                Assert.AreEqual(outside, BLAS.OpenBLASNumThreads, "OpenBLAS thread count not restored after parallel section");
                Assert.LessOrEqual(changes, 2);

                // provider switch is refused within a parallel section
                bool[] refused = new bool[2];
                ilPSP.Environment.ParallelFor(0, 2, delegate (int i) {
                    refused[i] = !BLAS.TryActivateProvider(BLASProvider.MKL, out _) || ilPSP.Environment.NumThreads <= 1;
                });
                Assert.IsTrue(refused.All(b => b));
                Assert.AreEqual(BLASProvider.OpenBLAS, BLAS.Provider);
            });
        }

        static int Ntr = 'n', Ntr2 = 'n';
        static double One = 1.0, Zero = 0.0;

        /// <summary>
        /// calibration leaves provider, threading state and the ordinary BLAS instance unchanged
        /// </summary>
        [Test]
        public static void CalibrationRestoresState() {
            if(!BLAS.OpenBLASAvailable)
                Assert.Ignore("OpenBLAS unavailable: " + BLAS.OpenBLASUnavailableReason);
            foreach(var p in AvailableProviders()) {
                WithProvider(p, delegate () {
                    var before = BLAS.F77_BLAS;
                    int threadsBefore = BLAS.OpenBLASNumThreads;
                    var cal = BLAS.CalibrateProviders();
                    Console.WriteLine($"active provider {p}:");
                    Console.WriteLine(cal.ToString());
                    Assert.AreSame(before, BLAS.F77_BLAS);
                    Assert.AreEqual(p, BLAS.Provider);
                    if(p == BLASProvider.OpenBLAS)
                        Assert.AreEqual(threadsBefore, BLAS.OpenBLASNumThreads);
                    Assert.IsTrue(cal.Cases.Where(c => !c.Threaded).All(c => c.ErrorRatio <= 1.0), "correctness check");
                });
            }
        }
    }
}
