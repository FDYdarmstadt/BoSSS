using System;
using ilPSP.Utils;
using MPI.Wrappers.Utils;
using NUnit.Framework;

namespace ilPSP {
    /// <summary>Provider coupling and representative LAPACK operations through both imports and driver methods.</summary>
    [TestFixture]
    public unsafe class LAPACKProvider_Tests {
        [Test]
        public static void Operations([Values(BLASProvider.MKL, BLASProvider.OpenBLAS)] BLASProvider provider,
            [Values(false, true)] bool threaded) {
            if(provider == BLASProvider.OpenBLAS && !BLAS.OpenBLASAvailable)
                Assert.Ignore(BLAS.OpenBLASUnavailableReason);
            var previousProvider = BLAS.Provider;
            bool previousThreaded = previousProvider == BLASProvider.MKL
                ? ReferenceEquals(BLAS.F77_BLAS, BLAS.m_omp_BLAS) : BLAS.OpenBLASNumThreads > 1;
            try {
                Assert.IsTrue(BLAS.TryActivateProvider(provider, out string reason), reason);
                if(threaded) BLAS.ActivateOMP(); else BLAS.ActivateSEQ();
                var native = LAPACK.F77_LAPACK;
                Assert.AreEqual(provider, native.Provider);
                if(provider == BLASProvider.MKL) {
                    bool effectiveOMP = threaded && Environment.MaxNumOpenMPthreads > 1 && Environment.NumThreads > 1;
                    Assert.AreSame(effectiveOMP ? LAPACK.omp_F77_LAPACK : LAPACK.seq_F77_LAPACK, native);
                    Assert.AreSame(effectiveOMP ? BLAS.m_omp_BLAS : BLAS.m_seq_BLAS, BLAS.F77_BLAS);
                } else {
                    Assert.AreEqual(BLAS.F77_BLAS.LibHandle.val, native.LibHandle.val);
                    if(!threaded) Assert.AreEqual(1, BLAS.OpenBLASNumThreads);
                }

                // LU factorization, solve, condition estimate, and inverse of [[4,1],[1,3]].
                int n = 2, nrhs = 1, info;
                double[] a = { 4, 1, 1, 3 }, b = { 6, 7 }, work = new double[32];
                int[] pivots = new int[2], iwork = new int[2];
                LAPACK.DGETRFF(ref n, ref n, a, ref n, pivots, out info);
                Assert.AreEqual(0, info);
                fixed(double* pa = a, pb = b, pw = work) fixed(int* pp = pivots, pi = iwork) {
                    char trans = 'N';
                    native.DGETRS(ref trans, ref n, ref nrhs, pa, ref n, pp, pb, ref n, out info);
                    Assert.AreEqual(0, info);
                    Assert.AreEqual(1.0, b[0], 1e-13);
                    Assert.AreEqual(2.0, b[1], 1e-13);
                    int norm = '1';
                    double anorm = 5, reciprocal = 0;
                    native.DGECON_(ref norm, ref n, pa, ref n, ref anorm, ref reciprocal, pw, pi, ref info);
                    Assert.AreEqual(0, info);
                    Assert.AreEqual(11.0 / 25.0, reciprocal, 1e-13);
                }
                int lwork = work.Length;
                LAPACK.DGETRI(ref n, a, ref n, pivots, work, ref lwork, out info);
                Assert.AreEqual(0, info);
                Assert.AreEqual(3.0 / 11, a[0], 1e-13);
                Assert.AreEqual(-1.0 / 11, a[1], 1e-13);
                Assert.AreEqual(4.0 / 11, a[3], 1e-13);

                // Cholesky factorization/inversion, then an independent positive-definite solve.
                a = new double[] { 4, 1, 1, 3 };
                LAPACK.DPOTRF('L', n, a, n, out info);
                Assert.AreEqual(0, info);
                LAPACK.DPOTRI('L', n, a, n, out info);
                Assert.AreEqual(0, info);
                Assert.AreEqual(3.0 / 11, a[0], 1e-13);
                Assert.AreEqual(4.0 / 11, a[3], 1e-13);
                a = new double[] { 4, 1, 1, 3 }; b = new double[] { 6, 7 };
                fixed(double* pa = a, pb = b) {
                    int uplo = 'L';
                    native.DPOSV_(ref uplo, ref n, ref nrhs, pa, ref n, pb, ref n, out info);
                    Assert.AreEqual(0, info);
                }
                Assert.AreEqual(1, b[0], 1e-13);
                Assert.AreEqual(2, b[1], 1e-13);

                // Workspace-query convenience paths: least squares and QR followed by explicit Q.
                a = new double[] { 1, 0, 1, 0, 1, 1 }; b = new double[] { 2, 3, 5 };
                LAPACK.DGELSY(3, 2, a, b, 1, 1e-12);
                Assert.AreEqual(2, b[0], 1e-12);
                Assert.AreEqual(3, b[1], 1e-12);
                a = new double[] { 1, 0, 1, 0, 1, 1 };
                double[] tau = new double[2];
                fixed(double* pa = a, pt = tau) {
                    int m = 3, k = 2;
                    LAPACK.DGEQRF(ref m, ref n, pa, ref m, pt);
                    LAPACK.DORGQR(ref m, ref n, ref k, pa, ref m, pt);
                }
                double norm0 = 0, norm1 = 0, dot = 0;
                for(int i = 0; i < 3; i++) {
                    norm0 += a[i] * a[i]; norm1 += a[i + 3] * a[i + 3]; dot += a[i] * a[i + 3];
                }
                Assert.AreEqual(1, norm0, 1e-13);
                Assert.AreEqual(1, norm1, 1e-13);
                Assert.AreEqual(0, dot, 1e-13);
            } finally {
                Assert.IsTrue(BLAS.TryActivateProvider(previousProvider, out string reason), reason);
                if(previousThreaded) BLAS.ActivateOMP(); else BLAS.ActivateSEQ();
            }
        }

        sealed class MissingInterface : DynLibLoader {
#pragma warning disable 649
            public UnsafeLAPACK._DGETRF deliberately_missing_export;
#pragma warning restore 649
            internal MissingInterface(DynLibLoader owner) : base(owner, BoSSS_Prefix) { }
        }

        [Test]
        public static void BorrowedInterfacesDoNotUnloadOwner() {
            using(var owner = new UnsafeDBLAS(Parallelism.SEQ)) {
                var missing = DynLibLoader.TryLoad(() => new MissingInterface(owner), out string reason);
                Assert.IsNull(missing);
                StringAssert.Contains("deliberately_missing_export", reason);
                using(var lapack = new UnsafeLAPACK(owner)) {
                    Assert.AreEqual(owner.LibHandle.val, lapack.LibHandle.val);
                    Assert.IsNotNull(lapack.DGETRF);
                }
                // Both failed binding and disposal of a successful borrowed interface must leave the owner usable.
                int n = 2, inc = 1;
                double[] x = { 3, 4 };
                fixed(double* px = x)
                    Assert.AreEqual(5, owner.DNRM2(ref n, px, ref inc), 1e-14);
            }
        }

        [Test]
        public static void InvalidProviderPreservesSelection() {
            var before = BLAS.F77_BLAS;
            var lapack = LAPACK.F77_LAPACK;
            Assert.Throws<ArgumentOutOfRangeException>(() => BLAS.TryActivateProvider((BLASProvider)1234, out _));
            Assert.AreSame(before, BLAS.F77_BLAS);
            Assert.AreSame(lapack, LAPACK.F77_LAPACK);
        }
    }
}
