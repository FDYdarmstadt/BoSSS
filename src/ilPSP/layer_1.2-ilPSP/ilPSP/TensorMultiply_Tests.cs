using System;
using System.Linq;
using System.Reflection;
using NUnit.Framework;

namespace ilPSP {
    /// <summary>Stride, scaling and aliasing regressions for matrix-like Multiply programs.</summary>
    [TestFixture]
    public unsafe class TensorMultiply_Tests {
        delegate void Kernel(int r, int s, double* t, double* a, double* b,
            int* lr, int* ct, int* ca, int* cb, int* ls, int* sa, int* sb, double alpha, double beta);
        static readonly Kernel Reference = (Kernel)typeof(MultidimensionalArray).GetMethod(
            "Multiply_SumGOTO_GOTO", BindingFlags.NonPublic | BindingFlags.Static).CreateDelegate(typeof(Kernel));

        static MultidimensionalArray Make(int batch, int rows, int cols, bool padded, bool strided) {
            int[] shape = batch == 0 ? new[] { rows, cols } : new[] { batch, rows, cols };
            if(!padded && !strided) return MultidimensionalArray.Create(shape);
            // Fixing the last index of an extra dimension makes *every* logical
            // stride non-unit; ordinary subarray frames only pad the row stride.
            int rank = shape.Length;
            int[] parentShape = new int[rank + (strided ? 1 : 0)];
            int[] start = new int[parentShape.Length], end = new int[parentShape.Length];
            for(int d = 0; d < rank; d++) {
                parentShape[d] = shape[d] + 3; start[d] = 1; end[d] = shape[d];
            }
            if(strided) { parentShape[rank] = 2; start[rank] = 1; end[rank] = -1; }
            return MultidimensionalArray.Create(parentShape).ExtractSubArrayShallow(start, end);
        }

        static void Fill(MultidimensionalArray x, Random rng) {
            for(int i = 0; i < x.Storage.Length; i++) x.Storage[i] = rng.NextDouble() - .5;
        }

        static void Check(int m, int n, int k, int batches, bool transposeA, bool transposeB,
            bool padded, bool strided, double alpha, double beta, int alias = 0, bool nonfinite = false,
            bool broadcast = false, int strideMask = 7) {
            var a = Make(broadcast ? 0 : batches, transposeA ? k : m, transposeA ? m : k, padded, strided && (strideMask & 1) != 0);
            var b = Make(batches, transposeB ? n : k, transposeB ? k : n, padded, strided && (strideMask & 2) != 0);
            var t = Make(batches, m, n, padded, strided && (strideMask & 4) != 0);
            var rng = new Random(7919);
            Fill(a,rng); Fill(b,rng); Fill(t,rng);
            if(alias == 1) a = t.ResizeShallow(m,n);
            if(alias == 2) b = t.ResizeShallow(m,n);
            if(nonfinite) {
                t.Storage[t.Index(new int[t.Dimension])] = double.NaN;
                t.Storage[t.Index(new int[t.Dimension]) + t.GetCycle(t.Dimension-1)] = double.PositiveInfinity;
            }
            double[] initial = (double[])t.Storage.Clone();
            int rank = t.Dimension, row = rank-2, col = rank-1, ar = a.Dimension-2;
            int* lr = stackalloc int[rank]; int* ct = stackalloc int[rank];
            int* ca = stackalloc int[rank]; int* cb = stackalloc int[rank];
            int* ls = stackalloc int[1]; int* sa = stackalloc int[1]; int* sb = stackalloc int[1];
            for(int d = 0; d < rank; d++) { lr[d] = t.GetLength(d); ct[d] = t.GetCycle(d); ca[d] = cb[d] = 0; }
            if(batches > 0) { ca[0] = broadcast ? 0 : a.GetCycle(0); cb[0] = b.GetCycle(0); }
            ca[row] = a.GetCycle(ar + (transposeA ? 1 : 0));
            cb[col] = b.GetCycle(row + (transposeB ? 0 : 1));
            sa[0] = a.GetCycle(ar + (transposeA ? 0 : 1));
            sb[0] = b.GetCycle(row + (transposeB ? 1 : 0)); ls[0] = k;
            fixed(double* pt = t.Storage, pa = a.Storage, pb = b.Storage) {
                Reference(rank,1,pt+t.Index(new int[rank]),pa+a.Index(new int[a.Dimension]),
                    pb+b.Index(new int[rank]),lr,ct,ca,cb,ls,sa,sb,alpha,beta);
            }
            double[] expected = (double[])t.Storage.Clone();
            initial.CopyTo(t.Storage,0);
            string prefix = batches > 0 ? "b" : "";
            var program = MultidimensionalArray.MultiplyProgram.Compile(prefix+"ij",
                (broadcast ? "" : prefix)+(transposeA ? "ki" : "ik"), prefix+(transposeB ? "jk" : "kj"));
            t.Multiply(alpha,a,b,beta,ref program);
            for(int i = 0; i < expected.Length; i++) {
                // beta = 0: non-finite entries of T on input give an undefined result (owner decision 2026-10-22: no extra
                // pass for NaN propagation; BLAS paths do not read T); non-finite values with beta != 0 must propagate
                if(beta == 0.0 && !double.IsFinite(initial[i])) continue;
                if(double.IsNaN(expected[i])) Assert.IsTrue(double.IsNaN(t.Storage[i]), "NaN propagation");
                else if(double.IsInfinity(expected[i])) Assert.AreEqual(expected[i],t.Storage[i]);
                else Assert.That(t.Storage[i], Is.EqualTo(expected[i]).Within(2e-12*Math.Max(1,Math.Abs(expected[i]))),
                    "Storage entry " + i);
            }
        }

        /// <summary>All transpose combinations, frames, batches and alpha/beta values.</summary>
        [Test]
        public static void MatrixLayouts() {
            foreach(bool ta in new[] { false,true }) foreach(bool tb in new[] { false,true })
            foreach(bool padded in new[] { false,true }) foreach(int batches in new[] { 0,3 })
            foreach(double alpha in new[] { 0.0,-.75,1.0 }) foreach(double beta in new[] { 0.0,.25,1.0,-.5 })
                Check(17,19,23,batches,ta,tb,padded,false,alpha,beta);
        }

        /// <summary>
        /// As <see cref="MatrixLayouts"/>, but above the size limit below which the blocked SIMD path is
        /// preferred to direct DGEMM (MultidimensionalArray.GETT_PREFERRED_BELOW_MNK), i.e. on the BLAS paths.
        /// </summary>
        [Test]
        public static void MatrixLayoutsBlasSize() {
            Assert.Greater(29L * 31 * 37, MultidimensionalArray.GETT_PREFERRED_BELOW_MNK);
            foreach(bool ta in new[] { false,true }) foreach(bool tb in new[] { false,true })
            foreach(bool padded in new[] { false,true }) foreach(int batches in new[] { 0,3 })
            foreach(double alpha in new[] { -.75,1.0 }) foreach(double beta in new[] { 0.0,.25 })
                Check(29,31,37,batches,ta,tb,padded,false,alpha,beta);
        }

        /// <summary>Fallback for small, strided and aliasing layouts, plus broadcast batches.</summary>
        [Test]
        public static void FallbackAndScaling() {
            foreach(int size in new[] { 1,4,8,15,16,17 }) Check(size,size,size,0,false,false,false,false,.75,.25);
            foreach(int mask in new[] { 1,2,4,7 })
                Check(17,19,23,3,false,false,true,true,.75,.25,0,false,false,mask);
            Check(16,16,16,0,false,false,false,false,.75,.25,1);
            Check(16,16,16,0,false,false,false,false,.75,.25,2);
            Check(16,16,16,0,false,false,true,false,.75,0,0,true);
            Check(16,16,16,0,false,false,true,false,0,0,0,true);
            Check(17,19,23,3,true,true,true,false,.75,.25,0,false,true);
        }

        [Test]
        public static void ShortReductionSimdTailAndOuterStride() {
            // Note: since 2026-10-01, TryMultiplySimd is disabled, so this exercises the generated kernels;
            // keep it as regression test for a future SIMD short-reduction path (PLAN-4, Phase F).
            const int outer = 2, count = 10, sum = 3;
            var t = MultidimensionalArray.Create(outer, count);
            var a = MultidimensionalArray.Create(outer, sum, count);
            var b = MultidimensionalArray.Create(outer, sum);
            var rng = new Random(2718);
            for(int i = 0; i < t.Storage.Length; i++) t.Storage[i] = rng.NextDouble() - .5;
            for(int i = 0; i < a.Storage.Length; i++) a.Storage[i] = rng.NextDouble() - .5;
            for(int i = 0; i < b.Storage.Length; i++) b.Storage[i] = rng.NextDouble() - .5;

            double[] expected = (double[])t.Storage.Clone();
            for(int o = 0; o < outer; o++) for(int i = 0; i < count; i++) {
                double dot = 0.0;
                for(int l = 0; l < sum; l++)
                    dot += a.Storage[a.Index(o,l,i)] * b.Storage[b.Index(o,l)];
                expected[t.Index(o,i)] = .75 * dot + .25 * expected[t.Index(o,i)];
            }
            t.Multiply(.75, a, b, .25, "ij", "ilj", "il");
            for(int i = 0; i < expected.Length; i++)
                Assert.That(t.Storage[i], Is.EqualTo(expected[i]).Within(2e-14), "Storage entry " + i);
        }

        delegate bool TryKernelGeneral(int r, int s, double* t, double* a, double* b,
            int* lr, int* ct, int* ca, int* cb, int* ls, int* sa, int* sb, double alpha, double beta, long gettBelow, int minDim);

        /// <summary>
        /// Packed DGEMM for general contractions (multi-index M, N, K, batches, padded strides), called directly,
        /// against the generic reference kernel; e.g. T[i,j,k] = sum_{p,q} A[p,j,q] B[i,p,k,q] (LEC volume quadrature).
        /// </summary>
        [Test]
        public static void PackedGemmGeneral() {
            var tryGeneral = (TryKernelGeneral)typeof(MultidimensionalArray).GetMethod("TryPackedGemmGeneral",
                BindingFlags.NonPublic | BindingFlags.Static).CreateDelegate(typeof(TryKernelGeneral));
            var sw = typeof(MultidimensionalArray).GetField("s_PackedGemmGeneral", BindingFlags.NonPublic | BindingFlags.Static);
            bool swOld = (bool)sw.GetValue(null);
            sw.SetValue(null, true);
            try {
            var rng = new Random(4242);
            // (T, A, B index strings, extents of i,j,k (running), p,q (summed)), all above the size limit
            var cases = new (string T, string A, string B, int[] run, int[] sum)[] {
                ("ijk", "pjq", "ipkq", new[] { 8, 20, 16 }, new[] { 40, 3 }),   // M = j, N = (i,k), K = (p,q)
                ("ij", "ipq", "pjq", new[] { 20, 56 }, new[] { 128, 3 }),        // NEC volume quadrature
                ("ijk", "ijpq", "ipkq", new[] { 3, 48, 48 }, new[] { 60, 2 }),   // batch i
                ("ijk", "jp", "ipk", new[] { 9, 30, 25 }, new[] { 40 }),         // single sum, N = (i,k)
            };
            foreach(var c in cases)
            foreach(double beta in new[] { 0.0, 0.25 })
            foreach(bool padded in new[] { false, true }) {
                var ext = new System.Collections.Generic.Dictionary<char, int>();
                for(int d = 0; d < c.T.Length; d++) ext[c.T[d]] = c.run[d];
                ext['p'] = c.sum[0]; if(c.sum.Length > 1) ext['q'] = c.sum[1];
                MultidimensionalArray Mk(string idx) {
                    int[] shape = idx.Select(ch => ext[ch]).ToArray();
                    if(!padded) return MultidimensionalArray.Create(shape);
                    var parent = MultidimensionalArray.Create(shape.Select(l => l + 2).ToArray());
                    return parent.ExtractSubArrayShallow(shape.Select(_ => 1).ToArray(), shape.ToArray());
                }
                var T = Mk(c.T); var A = Mk(c.A); var B = Mk(c.B);
                foreach(var X in new[] { T, A, B }) for(int i = 0; i < X.Storage.Length; i++) X.Storage[i] = rng.NextDouble() - .5;

                int DT = c.T.Length, NS = c.sum.Length;
                int* lr = stackalloc int[DT]; int* ct = stackalloc int[DT]; int* ca = stackalloc int[DT]; int* cb = stackalloc int[DT];
                int* ls = stackalloc int[2]; int* sa = stackalloc int[2]; int* sb = stackalloc int[2];
                for(int d = 0; d < DT; d++) {
                    lr[d] = T.GetLength(d); ct[d] = T.GetCycle(d);
                    int ia = c.A.IndexOf(c.T[d]), ib = c.B.IndexOf(c.T[d]);
                    ca[d] = ia >= 0 ? A.GetCycle(ia) : 0; cb[d] = ib >= 0 ? B.GetCycle(ib) : 0;
                }
                for(int s = 0; s < NS; s++) {
                    char ch = "pq"[s];
                    ls[s] = ext[ch]; sa[s] = A.GetCycle(c.A.IndexOf(ch)); sb[s] = B.GetCycle(c.B.IndexOf(ch));
                }
                double[] initial = (double[])T.Storage.Clone();
                bool taken;
                fixed(double* pt = T.Storage, pa = A.Storage, pb = B.Storage) {
                    double* t0 = pt + T.Index(new int[DT]), a0 = pa + A.Index(new int[A.Dimension]), b0 = pb + B.Index(new int[B.Dimension]);
                    Reference(DT, NS, t0, a0, b0, lr, ct, ca, cb, ls, sa, sb, .75, beta);
                    double[] expected = (double[])T.Storage.Clone();
                    initial.CopyTo(T.Storage, 0);
                    taken = tryGeneral(DT, NS, t0, a0, b0, lr, ct, ca, cb, ls, sa, sb, .75, beta, 0, 0);
                    Assert.IsTrue(taken, $"{c.T},{c.A},{c.B}: general packed DGEMM path expected");
                    double maxRef = expected.Where(x => !double.IsNaN(x)).Max(x => Math.Abs(x));
                    for(int i = 0; i < expected.Length; i++) {
                        if(double.IsNaN(expected[i])) Assert.IsTrue(double.IsNaN(T.Storage[i]), "NaN propagation, entry " + i);
                        else Assert.That(T.Storage[i], Is.EqualTo(expected[i]).Within(1e-13 * maxRef), $"{c.T},{c.A},{c.B} beta={beta} padded={padded}: entry {i}");
                    }
                }
            }
            } finally {
                sw.SetValue(null, swOld);
            }
        }

        delegate bool TryKernel(int r, double* t, double* a, double* b,
            int* lr, int* ct, int* ca, int* cb, int k, int* sa, int* sb, double alpha, double beta);

        /// <summary>
        /// The direct DGEMV path must reject matrix layouts with a leading dimension smaller than the
        /// reduction length (overlapping rows; BLAS would reject the call and leave T unchanged), and it
        /// must not touch the destination before all eligibility checks are passed (beta == 0 pre-scaling).
        /// Such layouts cannot be created via the public API, so the selector is called directly.
        /// </summary>
        [Test]
        public static void GemvLeadingDimensionGuard() {
            var tryGemv = (TryKernel)typeof(MultidimensionalArray).GetMethod("TryDirectGemv",
                BindingFlags.NonPublic | BindingFlags.Static).CreateDelegate(typeof(TryKernel));
            const int m = 32, k = 24; // m*k >= 512, so that only the leading dimension can lead to rejection
            var rng = new Random(31337);
            double[] aStor = new double[m * k + k], bStor = new double[k], tStor = new double[m];
            for(int i = 0; i < aStor.Length; i++) aStor[i] = rng.NextDouble() - .5;
            for(int i = 0; i < bStor.Length; i++) bStor[i] = rng.NextDouble() - .5;
            for(int i = 0; i < tStor.Length; i++) tStor[i] = rng.NextDouble() - .5;

            int* lr = stackalloc int[1]; int* ct = stackalloc int[1]; int* ca = stackalloc int[1]; int* cb = stackalloc int[1];
            int* ls = stackalloc int[1]; int* sa = stackalloc int[1]; int* sb = stackalloc int[1];
            lr[0] = m; ct[0] = 1; cb[0] = 0; ls[0] = k; sa[0] = 1; sb[0] = 1;

            fixed(double* pa = aStor, pb = bStor, pt = tStor) {
                // invalid: row stride k-2 < k
                ca[0] = k - 2;
                double[] before = (double[])tStor.Clone();
                bool taken = tryGemv(1, pt, pa, pb, lr, ct, ca, cb, k, sa, sb, .75, 0.0);
                Assert.IsFalse(taken, "DGEMV path must reject lda < k");
                for(int i = 0; i < m; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(before[i]), BitConverter.DoubleToInt64Bits(tStor[i]), "destination modified by rejected path, entry " + i);

                // positive control: valid row stride k; must be taken and agree with the reference kernel
                ca[0] = k;
                double[] expected = (double[])tStor.Clone();
                fixed(double* pe = expected)
                    Reference(1, 1, pe, pa, pb, lr, ct, ca, cb, ls, sa, sb, .75, .25);
                Assert.IsTrue(tryGemv(1, pt, pa, pb, lr, ct, ca, cb, k, sa, sb, .75, .25), "DGEMV path expected for a valid layout");
                for(int i = 0; i < m; i++)
                    Assert.That(tStor[i], Is.EqualTo(expected[i]).Within(2e-14), "entry " + i);
            }
        }
    }
}
