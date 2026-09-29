using System;
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
    }
}
