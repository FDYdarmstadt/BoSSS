#if NET
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;

namespace ilPSP {

    /// <summary>
    /// Tests for the blocked AVX2 Multiply path (MultidimensionalArray_MultiplyGett.cs): the results must be
    /// bitwise identical to the generated reference kernels (same per-entry arithmetic and summation order),
    /// for randomized loop layouts (unit-stride, gathered and broadcast operands, 0/1/2 summation loops,
    /// vector tails, padding, outer loops, index transformations with skipped entries, non-finite values).
    /// </summary>
    [TestFixture]
    public unsafe class TensorMultiplyGett_Tests {

        delegate void Disp(int r, int s, double* t, double* a, double* b, int* len, int* ct, int* ca, int* cb,
            int* lenSum, int* sa, int* sb, double alpha, double beta);
        delegate void DispTrf(int r, int s, double* t, double* a, double* b, int* len, int* ct, int* ca, int* cb,
            int* lenSum, int* sa, int* sb, double alpha, double beta, int* trfT, int* trfA, int* trfB, int swT, int swA, int swB,
            int preT, int cycT, int postT, int preA, int cycA, int postA, int preB, int cycB, int postB);
        delegate bool Run(int r, int s, double* t, double* a, double* b, int* len, int* ct, int* ca, int* cb,
            int* lenSum, int* sa, int* sb, double alpha, double beta, bool profitableOnly);
        delegate bool RunTrf(int r, int s, double* t, double* a, double* b, int* len, int* ct, int* ca, int* cb,
            int* lenSum, int* sa, int* sb, double alpha, double beta, int* trfT, int* trfA, int* trfB, int swT, int swA, int swB,
            int preT, int cycT, int postT, int preA, int cycA, int postA, int preB, int cycB, int postB, bool profitableOnly);

        static T Bind<T>(string name) where T : Delegate =>
            (T)typeof(MultidimensionalArray).GetMethod(name, BindingFlags.NonPublic | BindingFlags.Static).CreateDelegate(typeof(T));

        static readonly Disp RefDisp = Bind<Disp>("Multiply_Dispatch");
        static readonly DispTrf RefDispTrf = Bind<DispTrf>("MultiplyWTrafo_Dispatch");
        static readonly Run GettRun = Bind<Run>("GettRun");
        static readonly RunTrf GettRunTrf = Bind<RunTrf>("GettRunWTrafo");

        /// <summary> random loop layout; operand storage layouts with random dimension order and padding </summary>
        class Layout {
            public int DT, NS;
            public int[] len, ct, ca, cb, lenSum, sa, sb;
            public int sizeT, sizeA, sizeB;
            public int swA, swB;     // transformed operands (loop 0)
            public int[] trf;        // index trafo for loop 0
        }

        // assigns strides to an operand which has the given running loops (and all summation loops)
        static int AssignStrides(Random rnd, Layout L, bool[] runMember, int[] cRun, int[] cSum, int trfExtent0) {
            var dims = new List<int>(); // >= 0: running loop d, < 0: summation loop -1-s
            for(int d = 0; d < L.DT; d++) if(runMember[d]) dims.Add(d);
            for(int s = 0; s < L.NS; s++) dims.Add(-1 - s);
            dims = dims.OrderBy(_ => rnd.Next()).ToList(); // innermost first
            int stride = 1;
            foreach(int d in dims) {
                int l = d >= 0 ? (d == 0 && trfExtent0 > 0 ? trfExtent0 : L.len[d]) : L.lenSum[-1 - d];
                if(d >= 0) cRun[d] = stride; else cSum[-1 - d] = stride;
                stride = stride * l + (rnd.NextDouble() < 0.3 ? rnd.Next(1, 4) : 0);
            }
            return stride + 2;
        }

        static Layout RandomLayout(Random rnd, bool trafo) {
            var L = new Layout();
            L.DT = rnd.Next(trafo ? 2 : 1, 5);
            L.NS = rnd.Next(0, 3);
            int[] lenChoices = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 16, 17 };
            L.len = Enumerable.Range(0, L.DT).Select(_ => lenChoices[rnd.Next(lenChoices.Length)]).ToArray();
            L.lenSum = Enumerable.Range(0, L.NS).Select(_ => rnd.Next(1, 8)).ToArray();
            L.ct = new int[L.DT]; L.ca = new int[L.DT]; L.cb = new int[L.DT];
            L.sa = new int[L.NS]; L.sb = new int[L.NS];

            // T: all running loops, random order (the innermost one has unit stride)
            var tMember = Enumerable.Repeat(true, L.DT).ToArray();
            L.sizeT = AssignStrides(rnd, L, tMember, L.ct, new int[L.NS], 0);

            // A, B: random subsets of the running loops
            var aMember = new bool[L.DT];
            var bMember = new bool[L.DT];
            for(int d = 0; d < L.DT; d++) {
                int c = rnd.Next(3); // 0: A, 1: B, 2: both
                aMember[d] = c != 1; bMember[d] = c != 0;
            }
            int ext0 = 0;
            if(trafo) {
                int sw = rnd.Next(3); // 0: A, 1: B, 2: both
                L.swA = sw != 1 ? 1 : 0; L.swB = sw != 0 ? 1 : 0;
                if(L.swA != 0) aMember[0] = true;
                if(L.swB != 0) bMember[0] = true;
                ext0 = rnd.Next(1, L.len[0] + 3);
                L.trf = Enumerable.Range(0, L.len[0]).Select(_ => rnd.Next(-1, ext0)).ToArray(); // -1: skipped
            }
            L.sizeA = AssignStrides(rnd, L, aMember, L.ca, L.sa, L.swA != 0 ? ext0 : 0);
            L.sizeB = AssignStrides(rnd, L, bMember, L.cb, L.sb, L.swB != 0 ? ext0 : 0);
            return L;
        }

        static bool SameBits(double x, double y) =>
            BitConverter.DoubleToInt64Bits(x) == BitConverter.DoubleToInt64Bits(y) || (double.IsNaN(x) && double.IsNaN(y));

        /// <summary>
        /// runs reference and blocked SIMD path on identical data; returns whether the SIMD path was taken
        /// </summary>
        static bool CheckOne(Random rnd, Layout L, double alpha, double beta, bool nonfinite, string info) {
            double[] T0 = new double[L.sizeT], A = new double[L.sizeA], B = new double[L.sizeB];
            for(int i = 0; i < T0.Length; i++) T0[i] = rnd.NextDouble() - .5;
            for(int i = 0; i < A.Length; i++) A[i] = rnd.NextDouble() - .5;
            for(int i = 0; i < B.Length; i++) B[i] = rnd.NextDouble() - .5;
            if(nonfinite) {
                T0[rnd.Next(T0.Length)] = double.NaN;
                T0[rnd.Next(T0.Length)] = double.PositiveInfinity;
                A[rnd.Next(A.Length)] = double.NegativeInfinity;
            }
            double[] Tr = (double[])T0.Clone(), Tg = (double[])T0.Clone();
            int[] ls = L.NS > 0 ? L.lenSum : new int[1], sa = L.NS > 0 ? L.sa : new int[1], sb = L.NS > 0 ? L.sb : new int[1];
            int[] trf = L.trf ?? new int[1];
            bool taken;
            fixed(double* pTr = Tr, pTg = Tg, pA = A, pB = B)
            fixed(int* len = L.len, ct = L.ct, ca = L.ca, cb = L.cb, pls = ls, psa = sa, psb = sb, ptrf = trf) {
                if(L.trf == null) {
                    RefDisp(L.DT, L.NS, pTr, pA, pB, len, ct, ca, cb, pls, psa, psb, alpha, beta);
                    taken = GettRun(L.DT, L.NS, pTg, pA, pB, len, ct, ca, cb, pls, psa, psb, alpha, beta, false);
                } else {
                    RefDispTrf(L.DT, L.NS, pTr, pA, pB, len, ct, ca, cb, pls, psa, psb, alpha, beta,
                        ptrf, ptrf, ptrf, 0, L.swA, L.swB, 0, 1, 0, 0, 1, 0, 0, 1, 0);
                    taken = GettRunTrf(L.DT, L.NS, pTg, pA, pB, len, ct, ca, cb, pls, psa, psb, alpha, beta,
                        ptrf, ptrf, ptrf, 0, L.swA, L.swB, 0, 1, 0, 0, 1, 0, 0, 1, 0, false);
                }
            }
            if(taken) {
                for(int i = 0; i < Tr.Length; i++)
                    if(!SameBits(Tr[i], Tg[i]))
                        Assert.Fail($"{info}: entry {i}: reference {Tr[i]:R}, SIMD path {Tg[i]:R}");
            } else {
                for(int i = 0; i < Tr.Length; i++)
                    if(!SameBits(T0[i], Tg[i]))
                        Assert.Fail($"{info}: path not taken, but destination modified at entry {i}");
            }
            return taken;
        }

        static string Describe(Layout L, double alpha, double beta) =>
            $"DT={L.DT} NS={L.NS} len=[{string.Join(" ", L.len)}] ct=[{string.Join(" ", L.ct)}] ca=[{string.Join(" ", L.ca)}] cb=[{string.Join(" ", L.cb)}] " +
            $"lenSum=[{string.Join(" ", L.lenSum)}] sa=[{string.Join(" ", L.sa)}] sb=[{string.Join(" ", L.sb)}] trf={(L.trf != null ? $"A{L.swA}B{L.swB}" : "-")} alpha={alpha} beta={beta}";

        static void RandomCases(bool trafo, int seed, int N) {
            if(!System.Runtime.Intrinsics.X86.Avx2.IsSupported) {
                Assert.Ignore("AVX2 not supported on this machine.");
                return;
            }
            var rnd = new Random(seed);
            double[] alphas = { 1.0, 0.75, -1.25 }, betas = { 0.0, 1.0, 0.25 };
            int taken = 0;
            for(int n = 0; n < N; n++) {
                var L = RandomLayout(rnd, trafo);
                double alpha = alphas[rnd.Next(alphas.Length)], beta = betas[rnd.Next(betas.Length)];
                bool nonfinite = rnd.NextDouble() < 0.1;
                if(CheckOne(rnd, L, alpha, beta, nonfinite, $"case {n}: " + Describe(L, alpha, beta)))
                    taken++;
            }
            Console.WriteLine($"{(trafo ? "transformed" : "plain")}: {taken} of {N} random cases on the SIMD path, all bitwise identical to the reference kernels");
            Assert.Greater(taken, N / 3, "too few random cases on the SIMD path");
        }

        /// <summary> random layouts without index transformation </summary>
        [Test]
        public static void RandomLayouts() {
            RandomCases(false, 4711, 4000);
        }

        /// <summary> random layouts with index transformation (incl. skipped entries) </summary>
        [Test]
        public static void RandomLayoutsWTrafo() {
            RandomCases(true, 815, 3000);
        }

        /// <summary>
        /// Integration via the public API, for program shapes from the census of PLAN-4 (all operand
        /// patterns of the micro-kernels); the result must be bitwise identical to the generated kernels.
        /// </summary>
        [Test]
        public static void CensusPrograms() {
            if(!System.Runtime.Intrinsics.X86.Avx2.IsSupported) {
                Assert.Ignore("AVX2 not supported on this machine.");
                return;
            }
            var rnd = new Random(1234);
            // (T, A, B, extents of T, extents of summation indices p, q)
            var cases = new (string T, string A, string B, int[] lenT, int[] lenS)[] {
                ("ijk", "T(i)pk", "ipj", new[] { 7, 6, 10 }, new[] { 5 }),        // edge quadrature, outer product (X re-used)
                ("ijk", "T(i)pj", "ipk", new[] { 7, 10, 6 }, new[] { 5 }),
                ("ijk", "pjq", "ipkq", new[] { 5, 6, 7 }, new[] { 9, 3 }),        // volume quadrature, gathered X, 2 sums
                ("ijk", "ijp", "T(i)jkp", new[] { 6, 7, 9 }, new[] { 3 }),        // gathered X, no re-use
                ("ijk", "ik", "ij", new[] { 9, 6, 13 }, new int[0]),              // no summation
                ("ijk", "ij", "ijk", new[] { 5, 9, 15 }, new int[0]),
                ("ijk", "ijk", "T(i)", new[] { 6, 5, 7 }, new int[0]),
                ("ij", "ip", "pj", new[] { 13, 10 }, new[] { 16 }),
                ("ij", "pqi", "pqj", new[] { 15, 15 }, new[] { 6, 4 }),
                ("ij", "ij", "ji", new[] { 9, 13 }, new int[0]),                  // unit-stride X, gathered Y
                ("ij", "ji", "ji", new[] { 9, 13 }, new int[0]),                  // both gathered
                ("ijk", "ijk", "ijk", new[] { 3, 5, 6 }, new int[0]),             // both unit stride
            };
            foreach(var c in cases) {
                foreach(double beta in new[] { 0.0, 1.0, 0.25 }) {
                    string Strip(string s) => s.Replace("T(", "").Replace(")", "");
                    var idx = new Dictionary<char, int>();
                    for(int d = 0; d < c.T.Length; d++) idx[c.T[d]] = c.lenT[d];
                    idx['p'] = c.lenS.Length > 0 ? c.lenS[0] : 1; idx['q'] = c.lenS.Length > 1 ? c.lenS[1] : 1;
                    int[] Shape(string s) => Strip(s).Select(ch => idx[ch]).ToArray();
                    var T = MultidimensionalArray.Create(Shape(c.T));
                    var A = MultidimensionalArray.Create(Shape(c.A));
                    var B = MultidimensionalArray.Create(Shape(c.B));
                    for(int i = 0; i < T.Storage.Length; i++) T.Storage[i] = rnd.NextDouble() - .5;
                    for(int i = 0; i < A.Storage.Length; i++) A.Storage[i] = rnd.NextDouble() - .5;
                    for(int i = 0; i < B.Storage.Length; i++) B.Storage[i] = rnd.NextDouble() - .5;
                    bool trafo = c.A.Contains("T(") || c.B.Contains("T(");
                    int[] trf = trafo ? Enumerable.Range(0, c.lenT[0]).Select(i => (i * 5 + 3) % c.lenT[0] - (i == 2 ? c.lenT[0] : 0)).ToArray() : null;

                    var mp = MultidimensionalArray.MultiplyProgram.Compile(c.T, c.A, c.B, trafo);
                    var Tg = T.CloneAs();
                    Tg.Multiply(0.75, A, B, beta, ref mp, trf);

                    // reference: generated kernels with the same loop description as the driver
                    var Tr = T.CloneAs();
                    ReferenceMultiply(Tr, A, B, 0.75, beta, ref mp, trf);
                    for(int i = 0; i < Tr.Storage.Length; i++)
                        if(!SameBits(Tr.Storage[i], Tg.Storage[i]))
                            Assert.Fail($"{c.T},{c.A},{c.B} beta={beta}: entry {i}: reference {Tr.Storage[i]:R}, Multiply {Tg.Storage[i]:R}");
                }
            }
        }

        /// <summary> the generated kernels, called with the loop description of the Multiply driver </summary>
        static void ReferenceMultiply(MultidimensionalArray T, MultidimensionalArray A, MultidimensionalArray B,
            double alpha, double beta, ref MultidimensionalArray.MultiplyProgram mp, int[] trf) {
            int DT = T.Dimension;
            int[] len = new int[DT], ct = new int[DT], ca = new int[DT], cb = new int[DT];
            int[] ls = new int[4], sa = new int[4], sb = new int[4];
            for(int d = 0; d < DT; d++) { len[d] = T.GetLength(d); ct[d] = T.GetCycle(d); }
            var getCycles = typeof(MultidimensionalArray.MultiplyProgram).GetMethod("GetCyclesAndRanks", BindingFlags.NonPublic | BindingFlags.Instance);
            int NS = (int)typeof(MultidimensionalArray.MultiplyProgram).GetField("NoOfSumCycles", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(mp);
            int iTrafo = (int)typeof(MultidimensionalArray.MultiplyProgram).GetField("iTrafoIdx", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(mp);
            int swA = (int)typeof(MultidimensionalArray.MultiplyProgram).GetField("TrfA0Sw", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(mp);
            int swB = (int)typeof(MultidimensionalArray.MultiplyProgram).GetField("TrfB0Sw", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(mp);
            fixed(int* pca = ca, pcb = cb, pls = ls, psa = sa, psb = sb) {
                object boxed = mp;
                getCycles.Invoke(boxed, new object[] { Pointer.Box(pca, typeof(int*)), Pointer.Box(pcb, typeof(int*)), Pointer.Box(pls, typeof(int*)),
                    Pointer.Box(psa, typeof(int*)), Pointer.Box(psb, typeof(int*)), A, B });
            }
            if(NS == 2 && ls[1] > ls[0]) { // as in the driver: inner summation loop is the smaller one
                (ls[0], ls[1]) = (ls[1], ls[0]); (sa[0], sa[1]) = (sa[1], sa[0]); (sb[0], sb[1]) = (sb[1], sb[0]);
            }
            if(iTrafo > 0) {
                int k = iTrafo;
                (len[0], len[k]) = (len[k], len[0]); (ct[0], ct[k]) = (ct[k], ct[0]);
                (ca[0], ca[k]) = (ca[k], ca[0]); (cb[0], cb[k]) = (cb[k], cb[0]);
            }
            fixed(double* pT = T.Storage, pA = A.Storage, pB = B.Storage)
            fixed(int* plen = len, pct = ct, pca = ca, pcb = cb, pls = ls, psa = sa, psb = sb, ptrf = trf ?? new int[1]) {
                if(iTrafo >= 0)
                    RefDispTrf(DT, NS, pT, pA, pB, plen, pct, pca, pcb, pls, psa, psb, alpha, beta,
                        ptrf, ptrf, ptrf, 0, swA, swB, 0, 1, 0, 0, 1, 0, 0, 1, 0);
                else
                    RefDisp(DT, NS, pT, pA, pB, plen, pct, pca, pcb, pls, psa, psb, alpha, beta);
            }
        }
    }
}
#endif
