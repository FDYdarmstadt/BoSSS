/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using MPI.Wrappers;
using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using System.IO;
using ilPSP.LinSolvers;
using System.Collections;
using System.Diagnostics;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Multigrid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.MultigridTest {

    [TestFixture]
    class TestProgram {

        [TestFixtureSetUp]
        public static void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);

            //GridCommons grd = Grid2D.Cartesian2DGrid(RandomSpacing(), RandomSpacing());
            //grid = new GridData(Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-7, 7, 8), GenericBlas.Linspace(-1, 1, 2)));
            //grid = new GridData(Grid2D.Cartesian2DGrid(new double[] { -6, -4, -2, 2, 4, 6 }, GenericBlas.Linspace(-1, 1, 2)));
            grid = new GridData(Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1.5, 1.5, 17), GenericBlas.Linspace(-1.5, 1.5, 17)));
            MgSeq = CoarseningAlgorithms.CreateSequence(grid);

            for (int p = 0; p <= 3; p++) { // loop over polynomial degrees...
                var uMapping = new UnsetteledCoordinateMapping(new Basis(grid, p));

                var MgMapSeq = new MultigridMapping[MgSeq.Length];
                var BasisSeq = AggregationGridBasis.CreateSequence(MgSeq, uMapping.BasisS);
                for (int iLevel = 0; iLevel < MgSeq.Length; iLevel++) {
                    MgMapSeq[iLevel] = new MultigridMapping(uMapping, BasisSeq[iLevel], new int[] { p }); 
                }
                MultigrigMap.Add(p, MgMapSeq);
            }

        }

        /// <summary>
        /// Sequence of aggregation grids;
        /// index: grid level, coarser grids are found at higher indices.
        /// </summary>
        internal static AggregationGrid[] MgSeq;


        internal static GridData grid;
        
        /// <summary>
        /// Multigrid mappings containing one DG field, for different polynomial degrees (to test), for each multigrid level. 
        /// key: polynomial degree 
        /// value: sequence of multigrid mappings, corresponding with the grid sequence <see cref="MgSeq"/>.
        /// </summary>
        static Dictionary<int, MultigridMapping[]> MultigrigMap = new Dictionary<int, MultigridMapping[]>();

        [TestFixtureTearDown]
        public static void Cleanup() {
            //Console.Out.Dispose();
            csMPI.Raw.mpiFinalize();
        }

        
        static void SetTestValue(ConventionalDGField f) {
            switch (f.Basis.Degree) {

                case 0:
                    f.ProjectField(((x, y) => 1.0));
                    break;
                case 1:
                    f.ProjectField(((x, y) => 1.0 + x - y));
                    break;
                case 2:
                    f.ProjectField(((x, y) => 1.0 + x - y + x * y));
                    break;
                case 3:
                    f.ProjectField(((x, y) => x * x * x + 3.2 * x * y * y - 2.4 * y * y * y));
                    break;
                default:
                    throw new NotImplementedException();
            }
        }

        static void SetTestValue(XDGField Xf, IDictionary<SpeciesId, double> mod) {
            var Tracker = Xf.Basis.Tracker;

            foreach (var S in Tracker.SpeciesIdS) {
                var f = Xf.GetSpeciesShadowField(S);

                SetTestValue(f);
                f.Scale(mod[S]);
                //f.AccConstant(accVal[S]);
            }
        }


        /// <summary>
        /// Tests if the prolongation of a random vector on a coarse grid has jumps (which it should not have).
        /// </summary>
        [Test]
        public static void ProlongationTest([Values(1, 2, 3)] int p) {
            var MgMap = MultigrigMap;

            ProlongationTestRec(p, MgMap[p]);

        }

        private static void ProlongationTestRec(int p, MultigridMapping[] MgMap) {
            var AggBasis = MgMap[0].AggBasis;
            //var AggBasis = aggBasisEs[1]; 
            //{
            var Test = new SinglePhaseField(new Basis(grid, p));

            Random rnd = new Random();
            double[] RestVec = AggBasis[0].LocalDim.ForLoop(i => rnd.NextDouble());
            var PrlgVec = Test.CoordinateVector;

            AggBasis[0].ProlongateToFullGrid(PrlgVec, RestVec);

            int J = grid.Cells.NoOfLocalUpdatedCells;
            var aggGrd = AggBasis[0].AggGrid;
            double Err = 0;
            for (int jagg = 0; jagg < aggGrd.iLogicalCells.NoOfLocalUpdatedCells; jagg++) {
                BitArray CompCellMask = new BitArray(J);
                foreach (int jCell in aggGrd.iLogicalCells.AggregateCellToParts[jagg])
                    CompCellMask[jCell] = true;


                SubGrid CompCellSubGrid = new SubGrid(new CellMask(grid, CompCellMask));

                Err += JumpNorm(Test, CompCellSubGrid.InnerEdgesMask).Pow2();
            }

            //Console.WriteLine("Jump norm of random prolongation: {0}", Err);
            //Debug.Assert(Err < 1.0e-10);
            Assert.Less(Err, 1.0e-10);

            if (MgMap.Length > 1)
                ProlongationTestRec(p, MgMap.Skip(1).ToArray());
        }

        static double JumpNorm(DGField f, EdgeMask em) {
            GridData grd = grid;
            int D = grd.SpatialDimension;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;

            double Unorm = 0;


            EdgeQuadrature.GetQuadrature(
                new int[] { D + 1 }, grd,
                (new EdgeQuadratureScheme(true, em)).Compile(grd, f.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++) {
                        int iEdge = j + i0;
                        int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];

                        MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                        MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);
                        MultidimensionalArray Grad_uIN = MultidimensionalArray.Create(1, NoOfNodes, D);
                        MultidimensionalArray Grad_uOT = MultidimensionalArray.Create(1, NoOfNodes, D);

                        NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN);
                        NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT);

                        f.Evaluate(jCell_IN, 1, NS_IN, uIN);
                        f.Evaluate(jCell_OT, 1, NS_OT, uOT);
                        f.EvaluateGradient(jCell_IN, 1, NS_IN, Grad_uIN);
                        f.EvaluateGradient(jCell_OT, 1, NS_OT, Grad_uOT);

                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });
                        var Grad_uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 1 }, new int[] { j, NoOfNodes - 1, D });
                        uDiff.Acc(+1.0, uIN);
                        uDiff.Acc(-1.0, uOT);
                        Grad_uDiff.Acc(+1.0, Grad_uIN);
                        Grad_uDiff.Acc(-1.0, Grad_uOT);
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            double totNorm = Unorm.MPISum().Sqrt();
            return totNorm;
        }


        /// <summary>
        /// Tests whether restriction and prolongation on a polynomial field 
        /// results in the original data.
        /// </summary>
        [Test]
        public static void PolynomialRestAndPrlgTest([Values(1, 2, 3)] int p) {
            PolynomialRestAndPrlgTestRec(p, MultigrigMap[p]);
        }

        static void PolynomialRestAndPrlgTestRec(int p, MultigridMapping[] MgMapSeq) {
            var AggBasis = MgMapSeq.First().AggBasis[0];
            //var AggBasis = aggBasisEs[1]; 
            //{
            var Orig = new SinglePhaseField(new Basis(grid, p));
            SetTestValue(Orig);
            var Test = Orig.CloneAs();

            double[] OrigVec = Orig.CoordinateVector.ToArray();
            double[] RestVec = new double[AggBasis.LocalDim];
            double[] PrlgVec = new double[OrigVec.Length];


            AggBasis.RestictFromFullGrid(OrigVec, RestVec);
            Test.Clear();
            AggBasis.ProlongateToFullGrid(PrlgVec, RestVec);

            BlockMsrMatrix RestOp = new BlockMsrMatrix(MgMapSeq[0], MgMapSeq[0].ProblemMapping);
            AggBasis.GetRestrictionMatrix(RestOp, MgMapSeq[0], 0);
            double[] RestVec2 = new double[RestVec.Length];
            RestOp.SpMV(1.0, OrigVec, 0.0, RestVec2);


            Test.Clear();
            Test.CoordinateVector.Acc(1.0, PrlgVec);


            Test.Acc(-1.0, Orig);
            double ErrNorm = Test.L2Norm();
            Console.WriteLine("Polynomial Restriction/Prolongation test (p={1}, level={2}): {0}", ErrNorm, p, -MgMapSeq.Length);
            Assert.Less(ErrNorm, 1.0e-8);

            if (MgMapSeq.Length > 1)
                PolynomialRestAndPrlgTestRec(p, MgMapSeq.Skip(1).ToArray());

        }

        /// <summary>
        /// compares the matrix-version of the restriction operator
        /// (<see cref="AggregationGridBasis.GetRestrictionMatrix"/>)
        /// with the direct implementation
        /// (<see cref="AggregationGridBasis.RestictFromFullGrid"/>).
        /// </summary>
        [Test]
        public static void RestictionMatrixTest([Values(1, 2, 3)] int p) {
            RestictionMatrixTestRec(p, MultigrigMap[p]);
        }

        static void RestictionMatrixTestRec(int p, IEnumerable<MultigridMapping> MgMapSeq) {
            var currentLevelMap = MgMapSeq.First();
            AggregationGridBasis AggBasis = currentLevelMap.AggBasis[0];

            var map = new UnsetteledCoordinateMapping(new Basis(grid, p));

            Random rnd = new Random();
            double[] OrigVec = map.LocalLength.ForLoop(i => rnd.NextDouble());
            double[] RestVec = new double[AggBasis.LocalDim];
            double[] PrlgVec = new double[OrigVec.Length];
            double[] RestVec2 = new double[RestVec.Length];
            double[] PrlgVec2 = new double[OrigVec.Length];


            AggBasis.RestictFromFullGrid(OrigVec, RestVec);
            AggBasis.ProlongateToFullGrid(PrlgVec, RestVec);

            BlockMsrMatrix RestOp = new BlockMsrMatrix(MgMapSeq.First(), MgMapSeq.First().ProblemMapping);
            AggBasis.GetRestrictionMatrix(RestOp, MgMapSeq.First(), 0);
            RestOp.SpMV(1.0, OrigVec, 0.0, RestVec2);

            BlockMsrMatrix PrlgOp = RestOp.Transpose();
            PrlgOp.SpMV(1.0, RestVec2, 0.0, PrlgVec2);

            double RestErrNorm = GenericBlas.L2Dist(RestVec2, RestVec);
            double PrlgErrNorm = GenericBlas.L2Dist(PrlgVec2, PrlgVec);
            double LostInfNorm = GenericBlas.L2Dist(OrigVec, PrlgVec2);
            //Console.WriteLine("Rest. matrix test: {0}, Prolong. matrix test {1}, Lost info {2}", RestErrNorm, PrlgErrNorm, LostInfNorm);
            Assert.IsTrue(RestErrNorm < 1.0e-10);
            Assert.IsTrue(PrlgErrNorm < 1.0e-10);
            
            // restriction onto level itself
            BlockMsrMatrix RestMtx = currentLevelMap.FromOtherLevelMatrix(currentLevelMap);
            BlockMsrMatrix ShldBeEye = BlockMsrMatrix.Multiply(RestMtx, RestMtx.Transpose());
            ShldBeEye.AccEyeSp(-1.0);
            double errNorm = ShldBeEye.InfNorm();
            Console.WriteLine("Id norm {0} \t (level {1})", errNorm, currentLevelMap.AggGrid.MgLevel);
            Assert.IsTrue(errNorm < 1.0e-8);


            // recursion
            if (MgMapSeq.Count() > 1)
                RestictionMatrixTestRec(p, MgMapSeq.Skip(1));
        }

        /// <summary>
        /// basic test for the restriction operator for a system.
        /// </summary>
        [Test]
        public static void RestrictionOfSystemOpTest() {

            Basis B1 = new Basis(grid, 0), B2 = new Basis(grid, 2);
            var Map = new UnsetteledCoordinateMapping(B1, B2);

            AggregationGridBasis[][] aB = AggregationGridBasis.CreateSequence(TestProgram.MgSeq.Take(2), new Basis[] { B1, B2 });

            var Lev0 = new MultigridMapping(Map, aB[0], new int[] { B1.Degree, B2.Degree });
            var Lev1 = new MultigridMapping(Map, aB[1], new int[] { B1.Degree, B2.Degree });


            int[] I0col = Lev0.GetSubvectorIndices(new int[] { 0 });
            int[] I1col = Lev0.GetSubvectorIndices(new int[] { 1 });
            int[] I0row = Lev1.GetSubvectorIndices(new int[] { 0 });
            int[] I1row = Lev1.GetSubvectorIndices(new int[] { 1 });

            var RestMtx = Lev1.FromOtherLevelMatrix(Lev0);

            MsrMatrix Rest00 = new MsrMatrix(I0row.Length, I0col.Length, 1, 1);
            RestMtx.WriteSubMatrixTo(Rest00, I0row, default(int[]), I0col, default(int[]));
            MsrMatrix Rest01 = new MsrMatrix(I0row.Length, I1col.Length, 1, 1);
            RestMtx.WriteSubMatrixTo(Rest01, I0row, default(int[]), I1col, default(int[]));
            MsrMatrix Rest10 = new MsrMatrix(I1row.Length, I0col.Length, 1, 1);
            RestMtx.WriteSubMatrixTo(Rest10, I1row, default(int[]), I0col, default(int[]));
            MsrMatrix Rest11 = new MsrMatrix(I1row.Length, I1col.Length, 1, 1);
            RestMtx.WriteSubMatrixTo(Rest11, I1row, default(int[]), I1col, default(int[]));

            Assert.IsTrue(Rest10.InfNorm() == 0.0);
            Assert.IsTrue(Rest01.InfNorm() == 0.0);
            Assert.IsTrue(Rest00.InfNorm() != 0.0);
            Assert.IsTrue(Rest11.InfNorm() != 0.0);


        }


        internal class XDGTestSetup {

            public XDGTestSetup(
                int p,
                double AggregationThreshold,
                int TrackerWidth,
                MultigridOperator.Mode mumo,
                XQuadFactoryHelper.MomentFittingVariants momentFittingVariant,
                ScalarFunction LevSetFunc = null) {

                // Level set, tracker and XDG basis
                // ================================

                if (LevSetFunc == null) {
                    LevSetFunc = ((_2D)((x, y) => 0.8 * 0.8 - x * x - y * y)).Vectorize();
                }
                LevSet = new LevelSet(new Basis(grid, 2), "LevelSet");
                LevSet.Clear();
                LevSet.ProjectField(LevSetFunc);
                LsTrk = new LevelSetTracker(grid, XQuadFactoryHelper.MomentFittingVariants.Classic, TrackerWidth, new string[] { "A", "B" }, LevSet);
                LsTrk.UpdateTracker();

                XB = new XDGBasis(LsTrk, p);

                XSpatialOperator Dummy = new XSpatialOperator(1, 0, 1, QuadOrderFunc.SumOfMaxDegrees(RoundUp: true), "C1", "u");
                //Dummy.EquationComponents["c1"].Add(new 
                Dummy.Commit();

                //Tecplot.PlotFields(new DGField[] { LevSet }, "agglo", 0.0, 3);


                // operator
                // ========

                Debug.Assert(p <= 4);
                XDGBasis opXB = new XDGBasis(LsTrk, 4); // we want to have a very precise quad rule
                var map = new UnsetteledCoordinateMapping(opXB);

                int quadOrder = Dummy.QuadOrderFunction(map.BasisS.Select(bs => bs.Degree).ToArray(), new int[0], map.BasisS.Select(bs => bs.Degree).ToArray());
                //agg = new MultiphaseCellAgglomerator(new CutCellMetrics(momentFittingVariant, quadOrder, LsTrk, LsTrk.SpeciesIdS.ToArray()), AggregationThreshold, false);
                agg = LsTrk.GetAgglomerator(LsTrk.SpeciesIdS.ToArray(), quadOrder, __AgglomerationTreshold: AggregationThreshold);
                               

                foreach (var S in LsTrk.SpeciesIdS)
                    Console.WriteLine("Species {0}, no. of agglomerated cells {1} ",
                        LsTrk.GetSpeciesName(S),
                        agg.GetAgglomerator(S).AggInfo.SourceCells.Count());

                // mass matrix factory
                // ===================

                // Basis maxB = map.BasisS.ElementAtMax(b => b.Degree);
                //MassFact = new MassMatrixFactory(maxB, agg);
                MassFact = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), quadOrder, 1).MassMatrixFactory;


                // Test field
                // ==========

                // set the test field: this is a polynomial function,
                // but different for each species; On this field, restriction followed by prolongation should be the identity
                this.Xdg_uTest = new XDGField(this.XB, "uTest");
                Dictionary<SpeciesId, double> dumia = new Dictionary<SpeciesId, double>();
                int i = 2;
                foreach (var Spc in LsTrk.SpeciesIdS) {
                    dumia.Add(Spc, i);
                    i -= 1;
                }
                SetTestValue(Xdg_uTest, dumia);
                               

                // dummy operator matrix which fits polynomial degree p
                // ====================================================

                Xdg_opMtx = new BlockMsrMatrix(Xdg_uTest.Mapping, Xdg_uTest.Mapping);
                Xdg_opMtx.AccEyeSp(120.0);

                // XDG Aggregation BasiseS
                // =======================

                //XAggB = MgSeq.Select(agGrd => new XdgAggregationBasis[] { new XdgAggregationBasis(uTest.Basis, agGrd) }).ToArray();
                XAggB = new XdgAggregationBasis[MgSeq.Length][];
                var _XAggB = AggregationGridBasis.CreateSequence(MgSeq, Xdg_uTest.Mapping.BasisS);
                for (int iLevel = 0; iLevel < XAggB.Length; iLevel++) {
                    XAggB[iLevel] = new[] { (XdgAggregationBasis)(_XAggB[iLevel][0]) };
                    XAggB[iLevel][0].Update(agg);
                }

                // Multigrid Operator
                // ==================



                Xdg_opMtx = new BlockMsrMatrix(Xdg_uTest.Mapping, Xdg_uTest.Mapping);
                Xdg_opMtx.AccEyeSp(120.0);

                XdgMultigridOp = new MultigridOperator( XAggB, Xdg_uTest.Mapping,
                    Xdg_opMtx,
                    MassFact.GetMassMatrix(Xdg_uTest.Mapping, false),
                    new MultigridOperator.ChangeOfBasisConfig[][] {
                        new MultigridOperator.ChangeOfBasisConfig[] {
                            new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] { 0 }, mode = mumo, Degree = p }
                        }
                    });

            }

            public XdgAggregationBasis[][] XAggB;
            public LevelSet LevSet;
            public LevelSetTracker LsTrk;
            public XDGField Xdg_uTest;
            //public SinglePhaseField Nonx_uTest;
            public XDGBasis XB;
            public MassMatrixFactory MassFact;
            public MultiphaseCellAgglomerator agg;
            public BlockMsrMatrix Xdg_opMtx;
            //public BlockMsrMatrix Nonx_opMtx;
            public MultigridOperator XdgMultigridOp;
            //public MultigridOperator NonxMultigridOp;
        }




        /// <summary>
        /// Tests whether restriction and prolongation on a polynomial field 
        /// results in the original data.
        /// </summary>
        [Test]
        public static void XDG_PolynomialRestAndPrlgTest(
            [Values(0, 1, 2, 3)] int p,
            [Values(0.0, 0.3)] double AggregationThreshold,
            [Values(0, 1)] int TrackerWidth) {

            XQuadFactoryHelper.MomentFittingVariants variant = XQuadFactoryHelper.MomentFittingVariants.OneStepGauss;

            var xt = new XDGTestSetup(p, AggregationThreshold, TrackerWidth, MultigridOperator.Mode.Eye, variant);


            // Basic Restriction & prolongation Test
            // -------------------------------------


            for (int iLevel = 0; iLevel < MgSeq.Length; iLevel++) {

                // create basis
                var XAggBasis = xt.XAggB[iLevel][0];

                // do restriction/prolongation
                double[] RestVec = new double[XAggBasis.LocalDim];
                XAggBasis.RestictFromFullGrid(xt.Xdg_uTest.CoordinateVector, RestVec);
                var Test = xt.Xdg_uTest.CloneAs();
                Test.Clear();
                XAggBasis.ProlongateToFullGrid(Test.CoordinateVector, RestVec);
                xt.agg.Extrapolate(Test.Mapping);

                // compare/test
                var ERR = xt.Xdg_uTest.CloneAs();
                ERR.Acc(-1.0, Test);
                double ERR_NORM = ERR.L2Norm();

                Console.WriteLine("Restriction/Prolongation err (p={0}, level={1}, width={2}, agg={3}): {4}",
                    p, iLevel, TrackerWidth, AggregationThreshold, ERR_NORM);
                Assert.LessOrEqual(ERR_NORM, 1.0e-6);


            }
        }


        /// <summary>
        /// tests the matrix
        /// version of the restriction and prolongation operator.
        /// </summary>
        [Test]
        public static void XDG_MatrixPolynomialRestAndPrlgTest(
            [Values(0, 1, 2, 3)] int p,
            [Values(0.0, 0.3)] double AggregationThreshold,
            [Values(0, 1)] int TrackerWidth) {

            XQuadFactoryHelper.MomentFittingVariants variant = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            var xt = new XDGTestSetup(p, AggregationThreshold, TrackerWidth, MultigridOperator.Mode.Eye, variant);


            // test matrix version of the restriction operator
            // -----------------------------------------------

            List<MultigridMapping> MultigridMaps = new List<MultigridMapping>();
            for (var mgop = xt.XdgMultigridOp; mgop != null; mgop = mgop.CoarserLevel) {
                MultigridMaps.Add(mgop.Mapping);
            }

            for (int iLevel = 0; iLevel < MgSeq.Length; iLevel++) {
                MultigridMapping mgMap = MultigridMaps[iLevel];
                var XAggBasis = mgMap.AggBasis[0];

                // set the test field: 
                XDGField Test = new XDGField(xt.XB, "Test");
                Random rand = new Random();
                for (int i = 0; i < Test.CoordinateVector.Count; i++) {
                    Test.CoordinateVector[i] = rand.NextDouble();
                }
                xt.agg.ClearAgglomerated(Test.CoordinateVector, Test.Mapping);

                // do restriction/prolongation (Reference)
                double[] RestVecRef = new double[XAggBasis.LocalDim];
                XAggBasis.RestictFromFullGrid(Test.CoordinateVector, RestVecRef);

                // and now with the matrix:
                BlockMsrMatrix RestMtx = new BlockMsrMatrix(mgMap, mgMap.ProblemMapping);
                XAggBasis.GetRestrictionMatrix(RestMtx, mgMap, 0);
                double[] RestVec = new double[mgMap.LocalLength];
                RestMtx.SpMV(1.0, Test.CoordinateVector, 0.0, RestVec);
                
                double[] X1 = new double[xt.XdgMultigridOp.Mapping.LocalLength];
                XDGField X2 = new XDGField(Test.Basis);
                xt.XdgMultigridOp.TransformSolInto(Test.CoordinateVector, X1);
                xt.XdgMultigridOp.TransformSolFrom(X2.CoordinateVector, X1);
                //xt.agg.Extrapolate(X2.CoordinatesAsVector, X2.Mapping);
                var ERR2 = Test.CloneAs();
                ERR2.Acc(-1.0, X2);
                double ERR2Norm = ERR2.L2Norm();
                //Console.WriteLine("MultigridOperator TranformInto/FransformFrom mismatch: " + ERR2Norm);
                Assert.LessOrEqual(ERR2Norm, 1.0e-8);

                // compare
                double ERR = 0.0;
                int Nmax = XAggBasis.MaximalLength;
                for (int jAgg = 0; jAgg < mgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells; jAgg++) {
                    int i0Ref = jAgg * Nmax;
                    int i0Tst = mgMap.LocalUniqueIndex(0, jAgg, 0);
                    int N = mgMap.GetLength(jAgg);

                    for (int n = 0; n < N; n++) {
                        double dist = RestVecRef[i0Ref + n] - RestVec[i0Tst + n];
                        ERR += dist.Pow2();
                    }
                }
                Console.WriteLine("Restriction matrix test (iLevel = {0}): {1}", iLevel, ERR);
                Assert.LessOrEqual(ERR, 1.0e-8);

                //
                double[] PrlgVecA = new double[XAggBasis.LocalDim];
                double[] PrlgVecB = new double[mgMap.LocalLength];
                for (int jAgg = 0; jAgg < mgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells; jAgg++) {
                    int i0Ref = jAgg * Nmax;
                    int i0Tst = mgMap.LocalUniqueIndex(0, jAgg, 0);
                    int N = mgMap.GetLength(jAgg);

                    for (int n = 0; n < N; n++) {
                        double rndVal = rand.NextDouble();
                        PrlgVecA[i0Ref + n] = rndVal;
                        PrlgVecB[i0Tst + n] = rndVal;
                    }
                }

                XDGField QA = new XDGField(Test.Basis);
                XDGField QB = new XDGField(Test.Basis);

                XAggBasis.ProlongateToFullGrid(QA.CoordinateVector, PrlgVecA);
                var PrlgMtx = RestMtx.Transpose();
                PrlgMtx.SpMV(1.0, PrlgVecB, 0.0, QB.CoordinateVector);

                XDGField ERR5 = QA.CloneAs();
                ERR5.Acc(-1.0, QB);
                double ERR5_Norm = ERR5.L2Norm();
                Console.WriteLine("Prolongation matrix test (iLevel = {0}): {1}", iLevel, ERR5_Norm);
                Assert.LessOrEqual(ERR5_Norm, 1.0e-8);

            }
        }

        /// <summary>
        /// tests the matrix
        /// version of the restriction and prolongation operator.
        /// </summary>
        [Test]
        public static void XDG_MatrixPolynomialRestAndPrlgTest_2(
            [Values(0, 1, 2, 3)] int p,
            [Values(0.0, 0.3)] double AggregationThreshold,
            [Values(0, 1)] int TrackerWidth,
            [Values(MultigridOperator.Mode.Eye, MultigridOperator.Mode.IdMass)] MultigridOperator.Mode mode) {

            if (AggregationThreshold < 0.1 && p >= 3 && mode == MultigridOperator.Mode.IdMass)
                // this test combination is not supposed to work:
                // without agglomeration, for high p, the mass matrix may be indefinite in small cut-cells
                // => Cholesky decomposition on mass matrix fails, i.e. 'mode == IdMass' cannot succeed.
                return;


            XQuadFactoryHelper.MomentFittingVariants variant = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            var xt = new XDGTestSetup(p, AggregationThreshold, TrackerWidth, mode, variant);


            // Restriction & prolongation together with orthonormalization
            // -----------------------------------------------------------


            for (var mgop = xt.XdgMultigridOp; mgop != null; mgop = mgop.CoarserLevel) {
                var Itself = mgop.Mapping.FromOtherLevelMatrix(mgop.Mapping);
                Itself.AccEyeSp(-1.0);
                double Itslef_Norm = Itself.InfNorm();
                Console.WriteLine("Level {0}, Restriction onto itself {1}", mgop.Mapping.AggGrid.MgLevel, Itslef_Norm);
                Assert.LessOrEqual(Itslef_Norm, 1.0e-8);
            }

            {
                // test change of basis on top level

                XDGField uTestRnd = new XDGField(xt.XB);
                Random rnd = new Random();
                for (int i = 0; i < uTestRnd.CoordinateVector.Count; i++) {
                    uTestRnd.CoordinateVector[i] = rnd.NextDouble();
                }
                xt.agg.ClearAgglomerated(uTestRnd.CoordinateVector, uTestRnd.Mapping);

                // perform change of basis on top level ...
                int Ltop = xt.XdgMultigridOp.Mapping.LocalLength;
                double[] uTest_Fine = new double[Ltop];
                xt.XdgMultigridOp.TransformSolInto(uTestRnd.CoordinateVector, uTest_Fine);

                // .. and back
                XDGField uError2 = uTestRnd.CloneAs();
                uError2.Clear();
                xt.XdgMultigridOp.TransformSolFrom(uError2.CoordinateVector, uTest_Fine);

                // compare: 
                uError2.Acc(-1.0, uTestRnd);
                double NORM_uError = uError2.L2Norm();

                // output
                Console.WriteLine("Top level change of basis error: {0}", NORM_uError);
                Assert.LessOrEqual(NORM_uError, 1.0e-8);

            }

            {

                // perform change of basis on top level
                int Ltop = xt.XdgMultigridOp.Mapping.LocalLength;
                double[] uTest_Fine = new double[Ltop];
                xt.XdgMultigridOp.TransformSolInto(xt.Xdg_uTest.CoordinateVector, uTest_Fine);


                // check for each level of the multigrid operator...
                for (int iLevel = 0; iLevel < MgSeq.Count() - 1; iLevel++) {
                    double[] uTest_Prolonged = new double[Ltop];

                    XDG_Recursive(0, iLevel, xt.XdgMultigridOp, uTest_Fine, uTest_Prolonged);

                    XDGField uError = xt.Xdg_uTest.CloneAs();
                    uError.Clear();
                    xt.XdgMultigridOp.TransformSolFrom(uError.CoordinateVector, uTest_Prolonged);
                    xt.agg.Extrapolate(uError.Mapping);

                    uError.Acc(-1.0, xt.Xdg_uTest);
                    double NORM_uError = uError.L2Norm();

                    Console.WriteLine("Rest/Prlg error, level {0}: {1}", iLevel, NORM_uError);
                    Assert.LessOrEqual(NORM_uError, 1.0e-8);
                }
            }
        }


        /// <summary>
        /// Restricts some solution vector down to a certain multigrid level and prolongates it back.
        /// </summary>
        /// <param name="i">
        /// Multigrid level index; usually 0 at start of recursion.
        /// </param>
        /// <param name="iEnd">
        /// End level index.
        /// </param>
        /// <param name="mgOp"></param>
        /// <param name="FineIn">
        /// Input; solution vector which will be restricted down to level <paramref name="iEnd"/>.
        /// </param>
        /// <param name="FineOut">
        /// Output; vector <paramref name="FineIn"/> restricted to level <paramref name="iEnd"/>, and prolongated back to level <paramref name="i"/>.
        /// </param>
        static void XDG_Recursive(int i, int iEnd, MultigridOperator mgOp, double[] FineIn, double[] FineOut) {
            MultigridMapping mgMap = mgOp.Mapping;

            int Lfin = mgMap.LocalLength;
            int Lcrs = mgOp.CoarserLevel.Mapping.LocalLength;

            Assert.IsTrue(FineIn.Length == Lfin);
            Assert.IsTrue(FineOut.Length == Lfin);

            double[] Coarse1 = new double[Lcrs];
            double[] Coarse2 = new double[Lcrs];

            mgOp.CoarserLevel.Restrict(FineIn, Coarse1);

            if (i == iEnd) {
                Coarse2.SetV(Coarse1);
            } else {
                XDG_Recursive(i + 1, iEnd, mgOp.CoarserLevel, Coarse1, Coarse2);
            }

            mgOp.CoarserLevel.Prolongate(1.0, FineOut, 0.0, Coarse2);
        }

        /// <summary>
        /// tests if the prolongation of an arbitrary restricted vector has jumps (which it should not have).
        /// </summary>
        [Test]
        public static void XDG_ProlongationTest(
            [Values(0, 1, 2, 3)] int p,
            [Values(0.0, 0.3)] double AggregationThreshold,
            [Values(0, 1)] int TrackerWidth,
            [Values(MultigridOperator.Mode.Eye, MultigridOperator.Mode.IdMass)] MultigridOperator.Mode mode) {

            XQuadFactoryHelper.MomentFittingVariants variant = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            var xt = new XDGTestSetup(p, AggregationThreshold, TrackerWidth, MultigridOperator.Mode.Eye, variant
                //, ((Func<double[], double>)(X => X[0] + 0.75)).Vectorize()
                );
            int Jup = grid.Cells.NoOfLocalUpdatedCells;

            
            Random rnd = new Random();
            int Ltop = xt.XdgMultigridOp.Mapping.LocalLength; // Number of DOF's on top multigrid level.


            double[] RndVec = Ltop.ForLoop(i => rnd.NextDouble());
            double[] NoJmpVec = new double[Ltop];

            for (int iLevel = 0; iLevel < MgSeq.Length - 1; iLevel++) {
                XDG_Recursive(0, iLevel, xt.XdgMultigridOp, RndVec, NoJmpVec); // restrict RndVec downt to level 'iLevel', and back up

                // right now, the XDG field defined by 'NoJmpVec' should be a member 
                // of the aggregated XDG space on level 'iLevel';
                // so, there should be no inter-element jumps on the fine level, for each aggregated cell.
                // Let's test that!

                XDGField Test = new XDGField(xt.XB, "Test");
                xt.XdgMultigridOp.TransformSolFrom(Test.CoordinateVector, NoJmpVec);
                //xt.agg.Extrapolate(Test.Mapping);
                var aggGrd = MgSeq[iLevel];

                foreach (var spc in xt.LsTrk.SpeciesIdS) {
                    var Test_spc = Test.GetSpeciesShadowField(spc);
                    var SpcMask = xt.LsTrk.Regions.GetSpeciesMask(spc);

                    BitArray AggSourceBitmask = xt.agg.GetAgglomerator(spc).AggInfo.SourceCells.GetBitMask();

                    double Err = 0;
                    for (int jagg = 0; jagg < aggGrd.iLogicalCells.NoOfLocalUpdatedCells; jagg++) {
                        BitArray CompCellMask = new BitArray(Jup);
                        foreach (int jCell in aggGrd.iLogicalCells.AggregateCellToParts[jagg]) {
                            if(!AggSourceBitmask[jCell])
                                CompCellMask[jCell] = true;
                        }

                        SubGrid CompCellSubGrid = new SubGrid((new CellMask(grid, CompCellMask)).Intersect(SpcMask));

                        Err += JumpNorm(Test_spc, CompCellSubGrid.InnerEdgesMask).Pow2();
                    }

                    
                    Console.WriteLine("prolongation jump test (level {0}, species {2}): {1}", iLevel, Err, xt.LsTrk.GetSpeciesName(spc));
                    Assert.LessOrEqual(Err, 1.0e-8);
                }
            }
        }
    }
}
