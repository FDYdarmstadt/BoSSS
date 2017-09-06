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
using System.IO;
using System.Collections.Generic;
using System.Linq;

using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.Tecplot;
using MPI.Wrappers;
using ilPSP.Tracing;
using System.Globalization;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution;

namespace BoSSS.Application.XdgNastyLevsetLocationTest {

    /// <summary>
    /// A level-set integration test, when the level-set is located at nasty 
    /// positions (passing through corners, parallel to edges, etc.);
    /// In this test, we are more concerned on the correct handling of the 'pathologic' cases, so the integrand for the test is just 1.
    /// This choice makes the exact computation of the exact solution reasonable easy.
    /// </summary>
    class XdgNastyLevsetLocationTest : BoSSS.Solution.Application {

        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            BoSSS.Solution.Application._Main(args, true, "", delegate() {
                var p = new XdgNastyLevsetLocationTest();
                p.m_GridPartitioningType = GridPartType.METIS;
                return p;
            });
        }

        //internal ITest test = new Parallel(GetTestRange(),GetTestRange());
        //internal ITest test = new Parallel(GetTestRange(), new double[] { 0.0 });
        //internal ITest test = new Parallel(0.0, 1.0e-12);
        //internal ITest test = new Parallel(1.0e-12, 0.0);


        internal ITest test = new Schraeg(GetTestRange(), GetTestRange());
        //internal ITest test = new Schraeg( intercept: 5e-13);

        internal XQuadFactoryHelper.MomentFittingVariants momentFittingVariant =
            XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

        internal int QUAD_ORDER = 4;

        protected override GridCommons CreateOrLoadGrid() {
            return test.GetGrid();
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
        }

        LevelSet Phi;
        protected override void CreateFields() {
            Phi = new LevelSet(new Basis(this.GridData, 2), "LevelSet");
            base.LsTrk = new LevelSetTracker(this.GridData, 1, new string[] { "A", "B" }, Phi);
        }


        protected override void SetInitial() {
            Phi.ProjectField(test.GetLevelSet);
            this.LsTrk.UpdateTracker();
        }

        public static double[] GetTestRange() {
            List<double> ret = new List<double>();

            ret.AddRange(new double[] { 1e-64, 1e-32, 1e-24, 1e-16 });
            for (int i = 15; i >= 7; i--) {
                double[] b = new double[] { 1, 2, 5 };
                b.ScaleV(Math.Pow(10.0, -i));
                ret.AddRange(b);
            }

            var Neg = ret.ToArray();
            Neg.ScaleV(-1.0);
            ret.AddRange(Neg);

            ret.Insert(0, 0.0);

            return ret.ToArray();
        }

        internal bool IsPassed;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            
            IsPassed = true;
            try {
                int testcnt = 0;


                XQuadFactoryHelper.CheckQuadRules = true;

                while (this.test.NextTestCase()) {
                    testcnt++;
                    if (testcnt % 100 == 0)
                        Console.WriteLine("Test {0}", testcnt);

                    this.Phi.ProjectField(this.test.GetLevelSet);
                    this.LsTrk.UpdateTracker();

                    var schemes = new XQuadSchemeHelper(LsTrk, this.momentFittingVariant, LsTrk.SpeciesIdS.ToArray());
                    
                    var cutCells = LsTrk._Regions.GetCutCellSubGrid().VolumeMask;

                    var volSchemeA = schemes.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("A"), IntegrationDomain: cutCells);
                    var volSchemeB = schemes.GetVolumeQuadScheme(this.LsTrk.GetSpeciesId("B"), IntegrationDomain: cutCells);
                    CellQuadratureScheme surfScheme = schemes.GetLevelSetquadScheme(0, cutCells);
                    var edgeSchemeA = schemes.GetEdgeQuadScheme(this.LsTrk.GetSpeciesId("A"));
                    var edgeSchemeB = schemes.GetEdgeQuadScheme(this.LsTrk.GetSpeciesId("B"));


                    if (test.EdgeTestSupported) {
                        /*
                        var EdgeA = EdgeQuadrature(edgeSchemeA.Compile(this.GridDat, this.QUAD_ORDER));
                        var EdgeAref = EdgeAreaRef("A");
                        var EdgeAErr = EdgeA.CloneAs();
                        EdgeAErr.Acc(-1.0, EdgeAref);
                        double EdgeAL2err = EdgeAErr.L2Norm();
                        Console.WriteLine("Edge error for species A: " + EdgeAL2err);
                        */

                        var EdgeB = EdgeQuadrature(edgeSchemeB.Compile(this.GridData, this.QUAD_ORDER));
                        var EdgeBref = EdgeAreaRef("B");
                        var EdgeBErr = EdgeB.CloneAs();
                        EdgeBErr.Acc(-1.0, EdgeBref);
                        double EdgeBL2err = EdgeBErr.L2Norm();
                        if (EdgeBL2err >= 1.0e-5)
                            throw new ApplicationException("Edge error for species B: " + EdgeBL2err);
                        //Console.WriteLine("Edge error for species B: " + EdgeBL2err);

                    }

                    if (test.LevelsetTestSupported) {
                        var Surf = CellQuadrature(surfScheme.Compile(this.GridData, this.QUAD_ORDER));
                        var SurfRef = SurfAreaRef();
                        var surfErr = Surf.CloneAs();
                        surfErr.Acc(-1.0, SurfRef);
                        double surfL2err = surfErr.L2Norm();

                        if (surfL2err >= 1.0e-5)
                            //Console.WriteLine("Level-Set surface error " + surfL2err);
                            throw new ApplicationException("Level-Set surface error " + surfL2err);
                    }


                    if (test.VolumeTestSupported) {
                        /*
                        var VolA = CellQuadrature(volSchemeA.Compile(this.GridDat, this.QUAD_ORDER));
                        var VolAref = CellVolumeRef(cutCells, "A");
                        var volAErr = VolA.CloneAs();
                        volAErr.Acc(-1.0, VolAref);
                        double volAL2err = volAErr.L2Norm();
                        Console.WriteLine("Volume error for species A: " + volAL2err);
                        */

                        var VolB = CellQuadrature(volSchemeB.Compile(this.GridData, this.QUAD_ORDER));
                        var VolBref = CellVolumeRef(cutCells, "B");
                        var volBErr = VolB.CloneAs();
                        volBErr.Acc(-1.0, VolBref);
                        double volBL2err = volBErr.L2Norm();
                        if (volBL2err >= 1.0e-5)
                            throw new ApplicationException(string.Format("Volume error for species B: {0}", volBL2err));
                        //Console.WriteLine("Volume error for species B: {0}", volBL2err);
                    }


                    if (test.BoundaryPlusLevelsetTestSupported) {
                        var LevelSetArea = CellQuadrature(surfScheme.Compile(this.GridData, this.QUAD_ORDER));
                        var BoundaryB = EdgeQuadrature2CellBoundary(edgeSchemeB.Compile(this.GridData, this.QUAD_ORDER));
                        var RefAreaB = CellBoundaryPlusLevelsetAreaRef(cutCells, "B");

                        var ErrB = RefAreaB.CloneAs();
                        ErrB.Acc(-1.0, BoundaryB);
                        ErrB.Acc(-1.0, LevelSetArea);
                        double areaBL2err = ErrB.L2Norm();
                        Console.WriteLine("levelset + cell-boundary error for species B: " + areaBL2err);
                    }

                    //break;
                }
            } catch (Exception e) {
                Console.WriteLine(e.GetType().Name + ": " + e.Message);
                Console.WriteLine(e.StackTrace);
                IsPassed = false;
            }



            base.TerminationKey = true;
            return 0.0;
        }




        private MultidimensionalArray EdgeAreaRef(string Species) {
            int E = this.GridData.Edges.Count;
            var ret = MultidimensionalArray.Create(E);

            for (int iEdge = 0; iEdge < E; iEdge++) {
                ret[iEdge] = this.test.EdgeArea(iEdge, Species, this.GridData);
            }

            return ret;
        }

        private MultidimensionalArray SurfAreaRef() {
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            for (int jCell = 0; jCell < J; jCell++) {
                ret[jCell] = this.test.LevelsetArea(jCell, this.GridData);
            }

            return ret;
        }

        private MultidimensionalArray CellVolumeRef(CellMask cutCells, string Species) {
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            for (int jCell = 0; jCell < J; jCell++) {
                ret[jCell] = this.test.CellVolume(jCell, Species, this.GridData);
            }

            return ret;
        }

        private MultidimensionalArray CellBoundaryPlusLevelsetAreaRef(CellMask cutCells, string Species) {
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            for (int jCell = 0; jCell < J; jCell++) {
                ret[jCell] = this.test.CellBoundaryPlusLevelsetArea(jCell, Species, this.GridData);
            }

            return ret;
        }


        private MultidimensionalArray CellQuadrature(ICompositeQuadRule<QuadRule> surfRule) {
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);

            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                new int[] { 1 }, this.GridData,
                surfRule,
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    //Integrand.Evaluate(i0, Length, 0, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    EvalResult.SetAll(1.0);
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    var A = ret.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 });
                    var B = ResultsOfIntegration.ExtractSubArrayShallow(-1, 0);
                    A.Set(B);

                }).Execute();

            return ret;
        }


        private MultidimensionalArray EdgeQuadrature(ICompositeQuadRule<QuadRule> surfRule) {
            int E = this.GridData.Edges.Count;
            var ret = MultidimensionalArray.Create(E);

            BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, this.GridData,
                surfRule,
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    EvalResult.SetAll(1.0);
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    ret.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 })
                        .Set(ResultsOfIntegration.ExtractSubArrayShallow(-1, 0));
                }).Execute();

            return ret;
        }

        private MultidimensionalArray EdgeQuadrature2CellBoundary(ICompositeQuadRule<QuadRule> surfRule) {
            int J = this.GridData.Cells.NoOfLocalUpdatedCells;
            var ret = MultidimensionalArray.Create(J);
            var E2C = this.GridData.Edges.CellIndices;

            BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, this.GridData,
                surfRule,
                delegate(int i0, int Length, QuadRule Qr, MultidimensionalArray EvalResult) { // Evaluate
                    EvalResult.SetAll(1.0);
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        int iEdge = i0 + i;
                        int jCell0 = E2C[iEdge, 0];
                        int jCell1 = E2C[iEdge, 1];

                        ret[jCell0] += ResultsOfIntegration[i, 0];
                        if (jCell1 >= 0)
                            ret[jCell1] += ResultsOfIntegration[i, 0];
                    }
                }).Execute();

            return ret;
        }


        void PlotEdgeRule(ICompositeQuadRule<QuadRule> edgeRule) {
            var E2C = this.GridData.Edges.CellIndices;
            //var Trafos = this.GridData.Edges.Edge2CellTrafos;
            var trfIdx = this.GridData.Edges.Edge2CellTrafoIndex;

            using (TextWriter tw = new StreamWriter("edges.csv")) {

                foreach (var ChunkRule in edgeRule) {
                    var Nodes = ChunkRule.Rule.Nodes;
                    var Weights = ChunkRule.Rule.Weights;

                    for (int iEdge = ChunkRule.Chunk.i0; iEdge < ChunkRule.Chunk.JE; iEdge++) {


                        int jCell = E2C[iEdge, 0];
                        int iTrf = trfIdx[iEdge, 0];

                        //var NodesCell = Trafos[iTrf].Transform(Nodes);
                        NodeSet NodesCell = Nodes.GetVolumeNodeSet(this.GridData, iTrf);
                        var NodesGlob = MultidimensionalArray.Create(NodesCell.Lengths);

                        this.GridData.TransformLocal2Global(NodesCell, NodesGlob, jCell);

                        int K = NodesGlob.GetLength(0);
                        int D = NodesGlob.GetLength(1);

                        for (int k = 0; k < K; k++) {
                            tw.Write(iEdge);
                            tw.Write("\t");

                            for (int d = 0; d < D; d++) {
                                tw.Write(NodesGlob[k, d]);
                                tw.Write("\t");
                            }

                            tw.WriteLine(Weights[k]);
                        }
                    }
                }

            }
        }


        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            Tecplot.PlotFields(new DGField[] { this.Phi }, "NastyLevset", phystime, superSampling);
        }
    }
}

