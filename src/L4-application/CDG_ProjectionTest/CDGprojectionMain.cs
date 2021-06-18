using System;

using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using System.Collections;
using NUnit.Framework;
using BoSSS.Foundation.XDG;

namespace BoSSS.Application.CDG_ProjectionTest {


    class CDGprojectionMain : BoSSS.Solution.Application {

        static void Main(string[] args) {

            //BoSSS.Solution.Application.InitMPI();
            //var AUT = new BoSSS.Application.CDG_ProjectionTest.AllUpTest();
            //AUT.AllUp(3, 2, 4);
            //BoSSS.Solution.Application.FinalizeMPI();
            //Assert.IsTrue(false, "Remove me");
            //return;

            BoSSS.Solution.Application._Main(args, true, delegate () {
                CDGprojectionMain p = new CDGprojectionMain();
                return p;
            });
        }

        internal int dimension = 3;
        internal int degree = 3;
        internal int gridResolution = 6;

        internal bool periodicX = false;
        internal bool periodicY = false;
        internal bool periodicZ = false;


        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if (MPISize <= 1)
                return null;

            // for debugging purposes
            //int i0 = this.Grid.CellPartitioning.i0;
            //int iE = this.Grid.CellPartitioning.iE;
            //int[] CellDist = new int[iE - i0];
            //for (int i = i0; i < iE; i++) {
            //    switch (MPIRank) {
            //        case 0: {
            //                if (i == i0) {
            //                    CellDist[i - i0] = 1;
            //                } else {
            //                    CellDist[i - i0] = 0;
            //                }
            //                break;
            //            }
            //        case 1: {
            //                CellDist[i - i0] = 1;
            //                break;
            //            }
            //    }
            //}
            //return CellDist;

            return base.ComputeNewCellDistribution(TimeStepNo, physTime);
        }


        protected override IGrid CreateOrLoadGrid() {

            double[] nodes = GenericBlas.Linspace(0, 1, (2 * gridResolution) + 1);
            double[] node2pi = GenericBlas.Linspace(0, 2 * Math.PI, (2 * gridResolution) + 1);
            double[] node2 = GenericBlas.Linspace(0, (1.0 / 2.0 * gridResolution), 3);
            double[] node3 = GenericBlas.Linspace(0, 1, (3 * gridResolution) + 1);

            double[] droplet_xy = GenericBlas.Linspace(0, 3, (2 * gridResolution) + 1);
            double[] droplet_z = GenericBlas.Linspace(-3, 3, (4 * gridResolution) + 1);

            GridCommons grid;
            if (dimension == 2)
                grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: periodicX, periodicY: periodicY);
                //grid = Grid2D.Cartesian2DGrid(node2pi, node2, periodicX: periodicX, periodicY: periodicY);
            else if (dimension == 3)
                //grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes, periodicX: periodicX, periodicY: periodicY, periodicZ: periodicZ);
                grid = Grid3D.Cartesian3DGrid(droplet_xy, droplet_xy, droplet_z, periodicX: periodicX, periodicY: periodicY, periodicZ: periodicZ);
            else
                throw new NotSupportedException("Not supported spatial dimension");

            return grid;

        }

        SinglePhaseField origin;
        SinglePhaseField result0;
        SinglePhaseField result1;

        ConstrainedDGField cdgField0;   // projection on the same DG basis
        ConstrainedDGField cdgField1;   // projection on a DG basis with degree + 1


        protected override void CreateFields() {

            Basis dgBasis = new Basis(this.GridData, degree);
            Basis cdgBasis0 = new Basis(this.GridData, degree);
            Basis cdgBasis1 = new Basis(this.GridData, degree + 1);

            origin = new SinglePhaseField(dgBasis, "origin");
            cdgField0 = new ConstrainedDGField(cdgBasis0);
            cdgField1 = new ConstrainedDGField(cdgBasis1);
            result0 = new SinglePhaseField(cdgBasis0, "result0");
            result1 = new SinglePhaseField(cdgBasis1, "result1");

        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {

        }

        public bool passed = true;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {

            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            BitArray oneBit = new BitArray(J);
            oneBit[0] = true;
            CellMask oneCell = new CellMask(this.GridData, oneBit);


            if (this.dimension == 2) {
                //// for debugging purposes
                //CellMask msk2D = CellMask.GetCellMask((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData),
                //    X => (X[0] > 0.0 && X[0] < 4.0 && X[1] > 0.0 && X[1] < 1.0));


                //Console.WriteLine("Test 2D projection function 1: sin(x)");
                //Func<double[], double> projFunc = X => Math.Sin(X[0]);
                //Console.WriteLine("project on full mask");
                //string name_disc = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-func1";
                //passed &= ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), null, name_disc);


                Console.WriteLine("Test 2D projection function 1: sin(x) + cos(x) + x - cos(y) - 1");
                Func<double[], double> projFunc = X => Math.Sin(X[0]) + Math.Cos(X[0]) + X[0] - (Math.Cos(X[1]) + 1);
                Console.WriteLine("project on full mask");
                string name_disc = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-func1";
                passed &= ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), null, name_disc);


                //// should be exact for p >= 3
                //Console.WriteLine("Test 2D projection function 2: x^2 + y^3 - xy");
                //projFunc = X => (X[0] * X[0]) + (X[1] * X[1] * X[1]) - (X[0] * X[1]);
                //Console.WriteLine("project on full mask");
                //name_disc = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-func2";
                //passed &= ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), null, name_disc);


            }

            if (this.dimension == 3) {

                //Console.WriteLine("Test 3D projection function 1: sin(x) + cos(x) + x - y + sin(z) - 1");
                //Func<double[], double> projFunc = X => (Math.Sin(X[0]) + Math.Cos(X[0]) + X[0] - (X[1] + 1) + Math.Sin(X[2]));
                //Console.WriteLine("project on full mask");
                //string name_disc = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-func1";
                //passed &= ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), null, name_disc);


                Console.WriteLine("Test 3D projection function: Legendre Polynomial");
                double r_0 = 1;
                double a_P = 0.5;  
                double a_0 = 0.94754;   
                Func<double[], double> projFunc = delegate (double[] X) {
                    double r = ((X[0]).Pow2() + (X[1]).Pow2() + (X[2]).Pow2()).Sqrt();
                    double r_xy = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                    double theta = Math.Atan2(r_xy, -X[2]);
                    double f = r_0 * (a_0 + a_P * 0.5 * (3.0 * (Math.Cos(theta)).Pow2() - 1.0));   
                    double phi = r - f;
                    return phi;
                };
                //Console.WriteLine("project on full mask");

                SinglePhaseField phiField = new SinglePhaseField(new Basis(this.GridData, degree));
                phiField.ProjectField(projFunc);
                var LevSet = new LevelSet(phiField.Basis, "LevelSet");
                LevSet.Acc(1.0, phiField);
                var LsTrk = new LevelSetTracker((GridData)phiField.GridDat, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, LevSet);
                LsTrk.UpdateTracker(0.0);
                CellMask near = LsTrk.Regions.GetNearFieldMask(1);
                //SinglePhaseField nearField = new SinglePhaseField(new Basis(this.GridData, 0));
                //nearField.AccConstant(1.0, near);
                //Tecplot.PlotFields(new DGField[]{ nearField }, "CDGproj_nearField", 0.0, 0);
                Console.WriteLine("project on near field mask; no of cells {0}", near.NoOfItemsLocally);

                string name_disc = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-funcL";
                passed &= ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), near, name_disc);


                //// should be exact for p >= 3
                //Console.WriteLine("Test 3D projection function 2: x^2 + y^3 - zx + x - z");
                //projFunc = X => (X[0] * X[0]) + (X[1] * X[1] * X[1]) - (X[2] * X[0]) + X[0] - X[2];
                //Console.WriteLine("project on full mask");
                //name_disc = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-func2";
                //passed &= ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), null, name_disc);


            }

            base.TerminationKey = true;
            return 0.0;

        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 3) {
            Tecplot.PlotFields(new DGField[] { this.origin, this.result0, this.result1 }, "CDGproj_" + timestepNo, physTime, 3);
        }


        bool ProjectFieldAndEvaluate(ScalarFunction func, CellMask domain, string name_disc) {

            bool _passed = true;

            origin.Clear();
            result0.Clear();
            result1.Clear();

            origin.ProjectField(func);
            double L2jumpOrigin = JumpNorm(origin, domain);
            Console.WriteLine("L2 jump origin field = {0}", L2jumpOrigin);


            // project and check cdgField0
            var returnFields = cdgField0.ProjectDGField(origin, domain);
            //Tecplot.PlotFields(returnFields, "CDGproj_patchField", 0.0, 3);
            cdgField0.AccToDGField(1.0, result0, domain);

            var errField = origin.CloneAs();
            errField.Acc(-1.0, result0, domain);

            double L2err0 = errField.L2Norm(domain);
            double L2jump = JumpNorm(result0, domain);

            bool checkL2err0 = (gridResolution > 2) ? (L2err0 < 1.0e-2) : true;
            if (checkL2err0 && L2jump < 1.0e-12) {
                Console.WriteLine("projection0 PASSED");
                _passed &= true;
            } else {
                Console.WriteLine("projection0 FAILED");
                _passed &= false;
                Console.WriteLine("L2err0 = {0}; L2jump = {1}", L2err0, L2jump);
                PlotCurrentState(0.0, new TimestepNumber(1, 0), 3);
            }
            Console.WriteLine("L2 jump result0 field = {0}; L2 error norm = {1}", L2jump, L2err0);

            //// project and check cdgField1
            //cdgField1.ProjectDGField(origin, domain);
            //cdgField1.AccToDGField(1.0, result1, domain);

            //errField = origin.CloneAs();
            //errField.AccLaidBack(-1.0, result1, domain);

            //double L2err1 = errField.L2Norm(domain);
            //L2jump = JumpNorm(result1, domain);

            //if (L2err1 <= L2err0 && L2jump < 1.0e-12) {
            //    Console.WriteLine("projection1 PASSED");
            //    _passed &= true;
            //} else {
            //    Console.WriteLine("projection1 FAILED");
            //    _passed &= false;
            //    Console.WriteLine("L2err1 = {0}; L2err0 = {1}; L2jump = {2}", L2err1, L2err0, L2jump);
            //    PlotCurrentState(0.0, new TimestepNumber(1, 1), 3);
            //}
            //Console.WriteLine("L2 jump result1 field = {0}", L2jump);


            //// check parallel simulations
            //// ==========================
            //int RefMPIsize = 1;
            //{
            //    var projCheck = new TestingIO(this.GridData, $"CDG_Projection-{name_disc}.csv", RefMPIsize);
            //    projCheck.AddDGField(this.result0);
            //    projCheck.AddDGField(this.result1);
            //    projCheck.DoIOnow();

            //    Assert.Less(projCheck.AbsError(this.result0), 1.0e-15, "Mismatch in projected result0 between single-core and parallel run.");
            //    Assert.Less(projCheck.AbsError(this.result1), 1.0e-15, "Mismatch in projected result1 between single-core and parallel run.");

            //}

            return _passed;
        }


        static double JumpNorm(DGField f, CellMask mask = null) {
            GridData grd = (GridData)f.GridDat;
            int D = grd.SpatialDimension;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;

            if (mask == null) {
                mask = CellMask.GetFullMask(grd);
            }
            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            f.MPIExchange();

            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { 1 }, grd,
                (new EdgeQuadratureScheme(true, innerEM)).Compile(grd, f.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++) {
                        int iEdge = j + i0;
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if (jCell_OT >= 0) {

                            int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT);

                            f.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            f.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);

                            //if (uDiff.L2Norm() > 1e-10)
                            //    Console.WriteLine("uDiff at edge {0} between cell {1} and cell {2}: {3}", iEdge, jCell_IN, jCell_OT, uDiff.L2Norm());
                        } else {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum();

            return Unorm.Sqrt();
        }


    }


}
