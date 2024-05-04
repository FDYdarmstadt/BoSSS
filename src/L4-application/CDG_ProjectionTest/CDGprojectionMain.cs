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

using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation;
using BoSSS.Foundation.ConstrainedDGprojection;
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
using ilPSP.Tracing;
using System.Collections.Generic;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System.Linq;

namespace BoSSS.Application.CDG_ProjectionTest {

    /// <summary>
    /// Tests for <see cref="BoSSS.Foundation.ConstrainedDGprojection.ConstrainedDGFieldMk3"/>
    /// </summary>
    public class CDGprojectionMain : BoSSS.Solution.Application {

        static void Main(string[] args) {

            BoSSS.Solution.Application.InitMPI();
            BoSSS.Application.CDG_ProjectionTest.AllUpTest.AllUp_patchwiseOnly_case2_sameDegree_p3(3, 4, 2);
            ////BoSSS.Application.CDG_ProjectionTest.AllUpTest.AllUp(2, 3, 2, 8, true, ProjectionStrategy.patchwiseOnly);
            ////BoSSS.Application.CDG_ProjectionTest.AllUpTest.AllUp(1, 3, 2, 4, false, ProjectionStrategy.patchwiseOnly);
            //BoSSS.Application.CDG_ProjectionTest.AllUpTest.AllUp(2, 2, 2, 2, 2, true, ProjectionStrategy.globalOnly);
            ////var AUT = new BoSSS.Application.CDG_ProjectionTest.AllUpTest();
            ////AUT.AllUp(4, 3, 4, 2, false, ProjectionStrategy.globalOnly);
            //////AUT.AllUp_Cube(3, 4, 4, true);
            //BoSSS.Solution.Application.FinalizeMPI();
            Assert.IsTrue(false, "Remove me");
            return;

            BoSSS.Solution.Application._Main(args, true, delegate () {
                CDGprojectionMain p = new CDGprojectionMain();
                return p;
            });
        }

        internal int dimension = 2;
        internal int gridResolution = 4;
        internal int AMRlevel = 0;
        internal int degree = 2;
        internal bool projectOnSameBasis = false;   // if false projection on a DG basis with degree + 1
        internal int projectionCase = 0;

        internal ProjectionStrategy projectStrategy = ProjectionStrategy.globalOnly;


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

        double[] refinementPoint;


        protected override IGrid CreateOrLoadGrid() {

            double[] nodes = GenericBlas.Linspace(0, 1, (2 * gridResolution) + 1);
            double[] node2pi = GenericBlas.Linspace(0, 2 * Math.PI, (2 * gridResolution) + 1);
            double[] node2piy = GenericBlas.Linspace(0, (2 * Math.PI / gridResolution), 3);

            double[] node2 = GenericBlas.Linspace(0, (1.0 / gridResolution), 3);    // for debug purposes
            //double[] node3 = GenericBlas.Linspace(0, 1, (3 * gridResolution) + 1);    

            double[] droplet_xy = GenericBlas.Linspace(0, 1.5, (2 * gridResolution) + 1);
            double[] droplet_z = GenericBlas.Linspace(-1.5, 1.5, (4 * gridResolution) + 1);

            double[] cube = GenericBlas.Linspace(-1.0, 1.0, (2 * gridResolution) + 1);

            GridCommons grid;
            if (dimension == 2) {
                switch (projectionCase) { 
                    case 0: {
                        grid = Grid2D.Cartesian2DGrid(node2pi, node2piy, periodicX: true, periodicY: false);
                        //grid = Grid2D.Cartesian2DGrid(node2pi, node2, periodicX: true, periodicY: false);
                        break;
                    }
                    case 1:
                    case 2: {
                        grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: false, periodicY: false);
                        //grid = Grid2D.Cartesian2DGrid(nodes, node2, periodicX: false, periodicY: false);
                        refinementPoint = new double[] { 0.55, 0.55 };
                        break;
                    }
                    case 3:
                    case 4:
                    case 5: {
                        grid = Grid2D.Cartesian2DGrid(droplet_z, droplet_z, periodicX: false, periodicY: false);
                        refinementPoint = new double[] { 0.6, 0.5 };
                        break;
                    }
                    case 6: {
                        grid = Grid2D.Cartesian2DGrid(cube, cube, periodicX: false, periodicY: false);
                        refinementPoint = new double[] { 0.25, 0.25 };
                        break;
                    }
                    default:
                        throw new ArgumentException("no such projection case available");
                }
            } else if (dimension == 3) {
                switch (projectionCase) {
                    case 1:
                    case 2: {
                        grid = Grid3D.Cartesian3DGrid(nodes, nodes, nodes, periodicX: false, periodicY: false, periodicZ: false);
                            refinementPoint = new double[] { 0.55, 0.55, 0.55 };
                            break;
                    }
                    case 3:
                    case 4:
                    case 5: {
                        grid = Grid3D.Cartesian3DGrid(droplet_xy, droplet_xy, droplet_z, periodicX: false, periodicY: false, periodicZ: false);
                            refinementPoint = new double[] { 0.6, 0.6, 0.5 };
                            break;
                    }
                    case 6: {
                        grid = Grid3D.Cartesian3DGrid(cube, cube, cube, periodicX: false, periodicY: false, periodicZ: false);
                            refinementPoint = new double[] { 0.25, 0.25, 0.25 };
                            break;
                    }
                    default:
                        throw new ArgumentException("no such projection case available");
                }
            } else
                throw new NotSupportedException("Not supported spatial dimension");

            return grid;

        }


        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            using (var tr = new FuncTrace()) {

                if (AMRlevel > 0) {

                    // Check grid changes
                    // ==================

                    int[] desiredLevels = new int[this.GridData.CellPartitioning.TotalLength];
                    this.GridData.LocatePoint(refinementPoint, out long GlobalID, out long GlobalIndex, out bool IsInside, out bool OnThisProc);
                    desiredLevels[GlobalIndex] = AMRlevel;

                    GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData, null);
                    bool AnyChange = gridRefinementController.ComputeGridChange(desiredLevels, out List<int> CellsToRefineList, out List<int[]> Coarsening);

                    int NoOfCellsToRefine = 0;
                    int NoOfCellsToCoarsen = 0;
                    if (AnyChange.MPIOr()) {
                        int[] glb = (new int[] { CellsToRefineList.Count, Coarsening.Sum(L => L.Length) }).MPISum();
                        NoOfCellsToRefine = glb[0];
                        NoOfCellsToCoarsen = glb[1];
                    }
                    long oldJ = this.GridData.CellPartitioning.TotalLength;

                    // Update Grid
                    // ===========
                    if (AnyChange.MPIOr()) {

                        Console.WriteLine(" Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                        Console.WriteLine(" Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                        newGrid = ((GridData)this.GridData).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                    } else {

                        newGrid = null;
                        old2NewGrid = null;
                    }

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }

            }
        }



        SinglePhaseField origin;
        SinglePhaseField result;


        Basis cdgBasis;

        protected override void CreateFields() {

            Basis dgBasis = new Basis(this.GridData, degree);
            int cdgDegree = projectOnSameBasis ? degree : degree + 1;
            cdgBasis = new Basis(this.GridData, cdgDegree);

            origin = new SinglePhaseField(dgBasis, "origin");
            
            result = new SinglePhaseField(cdgBasis, "result");

        }

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
        }


        protected override void SetInitial(double time) {

            this.Control.AdaptiveMeshRefinement = (this.AMRlevel > 0) ? true : false;
            this.Control.AMR_startUpSweeps = this.AMRlevel;

            base.SetInitial(time);
        }


        public bool passed = false;

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {


            if (this.dimension == 2) 
                passed = this.SetUp2DCasesAndRun();
 

            if (this.dimension == 3) 
                passed = this.SetUp3DCasesAndRun();


            base.TerminationKey = true;
            return 0.0;

        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 3) {
            Tecplot.PlotFields(new DGField[] { this.origin, this.result }, "CDGproj_" + timestepNo, physTime, 3);
        }


        bool SetUp2DCasesAndRun() {

            Func<double[], double> projFunc;
            CellMask mask = null;
            string caseName;

            switch (projectionCase) {
                case 0: {
                    Console.WriteLine("Test 2D projection function 1: sin(x)");
                    projFunc = X => Math.Sin(X[0]);
                    Console.WriteLine("project on full mask");
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case0";
                    break;
                }
                case 1: {
                    Console.WriteLine("Test 2D projection function 1: sin(x) + cos(x) + x - cos(y) - 1");
                    projFunc = X => Math.Sin(X[0]) + Math.Cos(X[0]) + X[0] - (Math.Cos(X[1]) + 1);
                    Console.WriteLine("project on full mask");
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case1";
                    break;
                }
                case 2: {
                    // should be exact for p >= 3
                    Console.WriteLine("Test 2D projection function 2: x^2 + y^3 - xy");
                    projFunc = X => (X[0] * X[0]) + (X[1] * X[1] * X[1]) - (X[0] * X[1]);
                    Console.WriteLine("project on full mask");
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case2";
                    break;
                }
                case 3: {
                    Console.WriteLine("Test 2D projection function: Legendre Polynomial second mode");
                    double r_0 = 1;
                    double a_P = 0.5;       //initial disturbance to corresponding Legendre polynomial P
                    double a_0 = 0.94754;   // for a_2 = 0.5 and r_0 = 1 (script available in maple)
                    projFunc = delegate (double[] X) {
                        double r = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                        double theta = Math.Atan2(X[0], -X[1]);
                        //double f = r_0;
                        double f = r_0 * (a_0 + a_P * 0.5 * (3.0 * (Math.Cos(theta)).Pow2() - 1.0));                                        // P_2
                        double phi = r - f;
                        return phi;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case3";
                    break;
                }
                case 4: {
                    Console.WriteLine("Test 2D projection function: Legendre Polynomial third mode");
                    double r_0 = 1;
                    double a_P = 0.5;       //initial disturbance to corresponding Legendre polynomial P
                    double a_0 = 0.9643;   // for a_3 = 0.5 and r_0 = 1 (script available in maple)
                    projFunc = delegate (double[] X) {
                        double r = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                        double theta = Math.Atan2(X[0], -X[1]);
                        double f = r_0 * (a_0 + a_P * 0.5 * (5.0 * (Math.Cos(theta)).Pow(3) - 3.0 * Math.Cos(theta)));                      // P_3
                        double phi = r - f;
                        return phi;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case4";
                    break;
                }
                case 5: {
                    Console.WriteLine("Test 2D projection function: Legendre Polynomial fourth mode");
                    double r_0 = 1;
                    double a_P = 0.5;       //initial disturbance to corresponding Legendre polynomial P
                    double a_0 = 0.97146;   // for a_4 = 0.5 and r_0 = 1 (script available in maple)
                    projFunc = delegate (double[] X) {
                        double r = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                        double theta = Math.Atan2(X[0], -X[1]);
                        double f = r_0 * (a_0 + a_P * 0.125 * (35.0 * (Math.Cos(theta)).Pow(4) - 30.0 * (Math.Cos(theta)).Pow(2) + 3.0));   // P_4
                        double phi = r - f;
                        return phi;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case5";
                    break;
                }
                case 6: {
                    Console.WriteLine("Test 2D projection function: Cube");
                    projFunc = delegate (double[] X) {
                        double[] pos = new double[2];
                        double anglev = 10;
                        double t = 0;
                        double angle = -(anglev * t) % (2 * Math.PI);
                        double particleRad = 0.261;
                        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                            Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                            + particleRad;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case6";
                    break;
                }
                default:
                    throw new ArgumentOutOfRangeException("no such projection case available");
            }

            return ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), mask, caseName);

        }


        bool SetUp3DCasesAndRun() {

            Func<double[], double> projFunc;
            CellMask mask = null;
            string caseName;

            switch (projectionCase) {
                case 1: {
                    Console.WriteLine("Test 3D projection function 1: sin(x) + cos(x) + x - y + sin(z) - 1");
                    projFunc = X => (Math.Sin(X[0]) + Math.Cos(X[0]) + X[0] - (X[1] + 1) + Math.Sin(X[2]));
                    Console.WriteLine("project on full mask");
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case1";
                    break;
                }
                case 2: {
                    // should be exact for p >= 3
                    Console.WriteLine("Test 3D projection function 2: x^2 + y^3 - zx + x - z");
                    projFunc = X => (X[0] * X[0]) + (X[1] * X[1] * X[1]) - (X[2] * X[0]) + X[0] - X[2];
                    Console.WriteLine("project on full mask");
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case2";
                    break;
                }
                case 3: {
                    Console.WriteLine("Test 3D projection function: Legendre Polynomial second mode");
                    double r_0 = 1;
                    double a_P = 0.5;       //initial disturbance to corresponding Legendre polynomial P
                    double a_0 = 0.94754;   // for a_2 = 0.5 and r_0 = 1 (script available in maple)
                    projFunc = delegate (double[] X) {
                        double r = ((X[0]).Pow2() + (X[1]).Pow2() + (X[2]).Pow2()).Sqrt();
                        double r_xy = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                        double theta = Math.Atan2(r_xy, -X[2]);
                        //double f = r_0;
                        double f = r_0 * (a_0 + a_P * 0.5 * (3.0 * (Math.Cos(theta)).Pow2() - 1.0));                                        // P_2
                        double phi = r - f;
                        return phi;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case3";
                    break;
                }
                case 4: {
                    Console.WriteLine("Test 3D projection function: Legendre Polynomial third mode");
                    double r_0 = 1;
                    double a_P = 0.5;       //initial disturbance to corresponding Legendre polynomial P
                    double a_0 = 0.9643;   // for a_3 = 0.5 and r_0 = 1 (script available in maple)
                    projFunc = delegate (double[] X) {
                        double r = ((X[0]).Pow2() + (X[1]).Pow2() + (X[2]).Pow2()).Sqrt();
                        double r_xy = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                        double theta = Math.Atan2(r_xy, -X[2]);
                        double f = r_0 * (a_0 + a_P * 0.5 * (5.0 * (Math.Cos(theta)).Pow(3) - 3.0 * Math.Cos(theta)));                      // P_3
                        double phi = r - f;
                        return phi;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case4";
                    break;
                }
                case 5: {
                    Console.WriteLine("Test 3D projection function: Legendre Polynomial fourth mode");
                    double r_0 = 1;
                    double a_P = 0.5;       //initial disturbance to corresponding Legendre polynomial P
                    double a_0 = 0.97146;   // for a_4 = 0.5 and r_0 = 1 (script available in maple)
                    projFunc = delegate (double[] X) {
                        double r = ((X[0]).Pow2() + (X[1]).Pow2() + (X[2]).Pow2()).Sqrt();
                        double r_xy = ((X[0]).Pow2() + (X[1]).Pow2()).Sqrt();
                        double theta = Math.Atan2(r_xy, -X[2]);
                        double f = r_0 * (a_0 + a_P * 0.125 * (35.0 * (Math.Cos(theta)).Pow(4) - 30.0 * (Math.Cos(theta)).Pow(2) + 3.0));   // P_4
                        double phi = r - f;
                        return phi;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case5";
                    break;
                }
                case 6: {
                    Console.WriteLine("Test 3D projection function: Cube");
                    projFunc = delegate (double[] X) {
                        double[] pos = new double[3];
                        double anglev = 10;
                        double t = 0;
                        double angle = -(anglev * t) % (2 * Math.PI);
                        double particleRad = 0.261;
                        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                                Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)),
                                                Math.Abs(X[2] - pos[2])))
                                                + particleRad;
                    };
                    mask = GetNearbandMask(NonVectorizedScalarFunction.Vectorize(projFunc));
                    caseName = $"dim{this.dimension}-deg{this.degree}-grdRes{this.gridResolution}-case6";
                    break;
                }
                default:
                    throw new ArgumentOutOfRangeException("no such projection case available");
            }

            return ProjectFieldAndEvaluate(NonVectorizedScalarFunction.Vectorize(projFunc), mask, caseName);

        }


        CellMask GetNearbandMask(ScalarFunction func) {

            SinglePhaseField phiField = new SinglePhaseField(new Basis(this.GridData, degree));
            phiField.ProjectField(func);
            var LevSet = new LevelSet(phiField.Basis, "LevelSet");
            LevSet.Acc(1.0, phiField);
            var LsTrk = new LevelSetTracker((GridData)phiField.GridDat, XQuadFactoryHelper.MomentFittingVariants.Classic, 1, new string[] { "A", "B" }, LevSet);
            LsTrk.UpdateTracker(0.0);
            CellMask near = LsTrk.Regions.GetNearFieldMask(1);
            //SinglePhaseField nearField = new SinglePhaseField(new Basis(this.GridData, 0));
            //nearField.AccConstant(1.0, near);
            //Tecplot.PlotFields(new DGField[]{ nearField }, "CDGproj_nearField", 0.0, 0);
            Console.WriteLine("project on near field mask; no of cells {0}", near.NoOfItemsLocally.MPISum());
            return near;
        }


        bool ProjectFieldAndEvaluate(ScalarFunction func, CellMask domain, string name_disc) {

            bool _passed = false;

            origin.Clear();
            result.Clear();

            origin.ProjectField(func);
            double L2jumpOrigin = JumpNorm(origin, domain);
            Console.WriteLine("L2 jump origin field = {0}", L2jumpOrigin);


            // project and check cdgField
            using(var cdgField = ConstrainedDGFieldMk3.Factory(cdgBasis, domain, this.projectStrategy)) {
                //var returnFields = cdgField.ProjectDGField(origin, domain);
                cdgField.ProjectDGField(origin);
                //if (!returnFields.IsNullOrEmpty())
                //    Tecplot.PlotFields(returnFields, "CDGproj_patchField", 0.0, 3);
                cdgField.AccToDGField(1.0, result, domain);
            }

            var errField = origin.CloneAs();
            errField.AccLaidBack(-1.0, result, domain);

            double L2err = errField.L2Norm(domain);
            double L2jump = JumpNorm(result, domain);

            bool checkL2err = (gridResolution > 4 && degree > 2) ? (L2err < 1.0e-2) : true;
            Console.WriteLine("========================");
            if (checkL2err && L2jump < 1.0e-8) {
                Console.WriteLine("// projection PASSED //");
                _passed = true;
            } else {
                Console.WriteLine("// projection FAILED //");
                _passed = false;
                //PlotCurrentState(0.0, new TimestepNumber(1, 0), 3);
            }
            Console.WriteLine("L2 jump result field = {0}; L2 error norm = {1}", L2jump, L2err);
            Console.WriteLine("========================");


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

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN, false);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT, false);

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
