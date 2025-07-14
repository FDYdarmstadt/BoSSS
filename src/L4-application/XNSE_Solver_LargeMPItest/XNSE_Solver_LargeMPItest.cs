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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using ilPSP;
using System.Diagnostics;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Linq;


namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Tests whether the XNSE solver (<see cref="XNSE_SolverMain"/>) also works MPI-parallel for larger cases (8 MPI cores)
    /// </summary>
    [TestFixture]
    public static class XNSE_Solver_LargeMPItest {

        [Test]
        static public void ParallelRotatingSphere() {
            var C = XNSE_Solver_MPItest.Rotating_Sphere(k: 1, Res: 20, SpaceDim: 3, useAMR: true, useLoadBal: true);
            //C.TracingNamespaces = "*";

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
            //Assert.IsTrue(false, "this should fail!");
        }

        [Test]
        public static void ParallelRotatingTilted3DTorus() {
            var C = HardcodedControl.RotatingTiltedXRigid(k: 1, Res: 20, SpaceDim: 3, AMR: true, AMRLevel: 1, TiltAngle: Math.PI / 4, SolverOn: true);

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        [Test]
        public static void RotatingTilted3DTorusAgg0() {
            var C = HardcodedControl.RotatingTilted3DTorusAgg0();

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        /// <summary>
        /// This test case checks the source to source cell agglomeration with METIS domain decomposition
        /// </summary>
        [Test]
        public static void SourceToSourceAgglomerationMETIS() {
            var C = HardcodedControl.TwoTorusesAggTestCase();
            C.GridPartType = GridPartType.METIS;
            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        /// <summary>
        /// This test case checks the source to source cell agglomeration with Hilbert curves
        /// </summary>
        [Test]
        public static void SourceToSourceAgglomerationHILBERT() {
            var C = HardcodedControl.TwoTorusesAggTestCase();
            C.GridPartType = GridPartType.Hilbert;
            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();
            }
        }

        static XNSE_Control DropletOnPlate_AMRtest(int D = 2, int p = 2, int kelem = 16) {

            XNSE_Control C = new XNSE_Control();




            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "AMR test";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.SetDGdegree(p);

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            //C.PhysicalParameters.betaS_A = 0.05;
            //C.PhysicalParameters.betaS_B = 0.05;

            //C.PhysicalParameters.betaL = 0;
            //C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            double xSize = 0.3;
            double ySize = 0.3;
            double zSize = 0.3;

            if(D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize * 0.5, xSize * 0.5, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-ySize * 0.5, ySize * 0.5, kelem + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.DefineEdgeTags(delegate (double[] X) {
                        if(Math.Abs(X[1] + ySize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_lower";
                        if(Math.Abs(X[1] - ySize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_upper";
                        if(Math.Abs(X[0] + xSize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_left";
                        if(Math.Abs(X[0] - xSize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_right";

                        throw new ArgumentOutOfRangeException("unable to determine edge name.");
                    });

                    return grd;
                };
            } else if(D == 3) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize * 0.5, xSize * 0.5, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-ySize * 0.5, ySize * 0.5, kelem + 1);
                    double[] Znodes = GenericBlas.Linspace(-zSize * 0.5, zSize * 0.5, kelem + 1);
                    var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    grd.EdgeTagNames.Add(5, "navierslip_linear_front");
                    grd.EdgeTagNames.Add(6, "navierslip_linear_back");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[2] + zSize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_lower";
                        if(Math.Abs(X[2] - zSize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_upper";
                        if(Math.Abs(X[0] + xSize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_left";
                        if(Math.Abs(X[0] - xSize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_right";
                        if(Math.Abs(X[1] + ySize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_front";
                        if(Math.Abs(X[1] - ySize * 0.5) <= 1.0e-8)
                            return "navierslip_linear_back";

                        throw new ArgumentOutOfRangeException("unable to determine edge name.");
                    });

                    return grd;
                };
            } else {
                throw new ArgumentOutOfRangeException("unsupported spatial dimension: " + D);
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.1;
            //double Theta_e = Math.PI / 2.0;
            //double s = 2 * R * Math.Sin(Theta_e);
            //double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            Func<double[], double, double> PhiFunc = ((X, t) => -1.0);
            double prescribedVel = 0.0;

            if(D == 2) {
                double[] center_0 = [0, 0];
                prescribedVel = ySize / 2.1;
                PhiFunc = ((X, t) => ((X[0] - (center_0[0] + (0.0 * prescribedVel))).Pow2() + (X[1] - (center_0[1] - (t * prescribedVel))).Pow2()).Sqrt() - R);
            } else if(D == 3) {
                double[] center_0 = [0, 0, 0];
                prescribedVel = zSize / 2.0;
                PhiFunc = ((X, t) => ((X[0] - center_0[0]).Pow2() + (X[1] - center_0[1]).Pow2() + (X[2] - (center_0[2] - (t * prescribedVel))).Pow2()).Sqrt() - R);
            } else {
                throw new ArgumentOutOfRangeException("unsupported spatial dimension: " + D);
            }

            C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0.0));

            C.InitialValues_Evaluators.Add("Pressure#A", X => 0.0);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            /*
            C.AddBoundaryValue("navierslip_linear_lower");
            C.AddBoundaryValue("navierslip_linear_upper");
            C.AddBoundaryValue("navierslip_linear_left");
            C.AddBoundaryValue("navierslip_linear_right");

            if(D == 3) {
                C.AddBoundaryValue("navierslip_linear_front");
                C.AddBoundaryValue("navierslip_linear_back");
            }
            */
            #endregion


            // exact solution
            // ==============

            C.Phi = PhiFunc;

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);



            // misc. solver options
            // ====================
            #region solver



            //C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.Prescribed;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            //C.RefineStrategy = XNSE_Control.RefinementStrategy.ContactLineRefined;
            //C.RefinementLevel = 1;

            #endregion

            // Adaptive Mesh Refinement (AMR)
            // ==============================

            #region amr_config

            C.AdaptiveMeshRefinement = true;
            C.AMR_startUpSweeps = 3;
            C.activeAMRlevelIndicators.Add(
                new AMRonNarrowband() { bandwidth = 0, maxRefinementLevel = 2 }
                );

            #endregion

            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.LieSplitting;

            C.TimesteppingMode = compMode;

            double dt = prescribedVel / 50.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 5;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }




        /// <summary>
        /// Testing adaptive mesh refinement in a MPI-parallel environment
        /// </summary>
        [Test]
        public static void AMRtest_2D([Values(GridPartType.METIS, GridPartType.Hilbert, GridPartType.clusterHilbert)] GridPartType gridPartType) {
            var C = DropletOnPlate_AMRtest(D: 2, p: 2, kelem: 8);
            C.NoOfTimesteps = 1;
            C.GridPartType = gridPartType;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 0;

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();

                BoSSS.Foundation.AMRtests.MeshSymmetryTest(solver.Grid);
            }
        }



        /// <summary>
        /// Testing adaptive mesh refinement in a MPI-parallel environment
        /// </summary>
        [Test]
        public static void AMRtest_3D([Values(GridPartType.METIS, GridPartType.Hilbert, GridPartType.clusterHilbert)] GridPartType gridPartType) {
            var C = DropletOnPlate_AMRtest(D: 3, p: 2, kelem: 8);
            C.NoOfTimesteps = 1;
            C.GridPartType = gridPartType;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 0;

            using(var solver = new XNSE()) {
                solver.Init(C);
                solver.RunSolverMode();

                BoSSS.Foundation.AMRtests.MeshSymmetryTest(solver.Grid);
            }
        }



        /// <summary>
        /// Initiates all the test cases
        /// </summary>
        static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI();
            AMRtest_2D(GridPartType.METIS);
            //ParallelRotatingSphere();
            //ParallelRotatingTilted3DTorus();
            BoSSS.Solution.Application.FinalizeMPI();
        }




    }
}
