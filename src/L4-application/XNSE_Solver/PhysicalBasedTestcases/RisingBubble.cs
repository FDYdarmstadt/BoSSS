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
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution.Control;
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the Rising Bubble testcase
    /// </summary>
    public static class RisingBubble {


        /// <summary>
        /// Control for various testing
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <param name="dt">
        /// Timestepsize, should be chosen as (1.0 / (double)kelem) / 16.0;
        /// </param>
        /// <returns></returns>
        public static XNSE_Control RB(int p = 2, int kelem = 20, double dt = 3.125e-3, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.SessionName = "RisingBubbleHeimann";
            C.ProjectDescription = "rising bubble";

            C.LogValues = XNSE_Control.LoggingValues.RisingBubble;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Testcase 1");
            C.PhysicalParameters.rho_A = 100;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 10;
            C.PhysicalParameters.Sigma = 24.5;


            //C.Tags.Add("Testcase 1 - higher parameters");
            //C.PhysicalParameters.rho_A = 1000;
            //C.PhysicalParameters.rho_B = 10000;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            //C.Tags.Add("Testcase 2");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 1.96;

            // Re = 3.5 ; Bo(Eo) = 1
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            //// Re = 35 ; Bo(Eo) = 100
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 2.45;

            //// Re = 70 ; Bo(Eo) = 10
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.05;
            //C.PhysicalParameters.mu_B = 5;
            //C.PhysicalParameters.Sigma = 24.5;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 2.0;

            //int kelem = 160;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                grd.AddPredefinedPartitioning("ZwoProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if (x < 0.5)
                        rank = 0;
                    else
                        rank = 1;

                    return rank;
                });

                grd.AddPredefinedPartitioning("VierProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if (x < 0.35)
                        rank = 0;
                    else if (x < 0.5)
                        rank = 1;
                    else if (x < 0.75)
                        rank = 2;
                    else
                        rank = 3;

                    return rank;
                });


                return grd;
            };

            C.GridPartType = GridPartType.Predefined;
            C.GridPartOptions = "VierProcSplit";


            #endregion


            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.5, 0.5 }; //0.5,0.5
            double radius = 0.25;

            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("58745416-3320-4e0c-a5fa-fc3a2c5203c7");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");

            //C.AddBoundaryCondition("wall_lower", VariableNames.LevelSet, PhiFunc);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            //C.EnforceLevelSetConservation = true;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_MinIterations = 1;
            C.Solver_ConvergenceCriterion = 1e-7;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;


            C.Option_LevelSetEvolution = LevelSetEvolution.ExtensionVelocity;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ContinuousDG;
            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 3;
            C.NoOfTimesteps = 10000; // (int)(3 / dt);
            C.saveperiod = 10;

            #endregion

            return C;

        }


        /// <summary>
        /// specified control for benchmark testing
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control RB_BenchmarkTest(int p = 2, int kelem = 20, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();


            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Bubble";
            C.ProjectDescription = "rising bubble";
            C.Tags.Add("benchmark setup");

            C.LogValues = XNSE_Control.LoggingValues.RisingBubble;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Testcase 1");
            C.PhysicalParameters.rho_A = 100;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 10;
            double sigma = 24.5;
            C.PhysicalParameters.Sigma = sigma;

            //C.Tags.Add("Testcase 2");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 1.96;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double h = 1.0 / (double)kelem;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 1.0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 2.0, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - 2.0) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - 1.0) <= 1.0e-8)
                        et = 4;
                    return et;
                });

                //grd.AddPredefinedPartitioning("ZwoProcSplit", delegate (double[] X) {
                //    int rank;
                //    double x = X[0];
                //    if (x < 0.5)
                //        rank = 0;
                //    else
                //        rank = 1;

                //    return rank;
                //});

                //grd.AddPredefinedPartitioning("VierProcSplit", delegate (double[] X) {
                //    int rank;
                //    double x = X[0];
                //    if (x < 0.35)
                //        rank = 0;
                //    else if (x < 0.5)
                //        rank = 1;
                //    else if (x < 0.75)
                //        rank = 2;
                //    else
                //        rank = 3;

                //    return rank;
                //});


                return grd;
            };

            //C.GridPartType = GridPartType.Predefined;
            //C.GridPartOptions = "VierProcSplit";

            #endregion


            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.5, 0.5 };
            double radius = 0.25;

            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("e356230c-57dd-4228-9986-0ec0d1335818");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 123);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");

            #endregion

            // Level-Set
            // =================
            #region Fourier

            int numSp = 1024;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                samplP[sp] = radius;
            }

            //double circum = 2.0 * Math.PI * radius;
            //double filter = (circum * 20.0) / ((double)numSp / 2.0);
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
            //    center = center,
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    centerMove = CenterMovement.Reconstructed,
            //};


            #endregion


            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            C.LSContiProjectionMethod = ContinuityProjectionOption.ContinuousDG;
            //C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpoint+levelset" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;


            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.FullBoussinesqScriven;
            C.PhysicalParameters.mu_I = 1 * sigma;
            C.PhysicalParameters.lambda_I = 2 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            //C.LS_TrackerWidth = 2;
            //C.AdaptiveMeshRefinement = true;
            //C.RefinementLevel = 2;


            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;

            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 3000;
            C.saveperiod = 10;


            #endregion

            return C;

        }


        /// <summary>
        /// Control for TestProgramm (not to be changed!!!)
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control RB_Test(int p = 2, int kelem = 20, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Bubble";
            C.ProjectDescription = "rising bubble";

            C.LogValues = XNSE_Control.LoggingValues.RisingBubble;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Testcase 1");
            C.PhysicalParameters.rho_A = 100;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 10;
            C.PhysicalParameters.Sigma = 24.5;


            //C.Tags.Add("Testcase 1 - higher parameters");
            //C.PhysicalParameters.rho_A = 1000;
            //C.PhysicalParameters.rho_B = 10000;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            //C.Tags.Add("Testcase 2");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 1.96;

            // Re = 3.5 ; Bo(Eo) = 1
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            //// Re = 35 ; Bo(Eo) = 100
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 2.45;

            //// Re = 70 ; Bo(Eo) = 10
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 0.05;
            //C.PhysicalParameters.mu_B = 5;
            //C.PhysicalParameters.Sigma = 24.5;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // grid generation
            // ===============
            #region grid


            double xSize = 1.0;
            double ySize = 2.0;

            //int kelem = 160;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);


                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                grd.AddPredefinedPartitioning("ZwoProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if (x < 0.5)
                        rank = 0;
                    else
                        rank = 1;

                    return rank;
                });

                grd.AddPredefinedPartitioning("VierProcSplit", delegate (double[] X) {
                    int rank;
                    double x = X[0];
                    if (x < 0.35)
                        rank = 0;
                    else if (x < 0.5)
                        rank = 1;
                    else if (x < 0.75)
                        rank = 2;
                    else
                        rank = 3;

                    return rank;
                });


                return grd;
            };

            C.GridPartType = GridPartType.Predefined;
            C.GridPartOptions = "VierProcSplit";


            #endregion


            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.5, 0.5 }; //0.5,0.5
            double radius = 0.25;

            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("58745416-3320-4e0c-a5fa-fc3a2c5203c7");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");

            //C.AddBoundaryCondition("wall_lower", VariableNames.LevelSet, PhiFunc);

            #endregion

            // Level-Set
            // =================
            #region Fourier

            int numSp = 640;
            double[] FourierP = new double[numSp];
            double[] samplP = new double[numSp];
            for (int sp = 0; sp < numSp; sp++) {
                FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                samplP[sp] = radius;
            }

            //double circum = 2.0 * Math.PI * radius;
            //double filter = (circum * 20.0) / ((double)numSp / 2.0);
            C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
                //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2.0*Math.PI, PeriodicFunc, radius, 1.0/(double)kelem) { 
                center = center,
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                centerMove = CenterMovement.Reconstructed,
                //curvComp_extended = false
            };

            C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
            C.FourierLevSetControl.Timestepper = FourierLevelSet_Timestepper.RungeKutta1901;

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 4;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;


            //C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.CompMode = AppControl._CompMode.Transient;
            //C.TimeStepper = XNSE_Control._Timestepper.BDF2;
            double dt = 3e-3; // (1.0 / (double)kelem) / 16.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 3; // (int)(3 / dt);
            C.saveperiod = 10;

            #endregion

            return C;

        }



        /// <summary>
        /// Control for two Rising Bubble 
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control BubbleMerger(int p = 2, int kelem = 40, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            _DbPath = @"D:\local\local_test_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Bubble";
            C.ProjectDescription = "bubble merger";
            C.Tags.Add("smolianski");

            #endregion

            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            #endregion

            // grid generation
            // ===============
            #region grid

            double h = 1.0 / (double)kelem;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 1.0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 2.0, 2 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - 2.0) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - 1.0) <= 1.0e-8)
                        et = 4;
                    return et;
                });

                return grd;
            };

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            // Bo = 250, Re = 35
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 100;
            C.PhysicalParameters.mu_A = 0.01;
            C.PhysicalParameters.mu_B = 0.1;
            C.PhysicalParameters.Sigma = 0.097;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion

            // Initial Values
            // ==============
            #region init

            // large bubble above
            double[] center_l = new double[] { 0.5, 1.0 };
            double radius_l = 0.25;
            Func<double[], double> bubble_l = (X => ((X[0] - center_l[0]).Pow2() + (X[1] - center_l[1]).Pow2()).Sqrt() - radius_l); // signed-distance form

            // small bubble under
            double[] center_s = new double[] { 0.5, 0.5 };
            double radius_s = 0.2;
            Func<double[], double> bubble_s = (X => ((X[0] - center_s[0]).Pow2() + (X[1] - center_s[1]).Pow2()).Sqrt() - radius_s); // signed-distance form

            Func<double[], double> PhiFunc = (X => Math.Min(bubble_l(X), bubble_s(X)));

            //C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => -0.981);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => -0.981);

            var database = new DatabaseInfo(_DbPath);
            Guid restartID = new Guid("9d1bbcd2-38d0-43d3-90d4-7ac7f535079c");
            C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);


            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.MUMPS };

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;
            //C.option_solver = C.PhysicalParameters.IncludeConvection ? "fixpoint+levelset" : "direct";
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.CompMode = AppControl._CompMode.Transient;

            double dt = 2e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 1500;
            C.saveperiod = 10;

            #endregion


            return C;
        }

    }


}
