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
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the Rising Bubble testcase
    /// </summary>
    public static class RisingBubble {


        /// <summary>
        /// Control for various testing/Florian; **If you want to modify, create your own copy!**
        /// - jan2021: used as a test for the implementation of Stokes-Extension-Velocity
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <param name="dt">
        /// size of time-step, should be chosen as <c>(1.0 / (double)kelem) / 16.0</c> 
        /// </param>
        /// <returns></returns>
        public static XNSE_Control RB_fk(int p = 2, int kelem = 20, double dt = 3.125e-3, string _DbPath = null) {
            //BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases.RisingBubble.RB(p:1,kelem:10)

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            //_DbPath = @"\\dc1\userspace\yotov\bosss-db\RisingBubble"";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.SessionName = "RisingBubbleHeimann";
            C.ProjectDescription = "rising bubble";

            //C.LogValues = XNSE_Control.LoggingValues.RisingBubble;

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.SetDGdegree(p);

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Testcase 1");
            //C.PhysicalParameters.rho_A = 100;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 24.5;

            C.Tags.Add("Testcase 2");
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 10;
            C.PhysicalParameters.Sigma = 1.96;


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

                grd.DefineEdgeTags(delegate (double[] X) {
                    string et = null;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = "wall_lower";
                    else if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = "wall_upper";
                    else if(Math.Abs(X[0]) <= 1.0e-8)
                        et = "freeslip_left";
                    else if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = "freeslip_right";
                    else
                        throw new ArgumentException("unknown boundary region");

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
            
            //C.GridPartType = GridPartType.Predefined;
            //C.GridPartOptions = "VierProcSplit";


            #endregion


            // Initial Values
            // ==============
            #region init

            double[] center = new double[] { 0.5, 0.5 }; //0.5,0.5
            double radius = 0.25;

            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius); // signed-distance form

            C.InitialValues_Evaluators.Add(VariableNames.LevelSetCG, PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            
            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            //Guid restartID = new Guid("322f07e1-2ac3-4ed4-af8b-1c46ab7e55a0");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 986);
            //C.ReInitControl.PrintIterations = true; 

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            //C.AddBoundaryCondition("wall_lower", VariableNames.LevelSet, PhiFunc);

            #endregion


            // misc. solver options
            // ====================
            #region solver



            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = true;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = dt/2.0;
            C.dtMin = dt/2.0;
            C.Endtime = 3;
            C.NoOfTimesteps = 20000;
            C.saveperiod = 1;

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
        public static XNSE_Control RB_BenchmarkTest(int p = 2, int kelem = 10, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();


            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\rieckmann\BoSSS_DB";

            // basic database options
            // ======================
            #region db

            C.DbPath = null;//_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Bubble";
            C.ProjectDescription = "rising bubble";
            C.Tags.Add("benchmark setup");

            //C.LogValues = XNSE_Control.LoggingValues.RisingBubble;
            //C.LogPeriod = 1;
            C.PostprocessingModules.Add(new RisingBubble2DBenchmarkQuantities());

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(p);


            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Testcase 1");
            //C.PhysicalParameters.rho_A = 100;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 10;
            //double sigma = 24.5;
            //C.PhysicalParameters.Sigma = sigma;

            C.Tags.Add("Testcase 2");
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 10;
            double sigma = 1.96;
            C.PhysicalParameters.Sigma = sigma;


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

                grd.AddPredefinedPartitioning("VierProcSplit", delegate (double[] X)
                {
                    int rank;
                    double x = X[0];
                    if (x < 0.3)
                        rank = 0;
                    else if (x < 0.5)
                        rank = 1;
                    else if (x < 0.7)
                        rank = 2;
                    else
                        rank = 3;

                    return rank;
                });


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

            double cahn = (1.0 / 4.164) * h * 8.0 / (p + 1);
            //Func<double[], double> PhiFunc = (X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2()); // quadratic form
            Func<double[], double> PhiFunc = (X => Math.Tanh((((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius)/(Math.Sqrt(2)*cahn))); // signed-distance form

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double> PeriodicFunc = x => radius;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            C.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("b90c5f79-9b82-47cd-b400-e9abbbd83e19");  //new Guid("2953cd96-ea27-4989-abd3-07e99d35de5f"); 
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 1140);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            #endregion

            // Level-Set
            // =================
            #region Fourier

            //int numSp = 1024;
            //double[] FourierP = new double[numSp];
            //double[] samplP = new double[numSp];
            //for (int sp = 0; sp < numSp; sp++) {
            //    FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
            //    samplP[sp] = radius;
            //}

            //double circum = 2.0 * Math.PI * radius;
            //double filter = (circum * 20.0) / ((double)numSp / 2.0);
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
            //    center = center,
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    centerMove = CenterMovement.Reconstructed,
            //};


            #endregion


            C.Option_LevelSetEvolution = LevelSetEvolution.Phasefield;

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;


            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;


            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.PhysicalParameters.mu_I = 1 * sigma;
            C.PhysicalParameters.lambda_I = 2 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            //C.LS_TrackerWidth = 2;
            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 2;
            C.BaseRefinementLevel = 2;
            C.LS_TrackerWidth = 2;            

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            double dt = 0.1 * 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = (int)Math.Ceiling(3.0/dt);
            C.Endtime = 300;
            C.SuperSampling = 1;
            C.saveperiod = 1;
            C.rollingSaves = true;


            #endregion

            return C;

        }

        /// <summary>
        /// predefined Control settings for the Rising Bubble Benchmark
        /// </summary>
        /// <param name="setup"></param>
        /// <returns></returns>
        public static XNSE_Control RB_ForWorksheet(int setup) {
            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = true;
            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.RisingBubble;
            C.PostprocessingModules.Add(new RisingBubble2DBenchmarkQuantities());
            #endregion

            // DG degrees
            // ==========
            #region degrees

            int p = 2;

            // need to be set by user via setDGdegree() in worksheet
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
                Degree = p + 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p + 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            #endregion


            // Physical Parameters
            // ===================
            #region physics  


            // Physical Parameters
            // ===================

            switch(setup) {
                case 1:
                C.Tags.Add("Testcase 1");
                C.PhysicalParameters.rho_A = 100;
                C.PhysicalParameters.rho_B = 1000;
                C.PhysicalParameters.mu_A = 1;
                C.PhysicalParameters.mu_B = 10;
                C.PhysicalParameters.Sigma = 24.5;
                break;
                case 2:
                C.Tags.Add("Testcase 2");
                C.PhysicalParameters.rho_A = 1;
                C.PhysicalParameters.rho_B = 1000;
                C.PhysicalParameters.mu_A = 0.1;
                C.PhysicalParameters.mu_B = 10;
                C.PhysicalParameters.Sigma = 1.96;
                break;
                default:
                throw new NotImplementedException();
            }

            // need to be set during job creation 

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid genration
            // ==============
            #region grid

            // need to be set by user via setGrid() in worksheet

            #endregion


            // boundary conditions
            // ===================
            #region BC

            // need to be set during job creation 

            #endregion


            // Initial Values
            // ==============
            #region init

            // need to be set during job creation 

            #endregion


            // additional parameters (evaluation)
            // ==================================

            // need to be set during job creation 


            // misc. solver options
            // ====================
            #region solver

            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion


            // Level-Set options (AMR)
            // =======================
            #region levset

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.Option_LevelSetEvolution = LevelSetEvolution.Phasefield;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            #endregion

            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;


            //C.dtMax = dt; // need to be set according to grid and DG degree
            //C.dtMin = dt;
            C.Endtime = 1000;
            //C.NoOfTimesteps = 0; 

            C.saveperiod = 1;

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
        public static XNSE_Control RB_Test(int p = 2, int kelem = 20, int method = 2) {

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            //string _DbPath = @"D:\local\local_Testcase_databases\Testcase_RisingBubble";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //string _DbPath = null; // @"D:\local\local_XNSE_StudyDB";
            string _DbPath = null; // "\\HPCCLUSTER\hpccluster-scratch\smuda\XNSE_studyDB";
            //string _DbPath = @"\\terminal03\Users\smuda\local\terminal03_XNSE_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; //C.DbPath != null;
            C.ProjectName = "RisingBubble";

            //C.LogValues = XNSE_Control.LoggingValues.RisingBubble;
            //C.LogPeriod = 3;
            C.PostprocessingModules.Add(new RisingBubble2DBenchmarkQuantities() { LogPeriod = 3 });

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
            //C.FieldOptions.Add("Curvature", new FieldOpts() {
            //    Degree = p,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            //C.FieldOptions.Add("DivergenceVelocity", new FieldOpts() {
            //    Degree = p,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Testcase 1");
            //C.PhysicalParameters.rho_A = 100;
            //C.PhysicalParameters.rho_B = 1000;
            //C.PhysicalParameters.mu_A = 1;
            //C.PhysicalParameters.mu_B = 10;
            //C.PhysicalParameters.Sigma = 24.5;


            //C.Tags.Add("Testcase 1 - higher parameters");
            //C.PhysicalParameters.rho_A = 1000;
            //C.PhysicalParameters.rho_B = 10000;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 100;
            //C.PhysicalParameters.Sigma = 245;

            C.Tags.Add("Testcase 2");
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 10;
            C.PhysicalParameters.Sigma = 1.96;

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


                //grd.EdgeTagNames.Add(1, "wall_lower");
                //grd.EdgeTagNames.Add(2, "wall_upper");
                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    string et = null;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = "wall_lower";
                    else if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = "wall_upper";
                    else if(Math.Abs(X[0]) <= 1.0e-8)
                        et = "freeslip_left";
                    else if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = "freeslip_right";
                    else
                        throw new ArgumentException("unknown boundary region");

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


            ////var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("28adc7a7-a336-435b-b560-d3c02f46c43d");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

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

            double circum = 2.0 * Math.PI * radius;
            double filter = (circum * 20.0) / ((double)numSp / 2.0);
            var Fouriercontrl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem) {
                //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2.0*Math.PI, PeriodicFunc, radius, 1.0/(double)kelem) { 
                center = center,
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                centerMove = CenterMovement.Reconstructed,
                //curvComp_extended = false
            };


            C.SetLevelSetMethod(method, Fouriercontrl);
            //C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.None;
            //C.SessionName = "RisingBubble_methodStudy_k2_" + C.methodTagLS + "_withVolumeCorrection";
            //C.EnforceLevelSetConservation = true;

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.BaseRefinementLevel = 1;
            C.AMR_startUpSweeps = 1;

            //C.ReInitOnRestart = true;
            //C.ReInitPeriod = 100;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            double dt = (1.0 / (double)kelem) / 16.0;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 3;
            C.NoOfTimesteps = (int)(3 / dt);
            C.saveperiod = 1;

            #endregion

            return C;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="setup"> physical setup </param>
        /// <param name="method"> method setup regarding the level set handling </param>
        /// <returns></returns>
        public static XNSE_Control RB_forWorksheet(int setup, bool restart = false) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = true;
            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.RisingBubble;
            C.PostprocessingModules.Add(new RisingBubble2DBenchmarkQuantities());
            #endregion


            // DG degrees
            // ==========
            #region degrees

            // need to be set by user via setDGdegree() in worksheet

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

            double rho_h = 0.0;                // heavy fluid
            double mu_h = 0.0;
            double rho_l = 0.0;                // light fluid
            double mu_l = 0.0;
            double sigma = 0.0;

            switch (setup) {
                case 0: {
                        // testcase 1:
                        rho_h = 1000;
                        mu_h = 10;
                        rho_l = 100;
                        mu_l = 1;
                        sigma = 24.5;
                        break;
                    }
                case 1: {
                        // testcase 2:
                        rho_h = 1000;
                        mu_h = 10;
                        rho_l = 1;
                        mu_l = 0.1;
                        sigma = 1.96;
                        break;
                    }
            }

            C.PhysicalParameters.rho_A = rho_l;
            C.PhysicalParameters.rho_B = rho_h;
            C.PhysicalParameters.mu_A = mu_l;
            C.PhysicalParameters.mu_B = mu_h;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid genration
            // ==============
            #region grid

            // need to be set by user via setGrid() in worksheet

            #endregion


            // boundary conditions
            // ===================
            #region BC

            // need to be set during job creation 

            #endregion


            // Initial Values
            // ==============
            #region init

            if (!restart) {
                C.AddInitialValue("Phi", "X => ((X[0] - 0.5).Pow2() + (X[1] - 0.5).Pow2()).Sqrt() - 0.25", false);

                C.AddInitialValue("GravityY#A", "X => -9.81e-1", false);
                C.AddInitialValue("GravityY#B", "X => -9.81e-1", false);
            }

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion


            // Level-Set options (AMR)
            // =======================
            #region levset

            // need to be set by user via setGrid() in worksheet

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            //C.dtMax = dt; // need to be set according to grid and DG degree
            //C.dtMin = dt;
            C.Endtime = 3;
            //C.NoOfTimesteps = 0; 

            C.saveperiod = 1;
            
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
            C.SetDGdegree(p);

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

            //var database = DatabaseInfo.Open(_DbPath);
            Guid restartID = new Guid("9d1bbcd2-38d0-43d3-90d4-7ac7f535079c");
            C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);


            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

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
