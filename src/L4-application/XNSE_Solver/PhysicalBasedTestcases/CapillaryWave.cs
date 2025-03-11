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
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {


    /// <summary>
    /// class providing Controls for the Capillary Wave Testcase
    /// </summary>
    public static class CapillaryWave {


        /// <summary>
        /// Control for various testing
        /// </summary>
        public static XNSE_Control CW(int p = 2, int xkelem = 16, string _DbPath = null) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            // basic database options
            // ======================
            #region db

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";
            //C.Tags.Add("Popinet");
            //C.Tags.Add("Testcase1");

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;

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
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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
            C.FieldOptions.Add("KineticEnergy", new FieldOpts() {
                Degree = 2 * p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

            double rho_l = 1e-2;
            double rho_h = 1e-2;
            double mu_l = 1e-3;
            double mu_h = 1e-3;
            double sigma = 3e-2;

            // Testcase1:
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-4;
            //double mu_h = 1e-4;
            //double sigma = 3e-2;

            // for spatial convergence
            //double rho_l = 1e-2;
            //double rho_h = 1e-2;
            //double mu_l = 1e-3;
            //double mu_h = 1e-3;
            //double sigma = 3e-4;


            C.PhysicalParameters.rho_A = rho_l;
            C.PhysicalParameters.rho_B = rho_h;
            C.PhysicalParameters.mu_A = mu_l;
            C.PhysicalParameters.mu_B = mu_h;
            C.PhysicalParameters.Sigma = sigma;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial disturbance
            // ===================

            double lambda = 1;
            double A0 = lambda / 10;
            Func<double, double> PeriodicFunc = x => A0 * Math.Sin(x * 2 * Math.PI / lambda);


            // grid genration
            // ==============
            #region grid

            double L = lambda;

            //int xkelem = 16;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-1.0 * L / 2.0, 1.0 * L / 2.0, (1 * xkelem) + 0);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (1.0 * L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (1.0 * L / 2.0)) <= 1.0e-8)
                        et = 2;
                    //if (Math.Abs(X[0]) <= 1.0e-8)
                    //    et = 3;
                    //if (Math.Abs(X[0] - L) <= 1.0e-8)
                    //    et = 4;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - PeriodicFunc(X[0])));

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // additional parameters
            // =====================

            //double[] param = new double[3];
            //param[0] = lambda;  // wavelength
            //param[1] = A0;      // initial disturbance
            //param[2] = 0.0;      // y-gravity
            //C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver

            //C.solveKineticEnergyEquation = true;
            //C.ComputeEnergyProperties = true;

            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.JumpGradJump2;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Phasefield;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // specialized Fourier Level-Set
            // =============================


            //int numSp = 640;
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    Timestepper = FourierLevelSet_Timestepper.ExplicitEuler,
            //    InterpolationType = Interpolationtype.CubicSplineInterpolation,
            //};


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            double rho = rho_l;         // Testcase1
            //double dt = Math.Sqrt(rho * Math.Pow((1 / (double)xkelem), 3) / (Math.PI * sigma));             // !!!
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            //double omega0 = Math.Sqrt(sigma * Math.Pow((2 * Math.PI / lambda), 3) / (2.0 * rho));
            //C.NoOfTimesteps = (int)Math.Ceiling((25 / omega0) / dt);                                        // !!!
            C.NoOfTimesteps = 10000;                                     // !!!

            C.saveperiod = 1;

            #endregion


            return C;

        }


        /// <summary>
        /// Control for various paramStudies
        /// </summary>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control[] CW_Paramstudy(string _DbPath) {

            List<XNSE_Control> R = new List<XNSE_Control>();

            // for spatial convergence study
            //int[] h_study = new int[] { 8, 16, 32, 64, 128 };
            //double dt = 1e-5;
            //double t_end = 1e-4;

            // for temporal convergence study
            double[] dt_study = new double[] { 2e-4, 1e-4, 5e-5, 2.5e-5, 1.25e-5 };
            //int h = 64;
            double t_end = dt_study[0];

            // settings for the paramStudy
            //int p = 2;                                                  // !!!!!!!!!!!!
            //string _DbPath = @"D:\local\CWp2_temporalConv_coupledBDF1";
            //int[] saveTs = new int[] { 1, 2, 4, 8, 16 };
            //string[] restartSess = new string[] { "107a29ae-0b83-4cc0-ac52-8734079d7015",
            //    "ff03a824-2485-4b56-8cc3-a8db763ef32d",
            //    "619c7689-7863-4ca9-81fb-f9e52e8e7fdd",
            //    "b400707f-fad8-440b-950c-6b79502b7636",
            //    "76843518-c343-4212-870a-7b072c655d52" };
            //int i = 0;

            //foreach (int h_k in h_study) {
            //    var C = CapillaryWave_Popinet(p, h_k, dt, t_end, _DbPath);
            //    C.Paramstudy_CaseIdentification = new Tuple<string, object>[] { new Tuple<string, object>("xkelem", h_k),
            foreach (double dt_h in dt_study) {
                //var C = CapillaryWave_Popinet(p, h, dt_h, t_end, _DbPath);
                //C.Paramstudy_CaseIdentification = new Tuple<string, object>[] { new Tuple<string, object>("dt", dt_h),
                //foreach (int save in saveTs) {
                //    var C = CapillaryWave_Popinet(p, h, dt_study[4], t_end, _DbPath, save);
                //    C.Paramstudy_CaseIdentification = new Tuple<string, object>[] { new Tuple<string, object>("saveDb", save),
                //};
                //R.Add(C);
            }

            return R.ToArray();
        }


        /// <summary>
        /// Control for the testcase according to Popinet
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xkelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control CW_Popinet(int p = 2, int xkelem = 16, string _DbPath = null) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";
            _DbPath = @"D:\rieckmann\BoSSS_DB";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";
            C.Tags.Add("Popinet");
            C.Tags.Add("Testcase1");

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;
            C.PostprocessingModules.Add(new WaveLikeLogging());

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
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics  


            // Testcase1:
            double rho_l = 1e-3;
            double rho_h = 1e-3;
            double mu_l = 1e-4;
            double mu_h = 1e-4;
            double sigma = 3e-2;

            // for spatial convergence
            //double rho_l = 1e-2;
            //double rho_h = 1e-2;
            //double mu_l = 1e-3;
            //double mu_h = 1e-3;
            //double sigma = 3e-4;


            // unstable configuration
            C.PhysicalParameters.rho_A = rho_l;
            C.PhysicalParameters.rho_B = rho_h;
            C.PhysicalParameters.mu_A = mu_l;
            C.PhysicalParameters.mu_B = mu_h;
            C.PhysicalParameters.Sigma = sigma;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial disturbance
            // ===================

            double lambda = 1;
            double A0 = lambda / 100;
            Func<double, double> PeriodicFunc = x => A0 * Math.Sin(x * 2 * Math.PI / lambda);


            // grid genration
            // ==============
            #region grid

            double L = lambda;

            //int xkelem = 16;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-3.0 * L / 2.0, 3.0 * L / 2.0, (3 * xkelem) + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (3.0 * L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (3.0 * L / 2.0)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", X => 0.0);

            #endregion


            // Initial Values
            // ==============
            #region init
            double cahn = (1.0 / 4.164) * L/(double)xkelem * 8.0 / 4.0;
            C.InitialValues_Evaluators.Add("Phi", (X => Math.Tanh((X[1] - PeriodicFunc(X[0]))/(Math.Sqrt(2) * cahn))));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // additional parameters
            // =====================

            double[] param = new double[4];
            param[0] = 1.0;
            param[1] = lambda;  // wavelength
            param[2] = A0;      // initial disturbance
            param[3] = 0.0;      // y-gravity
            C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.Phasefield;
            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint;

            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;

            C.SuperSampling = 2;
            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // specialized Fourier Level-Set
            // =============================


            //int numSp = 640;
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    Timestepper = FourierLevelSet_Timestepper.ExplicitEuler,
            //    InterpolationType = Interpolationtype.CubicSplineInterpolation,
            //};

            C.FourierLevSetControl = null;

            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            double rho = rho_l;         // Testcase1
            double dt = Math.Sqrt(rho * Math.Pow((1 / (double)xkelem), 3) / (Math.PI * sigma));             // !!!
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            double omega0 = Math.Sqrt(sigma * Math.Pow((2 * Math.PI / lambda), 3) / (2.0 * rho));
            C.NoOfTimesteps = (int)Math.Ceiling((25 / omega0) / dt);                                        // !!!
            //C.NoOfTimesteps = (int)Math.Ceiling(t_end / dt);                                     // !!!

            C.saveperiod = 1;

            #endregion


            return C;

        }

        /// <summary>
        /// Some Settings for the CW Test to read from Worksheet
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CW_ForWorksheet()
        {
            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = true;
            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;
            C.PostprocessingModules.Add(new WaveLikeLogging());

            #endregion

            // DG degrees
            // ==========
            #region degrees

            // need to be set by user via setDGdegree() in worksheet

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

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
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
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
        /// <param name="xkelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control CW_Test(int p = 3, int xkelem = 16, int method = 0) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            string _DbPath = null;
            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";
            //string _DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSE_studyDB";
            //string _DbPath = @"\\terminal03\Users\smuda\local\terminal03_XNSE_studyDB";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "CapillaryWave";
            C.SessionName = "CapillaryWave_Setup0_convStudy_k3_mesh1_AMR1";
            //C.SessionName = "CapillaryWave_Setup0_convStudy_k2_mesh3_restart"; //"6f816774-8c9c-44f6-afe3-98f77d1764f6"
            //Guid restart = new Guid("6f816774-8c9c-44f6-afe3-98f77d1764f6");
            //C.SessionName = "CapillaryWave_Setup0_convStudy_k3_mesh2_restart"; //"7f9130d3-eaab-4ac2-9844-fd91be6f1edf"
            //Guid restart = new Guid("7f9130d3-eaab-4ac2-9844-fd91be6f1edf");      

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;
            //C.LogPeriod = 10;
            C.PostprocessingModules.Add(new WaveLikeLogging() { LogPeriod = 10 });

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
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            //C.FieldOptions.Add("GravityY", new FieldOpts() {
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 2*p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics  


            // Testcase1:
            double rho_l = 1e-3;
            double rho_h = 1e-3;
            double mu_l = 1e-5;
            double mu_h = 1e-5;
            double sigma = 3e-2;

            // for spatial convergence
            //double rho_l = 1e-2;
            //double rho_h = 1e-2;
            //double mu_l = 1e-3;
            //double mu_h = 1e-3;
            //double sigma = 3e-4;


            // unstable configuration
            C.PhysicalParameters.rho_A = rho_l;
            C.PhysicalParameters.rho_B = rho_h;
            C.PhysicalParameters.mu_A = mu_l;
            C.PhysicalParameters.mu_B = mu_h;
            C.PhysicalParameters.Sigma = sigma;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial disturbance
            // ===================

            double lambda = 1;
            double A0 = lambda / 100;
            Func<double, double> PeriodicFunc = x => A0 * Math.Sin(x * 2 * Math.PI / lambda);


            // grid genration
            // ==============
            #region grid

            double L = lambda;

            //int xkelem = 16;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-3.0 * L / 2.0, 3.0 * L / 2.0, (3 * xkelem) + 0);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: !(method == 2 || method == 3));

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                if (method == 2 || method == 3) {
                    grd.EdgeTagNames.Add(3, "freeslip");
                }


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (3.0 * L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (3.0 * L / 2.0)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 3;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => 0.0);
            C.AddBoundaryValue("wall_upper", "VelocityX#B", X => 0.0);
            if (method == 2 || method == 3) {
                C.AddBoundaryValue("freeslip");
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - PeriodicFunc(X[0])));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            ////Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restart, null);

            #endregion

            // additional parameters
            // =====================

            double[] param = new double[4];
            param[0] = 1;        // wavenumber;
            param[1] = L;        // wavelength
            param[2] = A0;       // initial disturbance
            param[3] = 0.0;      // y-gravity
            C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver


            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.BaseRefinementLevel = 1;
            C.InitSignedDistance = false;


            #endregion


            // specialized Fourier Level-Set
            // =============================


            int numSp = 640;
            var FourierContrl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
                FourierEvolve = Fourier_Evolution.MaterialPoints,
            };

            C.SetLevelSetMethod(method, FourierContrl);
            //C.SessionName = "RisingBubble_methodStudy_k2_" + C.methodTagLS;


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            //double rho = rho_l;         // Testcase1
            //double dt = Math.Sqrt(rho * Math.Pow((1 / (double)xkelem), 3) / (Math.PI * sigma));             // !!!
            double dt = 2e-5;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 0.4;
            //double omega0 = Math.Sqrt(sigma * Math.Pow((2 * Math.PI / lambda), 3) / (2.0 * rho));
            //C.NoOfTimesteps = 10; // (int)Math.Ceiling((25 / omega0) / dt);                                        // !!!
            C.NoOfTimesteps = (int)Math.Ceiling(0.4 / dt);                                     // !!!

            C.saveperiod = 10;

            #endregion


            return C;

        }


        public static XNSE_Control CW_BrokenLevelSet(int p = 2, int xkelem = 8, string _DbPath = null) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            // basic database options
            // ======================
            #region db

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";
            //C.Tags.Add("Popinet");
            //C.Tags.Add("Testcase1");

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;

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
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityY", new FieldOpts() {
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
            C.FieldOptions.Add("KineticEnergy", new FieldOpts() {
                Degree = 2 * p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

            double rho = 10;
            double mu = 0.01;
            double sigma = 0.01;

            C.PhysicalParameters.rho_A = rho;
            C.PhysicalParameters.rho_B = rho;
            C.PhysicalParameters.mu_A = mu;
            C.PhysicalParameters.mu_B = mu;
            C.PhysicalParameters.Sigma = sigma;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ==============
            #region grid

            double L = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-2.0 * L / 2.0, 2.0 * L / 2.0, (2 * xkelem));
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (2.0 * L / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (2.0 * L / 2.0)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            #endregion


            // Initial Values
            // ==============
            #region init

            double A0 = L / 40;
            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - (A0 * Math.Abs(Math.Sin(X[0] * 2.0 * Math.PI / L) - (Math.Sqrt(2.0) / 2.0)) + 0.002 * Math.Sign(X[0] - L / 2))));
            //C.InitialValues_Evaluators.Add("Phi", (X => X[1] - (A0 * Math.Sin(X[0] * 2.0 * Math.PI / L) + 0.002 * Math.Sign(X[0] - L / 2))));
            //C.InitialValues_Evaluators.Add("Phi", (X => X[1] - (A0 * Math.Abs(Math.Sin(X[0] * 2.0 * Math.PI / L) - (Math.Sqrt(2.0) / 2.0)))));

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // additional parameters
            // =====================

            //double[] param = new double[3];
            //param[0] = lambda;  // wavelength
            //param[1] = A0;      // initial disturbance
            //param[2] = 0.0;      // y-gravity
            //C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver

            //C.solveKineticEnergyEquation = true;
            //C.ComputeEnergyProperties = true;
            //C.CheckInterfaceProps = true;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.surfaceDivergence;

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // specialized Fourier Level-Set
            // =============================


            int numSp = 640;
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    Timestepper = FourierLevelSet_Timestepper.ExplicitEuler,
            //    InterpolationType = Interpolationtype.CubicSplineInterpolation,
            //};


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            double dt = Math.Sqrt(rho * Math.Pow((1 / (double)xkelem), 3) / (Math.PI * sigma))/8.0;             // !!!
            //double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            //double omega0 = Math.Sqrt(sigma * Math.Pow((2 * Math.PI / lambda), 3) / (2.0 * rho));
            //C.NoOfTimesteps = (int)Math.Ceiling((25 / omega0) / dt);                                        // !!!
            C.NoOfTimesteps = 10000;                                     // !!!

            C.saveperiod = 1;

            #endregion


            return C;

        }


        public static XNSE_Control CW_OneCell(int p = 2) {


            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            // basic database options
            // ======================
            #region db

            C.DbPath = null; 
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";

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
                Degree = 2*p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            C.RegisterUtilitiesToIOFields = true;

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

            double rho = 10;
            double mu = 0.01;
            double sigma = 0.01;

            C.PhysicalParameters.rho_A = rho;
            C.PhysicalParameters.rho_B = rho;
            C.PhysicalParameters.mu_A = mu;
            C.PhysicalParameters.mu_B = mu;
            C.PhysicalParameters.Sigma = sigma;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ==============
            #region grid

            double L = 1.0;
            int cellFactor = 3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L * cellFactor, cellFactor + 1);
                double[] Ynodes = GenericBlas.Linspace(-(3.0 / 2.0) * L, (3.0 / 2.0) * L, 4);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall");
                grd.EdgeTagNames.Add(2, "pressure_outlet");
                //grd.EdgeTagNames.Add(3, "navierslip_linear");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (3.0 / 2.0) * L) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (3.0 / 2.0) * L) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] - L * cellFactor) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall");
            C.AddBoundaryValue("pressure_outlet");
            //C.AddBoundaryValue("navierslip_linear");

            #endregion


            // Initial Values
            // ==============
            #region init

            double A0 = L / 10;
            C.InitialValues_Evaluators.Add("Phi", X => X[1] - A0 * Math.Sin(X[0] * 2.0 * Math.PI / (L * cellFactor)) );

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.CheckJumpConditions = true;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            //C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.surfaceDivergence;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            double dt = Math.Sqrt(rho / (cellFactor * (Math.PI * sigma))) / 8.0;           
            //double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 20000;
            C.NoOfTimesteps = 10000;

            C.saveperiod = 1;

            #endregion


            return C;
        }


        /// <summary>
        /// control object for the use in BoSSSPad Worksheet mode
        /// additional parameters need to be set during job creation
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CW_forWorksheet() {


            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = true;
            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;
            C.PostprocessingModules.Add(new WaveLikeLogging());

            #endregion

            // DG degrees
            // ==========
            #region degrees

            // need to be set by user via setDGdegree() in worksheet

            #endregion


            // Physical Parameters
            // ===================
            #region physics  

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

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            #endregion
            
            
            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;


            //C.dtMax = dt; // need to be set according to grid and DG degree
            //C.dtMin = dt;
            C.Endtime = 1000;
            //C.NoOfTimesteps = 0; 

            C.saveperiod = 1;
            
            #endregion


            return C;

        }


    }


}
