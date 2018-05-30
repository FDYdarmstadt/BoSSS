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
using BoSSS.Solution.Timestepping;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {


    /// <summary>
    /// class providing Controls for the Capillary Wave Testcase
    /// </summary>
    public static class CapillaryWave {


        /// <summary>
        /// Control for various testing
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xkelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control CW(int p = 2, int xkelem = 32, string _DbPath = null) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            _DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";
            C.Tags.Add("Popinet");
            C.Tags.Add("Testcase1");

            C.LogValues = XNSE_Control.LoggingValues.Wavelike;

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

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", X => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", X => 0.0);

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - PeriodicFunc(X[0])));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // additional parameters
            // =====================

            double[] param = new double[3];
            param[0] = lambda;  // wavelength
            param[1] = A0;      // initial disturbance
            param[2] = 0.0;      // y-gravity
            C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // specialized Fourier Level-Set
            // =============================


            int numSp = 640;
            C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                Timestepper = FourierLevelSet_Timestepper.ExplicitEuler,
                InterpolationType = Interpolationtype.CubicSplineInterpolation,
            };


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

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
        public static XNSE_Control CW_Popinet(int p = 2, int xkelem = 32, string _DbPath = null) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            _DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";
            C.Tags.Add("Popinet");
            C.Tags.Add("Testcase1");

            C.LogValues = XNSE_Control.LoggingValues.Wavelike;

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

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", X => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", X => 0.0);

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - PeriodicFunc(X[0])));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // additional parameters
            // =====================

            double[] param = new double[3];
            param[0] = lambda;  // wavelength
            param[1] = A0;      // initial disturbance
            param[2] = 0.0;      // y-gravity
            C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // specialized Fourier Level-Set
            // =============================


            int numSp = 640;
            C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                Timestepper = FourierLevelSet_Timestepper.ExplicitEuler,
                InterpolationType = Interpolationtype.CubicSplineInterpolation,
            };


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

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
        /// Control for TestProgramm (not to be changed!!!)
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xkelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control CW_Test(int p = 2, int xkelem = 32, string _DbPath = null) {

            //int p = 2;
            //int xkelem = 32;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CWp3_spatialConv";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_CapillaryWave";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryWave";
            C.Tags.Add("Popinet");
            C.Tags.Add("Testcase1");

            C.LogValues = XNSE_Control.LoggingValues.Wavelike;

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

            C.AddBoundaryCondition("wall_lower", "VelocityX#A", X => 0.0);
            C.AddBoundaryCondition("wall_upper", "VelocityX#B", X => 0.0);

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - PeriodicFunc(X[0])));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid(restartSession);
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion

            // additional parameters
            // =====================

            double[] param = new double[3];
            param[0] = lambda;  // wavelength
            param[1] = A0;      // initial disturbance
            param[2] = 0.0;      // y-gravity
            C.AdditionalParameters = param;

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 100;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            //C.AdvancedDiscretizationOptions.surfTensionMode = SurfaceTensionMode.Curvature_Fourier;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // specialized Fourier Level-Set
            // =============================


            int numSp = 640;
            C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                Timestepper = FourierLevelSet_Timestepper.ExplicitEuler,
                InterpolationType = Interpolationtype.CubicSplineInterpolation,
            };


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            double rho = rho_l;         // Testcase1
            double dt = Math.Sqrt(rho * Math.Pow((1 / (double)xkelem), 3) / (Math.PI * sigma));             // !!!
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            double omega0 = Math.Sqrt(sigma * Math.Pow((2 * Math.PI / lambda), 3) / (2.0 * rho));
            C.NoOfTimesteps = 10; // (int)Math.Ceiling((25 / omega0) / dt);                                        // !!!
            //C.NoOfTimesteps = (int)Math.Ceiling(t_end / dt);                                     // !!!

            C.saveperiod = 1;

            #endregion


            return C;

        }

    }


}
