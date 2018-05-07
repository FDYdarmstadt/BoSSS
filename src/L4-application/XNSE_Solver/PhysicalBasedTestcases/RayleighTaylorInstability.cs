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
    /// class providing Controls for the Rayleigh Taylor Instability Testcase
    /// </summary>
    public static class RayleighTaylorInstability {

        /// <summary>
        /// Control for various testing
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xkelem"></param>
        /// <param name="dt"></param>
        /// <param name="t_end"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control RT(int p = 2, int xkelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            _DbPath = @"D:\local\local_Testcase_databases\Testcase_RTinstability";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/RT-Instability";

            C.LogValues = XNSE_Control.LoggingValues.Wavelike;
            C.WriteInterfaceP = true;

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

            #endregion


            // Physical Parameters
            // ===================
            #region physics  


            // Capillary wave: Laplace number La = 3e5: 
            //C.Tags.Add("La=3e5");
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-5;
            //double mu_h = 1e-5;
            //double sigma = 3e-2;

            // Capillary wave: Laplace number La = 3000: 
            //C.Tags.Add("La3000");
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-4;
            //double mu_h = 1e-4;
            //double sigma = 3e-2;

            // Capillary wave: Laplace number La = 120: 
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 5e-4;
            //double mu_h = 5e-4;
            //double sigma = 3e-2;

            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-6;
            //double mu_h = 1e-4;
            //double sigma = 3e-2;

            //double rho_l = 1;
            //double rho_h = 1;
            //double mu_l = 1e-3;
            //double mu_h = 1e-3;
            //double sigma = 1;

            // Air(light)-Water(heavy)
            //double rho_h = 1e-3;          // kg / cm^3
            //double rho_l = 1.2e-6;        // kg / cm^3
            //double mu_h = 1e-5;           // kg / cm * s
            //double mu_l = 17.1e-8;        // kg / cm * s
            //double sigma = 72.75e-3;      // kg / s^2         // lambda_crit = 1.7121 ; lambda < lambda_c: stable

            // same kinematic viscosities
            //double rho_l = 1e-5;          // kg / cm^3
            //double mu_l = 1e-8;           // kg / cm * s
            //double rho_h = 1e-3;        // kg / cm^3
            //double mu_h = 1e-6;        // kg / cm * s         // Atwood number = 0.98

            //double rho_l = 7e-4;          // kg / cm^3
            //double mu_l = 7e-5;           // kg / cm * s
            //double rho_h = 1e-3;        // kg / cm^3
            //double mu_h = 1e-4;        // kg / cm * s           // Atwood number = 0.1765

            //double sigma = 72.75e-3;      // kg / s^2          


            //double rho_l = 1e-3;
            //double mu_l = 1e-4;
            //double rho_h = 1;
            //double mu_h = 1e-1;

            //double sigma = 100;      // kg / s^2          

            double rho_l = 1e-1;
            double mu_l = 1e-2;
            double rho_h = 10;
            double mu_h = 1;

            double sigma = 50;      // kg / s^2          

            // Water-Oil 
            //double rho_ = 8.63e-4;
            //double rho_B = 1.2e-6;
            //double mu_A = 2e-4;
            //double mu_B = 17.1e-8;


            // stable configuration
            //C.PhysicalParameters.rho_A = rho_h;
            //C.PhysicalParameters.rho_B = rho_l;
            //C.PhysicalParameters.mu_A = mu_h;
            //C.PhysicalParameters.mu_B = mu_l;
            //C.PhysicalParameters.Sigma = sigma;


            // unstable configuration
            C.PhysicalParameters.rho_A = rho_l;
            C.PhysicalParameters.rho_B = rho_h;
            C.PhysicalParameters.mu_A = mu_l;
            C.PhysicalParameters.mu_B = mu_h;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.Sigma = 0.0;   // free surface boundary condition

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial displacement
            // ===================

            int kmode = 1;
            double H0 = 0.01;

            //double lambda_stable = 1;
            //double lambda_instable = 4;
            //Func<double, double> displacement = x => ((lambda_stable / 100) * Math.Sin(x * 2 * Math.PI / lambda_stable )) + ((lambda_instable / 100) * Math.Sin(x * 2 * Math.PI / lambda_instable));
            Func<double, double> displacement = x => (H0 * Math.Sin(x * 2 * Math.PI * (double)kmode));


            // grid genration
            // ==============
            #region grid

            double L = 1; //lambda_instable;
            double H = 1.0 * L;
            double H_interf = L / 4.0;

            //int xkelem = 20;
            int ykelem_Interface = 1 * xkelem;      // /2
            int ykelem_outer = xkelem / 2;              // /2

            C.GridFunc = delegate () {
                double[] xNodes = GenericBlas.Linspace(0, L, xkelem + 1);
                //double[] Ynodes_Interface = GenericBlas.Linspace(-(H_interf) / 2.0, (H_interf) / 2.0, ykelem_Interface);
                //Ynodes_Interface = Ynodes_Interface.GetSubVector(1, Ynodes_Interface.GetLength(0) - 2);
                //double[] Ynodes_lower = GenericBlas.Linspace(-(H + (H_interf / 2.0)), -(H_interf / 2.0), ykelem_outer + 1);
                //double[] Ynodes_upper = GenericBlas.Linspace((H_interf / 2.0), H + (H_interf / 2.0), ykelem_outer + 1);
                //double[] Ynodes = Ynodes_lower.Concat(Ynodes_Interface).ToArray().Concat(Ynodes_upper).ToArray();
                var _yNodes1 = Grid1D.TanhSpacing(-(H + (H_interf / 2.0)), -(H_interf / 2.0), ykelem_outer + 1, 1.5, false);
                _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                var _yNodes2 = GenericBlas.Linspace(-H_interf / 2.0, H_interf / 2.0, ykelem_Interface);
                _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
                var _yNodes3 = Grid1D.TanhSpacing((H_interf / 2.0), H + (H_interf / 2.0), ykelem_outer + 1, 1.5, true);
                var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);

                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + ((H_interf / 2.0) + H)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ((H_interf / 2.0) + H)) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            // nach Popinet
            //C.GridFunc = delegate () {
            //    double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
            //    double[] Ynodes = GenericBlas.Linspace(-3.0 * L / 2.0, 3.0 * L / 2.0, (3 * xkelem) + 1);
            //    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

            //    grd.EdgeTagNames.Add(1, "wall_lower");
            //    grd.EdgeTagNames.Add(2, "wall_upper");


            //    grd.DefineEdgeTags(delegate (double[] X) {
            //        byte et = 0;
            //        if (Math.Abs(X[1] + (3.0 * L / 2.0)) <= 1.0e-8)
            //            et = 1;
            //        if (Math.Abs(X[1] - (3.0 * L / 2.0)) <= 1.0e-8)
            //            et = 2;
            //        if (Math.Abs(X[0]) <= 1.0e-8)
            //            et = 3;
            //        if (Math.Abs(X[0] - L) <= 1.0e-8)
            //            et = 4;

            //        return et;
            //    });

            //    return grd;
            //};


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double, double> PeriodicFunc = x => displacement(x);

            //superposed higher frequent disturbance
            //double lambda_dist = lambda / (double)(xkelem / 2);
            //if (lambda_dist > 0.0) {
            //    double A0_dist = A0 / 5.0;
            //    Func<double, double> disturbance = x => A0_dist * Math.Sin(x * 2 * Math.PI / lambda_dist);
            //    PeriodicFunc = x => (displacement(x) + disturbance(x));
            //}


            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - (PeriodicFunc(X[0]))));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            double g = 9.81; // e2;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("0141eba2-8d7b-4593-8595-bfff12dbfc40");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LinearSolver = DirectSolver._whichSolver.MUMPS;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            //C.AdvancedDiscretizationOptions.surfTensionMode = SurfaceTensionMode.Curvature_Fourier;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // additional parameters
            // =====================

            double[] param = new double[4];
            param[0] = kmode;   // wavenumber
            param[1] = L;  // wavelength
            param[2] = H0;      // initial disturbance
            param[3] = g;      // y-gravity
            C.AdditionalParameters = param;

            // specialized Fourier Level-Set
            // =============================

            int numSp = 640;
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
            //    FourierEvolve = Fourier_Evolution.MaterialPoints,
            //    Timestepper = FourierLevelSet_Timestepper.RungeKutta1901,
            //    InterpolationType = Interpolationtype.CubicSplineInterpolation
            //};


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            double dt = 1e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 4000;

            C.saveperiod = 4;

            #endregion


            return C;

        }

        /// <summary>
        /// Control for various paramStudies
        /// </summary>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control[] RT_Paramstudy(string _DbPath) {

            List<XNSE_Control> R = new List<XNSE_Control>();

            // for spatial convergence study
            //int[] h_study = new int[] { 16, 32, 64, 128, 256 };
            //double dt = 3e-6;
            //double t_end = dt;

            // for temporal convergence study
            double[] dt_study = new double[] { 8e-5, 4e-5, 2e-5, 1e-5, 5e-6 };
            int h = 64;
            double t_end = dt_study[0];

            // settings for the paramStudy
            int p = 2;                                                  // !!!!!!!!!!!!
            //string _DbPath = @"D:\local\RTp2_temporalConv_coupledCN";

            //foreach (int h_k in h_study) {
            //    var C = RayleighTaylorInstability(p, h_k, dt, t_end, _DbPath);
            //    C.Paramstudy_CaseIdentification = new Tuple<string, object>[] { new Tuple<string, object>("xkelem", h_k),
            foreach (double dt_h in dt_study) {
                var C = RT_Test(p, h, dt_h, t_end, _DbPath);
                C.Paramstudy_CaseIdentification = new Tuple<string, object>[] { new Tuple<string, object>("dt", dt_h),
                };
                R.Add(C);
            }

            return R.ToArray();
        }

        /// <summary>
        /// Control for the testcase according to Smolianksi
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control RT_Smolianski(int p = 2, int kelem = 32, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_RTinstability";
            _DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/RT-Instability";
            C.Tags.Add("unstable");
            //C.Tags.Add("smolianski");

            #endregion

            //C.LogValues = XNSE_Control.LoggingValues.Wavelike;
            //C.WriteInterfaceP = true;


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
                Degree = 4,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 8,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // grid genration
            // ==============
            #region grid

            double L = 1;
            double H = 4.0 * L;

            C.GridFunc = delegate () {
                // Smolianski
                //double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                //double[] Ynodes = GenericBlas.Linspace(0, H, (4 * kelem) + 1);
                //var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                double[] xNodes = GenericBlas.Linspace(0, L, kelem + 1);

                double[] _yNodes1 = GenericBlas.Linspace(0, 2.2, (int)(2.2 * kelem) + 1);
                _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
                var _yNodes2 = Grid1D.TanhSpacing(2.2, 4.0, (kelem / 2) + 1, 1.5, true);

                var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2);


                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };


            //double L = 1; //lambda_instable;
            //double H = 1.0 * L;
            //double H_interf = L / 4.0;

            ////int xkelem = 20;
            //int ykelem_Interface = 1 * kelem / 2;
            //int ykelem_outer = kelem / 2;

            //C.GridFunc = delegate () {
            //    double[] xNodes = GenericBlas.Linspace(0, L, kelem + 1);
            //    //double[] Ynodes_Interface = GenericBlas.Linspace(-(H_interf) / 2.0, (H_interf) / 2.0, ykelem_Interface);
            //    //Ynodes_Interface = Ynodes_Interface.GetSubVector(1, Ynodes_Interface.GetLength(0) - 2);
            //    //double[] Ynodes_lower = GenericBlas.Linspace(-(H + (H_interf / 2.0)), -(H_interf / 2.0), ykelem_outer + 1);
            //    //double[] Ynodes_upper = GenericBlas.Linspace((H_interf / 2.0), H + (H_interf / 2.0), ykelem_outer + 1);
            //    //double[] Ynodes = Ynodes_lower.Concat(Ynodes_Interface).ToArray().Concat(Ynodes_upper).ToArray();
            //    var _yNodes1 = Grid1D.TanhSpacing(-(H + (H_interf / 2.0)), -(H_interf / 2.0), ykelem_outer + 1, 1.5, false);
            //    _yNodes1 = _yNodes1.GetSubVector(0, (_yNodes1.Length - 1));
            //    var _yNodes2 = GenericBlas.Linspace(-H_interf / 2.0, H_interf / 2.0, ykelem_Interface);
            //    _yNodes2 = _yNodes2.GetSubVector(0, (_yNodes2.Length - 1));
            //    var _yNodes3 = Grid1D.TanhSpacing((H_interf / 2.0), H + (H_interf / 2.0), ykelem_outer + 1, 1.5, true);
            //    var yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);

            //    var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: true);

            //    grd.EdgeTagNames.Add(1, "wall_lower");
            //    grd.EdgeTagNames.Add(2, "wall_upper");


            //    grd.DefineEdgeTags(delegate (double[] X) {
            //        byte et = 0;
            //        if (Math.Abs(X[1] + ((H_interf / 2.0) + H)) <= 1.0e-8)
            //            et = 1;
            //        if (Math.Abs(X[1] - ((H_interf / 2.0) + H)) <= 1.0e-8)
            //            et = 2;
            //        if (Math.Abs(X[0]) <= 1.0e-8)
            //            et = 3;
            //        if (Math.Abs(X[0] - L) <= 1.0e-8)
            //            et = 4;

            //        return et;
            //    });

            //    return grd;
            //};

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");
            C.AddBoundaryCondition("freeslip_left");
            C.AddBoundaryCondition("freeslip_right");

            #endregion


            // Initial Values
            // ==============
            #region init

            int k = 1;
            double A0 = 0.05;
            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - (2.0 + A0 * Math.Cos(2.0 * Math.PI * X[0]))));

            C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);

            double g = 1.0;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("12b8c9f7-504c-40c4-a548-304a2e2aa14f");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 3420);

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            //C.PhysicalParameters.rho_A = 0.17;
            //C.PhysicalParameters.rho_B = 1.2;
            //C.PhysicalParameters.mu_A = 0.17e-2;
            //C.PhysicalParameters.mu_B = 1.2e-2;
            //C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.rho_A = 0.17;
            C.PhysicalParameters.rho_B = 1.2;
            C.PhysicalParameters.mu_A = 0.003;
            C.PhysicalParameters.mu_B = 0.003;
            C.PhysicalParameters.Sigma = 0.0;   // corresponding dt = 1e-3;
            //C.PhysicalParameters.Sigma = 0.015;   // corresponding dt = 4e-3;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion



            // additional parameters
            // =====================

            //double[] param = new double[4];
            //param[0] = k;   // wavenumber
            //param[1] = L;  // wavelength
            //param[2] = A0;      // initial disturbance
            //param[3] = g;      // y-gravity
            //C.AdditionalParameters = param;


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.LSContiProjectionMethod = ContinuityProjectionOption.SpecFEM;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LinearSolver = DirectSolver._whichSolver.MUMPS;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.Curvature_Projected;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            double dt = 1e-3;   
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 4000;

            C.saveperiod = 4;

            #endregion

            return C;

        }


        /// <summary>
        /// Control for TestProgramm (not to be changed!!!)
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xkelem"></param>
        /// <param name="dt"></param>
        /// <param name="t_end"></param>
        /// <param name="_DbPath"></param>
        /// <returns></returns>
        public static XNSE_Control RT_Test(int p = 2, int xkelem = 32, double dt = 8e-5, double t_end = 4e-4, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //_DbPath = @"D:\local\local_test_db";

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/RT-Instability";

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


            // Capillary wave: Laplace number La = 3e5: 
            //C.Tags.Add("La=3e5");
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-5;
            //double mu_h = 1e-5;
            //double sigma = 3e-2;

            // Capillary wave: Laplace number La = 3000: 
            //C.Tags.Add("La3000");
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-4;
            //double mu_h = 1e-4;
            //double sigma = 3e-2;

            // Capillary wave: Laplace number La = 120: 
            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 5e-4;
            //double mu_h = 5e-4;
            //double sigma = 3e-2;

            //double rho_l = 1e-3;
            //double rho_h = 1e-3;
            //double mu_l = 1e-6;
            //double mu_h = 1e-4;
            //double sigma = 3e-2;

            //double rho_l = 1;
            //double rho_h = 1;
            //double mu_l = 1e-3;
            //double mu_h = 1e-3;
            //double sigma = 1;

            // Air(light)-Water(heavy)
            //double rho_h = 1e-3;          // kg / cm^3
            //double rho_l = 1.2e-6;        // kg / cm^3
            //double mu_h = 1e-5;           // kg / cm * s
            //double mu_l = 17.1e-8;        // kg / cm * s
            //double sigma = 72.75e-3;      // kg / s^2         // lambda_crit = 1.7121 ; lambda < lambda_c: stable

            // same kinematic viscosities
            //double rho_l = 1e-5;          // kg / cm^3
            //double mu_l = 1e-8;           // kg / cm * s
            //double rho_h = 1e-3;        // kg / cm^3
            //double mu_h = 1e-6;        // kg / cm * s         // Atwood number = 0.98

            //double rho_l = 7e-4;          // kg / cm^3
            //double mu_l = 7e-5;           // kg / cm * s
            //double rho_h = 1e-3;        // kg / cm^3
            //double mu_h = 1e-4;        // kg / cm * s           // Atwood number = 0.1765

            //double sigma = 72.75e-3;      // kg / s^2          


            double rho_l = 1e-1;
            double mu_l = 1e-2;
            double rho_h = 1;
            double mu_h = 1e-1;

            double sigma = 100;      // kg / s^2          

            // Water-Oil 
            //double rho_ = 8.63e-4;
            //double rho_B = 1.2e-6;
            //double mu_A = 2e-4;
            //double mu_B = 17.1e-8;


            // stable configuration
            //C.PhysicalParameters.rho_A = rho_h;
            //C.PhysicalParameters.rho_B = rho_l;
            //C.PhysicalParameters.mu_A = mu_h;
            //C.PhysicalParameters.mu_B = mu_l;
            //C.PhysicalParameters.Sigma = sigma;


            // unstable configuration
            C.PhysicalParameters.rho_A = rho_l;
            C.PhysicalParameters.rho_B = rho_h;
            C.PhysicalParameters.mu_A = mu_l;
            C.PhysicalParameters.mu_B = mu_h;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.Sigma = 0.0;   // free surface boundary condition

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // Initial displacement
            // ===================

            int kmode = 1;
            double H0 = 0.01;

            //double lambda_stable = 1;
            //double lambda_instable = 4;
            //Func<double, double> displacement = x => ((lambda_stable / 100) * Math.Sin(x * 2 * Math.PI / lambda_stable )) + ((lambda_instable / 100) * Math.Sin(x * 2 * Math.PI / lambda_instable));
            Func<double, double> displacement = x => (H0 * Math.Sin(x * 2 * Math.PI * (double)kmode));

            // grid genration
            // ==============
            #region grid

            double L = 1; //lambda_instable;
            double H = 1.0 * L;
            double H_interf = L / 4.0;

            //int xkelem = 20;
            int ykelem_Interface = 1 * (xkelem / 4);
            int ykelem_outer = xkelem / 2;

            //C.GridFunc = delegate () {
            //    double[] Xnodes = GenericBlas.Linspace(0, L, xkelem + 1);
            //    double[] Ynodes_Interface = GenericBlas.Linspace(-(H_interf) / 2.0, (H_interf) / 2.0, ykelem_Interface + 1);
            //    Ynodes_Interface = Ynodes_Interface.GetSubVector(1, Ynodes_Interface.GetLength(0) - 2);
            //    double[] Ynodes_lower = GenericBlas.Linspace(-(H + (H_interf / 2.0)), -(H_interf / 2.0), ykelem_outer + 1);
            //    double[] Ynodes_upper = GenericBlas.Linspace((H_interf / 2.0), H + (H_interf / 2.0), ykelem_outer + 1);
            //    double[] Ynodes = Ynodes_lower.Concat(Ynodes_Interface).ToArray().Concat(Ynodes_upper).ToArray();


            //    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

            //    grd.EdgeTagNames.Add(1, "wall_lower");
            //    grd.EdgeTagNames.Add(2, "wall_upper");


            //    grd.DefineEdgeTags(delegate (double[] X) {
            //        byte et = 0;
            //        if (Math.Abs(X[1] + ((H_interf / 2.0) + H)) <= 1.0e-8)
            //            et = 1;
            //        if (Math.Abs(X[1] - ((H_interf / 2.0) + H)) <= 1.0e-8)
            //            et = 2;
            //        if (Math.Abs(X[0]) <= 1.0e-8)
            //            et = 3;
            //        if (Math.Abs(X[0] - L) <= 1.0e-8)
            //            et = 4;

            //        return et;
            //    });

            //    return grd;
            //};

            // nach Popinet
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

            C.AddBoundaryCondition("wall_lower");
            C.AddBoundaryCondition("wall_upper");

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double, double> PeriodicFunc = x => displacement(x);

            //superposed higher frequent disturbance
            //double lambda_dist = lambda / (double)(xkelem / 2);
            //if (lambda_dist > 0.0) {
            //    double A0_dist = A0 / 5.0;
            //    Func<double, double> disturbance = x => A0_dist * Math.Sin(x * 2 * Math.PI / lambda_dist);
            //    PeriodicFunc = x => (displacement(x) + disturbance(x));
            //}


            C.InitialValues_Evaluators.Add("Phi", (X => X[1] - (PeriodicFunc(X[0]))));

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            double g = 9.81e2;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("0141eba2-8d7b-4593-8595-bfff12dbfc40");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.Fourier;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            //C.AdvancedDiscretizationOptions.surfTensionMode = Solution.XNSECommon.SurfaceTensionMode.Curvature_Projected;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // additional parameters
            // =====================

            double[] param = new double[3];
            param[0] = L;  // wavelength
            param[1] = H0;      // initial disturbance
            param[2] = g;      // y-gravity
            C.AdditionalParameters = param;

            // specialized Fourier Level-Set
            // =============================

            int numSp = 640;
            C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, numSp, L, PeriodicFunc, 1.0 / (double)xkelem) {
                FourierEvolve = Fourier_Evolution.MaterialPoints,
                Timestepper = FourierLevelSet_Timestepper.RungeKutta1901,
                InterpolationType = Interpolationtype.CubicSplineInterpolation
            };


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.dt_increment = 20;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            //double dt = Math.Sqrt((rho_h+rho_l)*Math.Pow(1.0/(double)xkelem,3)/(2*Math.PI*sigma));
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = (int)Math.Ceiling(t_end / dt);

            C.saveperiod = 1;

            #endregion


            return C;

        }



    }


}

