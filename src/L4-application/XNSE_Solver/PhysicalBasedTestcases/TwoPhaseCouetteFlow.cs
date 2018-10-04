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
    /// class providing Controls for the Couette flow type testcases with two phases
    /// </summary>
    public static class TwoPhaseCouetteFlow {

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control Couette_GNBC(int tc = 1, int p = 2, int kelem = 16, double dt = 0.2, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            _DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Couette";
            C.ProjectDescription = "Couette flow by Gerbeau";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.LogPeriod = 10;

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


            // Physical Parameters
            // ===================
            #region physics

            switch (tc) {
                case 1: {
                        C.PhysicalParameters.rho_A = 0.81;
                        C.PhysicalParameters.rho_B = 0.81;
                        C.PhysicalParameters.mu_A = 1.95;
                        C.PhysicalParameters.mu_B = 1.95;
                        C.PhysicalParameters.Sigma = 5.5;

                        C.PhysicalParameters.betaS_A = 1.5;
                        C.PhysicalParameters.betaS_B = 1.5;

                        C.PhysicalParameters.betaL = 0.0;
                        C.PhysicalParameters.theta_e = Math.PI / 2.0;
                        break;
                    }
                case 2: {
                        C.PhysicalParameters.rho_A = 0.81;
                        C.PhysicalParameters.rho_B = 0.81;
                        C.PhysicalParameters.mu_A = 1.95;
                        C.PhysicalParameters.mu_B = 1.95;
                        C.PhysicalParameters.Sigma = 5.5;

                        C.PhysicalParameters.betaS_A = 0.591;
                        C.PhysicalParameters.betaS_B = 1.5;

                        C.PhysicalParameters.betaL = 0.0;
                        C.PhysicalParameters.theta_e = Math.Acos(0.38);
                        break;
                    }
                case 3: {
                        C.PhysicalParameters.rho_A = 8.1;
                        C.PhysicalParameters.rho_B = 8.1;
                        C.PhysicalParameters.mu_A = 1.95;
                        C.PhysicalParameters.mu_B = 1.95;
                        C.PhysicalParameters.Sigma = 0.55;

                        C.PhysicalParameters.betaS_A = 1.5;
                        C.PhysicalParameters.betaS_B = 1.5;

                        C.PhysicalParameters.betaL = 0.0;
                        C.PhysicalParameters.theta_e = Math.PI / 2.0;
                        break;
                    }
            }

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 27.2;
            double H = 13.6;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 4 * L, 8 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                grd.EdgeTagNames.Add(2, "navierslip_linear_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            //C.GridFunc = delegate () {
            //    double[] Xnodes = GenericBlas.Linspace(0, 2, 9 + 1);
            //    double[] Ynodes = GenericBlas.Linspace(0, 1, 5 + 1);
            //    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

            //    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
            //    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");

            //    grd.DefineEdgeTags(delegate (double[] X) {
            //        byte et = 0;
            //        if (Math.Abs(X[1]) <= 1.0e-8)
            //            et = 1;
            //        if (Math.Abs(X[1] - 1) <= 1.0e-8)
            //            et = 2;

            //        return et;
            //    });

            //    return grd;
            //};

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => Math.Abs(X[0] - 2 * L) - L);
            //Func<double[], double> PhiFunc = (X => Math.Abs(X[0] - 1) - 0.5);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //double U_slip = 0.16;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0); // -U_slip + 2 * U_slip * X[1] / H);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0); // -U_slip + 2 * U_slip * X[1] / H);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double U_wall = 0.0;
            switch(tc) {
                case 1:
                case 3: {
                        U_wall = 0.25;
                        break;
                    }
                case 2: {
                        U_wall = 0.2;
                        break;
                    }
            }

            C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", X => -U_wall);
            C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#B", X => -U_wall);
            C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#A", X => U_wall);
            C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#B", X => U_wall);

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 2;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.NonLinearSolver = NonlinearSolverMethod.Picard;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            //C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.FullBoussinesqScriven;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;

            
            //C.AdaptiveMeshRefinement = true;
            //C.RefinementLevel = 1;

            //C.LS_TrackerWidth = 2;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            //double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 200;
            C.NoOfTimesteps = (int)(200.0 / (dt));
            C.saveperiod = 2;

            #endregion


            return C;
        }

    }

}
