﻿/* =======================================================================
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
    /// class providing Controls for the Couette flow type testcases with two phases
    /// </summary>
    public static class TwoPhaseCouetteFlow {

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control Couette_GNBC(int tc = 2, int p = 3, int kelem = 12, double dt = 0.02, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
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
            //C.LogPeriod = 1;
            C.PostprocessingModules.Add(new MovingContactLineLogging());
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
                        // symmetric
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
                        // un-symmetric
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
                double[] Xnodes = GenericBlas.Linspace(0, 4 * L, 8 * kelem + 0);
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

            C.InitSignedDistance = false;

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

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            //C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.FullBoussinesqScriven;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.ContactLineRefined;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 0;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            //double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 200;
            C.NoOfTimesteps = (int)(200.0 / (dt));
            C.saveperiod = 2;

            #endregion


            return C;
        }


        public static XNSE_Control CouetteGNBC_forWorksheet(bool symmetric, bool restart) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = true;
            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            C.PostprocessingModules.Add(new MovingContactLineLogging());
            #endregion


            // DG degrees
            // ==========
            #region degrees

            // need to be set by user via setDGdegree() in worksheet

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 0.81;
            C.PhysicalParameters.rho_B = 0.81;
            C.PhysicalParameters.mu_A = 1.95;
            C.PhysicalParameters.mu_B = 1.95;
            C.PhysicalParameters.Sigma = 5.5;

            if(symmetric) {
                C.PhysicalParameters.betaS_A = 1.5;
                C.PhysicalParameters.betaS_B = 1.5;

                C.PhysicalParameters.betaL = 0.0;
                C.PhysicalParameters.theta_e = Math.PI / 2.0;

            } else {
                C.PhysicalParameters.betaS_A = 0.591;
                C.PhysicalParameters.betaS_B = 1.5;

                C.PhysicalParameters.betaL = 0.0;
                C.PhysicalParameters.theta_e = Math.Acos(0.38);

            }

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

            if (symmetric) {
                C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", "X => -0.25", false);
                C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#B", "X => -0.25", false);
                C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#A", "X => 0.25", false);
                C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#B", "X => 0.25", false);
            } else {
                C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", "X => -0.2", false);
                C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#B", "X => -0.2", false);
                C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#A", "X => 0.2", false);
                C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#B", "X => 0.2", false);

            }

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;


            #endregion


            // Initial Values
            // ==============
            #region init

            //double epsInit = 0.01;
            //Func<double[], double> PhiFunc = (X => (Math.Abs(X[0] - 2 * L) - L) + (X[1] - (H / 2.0)) * (epsInit / (H / 2.0)));
            if(!restart)
                C.AddInitialValue("Phi", "X => (Math.Abs(X[0] - 2 * 27.2) - 27.2) + (X[1] - (6.8)) * (0.01 / 6.8)", false);

            #endregion


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

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
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

    }

}
