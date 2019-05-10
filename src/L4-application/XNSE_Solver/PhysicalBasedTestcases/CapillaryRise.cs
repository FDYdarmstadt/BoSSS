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

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the Capillary rise testcases
    /// </summary>
    public static class CapillaryRise {


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CapillaryRise_Tube(int p = 2, int kelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            _DbPath = @"\\fdyprime\userspace\smuda\cluster\CapillaryRise\CapillaryRise_db";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryRise";
            C.ProjectDescription = "Capillary rise in tube";

            C.ContinueOnIoError = false;

            C.LogValues = XNSE_Control.LoggingValues.CapillaryHeight;
            C.LogPeriod = 10;

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

            C.PhysicalParameters.rho_A = 9.97e-4;    
            C.PhysicalParameters.rho_B = 1.204e-6;     
            C.PhysicalParameters.mu_A = 1e-5;        
            C.PhysicalParameters.mu_B = 17.1e-8;     
            C.PhysicalParameters.Sigma = 72.75e-3;

            C.PhysicalParameters.betaS_A = 1e-3;
            C.PhysicalParameters.betaS_B = 1e-5;

            C.PhysicalParameters.betaL = 1e-4;
            C.PhysicalParameters.theta_e = Math.PI / 3.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 0.75e-1;
            double H = 10 * R;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-R, R, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 5 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "pressure_outlet_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                grd.EdgeTagNames.Add(4, "navierslip_linear_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + R) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - R) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double h0 = 1.0e-1;

            Func<double[], double> PhiFunc = (X => X[1] - h0);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("GravityY", X => -9.81e2);

            C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("ae0490a8-7ae3-474d-978e-ccdd1ed8726a"); 
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("pressure_outlet_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("navierslip_linear_left");
            C.AddBoundaryValue("navierslip_linear_right");


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 50000;
            C.saveperiod = 10;

            #endregion


            return C;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CapillaryRise_Tube_SFB1194(int p = 2, int kelemR = 2, int omegaTc = 3, bool startUp = false, bool symmetry = true, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";
            _DbPath = @"\\HPCCLUSTER\hpccluster-scratch\smuda\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/CapillaryRise";
            C.ProjectDescription = "A comparative study for SFB 1194";

            C.ContinueOnIoError = false;

            C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            C.LogPeriod = 100;

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

            // the gaseous phase should be close to air
            //C.PhysicalParameters.rho_B = 1.204e-6;  // kg / cm^3
            //C.PhysicalParameters.mu_B = 17.1e-8;    // kg / cm s

            double g = 0;
            double R = 0;
            double H = 0;
            double t_end = 0;
            double dt = 0;
            double t_startUp = 0;
            double dt_startUp = 0;
            Guid restartID = new Guid();
            int ts_restart = 0;
            switch(omegaTc) {
                case 1: {
                        C.Tags.Add("omega = 0.1");
                        R = 5e-3;           // cm
                        H = 4e-2;

                        C.PhysicalParameters.rho_A = 1663.8;
                        C.PhysicalParameters.mu_A = 0.01;
                        C.PhysicalParameters.Sigma = 0.2;       // kg / s^2

                        C.PhysicalParameters.rho_B = 1663.8 / 1000;
                        C.PhysicalParameters.mu_B = 0.01 / 1000;

                        C.PhysicalParameters.betaS_A = 0;
                        C.PhysicalParameters.betaS_B = 0;

                        C.PhysicalParameters.betaL = 0;
                        C.PhysicalParameters.theta_e = 3.0 * Math.PI / 18.0;

                        g = 1.04;         // cm / s^2

                        t_end = 13.86;
                        t_startUp = 0.5204; // 0.49032;
                        dt = 2.5e-5; // dt = 1e-4; // dt = 2e-5;
                        dt_startUp = 3e-5;

                        //restartID = new Guid("07ca4397-8eed-4769-b795-9725fe7d3cd7");
                        //restartID = new Guid("fa8454ce-c05a-4dea-a308-663b6be04ff7");
                        //restartID = new Guid("6380b408-e043-4ed3-8ae5-819d7566e241");
                        //restartID = new Guid("32d62e31-2b42-4a2c-850a-038befc43072");
                        //ts_restart = 13010;

                        //restartID = new Guid("868a08d9-5b44-4462-9a2c-121cfc07e5db");
                        //restartID = new Guid("c3cef42c-f9e5-4757-a618-e9a4002e445f"); // restart with ReInit
                        //restartID = new Guid("a191b8d0-1cca-431a-bf00-a1231c61dee6"); // restart2 with ReInit

                        //restartID = new Guid("02f7d2b1-7546-45b3-808f-d2b6f89d68bf"); // restart with Reinit highdt (1e-4)
                        //ts_restart = 186000;
                        //restartID = new Guid("34c3ee0e-9244-497b-9382-58b458f7b873");   // restart with Reinit highdt (5e-5)
                        //restartID = new Guid("1d31b648-d19f-4a5a-b374-cde08f8106a9");

                        //restartID = new Guid("a32df5ae-393d-416a-bdd0-24a4437136fb");   //restart with sigma dt (2.5E-5)

                        restartID = new Guid("1d679e3c-03f3-41ee-8f1e-7e80ec497926"); // restart with Reinit semi implicit
                        ts_restart = 270100;

                        //restartID = new Guid("2f9deaa6-fab9-48ac-9279-319a1efa5547");

                        //C.ClearVelocitiesOnRestart = true;
                        //C.ReInitPeriod = 250;

                        break;
                    }
                case 2: {
                        C.Tags.Add("omega = 0.5");
                        R = 5e-3;           // cm
                        H = 3e-2;

                        C.PhysicalParameters.rho_A = 133.0;
                        C.PhysicalParameters.mu_A = 0.01;
                        C.PhysicalParameters.Sigma = 0.1;       // kg / s^2

                        C.PhysicalParameters.betaS_A = 1e-3 / Math.Min(R / kelemR, H / (5 * kelemR));
                        C.PhysicalParameters.betaS_B = 1e-5 / Math.Min(R / kelemR, H / (5 * kelemR));

                        C.PhysicalParameters.betaL = 0;
                        C.PhysicalParameters.theta_e = 3.0 * Math.PI / 18.0;

                        g = 6.51;         // cm / s^2

                        t_end = 1.11;
                        t_startUp = 0.15;
                        dt = 2e-5;
                        dt_startUp = 4e-5;

                        restartID = new Guid("59e43164-1c32-4ece-b5ab-257844198fa4");   

                        break;
                    }
                case 3: {
                        C.Tags.Add("omega = 1");
                        R = 5e-3;           // cm
                        H = 4e-2;

                        C.PhysicalParameters.rho_A = 83.1;
                        C.PhysicalParameters.mu_A = 0.01;
                        C.PhysicalParameters.Sigma = 0.04;       // kg / s^2

                        C.PhysicalParameters.rho_B = 83.1 / 1000;
                        C.PhysicalParameters.mu_B = 0.01 / 1000;

                        C.PhysicalParameters.betaS_A = 0; // 8;
                        C.PhysicalParameters.betaS_B = 0; // 0.008;

                        C.PhysicalParameters.betaL = 0; // 0.04; // 4.004;
                        C.PhysicalParameters.theta_e = 3.0 * Math.PI / 18.0;

                        g = 4.17;         // cm / s^2

                        t_end = 0.7;
                        t_startUp = 0.1955;
                        dt = 1e-5;
                        dt_startUp = 5e-4;

                        //restartID = new Guid("e2a38f38-bcdb-4588-bd87-9914dc2989e4");   //startUp
                        //restartID = new Guid("3a1136f2-5363-43b0-8084-5c2ee6ce9d06");   //restart

                        //restartID = new Guid("f37c9194-1bfb-4250-8dc6-a4d1bbe01ed9");

                        restartID = new Guid("c572378f-edd9-4b96-9917-bb037ae4fdac");   //startUp2
                        C.ClearVelocitiesOnRestart = true;

                        break;
                    }
                case 4: {
                        C.Tags.Add("omega = 10");
                        R = 5e-3;           // cm
                        H = 3e-2;

                        C.PhysicalParameters.rho_A = 3.3255;
                        C.PhysicalParameters.mu_A = 0.01;
                        C.PhysicalParameters.Sigma = 0.01;       // kg / s^2

                        C.PhysicalParameters.betaS_A = 1e-3 / Math.Min(R / kelemR, H / (5 * kelemR));
                        C.PhysicalParameters.betaS_B = 1e-5 / Math.Min(R / kelemR, H / (5 * kelemR));

                        C.PhysicalParameters.betaL = 0;
                        C.PhysicalParameters.theta_e = 3.0 * Math.PI / 18.0;

                        g = 26.042;         // cm / s^2

                        t_end = 2.7713;
                        t_startUp = 0.098;
                        dt = 4e-4;
                        dt_startUp = 1e-5;

                        restartID = new Guid("26a16f96-a657-4834-8216-89d30a04c938");
                        break;
                    }
                case 5: {
                        C.Tags.Add("omega = 100");
                        R = 5e-3;           // cm
                        H = 3e-2;

                        C.PhysicalParameters.rho_A = 0.33255;
                        C.PhysicalParameters.mu_A = 0.01;
                        C.PhysicalParameters.Sigma = 0.001;       // kg / s^2

                        C.PhysicalParameters.betaS_A = 1e-3 / Math.Min(R / kelemR, H / (5 * kelemR));
                        C.PhysicalParameters.betaS_B = 1e-5 / Math.Min(R / kelemR, H / (5 * kelemR));

                        C.PhysicalParameters.betaL = 0;
                        C.PhysicalParameters.theta_e = 3.0 * Math.PI / 18.0;

                        g = 26.042;         // cm / s^2

                        t_end = 27.713;
                        t_startUp = 0.098;
                        dt = 1e-3;
                        dt_startUp = 1e-3;

                        restartID = new Guid("b59abf15-1358-48c1-8b79-6b10e1c1d729");
                        break;
                    }
            } 

            C.PhysicalParameters.IncludeConvection = !startUp;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            if(startUp) {

                C.GridFunc = delegate () {
                    double[] Xnodes;
                    if(symmetry)
                        Xnodes = GenericBlas.Linspace(0, R, kelemR + 1);
                    else
                        Xnodes = GenericBlas.Linspace(-R, R, 2 * kelemR + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, H, 8 * kelemR + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "wall_lower");
                    grd.EdgeTagNames.Add(2, "pressure_outlet_upper");

                    if(symmetry)
                        grd.EdgeTagNames.Add(3, "slipsymmetry_left");
                    else
                        grd.EdgeTagNames.Add(3, "navierslip_linear_left");

                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    //grd.EdgeTagNames.Add(4, "navierslip_localized_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[1] - H) <= 1.0e-8)
                            et = 2;
                        if(symmetry) {
                            if(Math.Abs(X[0]) <= 1.0e-8)
                                et = 3;
                        } else {
                            if(Math.Abs(X[0] + R) <= 1.0e-8)
                                et = 3;
                        }
                        if(Math.Abs(X[0] - R) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double h0 = 1e-2;

            Func<double[], double> PhiFunc = (X => X[1] - h0);

            if(startUp) {
                C.InitialValues_Evaluators.Add("Phi", PhiFunc);

                C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
                C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            } else {

                if(ts_restart > 0)
                    C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restartID, ts_restart);
                else
                    C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restartID, null);

            }

            #endregion


            // boundary conditions
            // ===================
            #region BC

            if(startUp) {
                C.AddBoundaryValue("wall_lower");
            } else {
                //C.AddBoundaryValue("wall_lower");
                C.ChangeBoundaryCondition("wall_lower", "pressure_outlet_lower");
                C.AddBoundaryValue("pressure_outlet_lower");
            }
            C.AddBoundaryValue("pressure_outlet_upper");

            if(symmetry)
                C.AddBoundaryValue("slipsymmetry_left");
            else
                C.AddBoundaryValue("navierslip_linear_left");

            //C.ChangeBoundaryCondition("navierslip_localized_right", "navierslip_linear_right");
            C.AddBoundaryValue("navierslip_linear_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            C.PhysicalParameters.sliplength = 1e-3;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 20;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = 1 * C.PhysicalParameters.Sigma;
            //C.PhysicalParameters.lambda_I = 2 * C.PhysicalParameters.Sigma;
            C.AdvancedDiscretizationOptions.UseLevelSetStabilization = false;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 1;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            //int Nsp = 256;
            //C.FourierLevSetControl = new FourierLevSetControl(FourierType.Planar, Nsp, R, X => h0, R/(double)Nsp);
            //C.FourierLevSetControl.Timestepper = FourierLevelSet_Timestepper.RungeKutta1901;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            C.dtMax = (startUp) ? dt_startUp : dt;
            C.dtMin = (startUp) ? dt_startUp : dt;
            C.Endtime = (startUp) ? Math.Max(Math.Sqrt(2 * R / g), 2 * R * C.PhysicalParameters.mu_A / C.PhysicalParameters.Sigma) * 4.0 : (t_startUp + t_end);
            C.NoOfTimesteps = (int)(C.Endtime / C.dtMin);
            C.saveperiod = 100;

            #endregion


            return C;
        }
    }
}
