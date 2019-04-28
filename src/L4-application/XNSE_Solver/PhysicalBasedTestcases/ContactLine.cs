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

    public static class ContactLine {


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control StaticDroplet_OnPlate(int p = 2, int kelem = 32, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            int D = 2;

            if(D == 3)
                C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            AppControl._CompMode compMode = AppControl._CompMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            //C.ProjectDescription = "Static droplet on plate";

            C.ContinueOnIoError = false;

            C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;

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
            if(D == 3) {
                C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }
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

            C.Tags.Add("Reusken");
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            double sigma = 0.5;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 3.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid


            if(D == 2) {

                double xSize = 0.25;
                double ySize = 0.25;

                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem / 2 + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            if(D == 3) {

                double xSize = 0.125;
                double ySize = 0.125;
                double zSize = 0.125;

                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-ySize, ySize, kelem + 1);
                    double[] Znodes = GenericBlas.Linspace(0, zSize, kelem / 2 + 1);
                    var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    grd.EdgeTagNames.Add(5, "navierslip_linear_front");
                    grd.EdgeTagNames.Add(6, "navierslip_linear_back");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[2]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[2] - zSize) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;
                        if(Math.Abs(X[1] + ySize) <= 1.0e-8)
                            et = 5;
                        if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                            et = 6;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.1;
            double Theta_e = Math.PI / 2.0;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0, -h };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);
            //Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());         // quadratic

            //C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double pJump = sigma / R;
            if(D == 3)
                pJump *= 2.0;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            if(D == 3) {
                center = new double[] { 0, 0, -h };
                PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() + (X[2] - center[2]).Pow2()).Sqrt() - R);
                //PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() + (X[2] - center[2]).Pow2()) - R.Pow2());    // quadratic
            }

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("navierslip_linear_lower");
            C.AddBoundaryValue("navierslip_linear_upper");
            C.AddBoundaryValue("navierslip_linear_left");
            C.AddBoundaryValue("navierslip_linear_right");

            if(D == 3) {
                C.AddBoundaryValue("navierslip_linear_front");
                C.AddBoundaryValue("navierslip_linear_back");
            }

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = 0.001;

            #endregion


            // exact solution
            // ==============

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            //C.ExactSolutionPressure.Add("B", (X, t) => 0.0);



            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.ComputeEnergy = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = (compMode == AppControl._CompMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.ContactLineRefined;
            C.BaseRefinementLevel = 1;
            C.RefinementLevel = 12;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CheckJumpConditions = true;

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = (compMode == AppControl._CompMode.Steady) ? LevelSetHandling.None : LevelSetHandling.LieSplitting;

            C.CompMode = compMode;
            double dt = 1e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 4000;
            C.saveperiod = 1;

            #endregion


            return C;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control SlidingDropletToEquilibriumState(int tc = 2, int p = 2, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            int D = 2;

            _DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Static droplet sliding to Equilibrium";

            C.ContinueOnIoError = false;

            C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;

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
            if(D == 3) {
                C.FieldOptions.Add("VelocityZ", new FieldOpts() {
                    Degree = p,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
            }
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

            double sigma = 0.0;
            switch(tc) {
                case 1: {

                        C.Tags.Add("Reusken");
                        C.PhysicalParameters.rho_A = 1;
                        C.PhysicalParameters.rho_B = 1;
                        C.PhysicalParameters.mu_A = 1;
                        C.PhysicalParameters.mu_B = 1;
                        sigma = 0.5;
                        C.PhysicalParameters.Sigma = sigma;

                        C.PhysicalParameters.betaS_A = 0.05;
                        C.PhysicalParameters.betaS_B = 0.05;

                        C.PhysicalParameters.betaL = 0;
                        C.PhysicalParameters.theta_e = Math.PI / 3.0;
                        break;
                    }
                case 2: {

                        C.Tags.Add("Reusken");
                        C.Tags.Add("dynamic testcase");
                        C.PhysicalParameters.rho_A = 2;
                        C.PhysicalParameters.rho_B = 1;
                        C.PhysicalParameters.mu_A = 2;
                        C.PhysicalParameters.mu_B = 1;
                        sigma = 0.4;
                        C.PhysicalParameters.Sigma = sigma;

                        C.PhysicalParameters.betaS_A = 0.5;
                        C.PhysicalParameters.betaS_B = 0.5;

                        C.PhysicalParameters.betaL = 0;
                        C.PhysicalParameters.theta_e = Math.PI / 3.0;
                        break;
                    }
                case 3: {

                        C.Tags.Add("Buscaglia");
                        C.PhysicalParameters.rho_A = 1;
                        C.PhysicalParameters.rho_B = 1;
                        C.PhysicalParameters.mu_A = 0.2e-6;
                        C.PhysicalParameters.mu_B = 1e-5;
                        sigma = 0.075;
                        C.PhysicalParameters.Sigma = sigma;

                        C.PhysicalParameters.betaS_A = 1e-5;
                        C.PhysicalParameters.betaS_B = 1e-5;

                        C.PhysicalParameters.betaL = 0.005;
                        C.PhysicalParameters.theta_e = Math.PI / 4.0;
                        break;
                    }
                default:
                    break;
            }

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 0.0;
            double ySize = 0.0;
            double zSize = 0.5;

            int kelem = 0;

            switch(tc) {
                case 1:
                case 2: {
                        xSize = 0.25;
                        ySize = 0.25;
                        //kelem = 4;      // mesh level 0
                        kelem = 16;
                        break;
                    }
                case 3: {
                        xSize = 0.2;
                        ySize = 0.2;
                        kelem = 150;    // 1.3e-3;
                        //kelem = 40;    // 5e-3;
                        break;
                    }
                default:
                    break;
            }

            if(D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, 2 * kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    //grd.EdgeTagNames.Add(2, "wall_upper");
                    //grd.EdgeTagNames.Add(3, "wall_left");
                    //grd.EdgeTagNames.Add(4, "wall_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            if(D == 3) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, xSize, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 1);
                    double[] Znodes = GenericBlas.Linspace(0, zSize, kelem + 1);
                    var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);

                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    grd.EdgeTagNames.Add(3, "navierslip_linear_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    grd.EdgeTagNames.Add(5, "navierslip_linear_front");
                    grd.EdgeTagNames.Add(6, "navierslip_linear_back");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[2]) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[2] - zSize) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                            et = 4;
                        if(Math.Abs(X[1]) <= 1.0e-8)
                            et = 5;
                        if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                            et = 6;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = (tc == 3) ? 0.125 : 0.1;
            double Theta_0 = Math.PI / 2.0;
            double s = 2 * R * Math.Sin(Theta_0);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0, -h };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double Pjump = sigma / R;
            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            if(D == 3) {
                center = new double[] { xSize / 2, ySize / 2, -h };
                PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() + (X[2] - center[2]).Pow2()).Sqrt() - R);
            }

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("64c51add-d166-4540-97a1-463dc8493fde");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("navierslip_linear_lower");
            C.AddBoundaryValue("navierslip_linear_upper");
            C.AddBoundaryValue("navierslip_linear_left");
            C.AddBoundaryValue("navierslip_linear_right");
            //C.AddBoundaryCondition("wall_upper");
            //C.AddBoundaryCondition("wall_left");
            //C.AddBoundaryCondition("wall_right");
            if(D == 3) {
                C.AddBoundaryValue("navierslip_linear_front");
                C.AddBoundaryValue("navierslip_linear_back");
            }

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

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

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 0.0;
            switch(tc) {
                case 1: {
                        //dt = 2e-4;
                        dt = 1e-3;
                        C.Endtime = 1000;
                        C.NoOfTimesteps = 4000;
                        C.saveperiod = 1;
                        break;
                    }
                case 2: {
                        //dt = 2e-4;
                        dt = 1e-3;
                        C.Endtime = 1000;
                        C.NoOfTimesteps = 4000;
                        C.saveperiod = 1;
                        break;
                    }
                case 3: {
                        //dt = 2e-7;
                        dt = 2e-5;
                        C.Endtime = 1000;
                        C.NoOfTimesteps = 400;
                        C.saveperiod = 1;
                        break;
                    }
                default:
                    break;
            }
            C.dtMax = dt;
            C.dtMin = dt;

            #endregion

            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control DropletOnMovingPlate(int kelem = 8, int p = 2, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Droplet";
            C.ProjectDescription = "Static droplet sliding to Equilibrium";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;

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


            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 0.1;
            double sigma = 0.5;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.betaS_A = 0.05;
            //C.PhysicalParameters.betaS_B = 0.05;

            C.PhysicalParameters.betaL = 0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 1.0;
            double ySize = 1.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-3 * xSize, xSize, 4 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + 3 * xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.4;
            double Theta_0 = Math.PI / 2.0;
            double s = 2 * R * Math.Sin(Theta_0);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0, -h };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double Pjump = sigma / R;
            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("64c51add-d166-4540-97a1-463dc8493fde");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double Uwall = -1.0;

            C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", (X, t) => Uwall);
            C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#B", (X, t) => Uwall);
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("pressure_outlet_left");
            C.AddBoundaryValue("pressure_outlet_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.ContactLine;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            C.PhysicalParameters.sliplength = 1e-5;

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

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

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-3;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            C.dtMax = dt;
            C.dtMin = dt;

            #endregion

            return C;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CornerFlowAtContactLine(int kelem = 24, int p = 3, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Classic;

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/ContactLine";
            C.ProjectDescription = "Static corner flow at contact line";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;

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
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1e-3;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1e-3;
            double sigma = 1.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.betaS_A = 0.05;
            //C.PhysicalParameters.betaS_B = 0.05;

            C.PhysicalParameters.betaL = 0;
            double theta_e = Math.PI * 1.0 / 4.0;
            C.PhysicalParameters.theta_e = theta_e;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 2e-1;
            double ySize = 1e-1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, 2 * (2 * kelem) + 2);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init


            Func<double[], double> PhiFunc = (X => X[1]/Math.Tan(theta_e) + X[0]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double Pjump = 0.0;
            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("64c51add-d166-4540-97a1-463dc8493fde");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double Uwall = -1.0;

            C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", (X, t) => Uwall);
            C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#B", (X, t) => Uwall);
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("pressure_outlet_left");
            C.AddBoundaryValue("pressure_outlet_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            C.PhysicalParameters.sliplength = 1e-5;

            #endregion


            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.ComputeEnergy = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;
            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib; // _DropIndefinite;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefinementLevel = 1;
            C.RefineNavierSlipBoundary = true;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-5;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            C.dtMax = dt;
            C.dtMin = dt;

            #endregion

            return C;
        }



    }
}
