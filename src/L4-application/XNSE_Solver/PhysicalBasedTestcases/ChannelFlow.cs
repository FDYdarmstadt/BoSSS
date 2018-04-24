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
    /// class providing Controls for the channel flow type testcases
    /// </summary>
    public static class ChannelFlow {


        /// <summary>
        /// control object for various testing
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control ChannelFlow_WithInterface(int p = 2, int kelem = 16, int wallBC = 1) {

            XNSE_Control C = new XNSE_Control();

            string _DbPath = @"D:\local\local_test_db";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/Channel";
            C.ProjectDescription = "Channel flow with vertical interface";

            C.ContinueOnIoError = false;

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

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            double sigma = 0.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.beta_S = 0.05;
            C.PhysicalParameters.Theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 2;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                switch (wallBC) {
                    case 0:
                        goto default;
                    case 1:
                        grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                        grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                        break;
                    case 2:
                        grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                        grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                        break;
                    default:
                        grd.EdgeTagNames.Add(1, "wall_lower");
                        grd.EdgeTagNames.Add(2, "wall_upper");
                        break;

                }
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

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

            #endregion


            // Initial Values
            // ==============
            #region init

            //Func<double[], double> PhiFunc = (X => X[0] - L / 2);

            double[] center = new double[] { H / 2.0, H / 2.0 };
            double radius = 0.25;


            C.InitialValues_Evaluators.Add("Phi",
                //(X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2())   // quadratic form
                (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius)  // signed-distance form
                );

            //C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 1.0;

            C.InitialValues_Evaluators.Add("VelocityX#A", X => U);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => U);

            double Pjump = sigma / radius;
            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);


            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("cf6bd7bf-a19f-409e-b8c2-0b89388daad6");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, 10);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            switch(wallBC) {
                case 0:
                    goto default;
                case 1:
                    C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#B", X => U);
                    C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#B", X => U);
                    break;
                case 2:
                    C.AddBoundaryCondition("navierslip_linear_lower");
                    C.AddBoundaryCondition("navierslip_linear_upper");
                    break;
                default:
                    C.AddBoundaryCondition("wall_lower");
                    C.AddBoundaryCondition("wall_upper");
                    break;

            }
            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => U);
            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", X => U);
            C.AddBoundaryCondition("pressure_outlet_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux;

            //C.LS_TrackerWidth = 2;
            //C.AdaptiveMeshRefinement = true;
            //C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 2e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 20;
            C.saveperiod = 10;

            #endregion


            return C;
        }



        // ==========================================
        // Control-objects for elementalTestProgramm
        // ==========================================


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CF_BoundaryTest(int bc = 3, bool Xperiodic = true, bool init_exact = false, bool slip = false) {

            int p = 2;
            int kelem = 16;


            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null; //_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Channel flow for BC testing";

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
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 2;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: Xperiodic);

                switch (bc) {
                    case 1: {
                            grd.EdgeTagNames.Add(1, "freeslip_lower");
                            grd.EdgeTagNames.Add(2, "freeslip_upper");
                            break;
                        }
                    case 2: {
                            grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                            grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                            break;
                        }
                    case 3: {
                            grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                            grd.EdgeTagNames.Add(2, "freeslip_upper");
                            break;
                        }
                    case 4: {
                            grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                            grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                            break;
                        }
                    default: {
                            throw new NotImplementedException("No such testcase available");
                        }
                }

                if (!Xperiodic) {
                    grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                    grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (!Xperiodic) {
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 1.0;

            if (init_exact) 
                C.InitialValues_Evaluators.Add("VelocityX#A", X => U);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            switch (bc) {
                case 1: {
                        C.AddBoundaryCondition("freeslip_lower");
                        if (slip)
                            C.AddBoundaryCondition("freeslip_upper", "VelocityX#A", X => -U);
                        else
                            C.AddBoundaryCondition("freeslip_upper");

                        if (!Xperiodic) {
                            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryCondition("pressure_outlet_right");
                        }
                        break;
                    }
                case 2: {
                        //C.AddBoundaryCondition("navierslip_linear_lower");
                        //if (slip)
                        //    C.AddBoundaryCondition("navierslip_linear_upper", "VelocityX#A", X => -U);
                        //else
                        //    C.AddBoundaryCondition("navierslip_linear_upper");

                        C.AddBoundaryCondition("navierslip_linear_lower", "VelocityX#A", X => -U);
                        C.AddBoundaryCondition("navierslip_linear_upper", "VelocityX#A", X => U);

                        if (!Xperiodic) {
                            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryCondition("pressure_outlet_right");
                        }
                        break;
                    }
                case 3: {
                        C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", X => U);
                        C.AddBoundaryCondition("freeslip_upper");
                        if (!Xperiodic) {
                            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryCondition("pressure_outlet_right");
                        }
                        break;
                    }
                case 4: {
                        C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", X => U);
                        C.AddBoundaryCondition("navierslip_linear_upper");
                        if (!Xperiodic) {
                            C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryCondition("pressure_outlet_right");
                        }
                        break;
                    }
                default: {
                        break;
                    }
            }

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.FilterConfiguration.FilterCurvatureCycles = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.CompMode = AppControl._CompMode.Steady;

            double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 100;
            C.saveperiod = 1;

            #endregion


            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CF_LevelSetMovementTest(int boundarySetup = 2, LevelSetEvolution lsEvo = LevelSetEvolution.FastMarching, 
            LevelSetHandling lsHandl = LevelSetHandling.Coupled_Once, XNSE_Control.TimesteppingScheme tsScheme = XNSE_Control.TimesteppingScheme.ImplicitEuler) {

            int p = 2;
            int kelem = 16;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null; //_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Two-phase Channel flow for testing the level set movement";

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
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 2;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);

                bool xPeriodic = (boundarySetup == 1) ? true : false;
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: xPeriodic);

                grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");

                switch (boundarySetup) {
                    case 1:
                        break;
                    case 2:
                        grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                        grd.EdgeTagNames.Add(4, "pressure_outlet_right");
                        break;
                    default:
                        throw new ArgumentException("invalid boundary setup");

                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (!xPeriodic) {
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 4;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc;

            switch (boundarySetup) {
                case 1: {
                        // horizontal interface
                        PhiFunc = (X => X[0] - L / 4.0);
                        break;
                    }
                case 2: {
                        // radial interface
                        double[] center = new double[] { L / 4.0, H / 2.0 };
                        double radius = 0.25;
                        PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius);
                        break;
                    }
                default:
                    PhiFunc = (X => -1);
                    break;
                
            }
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 1.0;

            switch (boundarySetup) {
                case 1:
                    C.InitialValues_Evaluators.Add("VelocityY#A", X => U);
                    C.InitialValues_Evaluators.Add("VelocityY#B", X => U);
                    break;
                case 2:
                    C.InitialValues_Evaluators.Add("VelocityX#A", X => U);
                    C.InitialValues_Evaluators.Add("VelocityX#B", X => U);
                    break;
                default:
                    throw new ArgumentException("invalid boundary setup");

            }


            #endregion


            // boundary conditions
            // ===================
            #region BC

            switch (boundarySetup) {
                case 1:
                    C.AddBoundaryCondition("velocity_inlet_lower", "VelocityY#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_lower", "VelocityY#B", X => U);
                    C.AddBoundaryCondition("velocity_inlet_upper", "VelocityY#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_upper", "VelocityY#B", X => U);
                    break;
                case 2:
                    C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_lower", "VelocityX#B", X => U);
                    C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_upper", "VelocityX#B", X => U);
                    C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#A", X => U);
                    C.AddBoundaryCondition("velocity_inlet_left", "VelocityX#B", X => U);
                    C.AddBoundaryCondition("pressure_outlet_right");
                    break;
                default:
                    break;
            }

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.NoOfMultigridLevels = 1;
            C.Solver_MaxIterations = 50;
            C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ContinuousDG;

            C.Option_LevelSetEvolution = lsEvo;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            #endregion


            // Timestepping
            // ============
            #region time

            C.CompMode = AppControl._CompMode.Transient;
            C.Timestepper_LevelSetHandling = lsHandl;
            C.Timestepper_Scheme = tsScheme;

            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10;
            C.saveperiod = 1;

            #endregion


            return C;
        }


    }
}
