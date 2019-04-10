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
using BoSSS.Solution.LevelSetTools.TestCases;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the channel flow type testcases
    /// </summary>
    public static class ChannelFlow {


        /// <summary>
        /// control object for various testing
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control ChannelFlow_WithInterface(int p = 2, int kelem = 8, int wallBC = 0) {

            XNSE_Control C = new XNSE_Control();

            string _DbPath = null; // @"D:\local\local_test_db";

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
            C.FieldOptions.Add("KineticEnergy", new FieldOpts() {
                Degree = 2*p,
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
            //C.PhysicalParameters.theta_e = Math.PI / 2.0;

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
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX:true);

                switch(wallBC) {
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
                //grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                ////grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                //grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    //if(Math.Abs(X[0]) <= 1.0e-8)
                    //    et = 3;
                    //if(Math.Abs(X[0] - L) <= 1.0e-8)
                    //    et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1); // X[0] - (H / 2.0)); // + (H/20)*Math.Cos(8 * Math.PI * X[0] / L));

            //double[] center = new double[] { H / 2.0, H / 2.0 };
            //double radius = 0.25;

            //C.InitialValues_Evaluators.Add("Phi",
            //    //(X => (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - radius.Pow2())   // quadratic form
            //    (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius)  // signed-distance form
            //    );

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 0.125;

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => (-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);
            //C.InitialValues_Evaluators.Add("VelocityX#B", X => (4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);

            //C.InitialValues_Evaluators.Add("Pressure#A", X => 2.0 - X[0]);

            //C.InitialValues_Evaluators.Add("KineticEnergy#A", X => 1.0 * ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U).Pow2() / 2.0);

            //double Pjump = sigma / radius;
            //C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);

            C.InitialValues_Evaluators.Add("GravityX#A", X => 5.0);

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
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#B", X => U);
                    break;
                case 2:
                    C.AddBoundaryValue("navierslip_linear_lower");
                    C.AddBoundaryValue("navierslip_linear_upper");
                    break;
                default:
                    C.AddBoundaryValue("wall_lower");
                    C.AddBoundaryValue("wall_upper");
                    break;

            }

            //C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => (- 4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);
            ////C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => U);
            //C.AddBoundaryValue("velocity_inlet_left", "KineticEnergy#A", X => 1.0 * ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U).Pow2() / 2.0);
            ////C.AddBoundaryValue("velocity_inlet_left", "KineticEnergy#B", X => U.Pow2() / 2); 
            ////C.AddBoundaryValue("pressure_outlet_left");
            //C.AddBoundaryValue("pressure_outlet_right");


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = true;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = false;
            C.RefinementLevel = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            C.Timestepper_Scheme = XNSE_Control.TimesteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.CompMode = AppControl._CompMode.Transient;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 300;
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
                        C.AddBoundaryValue("freeslip_lower");
                        if (slip)
                            C.AddBoundaryValue("freeslip_upper", "VelocityX#A", X => -U);
                        else
                            C.AddBoundaryValue("freeslip_upper");

                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                case 2: {
                        //C.AddBoundaryCondition("navierslip_linear_lower");
                        //if (slip)
                        //    C.AddBoundaryCondition("navierslip_linear_upper", "VelocityX#A", X => -U);
                        //else
                        //    C.AddBoundaryCondition("navierslip_linear_upper");

                        C.AddBoundaryValue("navierslip_linear_lower", "VelocityX#A", X => -U);
                        C.AddBoundaryValue("navierslip_linear_upper", "VelocityX#A", X => U);

                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                case 3: {
                        C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                        C.AddBoundaryValue("freeslip_upper");
                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
                        }
                        break;
                    }
                case 4: {
                        C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                        C.AddBoundaryValue("navierslip_linear_upper");
                        if (!Xperiodic) {
                            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                            C.AddBoundaryValue("pressure_outlet_right");
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
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations=50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
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
        public static XNSE_Control CF_LevelSetMovementTest(int boundarySetup = 2, double characteristicLength = 1.0, LevelSetEvolution lsEvo = LevelSetEvolution.FastMarching, 
            LevelSetHandling lsHandl = LevelSetHandling.Coupled_Once, XNSE_Control.TimesteppingScheme tsScheme = XNSE_Control.TimesteppingScheme.ImplicitEuler) {

            int p = 2;
            int kelem = 16;
            double cLength = characteristicLength;

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
                        grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                        grd.EdgeTagNames.Add(4, "pressure_outlet_right");
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
                        PhiFunc = (X => ((X[0]- cLength).Pow2()).Sqrt() - cLength/2);
                        break;
                    }
                case 2: {
                        // radial interface
                        double[] center = new double[] { L / 4.0, H / 2.0 };
                        double radius = cLength;
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
                    //C.InitialValues_Evaluators.Add("VelocityY#A", X => U);
                    //C.InitialValues_Evaluators.Add("VelocityY#B", X => U);
                    C.InitialValues_Evaluators.Add("VelocityX#A", X => U);
                    C.InitialValues_Evaluators.Add("VelocityX#B", X => U);
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
                    //C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", X => U);
                    //C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#B", X => U);
                    //C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#A", X => U);
                    //C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => U);
                    C.AddBoundaryValue("pressure_outlet_right");
                    break;
                case 2:
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                    C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => U);
                    C.AddBoundaryValue("pressure_outlet_right");
                    break;
                default:
                    break;
            }

            #endregion

            // advanced settings for Fourier-Level-Set
            // ======================
            #region Fourier level-set

            switch (lsEvo)
            {
                case LevelSetEvolution.Fourier:
                    {
                        switch (boundarySetup)
                        {
                            case 1:
                                {
                                    throw new ArgumentException("Fourier Level-Set not implemented in Line Movement Test");
                                    break;
                                }
                            case 2:
                                {
                                    int numSp = 640;
                                    double[] FourierP = new double[numSp];
                                    double[] samplP = new double[numSp];
                                    double[] center = new double[] { L / 4.0, H / 2.0 };
                                    double radius = cLength;
                                    for (int sp = 0; sp < numSp; sp++)
                                    {
                                        FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                                        samplP[sp] = radius;
                                    }

                                    C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem)
                                    {
                                        center = center,
                                        FourierEvolve = Fourier_Evolution.MaterialPoints,
                                        centerMove = CenterMovement.Reconstructed,
                                    };

                                    C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
                                    break;
                                }                               
                            default:
                                break;
                               
                        }
                        break;
                    }
                default:
                    break;
            }

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 4;
            C.LinearSolver.MinSolverIterations = 4;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
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

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control CF_LevelSetRotationTest(int boundarySetup = 1, double characteristicLength = 1.0, LevelSetEvolution lsEvo = LevelSetEvolution.FastMarching,
            LevelSetHandling lsHandl = LevelSetHandling.Coupled_Once, XNSE_Control.TimesteppingScheme tsScheme = XNSE_Control.TimesteppingScheme.ImplicitEuler)
        {

            int p = 2;
            int kelem = 16;
            double cLength = characteristicLength;

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null; //_DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Two-phase flow for testing the level set movement in solid body rotation";

            #endregion


            // DG degrees
            // ==========
            #region degrees

            C.FieldOptions.Add("VelocityX", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts()
            {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts()
            {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts()
            {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts()
            {
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

            double L = 1;
            double H = 1;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L/2, L/2, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-H/2, H/2, kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");                   
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "velocity_inlet_right");               

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + H/2) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H/2) <= 1.0e-8)
                        et = 2;                    
                    if (Math.Abs(X[0] + L/2) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L/2) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc;

            switch (boundarySetup)
            {
                case 1:
                    {
                        // elliptoid
                        double[] center = new double[] { 0.15, 0.0 };
                        double[] shape = new double[] { 1, 0.36 };
                        double radius = cLength;
                        PhiFunc = (X => ((X[0] - center[0]).Pow2() / shape[0] + (X[1] - center[1]).Pow2() / shape[1]).Sqrt() - radius);
                        break;
                    }
                case 2:
                    {
                        // slotted disk
                        double[] xCutout = new double[] { -0.1, 0.1 };
                        double yCutout = -0.1;
                        double radius = cLength;
                        ZalesaksDisk disk = new ZalesaksDisk(xCutout, yCutout, radius);
                        PhiFunc = (X => disk.SignedDistanceLevelSet(X));
                        break;
                    }
                default:
                    PhiFunc = (X => -1);
                    break;

            }
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("VelocityX#A", X => -X[1]);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => -X[1]);
            C.InitialValues_Evaluators.Add("VelocityY#A", X => X[0]);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => X[0]);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#A", X => -X[1]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#B", X => -X[1]);

            C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#B", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_upper", "VelocityY#B", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#B", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityY#A", X => X[0]);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityY#B", X => X[0]);

            #endregion

            // advanced settings for Fourier-Level-Set
            // ======================
            #region Fourier level-set

            switch (lsEvo)
            {
                case LevelSetEvolution.Fourier:
                    {
                        switch (boundarySetup)
                        {
                            case 1:
                                {
                                    int numSp = 640;
                                    double[] FourierP = new double[numSp];
                                    double[] samplP = new double[numSp];
                                    double[] center = new double[] { 0.15, 0.0 };
                                    double radius = cLength;
                                    for (int sp = 0; sp < numSp; sp++)
                                    {
                                        FourierP[sp] = sp * (2 * Math.PI / (double)numSp);
                                        samplP[sp] = radius / (Math.Cos(FourierP[sp]).Pow2() + Math.Sin(FourierP[sp]).Pow2() / 0.36).Sqrt();
                                    }

                                    C.FourierLevSetControl = new FourierLevSetControl(FourierType.Polar, 2 * Math.PI, FourierP, samplP, 1.0 / (double)kelem)
                                    {
                                        center = center,
                                        FourierEvolve = Fourier_Evolution.MaterialPoints,
                                        centerMove = CenterMovement.Reconstructed,
                                        PeriodicFunc = (X => radius / (Math.Cos(X).Pow2() + Math.Sin(X).Pow2() / 0.36).Sqrt())
                                    };

                                    C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.Curvature_Fourier;
                                    break;
                                }
                            case 2:
                                {
                                    throw new ArgumentException("Fourier Level-Set is not suitable for Slotted Disk");
                                }
                            default:
                                break;

                        }
                        break;
                    }
                default:
                    break;
            }

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergy = false;

            C.VelocityBlockPrecondMode = MultigridOperator.Mode.SymPart_DiagBlockEquilib;
            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 4;
            C.LinearSolver.MinSolverIterations = 4;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
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
