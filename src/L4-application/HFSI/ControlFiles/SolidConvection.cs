//using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSFE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Transactions;
using HFSISolver.SolidPhase;

namespace HFSISolver.ControlFiles {

    public static class SolidConvection {

        static class Constants {
            // careful order of declaration matters!!!
            // material parameters
            public static double rho_A = 200.0;
            public static double rho_B = 5.0;
            public static double rho_C = 400.0;

            public static double mu_A = 10.0;//1.0;//0.1;
            public static double mu_C = mu_A;

            public static double Sigma = 0.1;

            public static double c_A = 1.0;
            public static double c_C = 1.0;

            public static double k_A = 0.01;
            public static double k_C = 0.01;

            public static double hVap = 10000;
            public static double T_sat = 0; // 500K, shifted
            public static double T_wall = T_sat + 5.0;

            public static double alpha_A = k_A / (c_A * rho_A);
            //public static double alpha_B = k_B / (c_B * rho_B);
            public static double eps = 1.0 - rho_B / rho_A;

            public static double g = 9.81;

            public static double lambda = 2 * Math.PI * Math.Sqrt((3 * Sigma) / (g * Math.Abs(rho_A - rho_B)));
            public static double L = lambda / 2.0;


            public static double ST = 0.5 * L - 1e-4;//0.03;//0.01;
            public static double X0 = 0 * L;
            public static double Y0 = L * (1.0 / 4.0 + 1.0 / 16.0) + ST;
            public static double R = L / 4.0;


            // capillary timestep , for finest res + highest degree, use this, for comparability?!, Is very small 1e-7 => 1e5 - 1e6 timesteps necessary => artificial surface tension?!
            public static Func<int, int, double> dt = (res, p) => 0.5 * Math.Sqrt((rho_A + rho_B) * Math.Pow(L / ((double)res * ((double)p + 1)), 3) / (2 * Math.PI * Math.Abs(Sigma)));
        }

        public static class BoundaryAndInitialValueFactory {

            public static string GetPrefixCode() {
                using(var stw = new System.IO.StringWriter()) {

                    stw.WriteLine("static class BoundaryAndInitialValues {");
                    stw.WriteLine("     ");
                    stw.WriteLine("     public static double Phi(double[] X, double t){");
                    stw.WriteLine($"         return -(X[1] - 1.1 * {Constants.ST});");
                    //stw.WriteLine($"         return 1.0 * ((X[0]-{Constants.X0}).Pow2() + (X[1]-{Constants.Y0}).Pow2() - {Constants.R}.Pow2());");
                    stw.WriteLine("     }");
                    stw.WriteLine("     ");
                    stw.WriteLine("     public static double Phi1(double[] X, double t){");
                    stw.WriteLine($"         return -(X[1] - {Constants.ST});");
                    //stw.WriteLine($"         return (1.0);");
                    stw.WriteLine("     }");
                    stw.WriteLine("     ");
                    stw.WriteLine("     public static double TemperatureB(double[] X, double t){");
                    //stw.WriteLine("         return " + Constants.T_sat + ";");
                    stw.WriteLine("         return 1.0 * t ;");
                    stw.WriteLine("     }");
                    stw.WriteLine("     ");
                    stw.WriteLine("}");
                    return stw.ToString();
                }
            }

            static public Formula Get_TempB() {
                return new Formula("BoundaryAndInitialValues.TemperatureB", true, AdditionalPrefixCode: GetPrefixCode());
            }
            static public Formula Get_Phi() {
                return new Formula("BoundaryAndInitialValues.Phi", true, AdditionalPrefixCode: GetPrefixCode());
            }

            static public Formula Get_Phi1() {
                return new Formula("BoundaryAndInitialValues.Phi1", true, AdditionalPrefixCode: GetPrefixCode());
            }
        }

        public static HFSI_Control Convection(int p = 3, int AMRlvl = 0) {
            HFSI_Control C = new HFSI_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;

            //C.AgglomerationThreshold = 0.3;
            //C.NoOfMultigridLevels = 1;

            //int D = 3;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            //C.DbPath = @"C:\Users\miao\Documents\Database\Droplet-Leidenfrost";
            C.ProjectName = "Droplet";
            C.SessionName = "Droplet-check";
            C.ProjectDescription = "Droplet running on pc";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            //if(D == 3) {
            //    C.FieldOptions.Add("DisplacementZ", new FieldOpts() {
            //        Degree = p,
            //        SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            //    });
            //}

            // Physical Parameters
            // ===================
            #region physics

            // Physical Parameters
            // ===================
            C.solveCoupledHeatEquation = true;
            C.PhysicalParameters.Material = false;

            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;

            C.PhysicalParameters.rho_A = Constants.rho_A;

            C.PhysicalParameters.mu_A = Constants.mu_A;

            C.PhysicalParameters.Sigma = Constants.Sigma;

            C.ThermalParameters.rho_A = Constants.rho_A;
            C.ThermalParameters.rho_C = Constants.rho_C;

            C.ThermalParameters.c_A = Constants.c_A;
            C.ThermalParameters.c_C = Constants.c_C;

            C.ThermalParameters.k_A = Constants.k_A;
            C.ThermalParameters.k_C = Constants.k_C;

            C.ThermalParameters.hVap = Constants.hVap;
            C.ThermalParameters.T_sat = Constants.T_sat;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.sliplength = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.Material = new Solid() {
                Density = Constants.rho_C,
                Lame2 = 1,
                Viscosity = Constants.mu_C
            };
            #endregion


            // grid generation
            // ===============
            #region grid
            int res = 4;



            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-Constants.L, Constants.L, 2 * res + 1);
                double[] Ynodes = GenericBlas.Linspace(0, Constants.L, res + 1);
                /*
                double trsh = 0.05 * Constants.ST;
                double[] Ynodes = GenericBlas.Linspace(0, Constants.ST - 2 * trsh, res + 1);
                Ynodes = Ynodes.Cat(GenericBlas.Linspace(Constants.ST + trsh, 2 * Constants.L, 2 * res + 1));
                */
                //double trsh = 0.05 * Constants.R;
                //double[] Ynodes = GenericBlas.Linspace(0, Constants.ST - trsh, res + 1);
                //Ynodes = Ynodes.Cat(GenericBlas.Linspace(Constants.ST - trsh + trsh * 3 / 5, Constants.ST + 2 * trsh - trsh * 3 / 5, 3 + 1));
                //Ynodes = Ynodes.Cat(GenericBlas.Linspace(Constants.ST + 2 * trsh, 2 * Constants.L, res + 1));
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_top_ZeroGradient");
                grd.EdgeTagNames.Add(2, "wall_ConstantTemperature");
                //grd.EdgeTagNames.Add(3, "wall_ConstantTemperature");
                //grd.EdgeTagNames.Add(4, "wall_ConstantTemperature");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - Ynodes.Last()) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - Ynodes.First()) <= 1.0e-8)
                        et = 2;
                    //if(Math.Abs(X[0] - Xnodes.Last()) <= 1.0e-8)
                    //    et = 2;
                    //if(Math.Abs(X[0] - Xnodes.First()) <= 1.0e-8)
                    //    et = 2;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            //C.AddInitialValue("Phi", BoundaryAndInitialValueFactory.Get_Phi());
            C.AddInitialValue("Phi2", BoundaryAndInitialValueFactory.Get_Phi1());
            //C.AddInitialValue("Temperature#A", "X => " + Constants.T_sat + 1.0, false);
            //C.AddInitialValue("Temperature#C", "X => " + Constants.T_sat + 3.0, false);
            //C.AddInitialValue("Temperature#A", $"X => + Math.Sin(Math.PI * X[0] / {Constants.L})", false);
            //C.AddInitialValue("Temperature#C", $"X => + Math.Sin(Math.PI * X[0] / {Constants.L})", false);
            //C.AddInitialValue("Temperature#B", BoundaryAndInitialValueFactory.Get_TempB());
            //C.AddInitialValue("Temperature#C", $"X => + {Constants.T_sat} + ({Constants.T_wall}-{Constants.T_sat}) * (1.0- X[1]/{Constants.ST})", false);
            //C.AddInitialValue("GravityY#A", new Formula($"X => -{Constants.g}", false));
            //C.AddInitialValue("GravityY#B", new Formula($"X => -{Constants.g}", false));
            //C.AddInitialValue("GravityY#C", new Formula($"X => -{Constants.g}", false));

            #endregion


            // boundary conditions
            // ===================
            #region BC

            Func<double[], double, double> Get_Temp = delegate (double[] X, double t) {
                return Math.Sin(Math.PI * X[0] / Constants.L);
            };
            Func<double[], double, double> Get_Temp1 = delegate (double[] X, double t) {
                //return 5.0 * Math.Sin(Math.PI * (X[0] - 0.01 * t) / Constants.L);
                return 2.5 * Math.Sin(Math.PI * (5.0 * t));
            };
            //if(!periodic) C.AddBoundaryValue("freeslip_ZeroGradient");
            //C.AddBoundaryValue("wall_top_ZeroGradient", "Temperature#A", Get_Temp1);
            //C.AddBoundaryValue("wall_top_ZeroGradient", "Temperature#B", Get_Temp1);
            //C.AddBoundaryValue("wall_top_ZeroGradient", "Temperature#C", Get_Temp1);
            //C.AddBoundaryValue("wall_top_ZeroGradient", "VelocityX#A", "X => " + 1.0, false);
            //C.AddBoundaryValue("wall_top_ZeroGradient", "VelocityY#A", "X => " + 0.0, false);
            C.AddBoundaryValue("wall_top_ZeroGradient", "VelocityX#A", Get_Temp1);
            C.AddBoundaryValue("wall_top_ZeroGradient", "VelocityY#A", "X => " + 0.0, false);
            //C.AddBoundaryValue("wall_top_ConstantTemperature", "DisplacementX#C", Get_DisX);
            //C.AddBoundaryValue("wall_top_ConstantTemperature", "DisplacementY#C", "X => " + 0.0, false);
            C.AddBoundaryValue("wall_ConstantTemperature", "Temperature#A", Get_Temp);
            C.AddBoundaryValue("wall_ConstantTemperature", "Temperature#B", Get_Temp);
            C.AddBoundaryValue("wall_ConstantTemperature", "Temperature#C", Get_Temp);
            //C.AddBoundaryValue("wall_ConstantTemperature", "Temperature#C", "X => " + Constants.T_wall, false);

            /*
            C.AddBoundaryValue("wall_ConstantHeatFlux_lower");
            C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper");
            C.AddBoundaryValue("wall_ConstantHeatFlux_left");
            C.AddBoundaryValue("wall_ConstantHeatFlux_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = 0.001;
            */
            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 1;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 15;
            C.NonLinearSolver.MinSolverIterations = 3;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.FailOnSolverFail = false;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;
            C.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;


            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            //C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = false;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = AMRlvl });
            //C.AMR_startUpSweeps = AMRlvl;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 0, maxRefinementLevel = AMRlvl });
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 1, maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;

            //C.ReInitPeriod = 5;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;

            //C.ReInitPeriod = 5;


            C.TimesteppingMode = compMode;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtFixed = 1.0 * Constants.dt((int)(res * Math.Pow(2, AMRlvl)), p); // Use same timestep for all, better comparability
            //C.dtFixed = 1;
            C.Endtime = 1.6;//0.5;
            C.NoOfTimesteps = (int)Math.Ceiling(C.Endtime / C.dtFixed);
            C.saveperiod = 1;
            C.rollingSaves = false;
            /*
            double dt = 5e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            //C.Endtime = 10;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;
            */
            #endregion

            return C;
        }

    }
}
