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

//using MathNet.Numerics;

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
using BoSSS.Foundation.XDG;
using MathNet.Numerics;
using BoSSS.Solution.LevelSetTools;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    // ====================================
    // utility functions for initial values
    // ====================================
    #region Init-Utils

    /// <summary>
    /// implementation of the error function acfcording to numerical recipes third edition
    /// </summary>
    public static class ErrorFuncUtil {

        static readonly double[] cof = new double[] { -1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2, -9.561514786808631e-3, -9.46595344482036e-4,
                                                        3.66839497852761e-4, 4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8,
                                                        -8.5238095915e-8, 6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11, 2.394038e-12,
                                                        -6.886027e-12, 8.94487e-13, 3.13092e-13, -1.12708e-13, 3.81e-16, 7.106e-15, -1.523e-15, -9.4e-17, 1.21e-16, -2.8e-17};

        public static double erf(double x) {
            if (x >= 0.0)
                return 1.0 - erfccheb(x);
            else
                return erfccheb(-x) - 1.0;
        }

        public static double erfc(double x) {
            if (x >= 0.0)
                return erfccheb(x);
            else
                return 2.0 - erfccheb(-x);
        }

        static double erfccheb(double z) {

            if (z < 0.0)
                throw new ArgumentOutOfRangeException("erfccheb requires nonegative argument");

            double t = 2.0 / (2.0 + z);
            double ty = 4.0 * t - 2.0;

            double tmp;
            double d = 0.0;
            double dd = 0.0;
            for (int j = cof.Length - 1; j > 0; j--) {
                tmp = d;
                d = ty * d - dd + cof[j];
                dd = tmp;
            }

            return t * Math.Exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
        }

        public static double inverf(double p) {
            return inverfc(1.0 - p);
        }

        static double inverfc(double p) {

            if (p >= 2.0)
                return -100.0;
            if (p <= 2.0)
                return 100.0;

            double pp = (p < 1.0) ? p : 2.0 - p;
            double t = Math.Sqrt(-2.0*Math.Log(pp/2.0));
            double x = -0.70711 * ((2.30753 + t * 0.27061) / (1.0+t*(0.99229+t*0.04481)) - t);

            double err;
            for (int j = 0; j < 2; j++) {
                err = erfc(x) - pp;
                x += err / (1.12837916709551257 * Math.Exp(-Math.Sqrt(x)) - x * err);
            }

            return (p < 1.0) ? x : -x;
        }

    }

    #endregion




    /// <summary>
    /// class providing Controls for the Capillary rise testcases
    /// </summary>
    public static class HeatedWall {

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control HeatedWall_StartUp(int p = 2, int kelemR = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            bool solveHeat = false;
            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; // C.DbPath != null;
            C.ProjectName = "XNSE/HeatedWall";
            C.ProjectDescription = "Leikonfiguration for SFB 1194";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.LogPeriod = 100;
            C.PostprocessingModules.Add(new MovingContactLineLogging() { LogPeriod = 100 });

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 0.1;
            C.PhysicalParameters.Sigma = 1.0;

            C.solveCoupledHeatEquation = solveHeat;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 1.0;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 1.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 3.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 4.0;
            double H = 12.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, R, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 3 * kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                if(solveHeat) {
                    grd.EdgeTagNames.Add(1, "wall_Dirichlet_lower");
                    grd.EdgeTagNames.Add(2, "pressure_outlet_ZeroGradient_upper");
                    grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_ZeroGradient_right");
                } else {
                    grd.EdgeTagNames.Add(1, "wall_lower");
                    grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                    grd.EdgeTagNames.Add(3, "slipsymmetry_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
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

            double h0 = 3.0;

            Func<double[], double> PhiFunc = (X => X[1] - h0);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double g = 1.5;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            //C.InitialValues_Evaluators.Add("Temperature#B", X => 10);

            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            if(solveHeat) {
                C.AddBoundaryValue("wall_Dirichlet_lower", "Temperature#A", (X, t) => 0.0);
                C.AddBoundaryValue("wall_Dirichlet_lower", "Temperature#B", (X, t) => 0.0);
                C.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");

                C.AddBoundaryValue("slipsymmetry_ZeroGradient_left");
                C.AddBoundaryValue("navierslip_linear_ZeroGradient_right");
            } else {

                C.AddBoundaryValue("wall_lower");
                C.AddBoundaryValue("pressure_outlet_upper");

                C.AddBoundaryValue("slipsymmetry_left");
                C.AddBoundaryValue("navierslip_linear_right");
            }

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            C.PhysicalParameters.sliplength = 0.1;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.RefinementLevel = 1;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-1;
            C.dtMin = 1e-1;
            C.Endtime = 10000;
            C.NoOfTimesteps = 100;
            C.saveperiod = 1;

            #endregion


            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control HeatedWall_Run(int p = 2, int kelemR = 1, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            bool solveHeat = false;
            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";
            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/HeatedWall";
            C.ProjectDescription = "Leikonfiguration for SFB 1194";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.LogPeriod = 100;

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // numerical values for various testing (A:liquid and B:vapor state)
            //double rho_l = 1.0;
            //C.PhysicalParameters.rho_A = rho_l;
            //C.PhysicalParameters.rho_B = 1.0e-1;
            //double mu_l = 0.5;
            //C.PhysicalParameters.mu_A = mu_l;
            //C.PhysicalParameters.mu_B = 0.25e-1;
            //C.PhysicalParameters.Sigma = 7.5e-1;

            //double g = 29.81;

            //C.solveCoupledHeatEquation = solveHeat;
            //C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            //C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            //double c_l = 1.0e+4;
            //C.ThermalParameters.c_A = c_l;
            //C.ThermalParameters.c_B = 1.0e+3;
            //double k_l = 1.0e-1;
            //C.ThermalParameters.k_A = k_l;
            //C.ThermalParameters.k_B = 1e-2;

            //if (C.solveCoupledHeatEquation) {
            //    C.ThermalParameters.hVap = 1.0e6;
            //    //C.ThermalParameters.hVap_B = -1.0e6; 
            //}
            //double Tsat = 329.75;
            //C.ThermalParameters.T_sat = Tsat;

            ////C.ThermalParameters.pc = 0.0;
            ////C.ThermalParameters.fc = 1.0;
            ////C.ThermalParameters.Rc = 1e7;

            //double Lslip = 1e-2;
            //C.PhysicalParameters.betaL = 1.1 / (2 * Lslip);
            //C.PhysicalParameters.theta_e = Math.PI * (1.0 / 3.0);

            //double L = 1.0;


            // FC-72 (A:liquid and B:vapor state)
            double rho_l = 1619.82;  
            C.PhysicalParameters.rho_A = rho_l;    
            double rho_v = 13.36; 
            C.PhysicalParameters.rho_B = rho_v;  
            double mu_l = 4.5306e-4;      
            C.PhysicalParameters.mu_A = mu_l;
            double mu_v = 9.4602e-6;      
            C.PhysicalParameters.mu_B = mu_v;

            C.PhysicalParameters.Sigma = 8.273e-3;
            double g = 9.81;

            C.solveCoupledHeatEquation = solveHeat;
            C.ThermalParameters.rho_A = rho_l;
            C.ThermalParameters.rho_B = rho_v;
            double c_l = 1098.41;   
            C.ThermalParameters.c_A = c_l;
            double c_v = 885.04;    
            C.ThermalParameters.c_B = c_v;
            double k_l = 5.216e-2;
            C.ThermalParameters.k_A = k_l;
            double k_v = 8.64e-3;
            C.ThermalParameters.k_B = k_v;

            if (C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = 8.4515e+8;
                //C.ThermalParameters.hVap_B = -8.4515e+8;
            }

            double Tsat = 329.75;   
            C.ThermalParameters.T_sat = Tsat;
            C.ThermalParameters.p_sat = 1000;   // 1bar
            C.ThermalParameters.Rc = 24.45;
            C.ThermalParameters.fc = 0.5;

            //double A = 4.37e-17; //4.37e-17;    // dispersion constant

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI * (5.0 / 36.0);
            double Lslip = 1e-4;

            double L = 1e-2;    //cm


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = !solveHeat;

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 1.0 * L;
            double H = 5.0 * L;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, R, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, (5*kelemR) + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                if(solveHeat) {
                    grd.EdgeTagNames.Add(1, "wall_ZeroGradient_lower");
                    grd.EdgeTagNames.Add(2, "pressure_outlet_ZeroGradient_upper");
                    grd.EdgeTagNames.Add(3, "slipsymmetry_ConstantTemperature_left");
                    //grd.EdgeTagNames.Add(4, "navierslip_linear_ConstantTemperature_right");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_ConstantHeatFlux_right");
                } else {
                    grd.EdgeTagNames.Add(1, "wall_lower");
                    //grd.EdgeTagNames.Add(1, "pressure_outlet_lower");
                    grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                    grd.EdgeTagNames.Add(3, "slipsymmetry_left");
                    //grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                    grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - R) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = true;
            C.BaseRefinementLevel = 4;
            C.RefinementLevel = 4;
            C.AMR_startUpSweeps = 4;

            double hmicro = (R / kelemR) / 2.0;


            #endregion


            // Initial Values
            // ==============
            #region init

            double h0 = 3.2*L;

            Func<double[], double> PhiFunc = (X => X[1] - h0);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            // disjoining pressure field
            //double A = 1e-2;
            //double delta = 1e-2;
            //C.DisjoiningPressureFunc = (X => A / (Math.Abs(X[0] - (R+delta))).Pow(3));

            double U = 5e-3;
            double HeatFluxA = 500;
            double HeatFluxB = 0.0; // HeatFluxA * (k_v / k_l);

            if (C.solveCoupledHeatEquation) {

                double Twall = Tsat + (HeatFluxA / k_l) * hmicro;

                Func<double[], double> TempAFunc = delegate (double[] X) {

                    double Temp = Tsat;

                    if (X[0] > (R-hmicro)) {
                        Temp += (HeatFluxA / k_l) * (hmicro + (X[0] - R));
                    }

                    return Temp;
                };
                C.InitialValues_Evaluators.Add("Temperature#A", TempAFunc);

                Func<double[], double> TempBFunc = delegate (double[] X) {

                    //double Twall = Tsat + (HeatFluxB / k_v) * hmicro;
                    double Temp = Tsat;

                    //if (X[0] > (R - hmicro)) {
                    //    Temp += (HeatFluxB / k_v) * (hmicro + (X[0] - R));
                    //}

                    return Temp;
                };
                C.InitialValues_Evaluators.Add("Temperature#B", TempBFunc);
            }

            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC



            if(solveHeat) {
                C.AddBoundaryValue("wall_ZeroGradient_lower");
                C.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");

                C.AddBoundaryValue("slipsymmetry_ConstantTemperature_left", "Temperature#A", (X, t) => Tsat);
                C.AddBoundaryValue("slipsymmetry_ConstantTemperature_left", "Temperature#B", (X, t) => Tsat);

                //C.AddBoundaryValue("navierslip_linear_ConstantTemperature_right", "VelocityY#A", (X, t) => U);
                //C.AddBoundaryValue("navierslip_linear_ConstantTemperature_right", "VelocityY#B", (X, t) => U);
                //C.AddBoundaryValue("navierslip_linear_ConstantTemperature_right", "Temperature#A", (X, t) => Twall);
                //C.AddBoundaryValue("navierslip_linear_ConstantTemperature_right", "Temperature#B", (X, t) => Tsat);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#A", (X, t) => U);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#B", (X, t) => U);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "HeatFluxX#A", (X, t) => HeatFluxA);
                //C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "HeatFluxX#B", (X, t) => HeatFluxB);
            } else {

                C.AddBoundaryValue("wall_lower");
                //C.AddBoundaryValue("pressure_outlet_lower");
                C.AddBoundaryValue("pressure_outlet_upper");

                C.AddBoundaryValue("slipsymmetry_left");
                //C.AddBoundaryValue("pressure_outlet_left");

                C.AddBoundaryValue("navierslip_linear_right", "VelocityY#A", (X, t) => U);
                C.AddBoundaryValue("navierslip_linear_right", "VelocityY#B", (X, t) => U);
            }

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            C.PhysicalParameters.sliplength = Lslip;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;


            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.CurvatureNeeded = false;


            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;


            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-4;
            C.dtMin = 1e-4;
            C.Endtime = 1;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 100;

            #endregion


            return C;
        }


        public static XNSE_Control HeatedWall_Run2(int p = 1, int kelemR = 8, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            bool Water = true;
            bool FC72 = false;
            bool contactAngle90 = true;

            bool startUp_Interface = true;
            bool pressureCondition = false;
            bool initmeniscus = false; 
            bool startUp_Heat = false;
            bool startUp_withEvap = false;
            bool run = false;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";
            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";

            //Guid restartID = new Guid("53b135e3-366e-48a0-adb8-6276423f9934");
            //TimestepNumber RestartTime = new TimestepNumber("1000");



            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            if(startUp_Interface)
                C.ProjectName = "HeatWall_startUp_Interface";
            if (startUp_Heat)
                C.ProjectName = "HeatWall_startUp_Heat";
            if (run)
                C.ProjectName = "HeatWall_run";
            C.ProjectDescription = "Leikonfiguration for SFB 1194";

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
                Degree = Math.Max(p,2),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 2*p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            double g = 9.81;
            double Tsat = 0.0;
            double tscale_phys = 0.0;

            if (Water) {

                double mu_scl = 1;

                //physical values for water(A: liquid and B: vapor state)
                double rho_l = 997; // kg/m3     //9.97e-7;   // kg / mm3
                C.PhysicalParameters.rho_A = rho_l;
                double rho_v = 2.31e-2; //kg/m3      //2.31e-11;  // kg / mm3
                C.PhysicalParameters.rho_B = rho_v;
                double mu_l = mu_scl * 0.891e-3; // kg/m*sec    //8.91e-7;    // kg / mm*sec
                C.PhysicalParameters.mu_A = mu_l;
                double mu_v = mu_scl * 0.987e-5; //kg/m*sec      //9.87e-9;    // kg / mm*sec
                C.PhysicalParameters.mu_B = mu_v;
                double sigma = 71.97e-3; //N/m    //7.178;     // kg / sec^2
                C.PhysicalParameters.Sigma = 0.0; // sigma;

                //C.AdvancedDiscretizationOptions.LFFB = 1.0;

                C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
                C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
                double c_l = 4180; //J/kgK   //4.18e+9;    // mm^2 / sec^2 * K
                C.ThermalParameters.c_A = c_l;
                double c_v = 1870; //J/kgK    //1.87e+9;    // mm^2 / sec^2 * K
                C.ThermalParameters.c_B = c_v;
                double k_l = 0.607; //W/mK   //607; // kg / sec^3 * K
                C.ThermalParameters.k_A = k_l;
                double k_v = 0.0186; //W/mK    //18.6;       // kg / sec^3 * K
                C.ThermalParameters.k_B = k_v;

                double hVap = 2442e3; //2.442e12; // mm^2 * kg / sec^2
                C.ThermalParameters.hVap = hVap;

                Tsat = 298.15;
                C.ThermalParameters.T_sat = Tsat;

                C.ThermalParameters.p_sat = 31.68e2;//N/m2  //3.168;   // kg / sec^2 * mm
                //C.ThermalParameters.Rc = 461.52; //J/kgK   //4.6152e8;   // mm^2 / sec * K    
                //C.ThermalParameters.fc = 0.5;
                //double A = 4.37e-17; //4.37e-17;    // dispersion constant

                C.PhysicalParameters.theta_e = (contactAngle90) ? Math.PI / 2.0 : Math.PI * (2.0 / 9.0);    // 40^


                tscale_phys = 47;
            }

            if (FC72) {

                // FC-72 (A:liquid and B:vapor state)
                double rho_l = 1619.82; 
                C.PhysicalParameters.rho_A = rho_l;
                double rho_v = 13.36;  
                C.PhysicalParameters.rho_B = rho_v;
                double mu_l = 4.5306e-4; 
                C.PhysicalParameters.mu_A = mu_l;
                double mu_v = 9.4602e-6;
                C.PhysicalParameters.mu_B = mu_v;
                double sigma = 8.273e-3; 
                C.PhysicalParameters.Sigma = sigma;


                C.ThermalParameters.rho_A = rho_l;
                C.ThermalParameters.rho_B = rho_v;
                double c_l = 1098.41; 
                C.ThermalParameters.c_A = c_l;
                double c_v = 885.04; 
                C.ThermalParameters.c_B = c_v;
                double k_l = 5.216e-2;
                C.ThermalParameters.k_A = k_l;
                double k_v = 0.864e-2;
                C.ThermalParameters.k_B = k_v;

                double hvap = 84515;
                C.ThermalParameters.hVap = hvap;

                Tsat = 329.75;
                C.ThermalParameters.T_sat = Tsat;
                C.ThermalParameters.p_sat = 1000;   // 1bar
                //C.ThermalParameters.Rc = 2.445e+5;
                //C.ThermalParameters.fc = 0.5;

                //double A = 4.37e-17; //4.37e-17;    // dispersion constant

                C.PhysicalParameters.theta_e = (contactAngle90) ? Math.PI / 2.0 : Math.PI * (5.0 / 36.0);   // 25^


                tscale_phys = 177;
            }
       
             
            C.solveCoupledHeatEquation = true;
            //C.useSolutionParamUpdate = true;

            if (startUp_Interface || (startUp_Heat && !startUp_withEvap)) {
                C.ThermalParameters.hVap = 0.0;
                C.PhysicalParameters.Material = true;
                C.PhysicalParameters.IncludeConvection = false;
                C.ThermalParameters.IncludeConvection = false;
            } else {
                C.PhysicalParameters.Material = false;
                C.PhysicalParameters.IncludeConvection = true;
                C.ThermalParameters.IncludeConvection = true;
            }

            

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 30e-3;
            double H = 2 * R;

            if (startUp_Interface || startUp_Heat) {

                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, R, 2 + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, H, (2 * 2) + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    if (pressureCondition) {
                        grd.EdgeTagNames.Add(1, "pressure_Dirichlet_ZeroGradient_lower");
                        grd.EdgeTagNames.Add(3, "pressure_Dirichlet_ConstantTemperature_left");
                    } else {
                        grd.EdgeTagNames.Add(1, "wall_ZeroGradient_lower");
                        grd.EdgeTagNames.Add(3, "freeslip_ConstantTemperature_left");
                    }

                    grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ZeroGradient_upper");
                    //grd.EdgeTagNames.Add(4, "navierslip_linear_ConstantHeatFlux_right");
                    grd.EdgeTagNames.Add(4, "freeslip_ConstantHeatFlux_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - H) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - R) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            if (run) {
                throw new NotImplementedException("LoadSession");
            }


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.CurvatureRefined;
            C.RefineNavierSlipBoundary = !startUp_Interface;
            C.BaseRefinementLevel = 2;
            C.RefinementLevel = 2;
            C.AMR_startUpSweeps = 2;
            //C.ReInitPeriod = 1;

            double hmin = (R / Math.Pow(2.0, C.RefinementLevel + 1));
            double tscale_grid = 2e-6;

            #endregion


            // Initial Values
            // ==============
            #region init

            if (startUp_Interface || startUp_Heat) {

                double h0 = 0.4 * H;
                Func<double[], double> PhiFunc;
                if (initmeniscus) {
                    double rho_l = C.PhysicalParameters.rho_A;
                    double sigma = C.PhysicalParameters.Sigma;
                    double thetaE = C.PhysicalParameters.theta_e;
                    Func<double, double> meniscus = x => Math.Sqrt(sigma / (rho_l * g)) * 1 / Math.Tan(thetaE) * Math.Exp(-Math.Pow(((rho_l * g) / sigma), 0.5) * x);

                    //double a = Math.Sqrt(sigma / (rho_l * g));
                    //Console.WriteLine("a = {0}", a);
                    //double x0 = Math.Sqrt(2.0) * a;
                    //Console.WriteLine("x0 = {0}", x0);
                    //Func<double, double> meniscus = x => a * (Math.Log((2.0 * a / x) + Math.Sqrt((2.0 * a / x).Pow2() - 1))
                    //- Math.Log((2.0 * a / x0) + Math.Sqrt((2.0 * a / x0).Pow2() - 1)) + Math.Sqrt(4.0 - (x0 / a).Pow2()) - Math.Sqrt(4.0 - (x / a).Pow2()));

                    PhiFunc = (X => X[1] - (h0 + meniscus(R - X[0])));
                } else {
                    PhiFunc = (X => X[1] - h0);
                }
                C.InitialValues_Evaluators.Add("Phi", PhiFunc);

                C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
                C.InitialValues_Evaluators.Add("GravityY#B", X => -g);
                C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);
                C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat);
                C.InitialValues_Evaluators.Add("Pressure#A", X => C.ThermalParameters.p_sat); // -C.PhysicalParameters.rho_A * g * (X[1] - h0) + C.ThermalParameters.p_sat);
                C.InitialValues_Evaluators.Add("Pressure#B", X => C.ThermalParameters.p_sat);
            } else {
                //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);
            }

            #endregion


            // boundary conditions
            // ===================
            #region BC

            if (startUp_Interface) {
                if (pressureCondition) {
                    C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_lower", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Temperature#A", (X, t) => Tsat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Temperature#B", (X, t) => Tsat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                } else {
                    C.AddBoundaryValue("wall_ZeroGradient_lower");
                    C.AddBoundaryValue("freeslip_ConstantTemperature_left", "Temperature#A", (X, t) => Tsat);
                    C.AddBoundaryValue("freeslip_ConstantTemperature_left", "Temperature#B", (X, t) => Tsat);
                }

                C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                //C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#A", (X, t) => 0.0);
                //C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#B", (X, t) => 0.0);
                //C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "HeatFluxX#A", (X, t) => 0.0);
                C.AddBoundaryValue("freeslip_ConstantHeatFlux_right");
            }

            if (startUp_Heat) {
                double HeatFluxA = 0.0;
                if (Water) {
                    HeatFluxA = 4500.0;
                } else if (FC72) {
                    HeatFluxA = 3000.0;
                }
                if (pressureCondition) {
                    C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_lower", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Temperature#A", (X, t) => Tsat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Temperature#B", (X, t) => Tsat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_left", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                } else {
                    C.AddBoundaryValue("wall_ZeroGradient_lower");
                    C.AddBoundaryValue("slipsymmetry_ConstantTemperature_left", "Temperature#A", (X, t) => Tsat);
                    C.AddBoundaryValue("slipsymmetry_ConstantTemperature_left", "Temperature#B", (X, t) => Tsat);
                }

                C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#A", (X, t) => 0.0);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#B", (X, t) => 0.0);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "HeatFluxX#A", (X, t) => HeatFluxA);
            }


            double Uwall = 0.1;

            if (run) {
                double HeatFluxA = 0.0;
                if (Water) {
                    HeatFluxA = 4500;
                } else if (FC72) {
                    HeatFluxA = 3000.0;
                }

                double rampUpTime = 1e-4;

                C.AddBoundaryValue("wall_ZeroGradient_lower");
                C.AddBoundaryValue("pressure_outlet_ZeroGradient_upper", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                C.AddBoundaryValue("slipsymmetry_ConstantTemperature_left", "Temperature#A", (X, t) => Tsat);
                C.AddBoundaryValue("slipsymmetry_ConstantTemperature_left", "Temperature#B", (X, t) => Tsat);
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#A", (X, t) => Uwall * Math.Tanh((Math.PI / rampUpTime) * t));
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "VelocityY#B", (X, t) => Uwall * Math.Tanh((Math.PI / rampUpTime) * t));
                C.AddBoundaryValue("navierslip_linear_ConstantHeatFlux_right", "HeatFluxX#A", (X, t) => HeatFluxA);

            }


            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = hmin;

            double beta_scl = 100;
            C.PhysicalParameters.betaS_A = 0.0; // beta_scl * C.PhysicalParameters.mu_A / (hmin / 10.0);
            C.PhysicalParameters.betaS_B = 0.0; // beta_scl * C.PhysicalParameters.mu_B / (hmin / 1000.0);
            //C.PhysicalParameters.betaL = 0.01;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.AgglomerationThreshold = 0.1;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            //C.AdvancedDiscretizationOptions.freeSurfaceFlow = true;
            //C.PhysicalParameters.pFree = C.ThermalParameters.p_sat;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.CurvatureNeeded = false;

            //C.InterVelocAverage = XNSE_Control.InterfaceVelocityAveraging.phaseA;

            #endregion


            // level-set
            // =========
            #region levelset

            if (startUp_Interface || run || (startUp_Heat && startUp_withEvap)) {
                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;
                C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            } else {
                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            }


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-4;
            C.dtMin = 1e-4;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 100000;
            C.saveperiod = 10;

            #endregion


            C.ClearVelocitiesOnRestart = false;


            // parameters
            // ==========

            double[] param = new double[1];
            param[0] = Uwall;   // wall velocity

            C.AdditionalParameters = param;


            return C;
        }


        public static XNSE_Control ColdWall(int p = 1, int kelemR = 8, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            bool Water = true;
            bool FC72 = false;
            bool contactAngle90 = true;

            bool startUp_Interface = true;
            bool pressureCondition = false;
            bool initmeniscus = false;
            bool run = false;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";
            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";

            //Guid restartID = new Guid("53b135e3-366e-48a0-adb8-6276423f9934");
            //TimestepNumber RestartTime = new TimestepNumber("1000");



            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            if (startUp_Interface)
                C.ProjectName = "COldWall_startUp_Interface";
            if (run)
                C.ProjectName = "ColdWall_run";
            C.ProjectDescription = "Leikonfiguration for SFB 1194";

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
                Degree = Math.Max(p, 2),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = 2 * p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });


            #endregion


            // Physical Parameters
            // ===================
            #region physics

            double g = 9.81;
            double Tsat = 0.0;
            double tscale_phys = 0.0;

            if (Water) {

                double mu_scl = 1;

                //physical values for water(A: liquid and B: vapor state)
                double rho_l = 997; // kg/m3     //9.97e-7;   // kg / mm3
                C.PhysicalParameters.rho_A = rho_l;
                double rho_v = 2.31e-2; //kg/m3      //2.31e-11;  // kg / mm3
                C.PhysicalParameters.rho_B = rho_v;
                double mu_l = mu_scl * 0.891e-3; // kg/m*sec    //8.91e-7;    // kg / mm*sec
                C.PhysicalParameters.mu_A = mu_l;
                double mu_v = mu_scl * 0.987e-5; //kg/m*sec      //9.87e-9;    // kg / mm*sec
                C.PhysicalParameters.mu_B = mu_v;
                double sigma = 71.97e-3; //N/m    //7.178;     // kg / sec^2
                C.PhysicalParameters.Sigma = 0.0; // sigma;

                //C.AdvancedDiscretizationOptions.LFFB = 1.0;

                C.PhysicalParameters.theta_e = (contactAngle90) ? Math.PI / 2.0 : Math.PI * (2.0 / 9.0);   


                tscale_phys = 47;
            }

            if (FC72) {

                // FC-72 (A:liquid and B:vapor state)
                double rho_l = 1619.82;
                C.PhysicalParameters.rho_A = rho_l;
                double rho_v = 13.36;
                C.PhysicalParameters.rho_B = rho_v;
                double mu_l = 4.5306e-4;
                C.PhysicalParameters.mu_A = mu_l;
                double mu_v = 9.4602e-6;
                C.PhysicalParameters.mu_B = mu_v;
                double sigma = 8.273e-3;
                C.PhysicalParameters.Sigma = sigma;


                C.PhysicalParameters.theta_e = (contactAngle90) ? Math.PI / 2.0 : Math.PI * (5.0 / 36.0);   // 25^


                tscale_phys = 177;
            }


            if (startUp_Interface) {
                C.ThermalParameters.hVap = 0.0;
                C.PhysicalParameters.Material = true;
                C.PhysicalParameters.IncludeConvection = false;
                C.ThermalParameters.IncludeConvection = false;
            } else {
                C.PhysicalParameters.Material = false;
                C.PhysicalParameters.IncludeConvection = true;
                C.ThermalParameters.IncludeConvection = true;
            }



            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 30e-3;
            double H = 2 * R;
            int cells = 32;

            if (startUp_Interface) {

                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, R, cells + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, H, (2 * cells) + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    if (pressureCondition) {
                        grd.EdgeTagNames.Add(1, "pressure_Dirichlet_lower");
                        grd.EdgeTagNames.Add(3, "pressure_Dirichlet_left");
                    } else {
                        grd.EdgeTagNames.Add(1, "wall_lower");
                        grd.EdgeTagNames.Add(3, "freeslip_left");
                    }

                    grd.EdgeTagNames.Add(2, "pressure_Dirichlet_upper");
                    //grd.EdgeTagNames.Add(4, "navierslip_linear_right");
                    grd.EdgeTagNames.Add(4, "freeslip_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - H) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 3;
                        if (Math.Abs(X[0] - R) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            if (run) {
                throw new NotImplementedException("LoadSession");
            }


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.CurvatureRefined;
            C.RefineNavierSlipBoundary = !startUp_Interface;
            C.BaseRefinementLevel = 2;
            C.RefinementLevel = 2;
            C.AMR_startUpSweeps = 2;
            //C.ReInitPeriod = 1;

            double hmin = (R / Math.Pow(2.0, C.RefinementLevel + 1));
            double tscale_grid = 2e-6;

            #endregion


            // Initial Values
            // ==============
            #region init

            if (startUp_Interface) {

                double h0 = 0.4 * H;
                Func<double[], double> PhiFunc;
                if (initmeniscus) {
                    double rho_l = C.PhysicalParameters.rho_A;
                    double sigma = C.PhysicalParameters.Sigma;
                    double thetaE = C.PhysicalParameters.theta_e;
                    Func<double, double> meniscus = x => Math.Sqrt(sigma / (rho_l * g)) * 1 / Math.Tan(thetaE) * Math.Exp(-Math.Pow(((rho_l * g) / sigma), 0.5) * x);

                    //double a = Math.Sqrt(sigma / (rho_l * g));
                    //Console.WriteLine("a = {0}", a);
                    //double x0 = Math.Sqrt(2.0) * a;
                    //Console.WriteLine("x0 = {0}", x0);
                    //Func<double, double> meniscus = x => a * (Math.Log((2.0 * a / x) + Math.Sqrt((2.0 * a / x).Pow2() - 1))
                    //- Math.Log((2.0 * a / x0) + Math.Sqrt((2.0 * a / x0).Pow2() - 1)) + Math.Sqrt(4.0 - (x0 / a).Pow2()) - Math.Sqrt(4.0 - (x / a).Pow2()));

                    PhiFunc = (X => X[1] - (h0 + meniscus(R - X[0])));
                } else {
                    PhiFunc = (X => X[1] - h0);
                }
                C.InitialValues_Evaluators.Add("Phi", PhiFunc);

                C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
                C.InitialValues_Evaluators.Add("GravityY#B", X => -g);
                //C.InitialValues_Evaluators.Add("Pressure#A", X => C.ThermalParameters.p_sat); // -C.PhysicalParameters.rho_A * g * (X[1] - h0) + C.ThermalParameters.p_sat);
                //C.InitialValues_Evaluators.Add("Pressure#B", X => C.ThermalParameters.p_sat);
            } else {
                //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);
            }

            #endregion


            // boundary conditions
            // ===================
            #region BC

            if (startUp_Interface) {
                if (pressureCondition) {
                    C.AddBoundaryValue("pressure_Dirichlet_lower", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                    C.AddBoundaryValue("pressure_Dirichlet_left", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                    C.AddBoundaryValue("pressure_Dirichlet_left", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                } else {
                    C.AddBoundaryValue("wall_lower");
                }

                C.AddBoundaryValue("pressure_Dirichlet_upper", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                C.AddBoundaryValue("pressure_Dirichlet_upper", "Pressure#A", (X, t) => C.ThermalParameters.p_sat);
                //C.AddBoundaryValue("navierslip_linear_right", "VelocityY#A", (X, t) => 0.0);
                //C.AddBoundaryValue("navierslip_linearx_right", "VelocityY#B", (X, t) => 0.0);
                C.AddBoundaryValue("freeslip_right");
            }

            double Uwall = 0.1;

            if (run) {
                double HeatFluxA = 0.0;
                if (Water) {
                    HeatFluxA = 4500;
                } else if (FC72) {
                    HeatFluxA = 3000.0;
                }

                double rampUpTime = 1e-4;

                C.AddBoundaryValue("wall_lower");
                C.AddBoundaryValue("pressure_outlet_upper", "Pressure#B", (X, t) => C.ThermalParameters.p_sat);
                C.AddBoundaryValue("navierslip_linear_right", "VelocityY#A", (X, t) => Uwall * Math.Tanh((Math.PI / rampUpTime) * t));
                C.AddBoundaryValue("navierslip_linear_right", "VelocityY#B", (X, t) => Uwall * Math.Tanh((Math.PI / rampUpTime) * t));

            }


            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = hmin;

            double beta_scl = 100;
            C.PhysicalParameters.betaS_A = 0.0; // beta_scl * C.PhysicalParameters.mu_A / (hmin / 10.0);
            C.PhysicalParameters.betaS_B = 0.0; // beta_scl * C.PhysicalParameters.mu_B / (hmin / 1000.0);
            //C.PhysicalParameters.betaL = 0.01;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.AgglomerationThreshold = 0.1;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 4e5;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.LinearSolver = LinearSolverCode.classic_mumps.GetConfig();
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.Default;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            //C.AdvancedDiscretizationOptions.freeSurfaceFlow = true;
            //C.PhysicalParameters.pFree = C.ThermalParameters.p_sat;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.CurvatureNeeded = false;

            //C.InterVelocAverage = XNSE_Control.InterfaceVelocityAveraging.phaseA;

            #endregion


            // level-set
            // =========
            #region levelset

            if (startUp_Interface || run) {
                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;
                C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            } else {
                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            }


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-4;
            C.dtMin = 1e-4;
            C.Endtime = 10.0;
            C.NoOfTimesteps = 100000;
            C.saveperiod = 10;

            #endregion


            C.ClearVelocitiesOnRestart = false;


            // parameters
            // ==========

            double[] param = new double[1];
            param[0] = Uwall;   // wall velocity

            C.AdditionalParameters = param;


            return C;
        }




        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control ThermodynamicEquilibrium_steadyStateTest(int p = 2, int kelemR = 17, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            bool steady = true;
            bool separated = false;

            //_DbPath = @"D:\local\local_XNSE_StudyDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/ThermodynamEquilib_steady";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            }); 

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Water (A: liquid, B: gaseous)
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1.0;    
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0; 

            C.solveCoupledHeatEquation = true;
            if (separated)
                C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.LDG;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 0.001;
            C.ThermalParameters.k_A = 1.0;
            double kv = 0.1;
            C.ThermalParameters.k_B = kv;

            if(C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = 100.0; 
                //C.ThermalParameters.hVap_B = -100.0;
            }

            double Tsat = 100.0;
            C.ThermalParameters.T_sat = Tsat;
            double pSat = 10;
            C.ThermalParameters.p_sat = pSat;


            bool includeConv = true;
            C.PhysicalParameters.IncludeConvection = includeConv; 
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid
            
            double L = 0.1;
            bool periodic = false;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: periodic);

                if (!steady) {
                    grd.EdgeTagNames.Add(1, "wall_ConstantHeatFlux_lower");
                    grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ZeroGradient_upper");
                } else {
                    grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature_lower");
                    grd.EdgeTagNames.Add(2, "velocity_inlet_ZeroGradient_upper");
                }
                if (!periodic)
                    grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)   
                        et = 1; 
                    if(Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;
                    if (!periodic) {
                        if (Math.Abs(X[0]) <= 1.0e-8 || Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 3;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double zi0 = 0.01;

            Func<double[], double> PhiFunc = (X => zi0 - X[1]);   

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double qv = 10.0;
            C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);
            C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat + (qv/kv)*(zi0 - X[1]));

            C.prescribedMassflux_Evaluator = (X, t) => -0.1;

            if (!steady) {
                C.InitialValues_Evaluators.Add("Pressure#A", X => pSat);
                C.InitialValues_Evaluators.Add("Pressure#B", X => pSat - (0.01) * (10 - 1));

                C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.9);
            } else {
                C.InitialValues_Evaluators.Add("VelocityY#A", X => -0.1);
                C.InitialValues_Evaluators.Add("VelocityY#B", X => -1.0);
            }

            //if (separated) {
            //    C.InitialValues_Evaluators.Add("HeatFluxY#B", X => qv);
            //}

            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            if (!steady) {
                //C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFlux#A", (X, t) => HeatFlux);
                if (separated) {
                    C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#B", (X, t) => qv);
                } else {
                    C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#B", (X, t) => -qv);
                }

                C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper", "Pressure#A", (X, t) => pSat);
            } else {
                C.AddBoundaryValue("pressure_outlet_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat + (qv / kv) * (zi0 - X[1]));
                C.AddBoundaryValue("velocity_inlet_ZeroGradient_upper", "VelocityY#A", (X, t) => -0.1);
            }

            if (!periodic)
                C.AddBoundaryValue("slipsymmetry_ZeroGradient");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.useSolutionParamUpdate = true;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = false; 
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 1;

            C.InitSignedDistance = false;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = steady? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = steady ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = steady ? AppControl._TimesteppingMode.Steady : AppControl._TimesteppingMode.Transient;
            C.dtMax = 5e-4;
            C.dtMin = 5e-4;
            C.Endtime = 10000;
            C.NoOfTimesteps = 200;
            C.saveperiod = 2;

            #endregion


            // additional parameters
            double[] param = new double[2];
            param[0] = L;           // domain height
            param[1] = L / 4.2;   // x probe

            C.AdditionalParameters = param;

            //C.LogValues = XNSE_Control.LoggingValues.EvaporationL;
            //C.LogPeriod = 2;
            //C.PostprocessingModules.Add(new EvaporationLogging() { LogPeriod = 2, mode = EvaporationLogging.Mode.LineInterface });

            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control ThermodynamicEquilibrium_unsteadyTest(int p = 3, int kelemR = 17, string _DbPath = null, int setUp = 1) {

            XNSE_Control C = new XNSE_Control();

            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            bool superheatedVapor = (setUp == 1) ? true : false;
            bool subcooledLiquid = (setUp == 2) ? true : false;
            bool superheatedLiquid = (setUp == 3) ? true : false;

            bool radial = false;

            //_DbPath = @"D:\local\local_XNSE_StudyDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";
            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/ThermodynamEquilib_unsteady";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";
            C.SessionName = "1Dheat_unsteadyTest_noHeatConv";

            //C.LogValues = XNSE_Control.LoggingValues.EvaporationL;
            //C.LogPeriod = 10;
            //C.PostprocessingModules.Add(new EvaporationLogging() { LogPeriod = 10, mode = EvaporationLogging.Mode.LineInterface });

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Water (A: liquid, B: gaseous)
            double rhol = 586.5;
            C.PhysicalParameters.rho_A = rhol;
            double rhov = 106.4;
            C.PhysicalParameters.rho_B = rhov;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            C.solveCoupledHeatEquation = true;
            //C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.LDG;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            double Cpl = 9.35e+3;
            C.ThermalParameters.c_A = Cpl;
            double Cpv = 15.4e+3;
            C.ThermalParameters.c_B = Cpv;
            double kl = 0.444;
            C.ThermalParameters.k_A = kl;
            double kv = 0.114;
            C.ThermalParameters.k_B = kv;

            double alpha_v = 0.6958e-7;
            double alpha_l = 0.8096e-7;

            //double hlg = 941.0;
            //C.ThermalParameters.hVap = hlg;

            double Tsat = 620.0;
            C.ThermalParameters.T_sat = Tsat;
            double P0 = 160e+5;
            C.ThermalParameters.p_sat = P0;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = false;

            //C.useSolutionParamUpdate = true;

            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid


            double Lv0 = 5e-3;
            double Ll0 = 5.3e-3;
            double L = Lv0 + Ll0;

            if (!radial) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, Lv0 / 2.0, (kelemR / 4) + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, L, kelemR + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                    grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                    //grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ConstantTemperature_upper"); 
                    grd.EdgeTagNames.Add(2, "pressure_outlet_ZeroGradient_upper"); 
                    grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - L) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0]) <= 1.0e-8 || Math.Abs(X[0] - (Lv0 / 2.0)) <= 1.0e-8)
                            et = 3;

                        return et;
                    });
                     
                    return grd;
                };

            } else {

                 C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0, L, kelemR + 1);
                    double[] Ynodes = GenericBlas.Linspace(0, L, kelemR + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "slipsymmetry_ZeroGradient");
                    //grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ConstantTemperature_upper"); 
                    grd.EdgeTagNames.Add(2, "pressure_outlet_ZeroGradient_upper");
                    //grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient");

                     grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1]) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - L) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0]) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 2;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init


            double deltaT = 0.0;

            double zi0 = Lv0;
            double t0 = 0.5;

            if (superheatedVapor) {

                double hlg = 941.0;
                C.ThermalParameters.hVap = hlg;

                double lmbdv = 1.71814;
                zi0 += 2.0 * lmbdv * Math.Sqrt(alpha_v * t0);

                deltaT = 5.0;

                Func<double, double, double> Tempv = 
                    //(zeta, t) => Tsat + deltaT * (ErrorFuncUtil.erf(lmbdv) - ErrorFuncUtil.erf(lmbdv + Math.Sqrt(1.0 / alpha_v) * (zeta / (2.0 * Math.Sqrt(t))))) / (1.0 + ErrorFuncUtil.erf(lmbdv));
                    (zeta, t) => Tsat + deltaT * (SpecialFunctions.Erf(lmbdv) - SpecialFunctions.Erf(lmbdv + Math.Sqrt(1.0 / alpha_v) * (zeta / (2.0 * Math.Sqrt(t))))) / (1.0 + SpecialFunctions.Erf(lmbdv));

                C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);

                C.InitialValues_Evaluators.Add("Temperature#B", X => Tempv((X[1] - zi0), t0));

                //if (C.conductMode != Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                //    Func<double, double, double> HeatFluxB =
                //          (zeta, t) => kv * (deltaT / (2.0 * Math.Sqrt(alpha_v * (t)))) * (2.0 / Math.Sqrt(Math.PI)) * Math.Exp(-(lmbdv + (zeta/(2.0*Math.Sqrt(alpha_v * (t))))).Pow2()) / (1 + SpecialFunctions.Erf(lmbdv));
                //    C.InitialValues_Evaluators.Add("HeatFluxY#B", X => HeatFluxB((X[1] - zi0 - 2.0 * lmbdv * Math.Sqrt(alpha_v * t0)), t0));
                //}

                //Func<double[], double, double> mdot = (X, t) => -(kv / hlg) * (deltaT / (2.0 * Math.Sqrt(alpha_v * (t0 + t)))) * (2.0 / Math.Sqrt(Math.PI)) * Math.Exp(-lmbdv.Pow2()) / (1 + SpecialFunctions.Erf(lmbdv));
                Func<double[], double, double> mdot = (X, t) => -rhov * Math.Sqrt(alpha_v) * lmbdv / Math.Sqrt(t0 + t);
                //C.prescribedMassflux_Evaluator = mdot;

                //Func<double[], double, double> Vl = (X, t) => -mdot(X, t) * ((1.0 / rhov) - (1.0 / rhol));
                Func<double[], double, double> Vl = (X, t) => -mdot(X, t) * ((1.0 / rhov) - (1.0 / rhol));

                C.InitialValues_Evaluators.Add("VelocityY#A", X => Vl(X, 0));
                C.InitialValues_Evaluators.Add("Pressure#B", X => P0 - mdot(X, 0).Pow2() * ((1.0 / rhov) - (1.0 / rhol)));

                C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);
                C.InitialValues_Evaluators.Add("Pressure#A", X => P0);


                Func<double, double> zi = t => Lv0 + 2.0 * lmbdv * Math.Sqrt(alpha_v * (t0 + t));
                //Func<double[], double> PhiFunc = (X => zi(t0) - X[1]);
                Func<double[], double, double> PhiFunc = ((X, t) => zi(t) - X[1]);

                C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0.0));
                C.Phi = PhiFunc;

            }


            //if (subcooledLiquid) {

            //    Func<double, double> zi = t => zi0 + 2.0 * lmbdl * Math.Sqrt(alpha_l * t);
            //    Func<double[], double> PhiFunc = (X => zi(t0) - X[1]);

            //    C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //    Func<double, double, double> Tempv =
            //        (zeta, t) => Tsat + -deltaT * (SpecialFunctions.Erf(lmbdl + Math.Sqrt(1.0 / alpha_l) * (zeta / (2.0 * Math.Sqrt(t)))) - SpecialFunctions.Erf(lmbdl)) / (1.0 - SpecialFunctions.Erf(lmbdl));

            //    C.InitialValues_Evaluators.Add("Temperature#A", X => Tempv((X[1] - zi0 - 2.0 * lmbdl * Math.Sqrt(alpha_l * t0)), t0));
            //    C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat);

            //    //if (C.conductMode != Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP) {
            //    //    Func<double, double, double> HeatFluxA =
            //    //          (zeta, t) => (kl * deltaT / (2.0 * Math.Sqrt(alpha_l * (t0 + t)))) * (2.0 / Math.Sqrt(Math.PI)) * Math.Exp(-lmbdl.Pow2()) / (1 - SpecialFunctions.Erf(lmbdl));
            //    //    C.InitialValues_Evaluators.Add("HeatFluxY#A", X => HeatFluxA((X[1] - zi0 - 2.0 * lmbdv * Math.Sqrt(alpha_v * t0)), t0));
            //    //}

            //    Func<double[], double, double> mdot = (X, t) => -(kl / hlg) * (-deltaT / (2.0 * Math.Sqrt(alpha_l * (t0 + t)))) * (2.0 / Math.Sqrt(Math.PI)) * Math.Exp(-lmbdl.Pow2()) / (1 - SpecialFunctions.Erf(lmbdl));
            //    //Func<double, double> mdot = t => rhol * Math.Sqrt(alpha_l) * lmbdl / Math.Sqrt(t0 + t);
            //    //C.prescribedMassflux_Evaluator = mdot;

            //    Func<double[], double, double> Vl = (X, t) => -mdot(X, t) * ((1.0 / rhov) - (1.0 / rhol));

            //    C.InitialValues_Evaluators.Add("VelocityY#A", X => Vl(X, 0));

            //    C.InitialValues_Evaluators.Add("Pressure#B", X => P0 - mdot(X, 0).Pow2() * ((1.0 / rhov) - (1.0 / rhol)));

            //    C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);
            //    C.InitialValues_Evaluators.Add("Pressure#A", X => P0);
            //}


            if (superheatedLiquid) {

                double hlg = 941.0e2;
                C.ThermalParameters.hVap = hlg;

                double lmbdl = 0.0872;

                deltaT = 1.0;


                Func<double, double, double> Tempv =
                    (zeta, t) => Tsat + deltaT * (SpecialFunctions.Erf(lmbdl + Math.Sqrt(1.0 / alpha_l) * (zeta / (2.0 * Math.Sqrt(t)))) - SpecialFunctions.Erf(lmbdl)) / (1.0 - SpecialFunctions.Erf(lmbdl));

                if (!radial) {
                    C.InitialValues_Evaluators.Add("Temperature#A", X => Tempv((X[1] - zi0 - 2.0 * (rhol / rhov) * lmbdl * Math.Sqrt(alpha_l * t0)), t0));
                    C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat);
                } else {
                    C.InitialValues_Evaluators.Add("Temperature#B", X => Tempv((Math.Sqrt(X[0].Pow2() + X[1].Pow2()) - zi0 - 2.0 * (rhol / rhov) * lmbdl * Math.Sqrt(alpha_l * t0)), t0));
                    C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);
                }

                //if (C.conductMode != Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP) {
                //    Func<double, double, double> HeatFluxA =
                //          (zeta, t) => (kl * deltaT / (2.0 * Math.Sqrt(alpha_l * (t0 + t)))) * (2.0 / Math.Sqrt(Math.PI)) * Math.Exp(-lmbdl.Pow2()) / (1 - SpecialFunctions.Erf(lmbdl));
                //    C.InitialValues_Evaluators.Add("HeatFluxY#A", X => HeatFluxA((X[1] - zi0 - 2.0 * lmbdv * Math.Sqrt(alpha_v * t0)), t0));
                //}

                //Func<double[], double, double> mdot = (X, t) => -(kl / hlg) * (deltaT / (2.0 * Math.Sqrt(alpha_l * (t0 + t)))) * (2.0/Math.Sqrt(Math.PI)) * Math.Exp(-lmbdl.Pow2()) / (1 - SpecialFunctions.Erf(lmbdl));
                Func<double[], double, double> mdot = (X, t) => -rhol * Math.Sqrt(alpha_l) * lmbdl / Math.Sqrt(t0 + t);
                //C.prescribedMassflux_Evaluator = mdot;

                Func<double[], double, double> Vl = (X, t) => -mdot(X, t) * ((1.0 / rhov) - (1.0 / rhol));

                if (!radial) {
                    C.InitialValues_Evaluators.Add("VelocityY#A", X => Vl(X, 0));
                    C.InitialValues_Evaluators.Add("Pressure#B", X => P0 - mdot(X, 0).Pow2() * ((1.0 / rhov) - (1.0 / rhol)));
                    C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);
                    C.InitialValues_Evaluators.Add("Pressure#A", X => P0);
                } else {
                    C.InitialValues_Evaluators.Add("VelocityX#B", X => Vl(X, 0) * (X[0] / Math.Sqrt(X[0].Pow2() + X[1].Pow2())));
                    C.InitialValues_Evaluators.Add("VelocityY#B", X => Vl(X, 0) * (X[1] / Math.Sqrt(X[0].Pow2() + X[1].Pow2())));

                    C.InitialValues_Evaluators.Add("Pressure#A", X => P0 - mdot(X, 0).Pow2() * ((1.0 / rhov) - (1.0 / rhol)));

                    C.InitialValues_Evaluators.Add("Pressure#B", X => P0);
                }

                Func<double, double> zi = t => zi0 + 2.0 * (rhol / rhov) * lmbdl * Math.Sqrt(alpha_l * (t0+t));
                Func<double[], double, double> PhiFunc;
                if (!radial)
                    PhiFunc = ((X, t) => zi(t) - X[1]);
                else
                    PhiFunc = ((X, t) => (X[0].Pow2() + X[1].Pow2()).Sqrt() - zi(t));

                C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0.0));
                C.Phi = PhiFunc;

            }


            //C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);
            //C.InitialValues_Evaluators.Add("Pressure#A", X => P0);


            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            if (superheatedVapor) {
                C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat + deltaT);
                //C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Temperature#A", (X, t) => Tsat);
                //C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Pressure#A", (X, t) => P0);
                C.AddBoundaryValue("pressure_outlet_ZeroGradient_upper", "Pressure#A", (X, t) => P0);
            }

            if (subcooledLiquid) {
                C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat);
                C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Temperature#A", (X, t) => Tsat - deltaT);
                C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Pressure#A", (X, t) => P0);
            }

            if (superheatedLiquid) {
                if (!radial) {
                    C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Temperature#A", (X, t) => Tsat + deltaT);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Pressure#A", (X, t) => P0);
                } else {
                    C.AddBoundaryValue("slipsymmetry_ZeroGradient");
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature", "Temperature#B", (X, t) => Tsat + deltaT);
                    C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature", "Pressure#B", (X, t) => P0);
                }
            }

            C.AddBoundaryValue("slipsymmetry_ZeroGradient");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.LinearSolver = LinearSolverCode.classic_mumps.GetConfig();
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 5;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 3;
            C.AMR_startUpSweeps = 3;

            //C.ReInitPeriod = 10;
            C.InitSignedDistance = false;

            #endregion


            // level-set
            // =========
            #region levelset

            //C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;
            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = (setUp == 1) ? 1e-2 : 1e-3;
            C.dtMin = (setUp == 1) ? 1e-2 : 1e-3;
            C.Endtime = 5;
            C.NoOfTimesteps = (setUp == 1) ? 1000 : 2000;
            C.saveperiod = 10;

            #endregion


            // additional parameters
            double[] param = new double[2];
            param[0] = L;           // domain height
            param[1] = Lv0 / 4.2;   // x probe

            C.AdditionalParameters = param;


            return C;
        }


        public static XNSE_Control StefanProblem_Sato(int p = 2, int kelem = 34, string _DbPath = null) {
            
            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";
            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/OneDimensionalVerification_Sato";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";

            //C.LogValues = XNSE_Control.LoggingValues.EvaporationL;
            //C.LogPeriod = 10;
            //C.PostprocessingModules.Add(new EvaporationLogging() { LogPeriod = 10, mode = EvaporationLogging.Mode.LineInterface });

            C.ContinueOnIoError = false;

            // additional parameters
            double[] param = new double[2];
            param[0] = 0.001;           // domain height L
            param[1] = 0.0001;   // x probe

            C.AdditionalParameters = param;

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // A: liquid (Water), B: vapor
            double rhol = 958.4;
            C.PhysicalParameters.rho_A = rhol;
            double rhov = 0.597;
            C.PhysicalParameters.rho_B = rhov;
            C.PhysicalParameters.mu_A = 2.8e-2;
            C.PhysicalParameters.mu_B = 1.26e-1;
            C.PhysicalParameters.Sigma = 0.0; // 0.059;

            C.solveCoupledHeatEquation = true;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            double Cpl = 4216;
            C.ThermalParameters.c_A = Cpl;
            double Cpv = 2030;
            C.ThermalParameters.c_B = Cpv;
            double kl = 0.679;
            C.ThermalParameters.k_A = kl;
            double kv = 0.025;
            C.ThermalParameters.k_B = kv;


            double hlg = 2.26e+6;
            C.ThermalParameters.hVap = hlg;
            //C.ThermalParameters.hVap_B = -hlg;

            double Tsat = 373.15;
            C.ThermalParameters.T_sat = Tsat;
            double P0 = 101.3e+3;
            C.ThermalParameters.p_sat = P0;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid


            double L = 0.001; ;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L / 8.0, (kelem / 8) + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ConstantTemperature_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double alpha_v = kv / (rhov * Cpv);
            double lmbd = 0.067;

            double t0 = 0.01;
            double deltaT = 10;

            Func<double, double> Xi = t => 2.0 * lmbd * Math.Sqrt(alpha_v * (t0 + t));
            Func<double[], double> PhiFunc = (X => Xi(0) - X[1]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            Func<double, double, double> TempV =
                (x, t) => (Tsat + deltaT) + (-deltaT / SpecialFunctions.Erf(lmbd)) * SpecialFunctions.Erf(x / (2.0 * Math.Sqrt(alpha_v * (t0 + t))));

            //if (C.conductMode != Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP) {
            //    Func<double, double, double> HeatFluxV =
            //          (x, t) => -kv * (-deltaT / SpecialFunctions.Erf(lmbd)) * (1.0 / Math.Sqrt(Math.PI * alpha_v * (t0 + t))) * Math.Exp(-(x / (2.0 * Math.Sqrt(alpha_v * (t0 + t)))).Pow2());
            //    //C.InitialValues_Evaluators.Add("HeatFluxY#B", X => HeatFluxV(X[1], 0));
            //}

            //C.prescribedMassflux_Evaluator = (X, t) => 0.0;

            C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);
            C.InitialValues_Evaluators.Add("Temperature#B", X => TempV(X[1], 0));

            C.InitialValues_Evaluators.Add("Pressure#A", X => P0);


            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat + deltaT);
            C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Temperature#A", (X, t) => Tsat);
            C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Pressure#A", (X, t) => P0);

            //C.AddBoundaryValue("slipsymmetry_ZeroGradient");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 40;
            C.NonLinearSolver.ConvergenceCriterion = 1e-7;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.useSolutionParamUpdate = true;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 2;
            C.AMR_startUpSweeps = 2;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 0.4;
            C.NoOfTimesteps = 390;
            C.saveperiod = 1;

            #endregion


            return C;
        }


        public static XNSE_Control SuckingProblem_Sato(int p = 2, int kelem = 34, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";
            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/OneDimensionalVerification_Sato";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";

            //C.LogValues = XNSE_Control.LoggingValues.EvaporationL;
            //C.LogPeriod = 10;
            //C.PostprocessingModules.Add(new EvaporationLogging() { LogPeriod = 10, mode = EvaporationLogging.Mode.LineInterface });

            C.ContinueOnIoError = false;

            // additional parameters
            double[] param = new double[2];
            param[0] = 0.006;           // domain height L
            param[1] = 0.0003;   // x probe

            C.AdditionalParameters = param;

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // A: liquid (Water), B: vapor
            double rhol = 958.4;
            C.PhysicalParameters.rho_A = rhol;
            double rhov = 0.597;
            C.PhysicalParameters.rho_B = rhov;
            C.PhysicalParameters.mu_A = 2.8e-2;
            C.PhysicalParameters.mu_B = 1.26e-1;
            C.PhysicalParameters.Sigma = 0.0; // 0.059;

            C.solveCoupledHeatEquation = true;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            double Cpl = 4216;
            C.ThermalParameters.c_A = Cpl;
            double Cpv = 2030;
            C.ThermalParameters.c_B = Cpv;
            double kl = 0.679;
            C.ThermalParameters.k_A = kl;
            double kv = 0.025;
            C.ThermalParameters.k_B = kv;


            double hlg = 2.26e+6;
            C.ThermalParameters.hVap = hlg;
            //C.ThermalParameters.hVap_B = -hlg;

            double Tsat = 373.15;
            C.ThermalParameters.T_sat = Tsat;
            double P0 = 101.3e+3;
            C.ThermalParameters.p_sat = P0;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid


            double L = 0.006;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L / 8.0, (kelem / 8) + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ConstantTemperature_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double alpha_l = kl / (rhol * Cpl);
            double lmbd = 0.0074;

            double t0 = 0.1;
            double deltaT = 5;

            double zi0 = 0.001;
            Func<double, double> zi = t => zi0 + 2.0 * lmbd * Math.Sqrt(alpha_l * (t0 + t));
            Func<double[], double> PhiFunc = (X => zi(0) - X[1]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double, double> TempL =
                //(zeta, t) => Tsat + deltaT * (SpecialFunctions.Erf(lmbd + Math.Sqrt(1.0 / alpha_l) * (zeta / (2.0 * Math.Sqrt(t0 + t)))) - SpecialFunctions.Erf(lmbd)) / (1.0 - SpecialFunctions.Erf(lmbd));
                (zeta, t) => Tsat + (deltaT / (L - zi(t)) )* zeta;

            //if (C.conductMode != Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP) {
            //    Func<double, double, double> HeatFluxL =
            //          (zeta, t) => -(kl * deltaT / (2.0 * Math.Sqrt(alpha_l * (t0 + t)))) * (2.0 / Math.Sqrt(Math.PI)) * Math.Exp(-lmbd.Pow2()) / (1 - SpecialFunctions.Erf(lmbd));
            //    C.InitialValues_Evaluators.Add("HeatFluxY#A", X => HeatFluxL((X[1] - zi0 - 2.0 * lmbd * Math.Sqrt(alpha_l * t0)), t0));
            //}

            Func<double[], double, double> mdot = (X, t) => -(kl / hlg) * (deltaT / (2.0 * Math.Sqrt(alpha_l * (t0 + t)))) * (2.0/Math.Sqrt(Math.PI)) * Math.Exp(-lmbd.Pow2()) / (1.0 - SpecialFunctions.Erf(lmbd));
            //C.prescribedMassflux_Evaluator = mdot;

            C.InitialValues_Evaluators.Add("Temperature#A", X => TempL((X[1] - zi0 - 2.0 * lmbd * Math.Sqrt(alpha_l * t0)), 0));
            C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat);

            C.InitialValues_Evaluators.Add("Pressure#A", X => P0);


            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat);
            C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Temperature#A", (X, t) => Tsat + deltaT);
            C.AddBoundaryValue("pressure_Dirichlet_ConstantTemperature_upper", "Pressure#A", (X, t) => P0);

            //C.AddBoundaryValue("slipsymmetry_ZeroGradient");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.useSolutionParamUpdate = true;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 2;
            C.AMR_startUpSweeps = 2;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 5e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 0.6;
            C.NoOfTimesteps = 100;
            C.saveperiod = 1;


            #endregion

            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control BubbleGrowth(int p = 2, int kelem = 17, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();


            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; // C.DbPath != null;
            C.ProjectName = "XNSE/BubbleGrowth";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Water (A: vapor, B: liquid)
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1000;
            C.PhysicalParameters.mu_A = 1.78e-5;
            C.PhysicalParameters.mu_B = 0.001;
            double sigma = 0.07;
            C.PhysicalParameters.Sigma = 0.07;

            C.solveCoupledHeatEquation = true;
            C.prescribedMassflux_Evaluator = (X, t) => 0.1;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1;
            C.ThermalParameters.c_B = 1;
            C.ThermalParameters.k_A = 1;
            C.ThermalParameters.k_B = 1;


            if (C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = -1;
                //C.ThermalParameters.hVap_B = 1;
            }

            double T_sat = 0.0;
            C.ThermalParameters.T_sat = T_sat;
            //double pSat = 10;
            //C.ThermalParameters.p_sat = pSat;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double l = 0.008;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-l / 2.0, l / 2.0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-l / 2.0, l / 2.0, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);


                grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (l / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - (l / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (l / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (l / 2.0)) <= 1.0e-8)
                        et = 1;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.001;
            double[] center = new double[] { 0, 0 };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / R);


            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("pressure_Outlet_ConstantTemperature", "Temperature#A", (X, t) => 0.0);
            C.AddBoundaryValue("pressure_Outlet_ConstantTemperature", "Temperature#B", (X, t) => 0.0);


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 1;
            C.AMR_startUpSweeps = 1;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-5;
            C.dtMin = 1e-5;
            C.Endtime = 0.01;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control BubbleGrowth_Scriven(int p = 3, int kelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            bool computeQuarter = true;


            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            //_DbPath = @"\\hpccluster\hpccluster-scratch\smuda\XNSFE_testDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";


            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/BubbleGrowth";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Water (A: vapor, B: liquid)
            double rho_vap = 0.59;
            C.PhysicalParameters.rho_A = rho_vap;
            double rho_liq = 958;
            C.PhysicalParameters.rho_B = rho_liq;
            C.PhysicalParameters.mu_A = 1.23e-1;
            C.PhysicalParameters.mu_B = 2.82e-4;
            double sigma = 0.059;
            C.PhysicalParameters.Sigma = 0.07; // sigma;

            C.solveCoupledHeatEquation = true;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            double Cp_vap = 2034;
            C.ThermalParameters.c_A = Cp_vap;
            double Cp_liq = 4216;
            C.ThermalParameters.c_B = Cp_liq;
            C.ThermalParameters.k_A = 0.026;
            double k_liq = 0.6;
            C.ThermalParameters.k_B = k_liq;

            double alpha = k_liq / (rho_liq * Cp_liq);
            double epsilon = 1 - (rho_vap / rho_liq);

            double L_vap = 2.257e+6;
            if (C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = L_vap;
                //C.ThermalParameters.hVap_B = L_vap;
            }

            double T_sat = 373.0;
            C.ThermalParameters.T_sat = T_sat;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = false; 
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double l = 12e-3;

            if (!computeQuarter) {

                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-l / 2.0, l / 2.0, kelem + 1);
                    double[] Ynodes = GenericBlas.Linspace(-l / 2.0, l / 2.0, kelem + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);


                    grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");


                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1] + (l / 2.0)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[0] - (l / 2.0)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - (l / 2.0)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[0] + (l / 2.0)) <= 1.0e-8)
                            et = 1;

                        return et;
                    });

                    return grd;
                };

            } else {

                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(0.0, l / 2.0, (kelem/2) + 1);
                    double[] Ynodes = GenericBlas.Linspace(0.0, l / 2.0, (kelem/2) + 1);
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);


                    grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");
                    grd.EdgeTagNames.Add(2, "slipsymmetry_ZeroGradient");


                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[1] + (0.0)) <= 1.0e-8)
                            et = 2;
                        if (Math.Abs(X[0] - (l / 2.0)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[1] - (l / 2.0)) <= 1.0e-8)
                            et = 1;
                        if (Math.Abs(X[0] + (0.0)) <= 1.0e-8)
                            et = 2;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 1e-3;
            double[] center = new double[] { 0, 0 };



            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / R);

            double Ja = 3.0;
            //double Ja = 10.0;
            double T_inf = (Ja * rho_vap * L_vap / (rho_liq * Cp_liq)) + T_sat;

            double beta = 3.3277; // 3.3293;
            //double beta = 10.231;

            double t_0 = (R / (2.0*beta)).Pow2() / alpha;

            //Func<double, double> intgrnd = (x => Math.Exp(-Math.Pow(x, 2) - ((2.0 * epsilon * beta.Pow(3)) / x)) / Math.Pow(x, 2));
            //Func<double[], double> Temp = (X => T_inf + (T_sat - T_inf) * (2.0 * beta.Pow(3) / Ja) * Math.Exp(beta.Pow2() + 2.0 * epsilon * beta.Pow2())
            //                                    * MathNet.Numerics.Integration.Integrate.OnClosedInterval(intgrnd, (Math.Sqrt(X[0].Pow2() + X[1].Pow2()) + R) / (2.0 * Math.Sqrt(alpha * t_0)), double.MaxValue));

            Func<double, double> intgrnd = (z => Math.Exp(-beta.Pow2() * (Math.Pow(1 - z, -2) - 2.0 * epsilon * z - 1)));
            Func<double[], double> Temp = (X => T_inf - 2.0 * beta.Pow2() * ((L_vap + (Cp_liq - Cp_vap) * (T_inf - T_sat)) / (Ja * L_vap)) * (T_inf - T_sat)
                                                * MathNet.Numerics.Integration.GaussKronrodRule.Integrate(intgrnd, 1 - (R / Math.Sqrt(X[0].Pow2() + X[1].Pow2())), 1, out double IntError, out double L1norm));

            C.InitialValues_Evaluators.Add("Temperature#B", Temp);
            C.InitialValues_Evaluators.Add("Temperature#A", X => T_sat);

            Func<double[], double, double> mdot = (X, t) => rho_vap * 2.0 * beta * Math.Sqrt(alpha / (t_0 + t));
            C.prescribedMassflux_Evaluator = mdot;

            //Func<double[], double, double> Vl = (X, t) => -mdot(X, t) * ((1.0 / rhov) - (1.0 / rhol));

            //C.InitialValues_Evaluators.Add("VelocityY#A", X => Vl(X, 0));
            //C.InitialValues_Evaluators.Add("Pressure#B", X => P0 - mdot(X, 0).Pow2() * ((1.0 / rhov) - (1.0 / rhol)));

            //C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);
            //C.InitialValues_Evaluators.Add("Pressure#A", X => P0);

            Func<double[], double, double> PhiFunc = ((X, t) => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - (2.0 * beta * Math.Sqrt(alpha * (t_0 + t))));

            C.InitialValues_Evaluators.Add("Phi", X => PhiFunc(X, 0.0));
            C.Phi = PhiFunc;


            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("pressure_Outlet_ConstantTemperature", "Temperature#B", (X, t) => T_inf);
            if (computeQuarter)
                C.AddBoundaryValue("slipsymmetry_ZeroGradient");


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 2;
            C.AMR_startUpSweeps = 2;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
             

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-5;
            C.dtMin = 1e-5;
            C.Endtime = 4.0*t_0;
            C.NoOfTimesteps = 60000;
            //C.NoOfTimesteps = 7000;
            C.saveperiod = 10;

            #endregion


            return C;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control DropVaporization(int p = 2, int kelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            bool OnWall = true;

            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; // C.DbPath != null;
            C.ProjectName = "XNSE/DropVaporization";

            //C.LogValues = XNSE_Control.LoggingValues.EvaporationC;
            //C.PostprocessingModules.Add(new EvaporationLogging() { mode = EvaporationLogging.Mode.CircleInterface });
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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Water (A: liquid, B: vapor)
            C.PhysicalParameters.rho_A = 200;
            C.PhysicalParameters.rho_B = 5;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 0.005;
            double sigma = 0.1;
            C.PhysicalParameters.Sigma = sigma;

            C.solveCoupledHeatEquation = true;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 400;
            C.ThermalParameters.c_B = 200;
            C.ThermalParameters.k_A = 40;
            C.ThermalParameters.k_B = 1;

            C.ThermalParameters.hVap = 1000;
            //C.ThermalParameters.hVap_B = -1000;

            double T_sat = 0.0;
            C.ThermalParameters.T_sat = T_sat;
            //double pSat = 10;
            //C.ThermalParameters.p_sat = pSat;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double g = 9.81;
            double lambda = 2 * Math.PI * Math.Sqrt(3*0.1/(-g*(5-200)));

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-lambda / 2.0, lambda / 2.0, kelem + 1);
                double[] Ynodes;
                if (OnWall)
                    Ynodes = GenericBlas.Linspace(0, lambda / 2.0, (kelem/2) + 1);
                else
                    Ynodes = GenericBlas.Linspace(-lambda / 2.0, lambda / 2.0, kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");
                if (OnWall)
                    grd.EdgeTagNames.Add(2, "freeslip_zeroGradient_bottom");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (!OnWall) {
                        if (Math.Abs(X[1] + (lambda / 2.0)) <= 1.0e-8)
                            et = 1;
                    } else {
                        if (Math.Abs(X[1] + 0) <= 1.0e-8)
                            et = 2;
                    }
                    if (Math.Abs(X[0] - (lambda / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (lambda / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + (lambda / 2.0)) <= 1.0e-8)
                        et = 1;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.02;
            double[] center = new double[] { 0, 0 };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            double T_wall = 5;
            Func<double[], double> TempFunc = delegate (double[] X) {

                double dR = 0.005;

                double r = Math.Sqrt(X[0].Pow2() + X[1].Pow2());
                if (r < R + dR) {
                    return T_sat;
                } else {
                    return T_wall;
                }
            };

            C.InitialValues_Evaluators.Add("Temperature#A", (X => T_sat));
            C.InitialValues_Evaluators.Add("Temperature#B", (X => TempFunc(X)));

            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / R);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("pressure_Outlet_ConstantTemperature", "Temperature#A", (X, t) => T_wall);
            C.AddBoundaryValue("pressure_Outlet_ConstantTemperature", "Temperature#B", (X, t) => T_wall);

            if (OnWall)
                C.AddBoundaryValue("freeslip_zeroGradient_bottom");


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 2;
            C.AMR_startUpSweeps = 2;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 5e-5;
            C.dtMin = 5e-5;
            C.Endtime = 6;
            C.NoOfTimesteps = 12000;
            C.saveperiod = 1;

            #endregion


            return C;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control FilmBoiling(int p = 3, int kelem = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/FilmBoiling";

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
            //C.FieldOptions.Add("GravityY", new FieldOpts() {
            //    Degree = p,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //});
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // A: liquid (Water), B: vapor
            double rhol = 958.4;
            C.PhysicalParameters.rho_A = rhol;
            double rhov = 0.597;
            C.PhysicalParameters.rho_B = rhov;
            C.PhysicalParameters.mu_A = 2.8e-2;
            C.PhysicalParameters.mu_B = 1.26e-1;
            C.PhysicalParameters.Sigma = 0.0; // 0.059;

            C.solveCoupledHeatEquation = true;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            double Cpl = 4216;
            C.ThermalParameters.c_A = Cpl;
            double Cpv = 2030;
            C.ThermalParameters.c_B = Cpv;
            double kl = 0.679;
            C.ThermalParameters.k_A = kl;
            double kv = 0.025;
            C.ThermalParameters.k_B = kv;


            double hlg = 2.26e+6;
            C.ThermalParameters.hVap = hlg;
            //C.ThermalParameters.hVap_B = -hlg;

            double Tsat = 373.15;
            C.ThermalParameters.T_sat = Tsat;
            double P0 = 101.3e+3;
            C.ThermalParameters.p_sat = P0;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 0.006;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L / 2.0, L / 2.0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0.0, L, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);


                grd.EdgeTagNames.Add(1, "navierslip_linear_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_ConstantTemperature_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + 0.0) <= 1.0e-8)
                        et = 1;
                    //if (Math.Abs(X[0] - (l / 2.0)) <= 1.0e-8)
                    //    et = 1;
                    if (Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;
                    //if (Math.Abs(X[0] + (l / 2.0)) <= 1.0e-8)
                    //    et = 1;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double alpha_l = kl / (rhol * Cpl);
            double lmbd = 0.0074;

            double t0 = 0.1;
            double deltaT = 5;

            double zi0 = 0.002;
            Func<double, double> zi = t => zi0 + 2.0 * lmbd * Math.Sqrt(alpha_l * (t0 + t));
            Func<double[], double> PhiFunc = (X => X[1] - zi(0));

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double, double, double> TempL =
                (zeta, t) => Tsat + deltaT * (SpecialFunctions.Erf(lmbd + Math.Sqrt(1.0 / alpha_l) * (zeta / (2.0 * Math.Sqrt(t0 + t)))) - SpecialFunctions.Erf(lmbd)) / (1.0 - SpecialFunctions.Erf(lmbd));


            C.InitialValues_Evaluators.Add("Temperature#A", X => TempL((zi(0) - X[1]), 0));
            C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat);

            C.InitialValues_Evaluators.Add("GravityX#A", X => -9.81);
            C.InitialValues_Evaluators.Add("GravityX#B", X => -9.81);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("navierslip_linear_ConstantTemperature_lower", "Temperature#A", (X, t) => Tsat + deltaT);
            C.AddBoundaryValue("navierslip_linear_ConstantTemperature_lower", "VelocityX#A", (X, t) => 0.5);
            C.AddBoundaryValue("pressure_Outlet_ConstantTemperature_upper", "Temperature#B", (X, t) => Tsat);

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 2;
            C.AMR_startUpSweeps = 2;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 5.0;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 10;

            #endregion


            return C;

        }


        public static XNSE_Control MicroLayerRegime_Test(int p = 2, int kelem = 32, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();


            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; // C.DbPath != null;
            C.ProjectName = "XNSE/MicroLayerRegime";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";

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
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics


            // Water (A: vapour, B: liquid Water)
            C.PhysicalParameters.rho_A = 0.5974;
            double rho_l = 958;
            C.PhysicalParameters.rho_B = rho_l;
            C.PhysicalParameters.mu_A = 1.228e-5;
            double mu_l = 2.82e-4;
            C.PhysicalParameters.mu_B = mu_l;
            double sigma = 0.058;
            C.PhysicalParameters.Sigma = sigma;

            C.solveCoupledHeatEquation = true;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 2034;
            double c_l = 4216;
            C.ThermalParameters.c_B = c_l;
            C.ThermalParameters.k_A = 0.024;
            double k_l = 0.677;
            C.ThermalParameters.k_B = k_l;


            if (C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = 2.256e+6;
                //C.ThermalParameters.hVap_B = 2.256e+6;
            }

            double T_sat = 373.12;
            C.ThermalParameters.T_sat = T_sat;
            //double pSat = 10;
            //C.ThermalParameters.p_sat = pSat;


            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 2e-3;
            double hmin = L / kelem;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);


                    grd.EdgeTagNames.Add(1, "navierslip_linear_ConstantTemperature_lower");
                    grd.EdgeTagNames.Add(2, "pressure_Outlet_ZeroGradient_right");
                    grd.EdgeTagNames.Add(3, "pressure_Outlet_ZeroGradient_upper");
                    grd.EdgeTagNames.Add(4, "slipsymmetry_ZeroGradient_left");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.06e-3;
            double Theta_e = Math.PI * (13.0 / 18.0);
            C.PhysicalParameters.theta_e = Theta_e;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0, h };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            double alpha_l = k_l / (rho_l * c_l);
            double g = 9.81;
            double beta_l = 7.52e-4;
            double T_w = 393.12;
            double KCval = mu_l * alpha_l / (rho_l * g * beta_l * (T_w - T_sat));

            double dKC = 7.14 * Math.Pow((KCval), (1.0 / 3.0));

            Func<double[], double> TempFunc = delegate (double[] X) {

                double Temp = T_sat;

                double dT = 20.0;
                if (X[1] < dKC) {
                    Temp += (dT / dKC) * (dKC - X[1]);
                }

                return Temp;
            };

            C.InitialValues_Evaluators.Add("Temperature#A", TempFunc);
            C.InitialValues_Evaluators.Add("Temperature#B", TempFunc);

            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / R);

            C.InitialValues_Evaluators.Add("GravityY#A", (X => -g));
            C.InitialValues_Evaluators.Add("GravityY#B", (X => -g));


            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("navierslip_linear_ConstantTemperature_lower", "Temperature#A", (X, t) => T_w);
            C.AddBoundaryValue("navierslip_linear_ConstantTemperature_lower", "Temperature#B", (X, t) => T_w);

            C.AddBoundaryValue("pressure_Outlet_ZeroGradient_right");
            C.AddBoundaryValue("pressure_Outlet_ZeroGradient_upper");
            C.AddBoundaryValue("slipsymmetry_ZeroGradient_left");


            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            C.PhysicalParameters.sliplength = hmin / 2.0;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = true;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = true;
            C.BaseRefinementLevel = 3;
            C.AMR_startUpSweeps = 3;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-7;
            C.dtMin = 1e-7;
            C.Endtime = 1e-3;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }

        
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control Evaporation_Test(int p = 2, int kelemR = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            bool solveHeat = true;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; // C.DbPath != null;
            C.ProjectName = "XNSE/HeatedWall";
            C.ProjectDescription = "Leikonfiguration for SFB 1194";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.LogPeriod = 100;
            C.PostprocessingModules.Add(new MovingContactLineLogging() { LogPeriod = 100 });

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            C.solveCoupledHeatEquation = solveHeat;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 1.0;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 1.0;

            //C.ThermalParameters.prescribedVolumeFlux = 0.1;
            //C.PhysicalParameters.prescribedVolumeFlux = C.ThermalParameters.prescribedVolumeFlux;
            C.ThermalParameters.hVap = 1.0;
            //C.ThermalParameters.hVap_B = -1.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 1.0;
            double H = 3.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, R, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 3 * kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);


                grd.EdgeTagNames.Add(1, "velocity_inlet_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_ZeroGradient_upper");
                grd.EdgeTagNames.Add(3, "freeslip_ZeroGradient_left");
                grd.EdgeTagNames.Add(4, "freeslip_ZeroGradient_right");

                //grd.EdgeTagNames.Add(1, "wall_lower");
                //grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");

                //grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                //grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
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

            double h0 = 1.05;

            Func<double[], double> PhiFunc = (X => X[1] - h0);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            //C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.1);
            //C.InitialValues_Evaluators.Add("VelocityY#B", X => 1.0);

            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("velocity_inlet_ConstantTemperature_lower", "Temperature#A", (X, t) => 5.0);
            C.AddBoundaryValue("velocity_inlet_ConstantTemperature_lower", "VelocityY#A", (X, t) => 0.1);
            C.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            C.AddBoundaryValue("freeslip_ZeroGradient_left");
            C.AddBoundaryValue("freeslip_ZeroGradient_right");


            //C.AddBoundaryValue("wall_lower");
            //C.AddBoundaryValue("pressure_outlet_upper");
            //C.AddBoundaryValue("freeslip_left");
            //C.AddBoundaryValue("freeslip_right");


            //C.AddBoundaryValue("velocity_inlet_lower", "VelocityY#A", (X, t) => 0.1);
            //C.AddBoundaryValue("pressure_outlet_upper");
            //C.AddBoundaryValue("freeslip_left");
            //C.AddBoundaryValue("freeslip_right");


            //C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            //C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.hmin_Grid;


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.RefinementLevel = 1;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-1;
            C.dtMin = 1e-1;
            C.Endtime = 10000;
            C.NoOfTimesteps = 100;
            C.saveperiod = 1;

            #endregion


            return C;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control XHeatEquation_Test(int p = 2, int kelemR = 16, string _DbPath = null) {

            XNSE_Control C = new XNSE_Control();

            bool Xheat = false;
            bool separated = false;

            //_DbPath = @"\\dc1\userspace\smuda\cluster\CapillaryRise\CapillaryRise_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = false; // C.DbPath != null;
            C.ProjectName = "XNSE/XHeatEquation_Test";

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
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            C.solveCoupledHeatEquation = true;
            if (separated)
                C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.LDG;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 1.0;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            //C.ThermalParameters.T_sat = 8.0;

            //C.ThermalParameters.hVap_A = 1.0;
            //C.ThermalParameters.hVap_B = -1.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.ThermalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double R = 1.0;
            double H = 3.0;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, R, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 3 * kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true); //, periodicY: true);


                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_ZeroGradient_upper");

                //grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                //grd.EdgeTagNames.Add(2, "wall_ConstantTemperature_upper");
                //grd.EdgeTagNames.Add(1, "wall_ConstantHeatFlux_lower");
                //grd.EdgeTagNames.Add(2, "wall_ConstantHeatFlux_upper");

                //grd.EdgeTagNames.Add(3, "freeslip_ZeroGradient_left");
                //grd.EdgeTagNames.Add(4, "freeslip_ZeroGradient_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    //if (Math.Abs(X[0]) <= 1.0e-8)
                    //    et = 3;
                    //if (Math.Abs(X[0] - R) <= 1.0e-8)
                    //    et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double h0 = 1.05;

            Func<double[], double> PhiFunc = (X => -1.0);
            if (Xheat)
                PhiFunc = (X => h0 - X[1]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => -g);


            if (!Xheat) {
                //C.InitialValues_Evaluators.Add("Temperature#A", X => 10 - 3 * X[1]);
                C.InitialValues_Evaluators.Add("Temperature#A", X => Math.Sin((2 * Math.PI / H) * X[1]) + 10);

                //C.InitialValues_Evaluators.Add("HeatFluxY#A", X => 3);

            } else {
                C.InitialValues_Evaluators.Add("Temperature#A", X => 10 - 3 * h0);
                C.InitialValues_Evaluators.Add("Temperature#B", X => 10 - 3 * X[1]);
            }


            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            if (!Xheat)
                C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#A", (X, t) => 10.0);
            else
                C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#B", (X, t) => 10.0);

            //C.AddBoundaryValue("wall_ConstantTemperature_upper", "Temperature#A", (X, t) => 10.0);
            C.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");

            //C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#A", (X, t) => 10.0);
            //C.AddBoundaryValue("wall_ConstantTemperature_upper", "Temperature#A", (X, t) => 1);
            //C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#A", (X, t) => 3);
            //C.AddBoundaryValue("wall_ConstantHeatFlux_upper", "HeatFluxY#A", (X, t) => 3);

            //C.AddBoundaryValue("freeslip_ZeroGradient_left");
            //C.AddBoundaryValue("freeslip_ZeroGradient_right");


            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.RefinementLevel = 1;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 5e-1;
            C.dtMin = 5e-1;
            C.Endtime = 10000;
            C.NoOfTimesteps = 200;
            C.saveperiod = 1;

            #endregion


            return C;
        }


    }
}
