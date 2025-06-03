using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Control;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Foundation;
using BoSSS.Solution.Statistic;
using NUnit.Framework;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using static BoSSS.Solution.AMRLevelIndicatorLibrary;

namespace BoSSS.Application.XNSFE_Solver.Tests {

    class StaticDropletTest {

        public static void StaticDropletScalingTest(int deg) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 12, 16 }) {
                var C = StaticDropletTestControl(deg, res, false, 0.0);
                C.SkipSolveAndEvaluateResidual = false;
                LaLa.Add(C);                
            }

            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = true, title = "ScalingStaticDropletTest-p" + deg });
        }

        public static void StaticDropletConvergenceTest(int deg, bool evap, double ls, double theta) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 16, 32 }) {
                var C = StaticDropletTestControl(deg, res, true, 0.1, evap, ls, theta);
                C.SkipSolveAndEvaluateResidual = false;
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }

            Console.WriteLine("Running XNSFE Static Droplet Convergence Tests (no evaporation)");
            List<IEnumerable<DGField>> solutionOnDifferentResolutions = new List<IEnumerable<DGField>>();
            foreach (var C in LaLa) {
                using (var solver = new XNSFE()) {
                    solver.Init(C);
                    solver.RunSolverMode();
                    solutionOnDifferentResolutions.Add(solver.CurrentState.Fields);
                }
            }

            DGFieldComparison.ComputeErrors(
                solutionOnDifferentResolutions, out var hS, out var DOFs, out var errorS, NormType.L2_embedded);

            Console.WriteLine("Study : " + deg + "-" + (evap ? "evap1" : "evap0") + "-Slip" + (ls.ToString("N2")) + "-Angle" + (theta.ToString("N2")));
            foreach (var fieldName in solutionOnDifferentResolutions.First().Select(f => f.Identification)) {
                Dictionary<string, double[][]> dataGroups = new Dictionary<string, double[][]>();
                double[] xValues = hS;
                double[] yValues = errorS[fieldName];
                dataGroups.Add(deg.ToString(), new double[2][] { xValues, yValues });
                var plt = new Plot2Ddata(dataGroups.ToArray()).WithLogX().WithLogY();
                var reg = plt.Regression();

                Console.WriteLine(fieldName + ":");
                foreach (var kvp in reg) {
                    Console.WriteLine("Slope: " + kvp.Value);
                }

                plt.ModFormat();
                using (var gp = plt.ToGnuplot()) {
                    // set terminal
                    int xRes = 1024;
                    int yRes = 768;
                    gp.Terminal = string.Format("pngcairo size {0},{1}", xRes, yRes);

                    string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                    gp.OutputFile = fieldName + "-" + deg + "-" + (evap ? "evap1" : "evap0") + "-Slip" + (ls.ToString("N2")) + "-Angle" + (theta.ToString("N2")) + "-" + DateNtime + ".png";

                    // call gnuplot
                    int exCode = gp.RunAndExit(); // run & close gnuplot
                }
            }

            //Assert.IsTrue(slope >= 1.7); //
        }

        public static void StaticDropletConvergenceTestAMR(int deg, bool evap, double ls, double theta) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 12, 16, 24, 32 }) {
                var C = StaticDropletTestControl(deg, res, true, 0.1, evap, ls, theta);
                C.AdaptiveMeshRefinement= true;
                C.AMR_startUpSweeps = 4;
                C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = 2 });
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
                C.SkipSolveAndEvaluateResidual = false;
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }        

            Console.WriteLine("Running XNSFE Static Droplet Convergence Tests (no evaporation)");
            List<IEnumerable<DGField>> solutionOnDifferentResolutions = new List<IEnumerable<DGField>>();
            foreach (var C in LaLa) {
                using (var solver = new XNSFE()) {
                    solver.Init(C);
                    solver.RunSolverMode();
                    solutionOnDifferentResolutions.Add(solver.CurrentState.Fields);
                }
            }

            DGFieldComparison.ComputeErrors(
                solutionOnDifferentResolutions, out var hS, out var DOFs, out var errorS, NormType.L2_embedded);

            Console.WriteLine("Study : " + deg + "-" + (evap ? "evap1" : "evap0") + "-Slip" + (ls.ToString("N2")) + "-Angle" + (theta.ToString("N2")));
            foreach (var fieldName in solutionOnDifferentResolutions.First().Select(f => f.Identification)) {
                Dictionary<string, double[][]> dataGroups = new Dictionary<string, double[][]>();
                double[] xValues = hS;
                double[] yValues = errorS[fieldName];
                dataGroups.Add(deg.ToString(), new double[2][] { xValues, yValues });
                var plt = new Plot2Ddata(dataGroups.ToArray()).WithLogX().WithLogY();
                var reg = plt.Regression();

                Console.WriteLine(fieldName + ":");
                foreach (var kvp in reg) {
                    Console.WriteLine("Slope: " + kvp.Value);
                }

                plt.ModFormat();
                using (var gp = plt.ToGnuplot()) {
                    // set terminal
                    int xRes = 1024;
                    int yRes = 768;
                    gp.Terminal = string.Format("pngcairo size {0},{1}", xRes, yRes);

                    string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                    gp.OutputFile = fieldName + "-" + deg + "-" + (evap? "evap1" : "evap0") + "-Slip" + (ls.ToString("N2")) + "-Angle" + (theta.ToString("N2")) + "-" + DateNtime + ".png";

                    // call gnuplot
                    int exCode = gp.RunAndExit(); // run & close gnuplot
                }
            }

            //Assert.IsTrue(slope >= 1.7); //
        }

        public static XNSFE_Control StaticDropletTestControl(int pDeg, int kelem, bool heat, double Agg = 0.1, bool evap = false, double ls = 0.0, double theta = 90.0) {

            XNSFE_Control C = new XNSFE_Control();

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = false;
            C.ContinueOnIoError = false;

            C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.SlipDropletLogging());
            C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.Dropletlike());
            C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MovingContactLineLogging());

            // misc. solver options
            // ====================
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            // Level-Set options (AMR)
            // =======================
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.None;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.SetDGdegree(pDeg);


            double a = 0.816;
            double b = 0.784;
            double y_off = b / (Math.Sqrt(1 + (Math.Pow(a, 2) / Math.Pow(b, 2)) * Math.Pow((Math.Tan(theta/90.0*Math.PI/2.0)), 2)));
            double x0 = a * Math.Cos(Math.Asin(y_off / b));

            Dictionary<byte, string> Boundaries = new Dictionary<byte, string>();
            int study = 2;
            switch (study) {
                case 1:
                default: {
                        Boundaries = new Dictionary<byte, string>() {
                            { 1, "NavierSlip_linear_ConstantHeatflux_upper"},
                            { 2, "NavierSlip_linear_ConstantHeatflux_lower"},
                            { 3, "NavierSlip_linear_ConstantTemperature_left"},
                            { 4, "NavierSlip_linear_ConstantTemperature_right"},
                        };

                        C.AddBoundaryValue(Boundaries[1]);
                        C.AddBoundaryValue(Boundaries[2]);
                        C.AddBoundaryValue(Boundaries[3], "Temperature#B", "X => X[0]", false);
                        C.AddBoundaryValue(Boundaries[4], "Temperature#B", "X => X[0]", false);
                        break;
                    }
                case 2: {
                        Boundaries = new Dictionary<byte, string>() {
                            { 1, "NavierSlip_linear_ConstantTemperature_upper"},
                            { 2, "NavierSlip_linear_ConstantTemperature_lower"},
                            { 3, "NavierSlip_linear_ConstantTemperature_left"},
                            { 4, "NavierSlip_linear_ConstantTemperature_right"},
                        };

                        C.AddBoundaryValue(Boundaries[1]);
                        C.AddBoundaryValue(Boundaries[2], "Temperature#A", $"X => 1.0/8.0 * Math.Sin(2*Math.PI*(X[0]+{x0})/(2*{x0}))", false);
                        C.AddBoundaryValue(Boundaries[3]);
                        C.AddBoundaryValue(Boundaries[4]);
                        break;
                    }
                case 3: { 
                        Boundaries = new Dictionary<byte, string>() {
                            { 1, "NavierSlip_linear_ConstantHeatflux_upper"},
                            { 2, "NavierSlip_linear_TemperatureSlip_lower"},
                            { 3, "NavierSlip_linear_ConstantTemperature_left"},
                            { 4, "NavierSlip_linear_ConstantTemperature_right"},
                        };

                        C.AddBoundaryValue(Boundaries[1]);
                        C.AddBoundaryValue(Boundaries[2]);
                        C.AddBoundaryValue(Boundaries[3], "Temperature#B", "X => X[0]*(3*(X[1]/1.5).Pow(2) - 2*(X[1]/1.5).Pow(3))", false);
                        C.AddBoundaryValue(Boundaries[4], "Temperature#B", "X => X[0]*(3*(X[1]/1.5).Pow(2) - 2*(X[1]/1.5).Pow(3))", false);
                        break;
                    }
                case 4: {
                        Boundaries = new Dictionary<byte, string>() {
                            { 1, "NavierSlip_linear_ConstantHeatflux_upper"},
                            { 2, "NavierSlip_linear_TemperatureSlip_lower"},
                            { 3, "NavierSlip_linear_ConstantTemperature_left"},
                            { 4, "NavierSlip_linear_ConstantTemperature_right"},
                        };

                        C.AddBoundaryValue(Boundaries[1]);
                        C.AddBoundaryValue(Boundaries[2], "Temperature#A", "X => Math.Sin(2*Math.PI*X[0]/3.0)", false);
                        C.AddBoundaryValue(Boundaries[2], "Temperature#B", "X => Math.Sin(2*Math.PI*X[0]/3.0)", false);
                        C.AddBoundaryValue(Boundaries[3]);
                        C.AddBoundaryValue(Boundaries[4]);
                        break;
                    }
                case 5: {
                        Boundaries = new Dictionary<byte, string>() {
                            { 1, "NavierSlip_linear_ConstantTemperature_upper"},
                            { 2, "NavierSlip_linear_ConstantTemperature_lower"},
                            { 3, "NavierSlip_linear_ConstantTemperature_left"},
                            { 4, "NavierSlip_linear_ConstantTemperature_right"},
                        };

                        C.AddBoundaryValue(Boundaries[1], "Temperature#B", "X => X[0]*X[1]", false);
                        C.AddBoundaryValue(Boundaries[2]);
                        C.AddBoundaryValue(Boundaries[3], "Temperature#B", "X => X[0]*X[1]", false);
                        C.AddBoundaryValue(Boundaries[4], "Temperature#B", "X => X[0]*X[1]", false);
                        break;
                    }
            }

            C.GridFunc = () => {
                GridCommons grd;
                double[] xNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);
                double[] yNodes = GenericBlas.Linspace(0.0, 3.0 / 2.0, kelem + 1);
                grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                // HIER ANPASSUNG DER RANDBEDINGUNGEN!
                Boundaries.ForEach(kvp => grd.EdgeTagNames.Add(kvp.Key, kvp.Value));

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - yNodes.Last()) <= 1.0e-8) // UPPER BOUNDARY
                        et = 1;
                    if (Math.Abs(X[1] - yNodes.First()) <= 1.0e-8) // LOWER BOUNDARY
                        et = 2;
                    if (Math.Abs(X[0] - xNodes.First()) <= 1.0e-8) // LEFT BOUNDARY
                        et = 3;
                    if (Math.Abs(X[0] - xNodes.Last()) <= 1.0e-8) // RIGHT BOUNDARY
                        et = 4;
                    return et;
                });

                return grd;
            };

            // Physical Parameters
            // ===================
            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1; // unequal density
            C.PhysicalParameters.mu_A = 0.5;
            C.PhysicalParameters.mu_B = 0.05;
            C.PhysicalParameters.Sigma = 0.1;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;
            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = theta / 90.0 * Math.PI/2.0;//Math.PI/2.0;

            C.ThermalParameters.IncludeConvection = false;
            C.IncludeRecoilPressure= false;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = evap ? 1000.0 : double.NegativeInfinity; // 0.0 means no evaporation -  for reference case!
            C.ThermalParameters.T_sat = 0.0;

            C.ThermalParameters.sliplength = 0.1;
            C.PhysicalParameters.slipI = ls;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            if(C.LinearSolver is OrthoMGSchwarzConfig) {
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseKickIn = 90000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).TargetBlockSize = 10000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseUsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).UsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).ConvergenceCriterion = 1e-12;
            }

            //VERDAMPFUNG IN MEDIUM A

            

            //C.AddBoundaryValue("NavierSlip_linear_ConstantHeatflux_1", "HeatFluxY#A", "X => -0.5", false);
            //C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_1", "Temperature#A", $"X => -(X[0]-{x0})*(X[0]+{x0})/{x0}.Pow2()", false);
            //C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_1", "Temperature#A", $"X => 1.0/8.0 * Math.Sin(2*Math.PI*(X[0]+{x0})/(2*{x0}))", false);
            //C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_3");
            //C.AddBoundaryValue("NavierSlip_linear_TemperatureSlip_3", "Temperature#A", "X => 0.0", false);
            //C.AddBoundaryValue("NavierSlip_linear_TemperatureSlip_3", "Temperature#B", "X => 0.0", false);
            //if (heat) {
            //    //C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_2", "Temperature#B", "X => X[0]*(3*(X[1]/1.5).Pow(2) - 2*(X[1]/1.5).Pow(3))", false);
            //    C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_2", "Temperature#B", "X => 0.0", false);
            //} else {
            //    C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_2");
            //}

            C.AddInitialValue("Phi", $"X => ((X[0]).Pow2() /({a}).Pow2() + (X[1]+({y_off})).Pow2() / ({b}).Pow2()) - 1.0", false);
            //C.AddInitialValue("Phi", $"X => X[0] < {-x0} ? X[1] - ({Math.Tan(C.PhysicalParameters.theta_e)}*(X[0]+{x0}) + {y_off}): (X[0] > {x0} ? X[1] - ({-Math.Tan(C.PhysicalParameters.theta_e)}*(X[0]-{x0}) + {y_off}) : ((X[0]).Pow2() /({a}).Pow2() + (X[1]).Pow2() / ({b}).Pow2()) - 1.0)", false);
            //C.AddInitialValue("Phi", $"X => X[0]*X[0]+X[1]-1.0", false);
            //double s = 0.0;
            //C.AddInitialValue("Phi", $"X => (0.747 - {s})*X[0]*X[0]*X[0]+{s}*X[0]+0.747-X[1]", false);
            C.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.LevelSetCGidx(1)].Degree = 4;
            C.AgglomerationThreshold = Agg;
            //C.PlotAgglomeration = true;

            // MAKE TIMESTEPPING SETTINGS MORE CLEAR!
            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = BoSSS.Solution.Timestepping.TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            return C;
        }

        public static XNSFE_Control EvaporatigDropletTestControl(int pDeg, int kelem, bool heat, double Agg = 0.1, bool evap = false, double ls = 0.0, double theta = 90.0) {

            XNSFE_Control C = new XNSFE_Control();

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = false;
            C.ContinueOnIoError = false;

            //C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.Dropletlike());
            //C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MovingContactLineLogging());
            //C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.SlipDropletLogging());

            // misc. solver options
            // ====================
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            // Level-Set options (AMR)
            // =======================
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.SetDGdegree(pDeg);

            double r = 0.8;
            double a = r;
            double b = r;
            double y_off = b / (Math.Sqrt(1 + (Math.Pow(a, 2) / Math.Pow(b, 2)) * Math.Pow((Math.Tan(theta / 90.0 * Math.PI / 2.0)), 2)));
            double x0 = a * Math.Cos(Math.Asin(y_off / b));

            Dictionary<byte, string> Boundaries = new Dictionary<byte, string>() {
                            { 1, "Pressure_outlet_ConstantTemperature_upper"},
                            { 2, "Freeslip_TemperatureSlip_lower"},
                            { 3, "NavierSlip_linear_ConstantHeatflux_left"},
                            { 4, "NavierSlip_linear_ConstantHeatflux_right"},
                        };

            C.AddBoundaryValue(Boundaries[1]);
            C.AddBoundaryValue(Boundaries[2], "Temperature#A", "X => 1.0", false);
            C.AddBoundaryValue(Boundaries[2], "Temperature#B", "X => 1.0", false);
            C.AddBoundaryValue(Boundaries[3]);
            C.AddBoundaryValue(Boundaries[4]);

            C.GridFunc = () => {
                GridCommons grd;
                double[] xNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);
                double[] yNodes = GenericBlas.Linspace(0.0, 3.0 / 2.0, kelem + 1);
                grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                // HIER ANPASSUNG DER RANDBEDINGUNGEN!
                Boundaries.ForEach(kvp => grd.EdgeTagNames.Add(kvp.Key, kvp.Value));

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - yNodes.Last()) <= 1.0e-8) // UPPER BOUNDARY
                        et = 1;
                    if (Math.Abs(X[1] - yNodes.First()) <= 1.0e-8) // LOWER BOUNDARY
                        et = 2;
                    if (Math.Abs(X[0] - xNodes.First()) <= 1.0e-8) // LEFT BOUNDARY
                        et = 3;
                    if (Math.Abs(X[0] - xNodes.Last()) <= 1.0e-8) // RIGHT BOUNDARY
                        et = 4;
                    return et;
                });

                return grd;
            };

            // Physical Parameters
            // ===================
            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1; // unequal density
            C.PhysicalParameters.mu_A = 5.0;
            C.PhysicalParameters.mu_B = 5.0;
            C.PhysicalParameters.Sigma = 0.01;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;
            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = theta / 90.0 * Math.PI / 2.0;//Math.PI/2.0;

            C.ThermalParameters.IncludeConvection = false;
            C.IncludeRecoilPressure = false;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = evap ? 1000.0 : 0.0; // 0.0 means no evaporation -  for reference case!
            C.ThermalParameters.T_sat = 0.0;

            C.ThermalParameters.sliplength = 0.1;
            C.PhysicalParameters.slipI = ls;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            if (C.LinearSolver is OrthoMGSchwarzConfig) {
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseKickIn = 90000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).TargetBlockSize = 10000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseUsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).UsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).ConvergenceCriterion = 1e-12;
            }

            C.AMR_startUpSweeps = 5;
            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonBoundary(2) { maxRefinementLevel = 3 });
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
            C.activeAMRlevelIndicators.Add(new AMRatContactLine() { maxRefinementLevel = 5 });
            C.ReInitPeriod = 10;

            C.AddInitialValue("Phi", $"X => ((X[0]).Pow2() /({a}).Pow2() + (X[1]+({y_off})).Pow2() / ({b}).Pow2()) - 1.0", false);
            C.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.LevelSetCGidx(1)].Degree = 4;
            C.AgglomerationThreshold = Agg;
            //C.PlotAgglomeration = true;

            // MAKE TIMESTEPPING SETTINGS MORE CLEAR!
            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = BoSSS.Solution.Timestepping.TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting;

            C.dtFixed = 0.001;
            C.NoOfTimesteps= 1000;
            C.Endtime = C.dtFixed * C.NoOfTimesteps;

            return C;
        }

    }
}
