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
using System.IO;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Application.XNSFE_Solver.Tests {
    class InterfaceConvergenceTests {

        static void ConvergenceStudy(List<XNSFE_Control> LaLa) {

            Console.WriteLine("Running XNSFE Convergence Study");
            List<IEnumerable<DGField>> solutionOnDifferentResolutions = new List<IEnumerable<DGField>>();
            foreach (var C in LaLa) {
                using (var solver = new XNSFE()) {
                    solver.Init(C);
                    solver.RunSolverMode();
                    var fields = solver.CurrentState.Fields.CloneNonshallow().ToList();
                    fields.AddRange(solver.RegisteredFields.Where(f => f.Identification == VariableNames.LevelSetDG));
                    solutionOnDifferentResolutions.Add(fields);
                }
            }

            DGFieldComparison.ComputeErrors(
                solutionOnDifferentResolutions, out var hS, out var DOFs, out var errorS, NormType.L2_embedded);

            foreach (var fieldName in solutionOnDifferentResolutions.First().Select(f => f.Identification)) {
                int deg = solutionOnDifferentResolutions.First().Where(f => f.Identification == fieldName).First().Basis.Degree;
                Dictionary<string, double[][]> dataGroups = new Dictionary<string, double[][]>();
                double[] xValues = hS;
                double[] yValues = errorS[fieldName];
                dataGroups.Add(deg.ToString(), new double[2][] { xValues, yValues });
                var plt = new Plot2Ddata(dataGroups.ToArray()).WithLogX().WithLogY();
                var reg = plt.Regression();

                Console.Write(fieldName + " (" + deg + ") " + ": ");
                foreach (var kvp in reg) {
                    Console.WriteLine("Slope: " + kvp.Value);
                }

                //plt.ModFormat();
                //using (var gp = plt.ToGnuplot()) {
                //    // set terminal
                //    int xRes = 1024;
                //    int yRes = 768;
                //    gp.Terminal = string.Format("pngcairo size {0},{1}", xRes, yRes);

                //    string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                //    gp.OutputFile = fieldName + "-" + deg + "-Setup" + (setup) + "-" + DateNtime + ".png";

                //    // call gnuplot
                //    int exCode = gp.RunAndExit(); // run & close gnuplot
                //}
            }
            Console.WriteLine();

            //Assert.IsTrue(slope >= 1.7); //
        }

        /// <summary>
        /// Test convergence of temperature equation in different setups.
        /// No convection.
        /// Full XNSFE, fixed interface temperature.
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="setup"></param>
        /// 0: Sin based dirichlet b.c.
        /// 1: b.c. based on exact solution, see "TemperatureConvergenceControl.m"
        /// <param name="theta">contact angle in deg. Only valid choices are 80 or 90</param>
        public static void TemperatureConvergence(int deg, int setup, double theta = 90.0) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 16, 32 }) {
                var C = TemperatureConvergenceControl(deg, res, theta, setup);
                C.SkipSolveAndEvaluateResidual = false;
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }

            ConvergenceStudy(LaLa);
        }

        
        /// <summary>
        /// Test convergence of parasitic currents.
        /// SDF circle
        /// Linear
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="lsdeg"></param>
        /// <param name="exponent"></param>
        public static void CurvatureConvergence(int deg, int lsdeg, double exponent) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 16, 32 }) {
                var C = CurvatureConvergenceControl(deg, lsdeg, res, exponent);
                C.SkipSolveAndEvaluateResidual = false;
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }

            ConvergenceStudy(LaLa);
        }


        
        /// <summary>
        /// Test convergence of evaporative flow.
        /// Quadratic free floating circle
        /// Linear (no convection)
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="tdeg">degree for temperature, seperate if specified > 0</param>
        public static void EvaporationConvergence(int deg, int tdeg) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 16, 32 }) {
                var C = EvaporationConvergenceControl(deg, tdeg, res);
                C.SkipSolveAndEvaluateResidual = false;
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }

            ConvergenceStudy(LaLa);
        }
        

        /// <summary>
        /// Used to compare the temperature convergence orders for circle with 80/90° contactangle.
        /// </summary>
        /// <param name="pDeg"></param>
        /// <param name="kelem"></param>
        /// <param name="theta"></param>
        /// <param name="setup"></param>
        /// 0 : temperature dirichlet value prescribed by a simple sine function on the lower boundary
        /// 1 : temperature dirichlet value prescribed by a spline based on an exact solution to the Poisson problem,
        /// using an alternating smooth kernel on the not resolved lower half of the circle, see TemperatureConvergenceControl.m for details
        public static XNSFE_Control TemperatureConvergenceControl(int pDeg, int kelem, double theta, int setup = 1) {

            XNSFE_Control C = new XNSFE_Control();

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = false;
            C.ContinueOnIoError = false;

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

            double r = 0.8;
            double y0 = r * Math.Cos(theta / 90.0 * Math.PI / 2.0);
            double x0 = r * Math.Sin(theta / 90.0 * Math.PI / 2.0);

            Dictionary<byte, string> Boundaries = new Dictionary<byte, string>();
            int study = setup;
            Boundaries = new Dictionary<byte, string>() {
                { 1, "NavierSlip_linear_ConstantTemperature_upper"},
                { 2, "NavierSlip_linear_ConstantTemperature_lower"},
                { 3, "NavierSlip_linear_ConstantTemperature_left"},
                { 4, "NavierSlip_linear_ConstantTemperature_right"},
            };

            C.AddBoundaryValue(Boundaries[1]);
            C.AddBoundaryValue(Boundaries[2]);
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
            C.PhysicalParameters.mu_A = 0.5;
            C.PhysicalParameters.mu_B = 0.05;
            C.PhysicalParameters.Sigma = 0.1;
            C.PhysicalParameters.theta_e = theta / 90.0 * Math.PI / 2.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;
            C.PhysicalParameters.betaL = 0.0;

            C.ThermalParameters.IncludeConvection = false;
            C.IncludeRecoilPressure = false;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = double.NegativeInfinity; // 0.0 means no evaporation -  for reference case!
            C.ThermalParameters.T_sat = 0.0;

            C.ThermalParameters.sliplength = 0.0;
            C.PhysicalParameters.slipI = 0.0;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NoOfMultigridLevels = 1;
            if (C.LinearSolver is OrthoMGSchwarzConfig) {
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseKickIn = 90000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).TargetBlockSize = 10000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseUsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).UsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).ConvergenceCriterion = 1e-12;
            }

            C.AddInitialValue("Phi", $"X => (X[0]).Pow2() /({r}).Pow2() + (X[1]+{y0}).Pow2() / ({r}).Pow2() - 1.0", false);

            string[] lines;
            if(theta == 90.0) {
                lines = File.ReadAllLines("../../../Tests/TemperatureConvergenceControl_boundarycond90.txt");
            } else if(theta == 80.0) {
                lines = File.ReadAllLines("../../../Tests/TemperatureConvergenceControl_boundarycond80.txt");
            } else {
                throw new NotImplementedException();
            }

            double[] nodes = new double[lines.Length];
            double[] values = new double[lines.Length];
            for (int i = 0; i < lines.Length; i++) {
                string[] xy = lines[i].Split('\t');
                nodes[i] = Convert.ToDouble(xy[0]);
                values[i] = Convert.ToDouble(xy[1]);
            }
            Spline1D TFunc = new Spline1D(nodes, values, 0, Spline1D.OutOfBoundsBehave.Clamp);

            switch (study) {
                case 0:
                    C.AddBoundaryValue(Boundaries[2], "Temperature#A", X => Math.Sin(2 * Math.PI * (X[0] - x0) / (2 * x0)));
                    break;
                case 1:
                    C.AddBoundaryValue(Boundaries[2], "Temperature#A", TFunc);
                    break;
                default:
                    throw new NotImplementedException();
            }

            // MAKE TIMESTEPPING SETTINGS MORE CLEAR!
            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = BoSSS.Solution.Timestepping.TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            return C;
        }

        /// <summary>
        /// Used to investigate convergence for pressure/velocity in dependence of the temperature field
        /// </summary>
        /// <param name="pDeg"></param>
        /// <param name="kelem"></param>
        public static XNSFE_Control EvaporationConvergenceControl(int pDeg, int tDeg, int kelem) {

            XNSFE_Control C = new XNSFE_Control();

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = false;
            C.ContinueOnIoError = false;

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
            C.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Temperature].Degree = tDeg;

            Dictionary<byte, string> Boundaries = new Dictionary<byte, string>();
            Boundaries = new Dictionary<byte, string>() {
                { 1, "NavierSlip_linear_ConstantTemperature_upper"},
                { 2, "NavierSlip_linear_ConstantTemperature_lower"},
                { 3, "NavierSlip_linear_ConstantTemperature_left"},
                { 4, "NavierSlip_linear_ConstantTemperature_right"},
            };

            //Func<double[], double> TFunc = X => X[0];
            //Func<double[], double> TFunc = X => X[0]*X[1];
            Func<double[], double> TFunc = X => Math.Sin(Math.PI * X[0] / 3.0) * Math.Sin(Math.PI * X[1] / 3.0);

            C.AddBoundaryValue(Boundaries[1], "Temperature#B", TFunc);
            C.AddBoundaryValue(Boundaries[2], "Temperature#B", TFunc);
            C.AddBoundaryValue(Boundaries[3], "Temperature#B", TFunc);
            C.AddBoundaryValue(Boundaries[4], "Temperature#B", TFunc);

            C.GridFunc = () => {
                GridCommons grd;
                double[] xNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);
                double[] yNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);

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
            C.PhysicalParameters.Sigma = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;
            C.PhysicalParameters.betaL = 0.0;

            C.ThermalParameters.IncludeConvection = false;
            C.IncludeRecoilPressure = false;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = 1000.0; // 0.0 means no evaporation -  for reference case!
            C.ThermalParameters.T_sat = 0.0;

            C.ThermalParameters.sliplength = 0.0;
            C.PhysicalParameters.slipI = 0.0;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NoOfMultigridLevels = 1;
            if (C.LinearSolver is OrthoMGSchwarzConfig) {
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseKickIn = 90000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).TargetBlockSize = 10000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseUsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).UsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).ConvergenceCriterion = 1e-12;
            }

            double r = 0.8;
            C.AddInitialValue("Phi", $"X => (X[0]).Pow2() /({r}).Pow2() + (X[1]).Pow2() / ({r}).Pow2() - 1.0", false);

            // MAKE TIMESTEPPING SETTINGS MORE CLEAR!
            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = BoSSS.Solution.Timestepping.TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            return C;
        }

        /// <summary>
        /// Used to investigate convergence for pressure/velocity for parasitic currents with a circular interface
        /// </summary>
        /// <param name="pDeg"></param>
        /// <param name="kelem"></param>
        public static XNSFE_Control CurvatureConvergenceControl(int pDeg, int lsDeg, int kelem, double exponent) {

            XNSFE_Control C = new XNSFE_Control();

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = false;
            C.ContinueOnIoError = false;

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


            Dictionary<byte, string> Boundaries = new Dictionary<byte, string>();
            Boundaries = new Dictionary<byte, string>() {
                { 1, "NavierSlip_linear_ConstantTemperature_upper"},
                { 2, "NavierSlip_linear_ConstantTemperature_lower"},
                { 3, "NavierSlip_linear_ConstantTemperature_left"},
                { 4, "NavierSlip_linear_ConstantTemperature_right"},
            };

            C.AddBoundaryValue(Boundaries[1]);
            C.AddBoundaryValue(Boundaries[2]);
            C.AddBoundaryValue(Boundaries[3]);
            C.AddBoundaryValue(Boundaries[4]);

            C.GridFunc = () => {
                GridCommons grd;
                double[] xNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);
                double[] yNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);

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
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;
            C.PhysicalParameters.betaL = 0.0;

            C.ThermalParameters.IncludeConvection = false;
            C.IncludeRecoilPressure = false;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = double.NegativeInfinity; // 0.0 means no evaporation -  for reference case!
            C.ThermalParameters.T_sat = 0.0;

            C.ThermalParameters.sliplength = 0.0;
            C.PhysicalParameters.slipI = 0.0;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NoOfMultigridLevels = 1;
            if (C.LinearSolver is OrthoMGSchwarzConfig) {
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseKickIn = 90000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).TargetBlockSize = 10000;
                ((OrthoMGSchwarzConfig)C.LinearSolver).CoarseUsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).UsepTG = false;
                ((OrthoMGSchwarzConfig)C.LinearSolver).ConvergenceCriterion = 1e-12;
            }

            double r = 0.8;
            C.AddInitialValue("Phi", $"X => Math.Pow((X[0]).Pow2() /({r}).Pow2() + (X[1]).Pow2() / ({r}).Pow2(), {exponent}) - 1.0", false);
            C.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.LevelSetCG].Degree = lsDeg;

            // MAKE TIMESTEPPING SETTINGS MORE CLEAR!
            // Timestepping
            // ============
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = BoSSS.Solution.Timestepping.TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            return C;
        }
    }
}
