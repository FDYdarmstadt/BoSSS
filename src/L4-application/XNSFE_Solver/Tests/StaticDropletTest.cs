using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
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

        public static void StaticDropletConvergenceTest(int deg) {

            List<XNSFE_Control> LaLa = new List<XNSFE_Control>();
            foreach (int res in new int[] { 4, 8, 16, 32, 64, 128 }) {
                var C = StaticDropletTestControl(deg, res, true, 0.1);
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
            }
            //Assert.IsTrue(slope >= 1.7); //
        }

        public static XNSFE_Control StaticDropletTestControl(int pDeg, int kelem, bool heat, double Agg = 0.1) {

            XNSFE_Control C = new XNSFE_Control();

            //C.DbPath = set by workflowMgm during job creation
            C.savetodb = false;
            C.ContinueOnIoError = false;

            //C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.Dropletlike());
            //C.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MovingContactLineLogging());

            // misc. solver options
            // ====================
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            // Level-Set options (AMR)
            // =======================
            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;


            C.SetDGdegree(pDeg);

            
            C.GridFunc = () => {
                GridCommons grd;
                double[] xNodes = GenericBlas.Linspace(-3.0 / 2.0, 3.0 / 2.0, 2 * kelem + 1);
                double[] yNodes = GenericBlas.Linspace(0.0, 3.0 / 2.0, kelem + 1);
                grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                // HIER ANPASSUNG DER RANDBEDINGUNGEN!
                grd.EdgeTagNames.Add(1, "NavierSlip_linear_ConstantHeatflux_1");
                grd.EdgeTagNames.Add(2, "NavierSlip_linear_ConstantTemperature_2");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - yNodes.Last()) <= 1.0e-8) // UPPER BOUNDARY
                        et = 1;
                    if (Math.Abs(X[1] - yNodes.First()) <= 1.0e-8) // LOWER BOUNDARY
                        et = 1;
                    if (Math.Abs(X[0] - xNodes.First()) <= 1.0e-8) // LEFT BOUNDARY
                        et = 2;
                    if (Math.Abs(X[0] - xNodes.Last()) <= 1.0e-8) // RIGHT BOUNDARY
                        et = 2;
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
            C.PhysicalParameters.theta_e = Math.PI/2.0;

            C.ThermalParameters.IncludeConvection = false;
            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = 0.0; // 0.0 means no evaporation -  for reference case!
            C.ThermalParameters.T_sat = 0.0;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            //VERDAMPFUNG IN MEDIUM A
            C.AddBoundaryValue("NavierSlip_linear_ConstantHeatflux_1");
            if (heat) {
                C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_2", "Temperature#B", "X => X[0]", false);
            } else {
                C.AddBoundaryValue("NavierSlip_linear_ConstantTemperature_2");
            }

            double a = 0.816;
            double b = 0.784;
            double y_off = 0.0; // no offset for theta = 90°
            C.AddInitialValue("Phi", $"X => ((X[0]).Pow2() /({a}).Pow2() + (X[1]+({y_off})).Pow2() / ({b}).Pow2()) - 1.0", false);
            //C.AddInitialValue("Phi", $"X => X[0]*X[0]+X[1]-1.0", false);
            //double s = 0.0;
            //C.AddInitialValue("Phi", $"X => (0.747 - {s})*X[0]*X[0]*X[0]+{s}*X[0]+0.747-X[1]", false);
            C.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.LevelSetCGidx(1)].Degree = 4;
            C.AgglomerationThreshold = Agg;
            C.PlotAgglomeration = true;

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
