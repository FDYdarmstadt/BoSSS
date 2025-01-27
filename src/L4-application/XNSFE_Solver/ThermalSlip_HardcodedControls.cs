using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.Control;
using ilPSP.Connectors.Matlab;
using ilPSP;

namespace BoSSS.Application.XNSFE_Solver {
    public static class ThermalSlip_HardcodedControls {

        public static XNSFE_Control SinglePhaseSlip() {
            // --control "cs:BoSSS.Application.XNSFE_Solver.ThermalSlip_HardcodedControls.SinglePhaseSlip()"

            var ctrl = new XNSFE_Control();

            ctrl.SetDGdegree(2);

            double L = 1;

            ctrl.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.3 * 1e-4,
                mu_B = 1.3 * 1e-6,

                Sigma = L,
                theta_e = 10.0 / 90.0 * Math.PI / 2.0,
                sliplength = 0.0
            };

            double Tr = 0.99;
            ctrl.ThermalParameters = new ThermalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,
                rho_C = 1.0,

                c_A = 1.0,
                c_B = 1.0,
                c_C = 1.0,

                k_A = 1.0,
                k_B = 1.0,
                k_C = 1.0,

                sliplength = -2.0 * (Tr - 1)/(Tr+1),

                hVap = 0.0
            };
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;

            ctrl.GridFunc = () => {
                var Nodes = GenericBlas.Linspace(-L, L, 5);
                var grd = Grid2D.Cartesian2DGrid(Nodes, Nodes);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_left");
                grd.EdgeTagNames.Add(2, "wall_TemperatureSlip_right");
                grd.EdgeTagNames.Add(3, "pressure_outlet_ZeroGradient_upper");
                grd.EdgeTagNames.Add(4, "pressure_outlet_ZeroGradient_lower");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Nodes.First()) <= 1.0e-12)
                        et = 1;
                    if (Math.Abs(X[0] - Nodes.Last()) <= 1.0e-12)
                        et = 2;
                    if (Math.Abs(X[1] - Nodes.First()) <= 1.0e-12)
                        et = 3;
                    if (Math.Abs(X[1] - Nodes.Last()) <= 1.0e-12)
                        et = 4;

                    return et;
                });

                return grd;
            };

            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", $"(X, t) => -1", true);
            ctrl.AddBoundaryValue("wall_TemperatureSlip_right", "Temperature#A", $"(X, t) => 1", true);
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.AddInitialValue("Temperature#A", $"(X, t) => (X[0])", true);

            ctrl.SkipSolveAndEvaluateResidual = false;

            ctrl.AddInitialValue("Phi", "X => -1.0", false);

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            return ctrl;
        }

        public static XNSFE_Control TwoPhaseSlip() {
            // --control "cs:BoSSS.Application.XNSFE_Solver.ThermalSlip_HardcodedControls.TwoPhaseSlip()"

            var ctrl = new XNSFE_Control();

            ctrl.SetDGdegree(2);

            double L = 1;

            ctrl.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.3 * 1e-4,
                mu_B = 1.3 * 1e-6,

                Sigma = L,
                theta_e = 10.0 / 90.0 * Math.PI / 2.0,
                sliplength = 0.0
            };

            ctrl.ThermalParameters = new ThermalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,
                rho_C = 1.0,

                c_A = 1.0,
                c_B = 1.0,
                c_C = 1.0,

                k_A = 1.0,
                k_B = 1.0,
                k_C = 1.0,

                sliplength = 0.001,
                T_sat = 0.0,

                hVap = 1.0
            };
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;
            ctrl.IncludeRecoilPressure = false;

            ctrl.GridFunc = () => {
                var Nodes = GenericBlas.Linspace(-L, L, 5);
                var grd = Grid2D.Cartesian2DGrid(Nodes, Nodes);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_left");
                grd.EdgeTagNames.Add(2, "wall_TemperatureSlip_right");
                grd.EdgeTagNames.Add(3, "pressure_outlet_ZeroGradient_upper");
                grd.EdgeTagNames.Add(4, "pressure_outlet_ZeroGradient_lower");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Nodes.First()) <= 1.0e-12)
                        et = 1;
                    if (Math.Abs(X[0] - Nodes.Last()) <= 1.0e-12)
                        et = 2;
                    if (Math.Abs(X[1] - Nodes.First()) <= 1.0e-12)
                        et = 3;
                    if (Math.Abs(X[1] - Nodes.Last()) <= 1.0e-12)
                        et = 4;

                    return et;
                });

                return grd;
            };

            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#B", $"(X, t) => 0", true);
            ctrl.AddBoundaryValue("wall_TemperatureSlip_right", "Temperature#A", $"(X, t) => 1", true);
            ctrl.AddBoundaryValue("wall_TemperatureSlip_right", "Temperature#B", $"(X, t) => 0", true);

            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.AddInitialValue("Temperature#A", $"(X, t) => (X[0])", true);

            ctrl.SkipSolveAndEvaluateResidual = false;

            double offset = 4.9 * L / 16.0;
            double ContactAngle = ctrl.PhysicalParameters.theta_e;
            ctrl.AddInitialValue("Phi", $"(X, t) => X[1] - {offset} - (X[0]-1)/Math.Tan({ContactAngle})", true);

            int level = 6;
            ctrl.AdaptiveMeshRefinement = level > 0;
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.AMRLevelIndicatorLibrary.AMRonBoundary(2) { maxRefinementLevel = level });
            ctrl.AMR_startUpSweeps = level;

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            return ctrl;
        }
        public static XNSFE_Control IBMSlip() {
            // --control "cs:BoSSS.Application.XNSFE_Solver.ThermalSlip_HardcodedControls.IBMSlip()"

            var ctrl = new XNSFE_Control();

            ctrl.SetDGdegree(2);

            double L = 1;

            ctrl.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.3 * 1e-4,
                mu_B = 1.3 * 1e-6,

                Sigma = L,
                theta_e = 10.0 / 90.0 * Math.PI / 2.0,
                sliplength = 0.0
            };

            ctrl.ThermalParameters = new ThermalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,
                rho_C = 1.0,

                c_A = 1.0,
                c_B = 1.0,
                c_C = 1.0,

                k_A = 1.0,
                k_B = 1.0,
                k_C = 1000.0,

                sliplength = 1.0,
                T_sat = 0.0,

                hVap = 1.0
            };
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;
            ctrl.IncludeRecoilPressure = false;
            ctrl.AdvancedDiscretizationOptions.IBM_ThermalBoundaryType = IBM_ThermalBoundaryType.ThermalSlip;

            ctrl.UseImmersedBoundary = true;

            ctrl.GridFunc = () => {
                var Nodes = GenericBlas.Linspace(-L, L, 5);
                var grd = Grid2D.Cartesian2DGrid(Nodes, Nodes);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_left");
                grd.EdgeTagNames.Add(2, "wall_ConstantTemperature_right");
                grd.EdgeTagNames.Add(3, "pressure_outlet_ZeroGradient_upper");
                grd.EdgeTagNames.Add(4, "pressure_outlet_ZeroGradient_lower");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Nodes.First()) <= 1.0e-12)
                        et = 1;
                    if (Math.Abs(X[0] - Nodes.Last()) <= 1.0e-12)
                        et = 2;
                    if (Math.Abs(X[1] - Nodes.First()) <= 1.0e-12)
                        et = 3;
                    if (Math.Abs(X[1] - Nodes.Last()) <= 1.0e-12)
                        et = 4;

                    return et;
                });

                return grd;
            };

            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", $"(X, t) => -1", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#C", $"(X, t) => 1", true);

            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.SkipSolveAndEvaluateResidual = false;

            ctrl.AddInitialValue("Phi", $"(X, t) => -1.0", true);
            ctrl.AddInitialValue("Phi2", $"(X, t) => X[0]", true);

            int level = 0;
            ctrl.AdaptiveMeshRefinement = level > 0;
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRonNarrowband() { levelSet = 1, maxRefinementLevel = level });
            ctrl.AMR_startUpSweeps = level;

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            return ctrl;
        }

        public static XNSFE_Control ThreePhaseSlip() {
            // --control "cs:BoSSS.Application.XNSFE_Solver.ThermalSlip_HardcodedControls.ThreePhaseSlip()"

            var ctrl = new XNSFE_Control();

            ctrl.DbPath = @"\\dc1\userspace\rieckmann\cluster\TemperatureSlip";
            ctrl.savetodb = false;

            ctrl.SessionName = "TemperatureSlipShowcase";
            ctrl.ProjectName = "TemperatureSlip";

            ctrl.SetDGdegree(2);

            double L = 1;

            ctrl.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.3 * 1e-4,
                mu_B = 1.3 * 1e-6,

                Sigma = L,
                theta_e = 10.0 / 90.0 * Math.PI / 2.0,
                sliplength = 0.0
            };

            double slip = 0.0;
            ctrl.ThermalParameters = new ThermalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,
                rho_C = 1.0,

                c_A = 1.0,
                c_B = 1.0,
                c_C = 1.0,

                k_A = 1.0,
                k_B = 1.0,
                k_C = 1000.0,

                sliplength = slip,
                T_sat = 0.0,

                hVap = 1.0
            };

            ctrl.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("TSlip", slip));

            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;
            ctrl.IncludeRecoilPressure = false;
            ctrl.AdvancedDiscretizationOptions.IBM_ThermalBoundaryType = IBM_ThermalBoundaryType.ThermalSlip;

            ctrl.UseImmersedBoundary = true;
            ctrl.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = true;

            ctrl.GridFunc = () => {
                var Nodes = GenericBlas.Linspace(-L, 1.5*L, 6);
                var grd = Grid2D.Cartesian2DGrid(Nodes, Nodes);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_left");
                grd.EdgeTagNames.Add(2, "wall_ConstantTemperature_right");
                grd.EdgeTagNames.Add(3, "pressure_outlet_ZeroGradient_upper");
                grd.EdgeTagNames.Add(4, "pressure_outlet_ZeroGradient_lower");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Nodes.First()) <= 1.0e-12)
                        et = 1;
                    if (Math.Abs(X[0] - Nodes.Last()) <= 1.0e-12)
                        et = 2;
                    if (Math.Abs(X[1] - Nodes.First()) <= 1.0e-12)
                        et = 3;
                    if (Math.Abs(X[1] - Nodes.Last()) <= 1.0e-12)
                        et = 4;

                    return et;
                });

                return grd;
            };

            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#B", $"(X, t) => 0", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#C", $"(X, t) => 1", true);

            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.SkipSolveAndEvaluateResidual = false;

            double offset = 4.9 * L / 16.0;
            double ContactAngle = ctrl.PhysicalParameters.theta_e;
            //ctrl.AddInitialValue("Phi", $"(X, t) => 1.0", true);
            ctrl.AddInitialValue("Phi", $"(X, t) => X[1] - {offset} - (X[0]-1)/Math.Tan({ContactAngle})", true);

            ctrl.AddInitialValue("Phi2", $"(X, t) => X[0] - 1.0", true);


            int level = 4;
            ctrl.AdaptiveMeshRefinement = level > 0;
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRonNarrowband() { levelSet = 1, maxRefinementLevel = level });
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRatInnerContactLine() { maxRefinementLevel = level*2 });
            ctrl.AMR_startUpSweeps = level*2;

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
            ctrl.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.RK_ImplicitEuler;

            return ctrl;
        }

        public static XNSFE_Control HeatedWallSimple_3Phase() {
            // --control "cs:BoSSS.Application.XNSFE_Solver.ThermalSlip_HardcodedControls.HeatedWallSimple_3Phase()"

            var ctrl = new XNSFE_Control();

            ctrl.DbPath = null;
            ctrl.SessionName = $"HeatedWall_Simple";
            ctrl.ProjectName = $"HeatedWall_Simple";
            ctrl.savetodb = false;

            ctrl.SetDGdegree(2);

            #region grid
            double L = 5.0;
            int kelemR = 1;
            string[] Bndy = new string[] {  "Inner",
                                        "wall_ZeroGradient_right",
                                        "pressure_outlet_ZeroGradient_top",
                                        "freeslip_ZeroGradient_left",
                                        "pressure_outlet_ZeroGradient_bottom"};

            ctrl.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L, 1, 6 * kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 3 * L, 5 * 3 * kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                for (byte i = 1; i < Bndy.Count(); i++) {
                    grd.EdgeTagNames.Add(i, Bndy[i]);
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Xnodes.Last()) < 1e-8)
                        return 1;
                    if (Math.Abs(X[0] - Xnodes.First()) < 1e-8)
                        return 3;
                    if (Math.Abs(X[1] - Ynodes.Last()) < 1e-8)
                        return 2;
                    if (Math.Abs(X[1] - Ynodes.First()) < 1e-8)
                        return 4;
                    return et;
                });

                return grd;
            };
            #endregion

            #region material
            ctrl.PhysicalParameters = new BoSSS.Solution.XNSECommon.PhysicalParameters() {
                rho_A = 1.0, // 958.0
                rho_B = 1.0, // 0.59,

                mu_A = 1, //2.82 * 1e-4,
                mu_B = 0.001, //1.23 * 1e-6,

                Sigma = 1.0,
                betaS_A = 1000, // sliplength is mu/beta
                betaS_B = 1000,
            };

            ctrl.ThermalParameters = new BoSSS.Solution.XheatCommon.ThermalParameters() {
                rho_A = 1.0, // 958.0
                rho_B = 1.0, //0.59,
                rho_C= 1.0,

                k_A = 1.0, // 0.6
                k_B = 1.0, // 0.026,
                k_C= 1.0,

                c_A = 0.0,
                c_B = 0.0,
                c_C= 0.0,

                hVap = 1,//2.257 * 1e6,
                T_sat = 0.0 // 373.0
            };

            ctrl.PhysicalParameters.IncludeConvection = true;
            ctrl.ThermalParameters.IncludeConvection = true;
            ctrl.PhysicalParameters.Material = false;
            #endregion

            #region Initial Condition - Exact Solution

            // solution for massflux and velocity at level set
            double y0 = 0.2 * L;

            // inital values
            double g = 4;
            ctrl.AddInitialValue("Phi", $"(X, t) => -{y0} + X[1]", true);
            ctrl.UseImmersedBoundary = true;
            ctrl.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = true;
            ctrl.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NavierSlip;
            ctrl.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Everywhere;
            ctrl.AddInitialValue("Phi2", $"(X, t) => X[0]", true);
            ctrl.AddInitialValue("Temperature#A", $"(X, t) => {ctrl.ThermalParameters.T_sat}", true);
            ctrl.AddInitialValue("Temperature#B", $"(X, t) => {ctrl.ThermalParameters.T_sat}", true);
            ctrl.AddInitialValue("Temperature#C", $"(X, t) => {ctrl.ThermalParameters.T_sat}", true);
            ctrl.AddInitialValue("VelocityX@Phi2", $"(X, t) => 0", true);
            ctrl.AddInitialValue("VelocityY@Phi2", $"(X, t) => 1", true);
            ctrl.AddInitialValue("GravityY#A", $"(X, t) => -{g}", true);

            ctrl.HeatSourceIBM = new Dictionary<string, string>();
            ctrl.HeatSourceIBM["AC"] = "(X,t) => 0.2";

            #endregion

            #region Boundary Conditions

            ctrl.AddBoundaryValue(Bndy[1]);
            ctrl.AddBoundaryValue(Bndy[3]);
            ctrl.AddBoundaryValue(Bndy[2]);
            ctrl.AddBoundaryValue(Bndy[4], "Pressure#A", $"(X, t) => {y0} * {ctrl.PhysicalParameters.rho_A} * {g}", true);

            #endregion

            #region AMR

            int level = 3;
            ctrl.AdaptiveMeshRefinement = level > 0;
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRonNarrowband() { levelSet=0,maxRefinementLevel = level });
            ctrl.AMR_startUpSweeps = level;

            #endregion

            #region Timestepping

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting;

            ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            ctrl.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1e-8;
            ctrl.NonLinearSolver.MaxSolverIterations = 10;

            ctrl.SkipSolveAndEvaluateResidual = false;

            ctrl.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 0.01;
            ctrl.Endtime = 15.0;
            ctrl.NoOfTimesteps = (int)(ctrl.Endtime / ctrl.dtFixed);

            #endregion
            //ctrl.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MassfluxLogging() { LogPeriod = 1 });
            ctrl.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MovingContactLineZwoLsLogging() { LogPeriod = 1 });

            ctrl.ImmediatePlotPeriod = 1;
            ctrl.SuperSampling = 2;

            return ctrl;
        }

        public static XNSFE_Control HeatedWall_3PhaseDemo(bool threephase = false) {
            // --control "cs:BoSSS.Application.XNSFE_Solver.ThermalSlip_HardcodedControls.HeatedWall_3PhaseDemo()"

            var ctrl = new XNSFE_Control();

            ctrl.DbPath = @"\\dc3\backup\rieckmann\cluster\databases\HeatedWall_Simple_3Phase";
            ctrl.SessionName = $"HeatedWall_"+(threephase ? "3" : "2")+"Phase_Demo";
            ctrl.ProjectName = $"HeatedWall_3Phase_Demo";
            ctrl.savetodb = false;

            ctrl.SetDGdegree(2);

            #region grid
            double L = 1.0;
            int kelemR = 1;
            string[] Bndy = new string[] {  "Inner",
                                        "NavierSlip_linear_ConstantHeatFlux_right",
                                        "pressure_outlet_ZeroGradient_top",
                                        "freeslip_ZeroGradient_left",
                                        "pressure_outlet_ZeroGradient_bottom"};

            ctrl.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-5*L, threephase ? L : 0, (threephase ? 6 : 5) * kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 5*L, 5 * kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                for (byte i = 1; i < Bndy.Count(); i++) {
                    grd.EdgeTagNames.Add(i, Bndy[i]);
                }

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Xnodes.Last()) < 1e-8)
                        return 1;
                    if (Math.Abs(X[0] - Xnodes.First()) < 1e-8)
                        return 3;
                    if (Math.Abs(X[1] - Ynodes.Last()) < 1e-8)
                        return 2;
                    if (Math.Abs(X[1] - Ynodes.First()) < 1e-8)
                        return 4;
                    return et;
                });

                return grd;
            };
            #endregion

            #region material
            ctrl.PhysicalParameters = new BoSSS.Solution.XNSECommon.PhysicalParameters() {
                rho_A = 1.0, // 958.0
                rho_B = 1.0, // 0.59,

                mu_A = 1, //2.82 * 1e-4,
                mu_B = 0.001, //1.23 * 1e-6,

                theta_e = Math.PI / 180.0 * 90.0,
                Sigma = 1.0,
                betaS_A = 1000, // sliplength is mu/beta
                betaS_B = 1000,
            };

            ctrl.ThermalParameters = new BoSSS.Solution.XheatCommon.ThermalParameters() {
                rho_A = 1.0, // 958.0
                rho_B = 1.0, //0.59,
                rho_C = 1.0,

                k_A = 1.0, // 0.6
                k_B = 1.0, // 0.026,
                k_C = 10.0,

                c_A = 0.0,
                c_B = 0.0,
                c_C = 0.0,

                sliplength = 0.001,

                hVap = 1,//2.257 * 1e6,
                T_sat = 0.0 // 373.0
            };
            
            ctrl.AdvancedDiscretizationOptions.IBM_ThermalBoundaryType = IBM_ThermalBoundaryType.ThermalSlip;
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;
            ctrl.IncludeRecoilPressure= false;
            ctrl.PhysicalParameters.Material = false;
            #endregion

            #region Initial Condition - Exact Solution

            // solution for massflux and velocity at level set
            double y0 = 2.2 * L;

            // inital values
            double g = 4;

            double a = Math.Sqrt(ctrl.PhysicalParameters.Sigma / (ctrl.PhysicalParameters.rho_A * g));
            double h = Math.Sqrt(2*a*a*(1-Math.Sin(ctrl.PhysicalParameters.theta_e)));

            //Func<double, double> MeniskusInitRat55(double dx, double h, double angle) {
            //    Console.WriteLine("starting Batchmode Connector");
            //    Func<double, double> Phi, PhiR;
            //    MultidimensionalArray coefficients = MultidimensionalArray.Create(6, 2);
            //    using (BatchmodeConnector matlab = new BatchmodeConnector()) {
            //        //note: BatchmodeCon maybe working on proc0 but savetotxt file, etc. (I/O) is full mpi parallel
            //        //so concider this as full mpi-parallel
            //        matlab.PutVector(new double[] { dx }, "a");
            //        matlab.PutVector(new double[] { h }, "h");
            //        matlab.Cmd(String.Format("fp = fimplicit(@(y,z) -y./a + acosh(2*a./z) - acosh(2*a/h)+sqrt(4-(h/a)^2) - sqrt(4-(z./a).^2), [0 5*a 0 h],'Meshdensity', 1000);"));
            //        matlab.Cmd(String.Format("x=fp.XData;"));
            //        matlab.Cmd(String.Format("y=fp.YData;"));
            //        matlab.Cmd(String.Format("[xData, yData] = prepareCurveData(x, y);"));
            //        matlab.Cmd(String.Format("ft = fittype('rat55');"));
            //        matlab.Cmd(String.Format("opts = fitoptions('Method', 'NonlinearLeastSquares');"));
            //        matlab.Cmd(String.Format("opts.Display = 'Off';"));
            //        matlab.Cmd(String.Format("opts.Algorithm = 'Levenberg-Marquardt';"));
            //        matlab.Cmd(String.Format("opts.StartPoint = [0.12553623135482 0.82239400675902 0.0251505014285022 0.414428880924032 0.731407467972937 0.781374000275963 0.367285915131369 0.744867856241673 0.892267188231107 0.24260338627967 0.12959697583754];"));
            //        matlab.Cmd(String.Format("[fitresult, gof] = fit(xData, yData, ft, opts);"));
            //        matlab.Cmd(String.Format("result=[fitresult.p1 0;fitresult.p2 fitresult.q1;fitresult.p3 fitresult.q2;fitresult.p4 fitresult.q3; fitresult.p5 fitresult.q4;fitresult.p6 fitresult.q5];"));
            //        matlab.GetMatrix(coefficients, "result");
            //        matlab.Execute();
            //    }
            //    Phi = X => Math.Max((coefficients[0, 0] * X.Pow(5) + coefficients[1, 0] * X.Pow(4) + coefficients[2, 0] * X.Pow(3) + coefficients[3, 0] * X.Pow(2) + coefficients[4, 0] * X + coefficients[5, 0]) /
            //                (X.Pow(5) + coefficients[1, 1] * X.Pow(4) + coefficients[2, 1] * X.Pow(3) + coefficients[3, 1] * X.Pow(2) + coefficients[4, 1] * X + coefficients[5, 1]), 0.0);

            //    PhiR = X => X >= 0 ? Phi(X) : Phi(0) - X / Math.Tan(angle);

            //    return PhiR;
            //}
            //var phi = MeniskusInitRat55(a, h, ctrl.PhysicalParameters.theta_e);
            //ctrl.InitialValues_Evaluators.Add("Phi", X => -y0 -phi(-X[0]) + X[1]);
            ctrl.InitialValues_Evaluators.Add("Phi", X => -y0 + X[1]);


            //ctrl.AddInitialValue("Phi", $"X => -{y0} + X[1] - MathNet.Numerics.RootFinding.Bisection.FindRoot(y => MathNet.Numerics.Trig.Acosh(2 * {a} / y) - MathNet.Numerics.Trig.Acosh(2 * {a} / {h}) + Math.Sqrt(4 - Math.Pow({h} / {a}, 2)) - Math.Sqrt(4 - Math.Pow(y / {a}, 2)) + X[0], 0.0, {h})", false);

            if (threephase) {
                ctrl.UseImmersedBoundary = true;
                ctrl.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = true;
                ctrl.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NavierSlip;
                ctrl.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Everywhere;
                ctrl.AddInitialValue("Phi2", $"(X, t) => X[0]", true);
                ctrl.AddInitialValue("Temperature#C", $"(X, t) => {ctrl.ThermalParameters.T_sat}", true);
                ctrl.AddInitialValue("VelocityX@Phi2", $"(X, t) => 0", true);
                ctrl.AddInitialValue("VelocityY@Phi2", $"(X, t) => 0", true);
                ctrl.HeatSourceIBM = new Dictionary<string, string>();
                ctrl.HeatSourceIBM["AC"] = "(X,t) => 0.2";
                ctrl.AddBoundaryValue(Bndy[1]);

            } else {
                ctrl.AddBoundaryValue(Bndy[1], "HeatFluxX#A", "X => 0.2", false);
            }

            ctrl.AddInitialValue("Temperature#A", $"(X, t) => {ctrl.ThermalParameters.T_sat}", true);
            ctrl.AddInitialValue("Temperature#B", $"(X, t) => {ctrl.ThermalParameters.T_sat}", true);
            ctrl.AddInitialValue("GravityY#A", $"(X, t) => -{g}", true);



            #endregion

            #region Boundary Conditions

            ctrl.AddBoundaryValue(Bndy[3]);
            ctrl.AddBoundaryValue(Bndy[2]);
            ctrl.AddBoundaryValue(Bndy[4], "Pressure#A", $"(X, t) => {y0} * {ctrl.PhysicalParameters.rho_A} * {g}", true);

            #endregion

            #region AMR

            int level = 5;
            ctrl.AdaptiveMeshRefinement = level > 0;
            ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRonNarrowband() { maxRefinementLevel = level });
            if (threephase) {
                ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRatInnerContactLine() { maxRefinementLevel = 2 * level });
            } else {
                ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater.AMRatContactLine() { edgeTags = new byte[] { 1 }, maxRefinementLevel = 2 * level });
            }

            ctrl.AMR_startUpSweeps = 2 * level;

            #endregion

            #region Timestepping

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            ctrl.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1e-8;
            ctrl.NonLinearSolver.MaxSolverIterations = 10;

            ctrl.SkipSolveAndEvaluateResidual = false;

            ctrl.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;
            //ctrl.dtFixed = 0.01;
            //ctrl.Endtime = 15.0;
            //ctrl.NoOfTimesteps = (int)(ctrl.Endtime / ctrl.dtFixed);

            #endregion
            //if (threephase) {
            //    ctrl.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MovingContactLineZwoLsLogging() { LogPeriod = 1 });
            //} else {
            //    ctrl.PostprocessingModules.Add(new BoSSS.Application.XNSFE_Solver.PhysicalBasedTestcases.MassfluxLogging() { LogPeriod = 1 });
            //}

            ctrl.ImmediatePlotPeriod = 1;
            ctrl.SuperSampling = 2;

            return ctrl;
        }
    }
}
