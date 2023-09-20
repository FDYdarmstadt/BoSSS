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
    }
}
