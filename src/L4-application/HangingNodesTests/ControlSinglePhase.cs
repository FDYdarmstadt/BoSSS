using BoSSS.Application.XNSFE_Solver;
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
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;

namespace HangingNodesTests {
    partial class Control {

        /// <summary>
        /// Test discrete maximum principle for heat equation
        /// true = acute corner 10°
        /// false = obtuse corner 170°
        /// </summary>
        /// <param name="setup"></param>
        /// <returns></returns>
        static public XNSFE_Control TestSkeleton1Phase(double size) {
            var ctrl = new XNSFE_Control();

            ctrl.SetDGdegree(2);

            //ctrl.FieldOptions["Temperature"] = new FieldOpts() {
            //    Degree = 2,
            //    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            //};

            double L = size;

            ctrl.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.0,
                mu_B = 1.0,

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

                hVap = 0.0
            };
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;

            ctrl.UseImmersedBoundary = false;
            ctrl.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = false;

            ctrl.AddInitialValue("Phi", $"(X, t) => -1.0", true);

            ctrl.GridFunc = () => {
                var Nodes = GenericBlas.Linspace(-L, L, 2);
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

                grd.AddPredefinedPartitioning("2Proc", X => X[0] < 0.0 ? 0 : 1);
                grd.AddPredefinedPartitioning("4Proc", X => X[0] < 0.0 ? (X[1] < 0.0 ? 0 : 1) : (X[1] < 0.0 ? 2 : 3));

                return grd;
            };

            //ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", $"(X, t) => -{L}", true);
            //ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#B", $"(X, t) => -{L}", true);
            //ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#A", $"(X, t) => {L}", true);
            //ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#B", $"(X, t) => {L}", true);
            //ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            //ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            //ctrl.AddInitialValue("Temperature#A", $"(X, t) => (X[0])", true);
            //ctrl.AddInitialValue("Temperature#B", $"(X, t) => (X[0])", true);

            //ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", $"(X, t) => -1", true);
            //ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#B", $"(X, t) => -1", true);
            //ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#A", $"(X, t) => 1", true);
            //ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#B", $"(X, t) => 1", true);
            //ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            //ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            //ctrl.AddInitialValue("Temperature#A", $"(X, t) => (X[0])/{L}", true);
            //ctrl.AddInitialValue("Temperature#B", $"(X, t) => (X[0])/{L}", true);

            double T = 1.0;
            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", $"(X, t) => {T}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#B", $"(X, t) => {T}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#A", $"(X, t) => {T}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#B", $"(X, t) => {T}", true);
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.AddInitialValue("Temperature#A", $"(X, t) => {T}", true);
            ctrl.AddInitialValue("Temperature#B", $"(X, t) => {T}", true);

            ctrl.SkipSolveAndEvaluateResidual = true;            

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.RK_ImplicitEuler;
            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            return ctrl;
        }

        /// <summary>
        /// </summary>
        /// <param name="setup"></param>
        /// <returns></returns>
        static public XNSFE_Control TestSkeletonShear1Phase(double size) {
            var ctrl = new XNSFE_Control();

            ctrl.SetDGdegree(2);

            double L = size;

            ctrl.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.0,
                mu_B = 1.0,

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

                hVap = 0.0
            };
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;

            ctrl.UseImmersedBoundary = false;
            ctrl.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = false;

            ctrl.AddInitialValue("Phi", $"(X, t) => -1.0", true);

            ctrl.GridFunc = () => {
                var Nodes = GenericBlas.Linspace(-L, L, 2);
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

                grd.AddPredefinedPartitioning("2Proc", X => X[0] < 0.0 ? 0 : 1);
                grd.AddPredefinedPartitioning("4Proc", X => X[0] < 0.0 ? (X[1] < 0.0 ? 0 : 1) : (X[1] < 0.0 ? 2 : 3));

                return grd;
            };

            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "VelocityY#A", $"(X, t) => -1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "VelocityY#B", $"(X, t) => -1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "VelocityY#A", $"(X, t) => 1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "VelocityY#B", $"(X, t) => 1000*{L}", true);
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.AddInitialValue("VelocityY#A", $"(X, t) => 1000*(X[0])", true);
            ctrl.AddInitialValue("VelocityY#B", $"(X, t) => 1000*(X[0])", true);

            ctrl.SkipSolveAndEvaluateResidual = true;

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.RK_ImplicitEuler;
            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;

            return ctrl;
        }
             
    }
}
