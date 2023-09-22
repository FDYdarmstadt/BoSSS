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
    public partial class Control {

        /// <summary>
        /// Test discrete maximum principle for heat equation
        /// true = acute corner 10°
        /// false = obtuse corner 170°
        /// </summary>
        /// <param name="setup"></param>
        /// <returns></returns>
        static public XNSFE_Control TestSkeleton(double size) {
            var ctrl = new XNSFE_Control();

            //ctrl.SuperSampling = 0;
            //ctrl.ImmediatePlotPeriod = 1;

            ctrl.SetDGdegree(2);           

            double L = size;

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

                hVap = 0.0
            };
            ctrl.PhysicalParameters.IncludeConvection = false;
            ctrl.ThermalParameters.IncludeConvection = false;

            ctrl.UseImmersedBoundary = true;
            ctrl.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = true;

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

                grd.AddPredefinedPartitioning("2Proc", X => X[0] < 0.0 ? 0 : 1);
                grd.AddPredefinedPartitioning("2ProcTranspose", X => X[1] < 0.0 ? 0 : 1);
                grd.AddPredefinedPartitioning("4Proc", X => X[0] < 0.0 ? (X[1] < 0.0 ? 0 : 1) : (X[1] < 0.0 ? 2 : 3));

                return grd;
            };

            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", $"(X, t) => -1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#B", $"(X, t) => -1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#C", $"(X, t) => -1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#A", $"(X, t) => 1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#B", $"(X, t) => 1000*{L}", true);
            ctrl.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#C", $"(X, t) => 1000*{L}", true);
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_upper");
            ctrl.AddBoundaryValue("pressure_outlet_ZeroGradient_lower");

            ctrl.AddInitialValue("Temperature#A", $"(X, t) => 1000*(X[0])", true);
            ctrl.AddInitialValue("Temperature#B", $"(X, t) => 1000*(X[0])", true);
            ctrl.AddInitialValue("Temperature#C", $"(X, t) => 1000*(X[0])", true);

            ctrl.SkipSolveAndEvaluateResidual = true;            

            ctrl.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            ctrl.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
            ctrl.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.None;

            ctrl.TimesteppingMode = BoSSS.Solution.Control.AppControl._TimesteppingMode.Steady;


            ctrl.CutCellQuadratureType = BoSSS.Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;

            return ctrl;
        }

        public static void SetAMR(XNSFE_Control ctrl, double size, byte setup) {

            if(setup != 0 && setup < 16)
                ctrl.AdaptiveMeshRefinement = true;
            
            double L = size;
            if((setup & 1) != 0) {
                ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.AMRLevelIndicatorLibrary.AMRInBoundingBox(new BoundingBox(new double[,] { { -L, -L }, { 0, 0 } })) { maxRefinementLevel = 1 });
            }
            if ((setup & 2) != 0) {
                ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.AMRLevelIndicatorLibrary.AMRInBoundingBox(new BoundingBox(new double[,] { { -L, 0 }, { 0, L } })) { maxRefinementLevel = 1 });
            }
            if ((setup & 4) != 0) {
                ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.AMRLevelIndicatorLibrary.AMRInBoundingBox(new BoundingBox(new double[,] { { 0, -L }, { L, 0 } })) { maxRefinementLevel = 1 });
            }
            if ((setup & 8) != 0) {
                ctrl.activeAMRlevelIndicators.Add(new BoSSS.Solution.AMRLevelIndicatorLibrary.AMRInBoundingBox(new BoundingBox(new double[,] { { 0, 0 }, { L, L } })) { maxRefinementLevel = 1 });
            }
            
            ctrl.AMR_startUpSweeps = (new int[] { 1 }).Max();
        }

        public static void SetLevelSet(XNSFE_Control ctrl, double size, int phases) {

            double L = size;
            double offset = 5.0 * L / 16.0;
            double ContactAngle = ctrl.PhysicalParameters.theta_e;
            switch (phases) {
                case 1:
                    ctrl.AddInitialValue("Phi", $"(X, t) => -1.0", true);
                    ctrl.AddInitialValue("Phi2", $"(X, t) => -1.0", true);
                    break;
                case 2:
                    ctrl.AddInitialValue("Phi", $"(X, t) => -1.0", true);
                    ctrl.AddInitialValue("Phi2", $"(X, t) => 1/{L} * (X[0])", true);
                    break;
                case 3:
                    ctrl.AddInitialValue("Phi", $"(X, t) => X[0] < 0.0 ? X[1] - {offset} - (X[0])/Math.Tan({ContactAngle}) : X[1] - {offset} ", true);
                    ctrl.AddInitialValue("Phi2", $"(X, t) => 1/{L} * (X[0])", true);
                    break;
                default:
                    throw new NotSupportedException();
            }
        }

        public static void SetParallel(XNSFE_Control ctrl, int procs) {
            if(procs == 1) {
                return;
            }

            ctrl.GridPartType = GridPartType.Predefined;
            switch (procs) {
                case 2:
                    ctrl.GridPartOptions = "2Proc";
                    break;
                case -2:
                    ctrl.GridPartOptions = "2ProcTranspose";
                    break;
                case 4:
                    ctrl.GridPartOptions = "4Proc";
                    break;
                default:
                    throw new NotSupportedException();
            }           
        }
    }
}
