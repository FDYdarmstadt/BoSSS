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
using BoSSS.Solution.NSECommon;
using ilPSP;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Utils;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.XDG;

namespace BoSSS.Application.XNSFE_Solver {

    /// <summary>
    /// A few example configurations.
    /// Test out interfacial slip.
    /// </summary>
    public static class HardcodedControl {

        /// <summary>
        /// Maintainer: rieckmann
        /// Simple shear flow with dividing interface, to test interfacial slip
        /// </summary>
        public static XNSFE_Control ShearFlowInterfaceSlip(int p = 2, int kelem = 5) {

            XNSFE_Control C = new XNSFE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Shear flow with slip interface";
            C.SuperSampling = 3;
            
            #endregion
            // ============================================

            // ============================================
            #region physics

            #region grid
            double H = 1;
            double L = 5;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-H, H, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-H, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + H) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 1;
                    return et;
                });

                return grd;
            };

            #endregion

            #region initial condition

            Func<double[], double> PhiFunc = (X => X[1]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            #endregion

            #region boundary condition

            C.AddBoundaryValue("wall_ConstantTemperature", "VelocityX#A", "X => X[1]", false);
            C.AddBoundaryValue("wall_ConstantTemperature", "VelocityX#B", "X => X[1]", false);

            #endregion

            #region material

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 0.001;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 0.1;
            C.PhysicalParameters.Sigma = 1.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.slipI = 1 * 2 * H; // 10% wider, i.e. at the interface should be +/- 0.1

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = false;

            C.ThermalParameters.IncludeConvection = false;
            C.ThermalParameters.hVap = 0.0;
            C.ThermalParameters.T_sat = 0.0;


            #endregion

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees

            C.SetDGdegree(p);

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.SkipSolveAndEvaluateResidual = false;
            
            // double dt = 1e-1;
            // C.dtMax = dt;
            // C.dtMin = dt;
            // C.Endtime = 1000;
            // C.NoOfTimesteps = 100;
            // C.saveperiod = 1;

            #endregion

            #endregion
            // ============================================

            return C;
        }

        /// <summary>
        /// Maintainer: rieckmann
        /// Artificial wedge, kept with steady-interface;
        /// Change the boundary conditions to investigate convergence
        /// For some B.C. there will be singularities, i.e. No-Slip
        /// in turn this should destroy optimal convergence 
        /// <param name="wall_slip">
        /// - 0: no-slip
        /// - 1: some-slip
        /// - 2: free-slip
        /// </param>
        /// <param name="interface_slip">
        /// - 0: no-slip
        /// - 1: some-slip
        /// - 2: free-slip
        /// </param>
        /// <param name="evap">
        /// - false: no-evaporation
        /// - true: evaporation
        /// </param>
        /// </summary>
        public static XNSFE_Control WedgeConvergence(int p = 2, int res = 0, byte wall_slip = 0, byte interface_slip = 0, bool evap = false) {

            XNSFE_Control C = new XNSFE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Wedge convergence test";
            C.SuperSampling = 3;
            
            #endregion
            // ============================================

            // ============================================
            #region physics

            #region grid
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 1, (int)Math.Pow(2, res) + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 1, (int)Math.Pow(2, res) + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");
                switch(wall_slip){
                    case 0:
                    default:{
                        grd.EdgeTagNames.Add(2, "wall_ConstantHeatflux");
                        break;
                    }
                    case 1:{
                        grd.EdgeTagNames.Add(2, "NavierSlip_linear_ConstantHeatflux");
                        break;
                    }
                    case 2:{
                        grd.EdgeTagNames.Add(2, "freeslip_ConstantHeatflux");
                        break;
                    }
                }
                grd.EdgeTagNames.Add(3, "pressure_outlet_ZeroGradient_bottom");
                grd.EdgeTagNames.Add(4, "wall_ZeroGradient_top");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0] - Xnodes.First()) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - Xnodes.Last()) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1] - Ynodes.First()) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - Ynodes.Last()) <= 1.0e-8)
                        et = 4;
                    return et;
                });

                return grd;
            };

            #endregion

            #region initial condition

            Func<double[], double> PhiFunc = (X => X[1]);
            double theta = 30.0 / 90.0 * Math.PI / 2.0;
            double dtheta = 5.0 / 90.0 * Math.PI / 2.0;
            double y0 = 0.6;
            C.AddInitialValue("Phi", $"X => (X[1]-{y0}) - (X[0] - 1.0)/Math.Tan({theta + dtheta})", false);

            #endregion

            #region boundary condition

            C.AddBoundaryValue("pressure_outlet_ConstantTemperature");
            switch(wall_slip){
                case 0:
                default:{
                    C.AddBoundaryValue("wall_ConstantHeatflux", "HeatFluxX#A", "X => 1.0", false);
                    break;
                }
                case 1:{
                    C.AddBoundaryValue("NavierSlip_linear_ConstantHeatflux", "HeatFluxX#A", "X => 1.0", false);
                    break;
                }
                case 2:{
                    C.AddBoundaryValue("freeslip_ConstantHeatflux", "HeatFluxX#A", "X => 1.0", false);
                    break;
                }
            }
            C.AddBoundaryValue("wall_ZeroGradient_top");
            C.AddBoundaryValue("pressure_outlet_ZeroGradient_bottom");
            #endregion

            #region material

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 0.001;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 0.01;
            C.PhysicalParameters.Sigma = 1.0;
            C.PhysicalParameters.theta_e = theta;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.sliplength = 0.1;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_SlipLength;
            switch(interface_slip){
                    case 0:
                    default:{
                        C.PhysicalParameters.slipI = 0.0;
                        break;
                    }
                    case 1:{
                        C.PhysicalParameters.slipI = 0.1;
                        break;
                    }
                    case 2:{
                        C.PhysicalParameters.slipI = double.PositiveInfinity;
                        break;
                    }
                }

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = !evap;

            C.ThermalParameters.IncludeConvection = false;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 0.1;
            C.ThermalParameters.hVap = evap ? 1.0 : 0.0;
            C.ThermalParameters.T_sat = 0.0;


            #endregion

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees

            C.SetDGdegree(p);

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.SkipSolveAndEvaluateResidual = false;
            
            // double dt = 1e-1;
            // C.dtMax = dt;
            // C.dtMin = dt;
            // C.Endtime = 1000;
            // C.NoOfTimesteps = 100;
            // C.saveperiod = 1;

            #endregion

            #endregion
            // ============================================

            return C;
        }

        /// <summary>
        /// Maintainer: rieckmann
        /// Static circular meniscus in a capillary,
        /// Investigation of the influence of slip, when evaporation is turned on.
        /// </summary>
        public static XNSFE_Control CapillaryEvaporationInterfaceSlip(int p = 2, int kelem = 5) {

            XNSFE_Control C = new XNSFE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Shear flow with slip interface";
            C.SuperSampling = 3;
            
            #endregion
            // ============================================

            // ============================================
            #region physics

            #region grid
            double H = 1;
            double L = 1;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 5*H, 5 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "freeslip_ZeroGradient");
                grd.EdgeTagNames.Add(2, "pressure_outlet_ConstantTemperature");
                grd.EdgeTagNames.Add(3, "wall_ConstantTemperature");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] - Ynodes.First()) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - Ynodes.Last()) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] - Xnodes.First()) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - Xnodes.Last()) <= 1.0e-8)
                        et = 1;
                    return et;
                });

                return grd;
            };

            #endregion

            #region initial condition
            double h0 = 2 * H;
            double alpha = Math.PI * 30 / 180;
            double R = L / Math.Sin(Math.PI/2.0 - alpha);
            Func<double[], double> PhiFunc = (X => X[1] - (h0 + R - Math.Sqrt(R*R - X[0]*X[0])));

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            #endregion

            #region boundary condition

            C.AddBoundaryValue("wall_ConstantTemperature", "Temperature#A", "X => 1.0", false);

            #endregion

            #region material

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 1.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.slipI = 0.0; // 0.1

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            C.ThermalParameters.IncludeConvection = false;
            C.ThermalParameters.hVap = 10.0;
            C.ThermalParameters.T_sat = 0.0;

            #endregion

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees

            C.SetDGdegree(p);

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 100;
            C.saveperiod = 1;

            #endregion

            #endregion
            // ============================================

            return C;
        }
    }
}