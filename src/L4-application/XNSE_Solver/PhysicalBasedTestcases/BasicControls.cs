using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// Collection of Simple Simulations, intended to verify XNSERefactor 
    /// </summary>
    public static class BasicControls {

        #region Basic XNSE Tests

        /// <summary>
        /// Simple Channelflow, Singlephase
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <returns></returns>
        public static XNSE_Control CFSinglePhase(int p = 2, int kelem = 4) {

            XNSE_Control C = new XNSE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Channel flow for BC testing";
            C.SuperSampling = 3;

            #endregion
            // ============================================

            // ============================================
            #region physics

            #region grid
            double H = 1;
            double L = 5;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 5 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            #region initial condition

            Func<double[], double> PhiFunc = (X => -1);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> UInit = X => -X[1] * (X[1] - H);

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);

            #endregion

            #region boundary condition

            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", UInit);
            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("pressure_outlet_right");

            #endregion

            #region material
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
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

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

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

        /// <summary>
        /// Simple Channelflow, Dualphase with interface along channel direction, same velocity profile as Singlephase
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <returns></returns>
        public static XNSE_Control CFDualPhase(int p = 2, int kelem = 4) {

            XNSE_Control C = new XNSE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Channel flow for BC testing";
            C.SuperSampling = 3;

            #endregion
            // ============================================

            // ============================================
            #region physics

            #region grid
            double H = 1;
            double L = 5;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, 5 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            #region initial condition

            Func<double[], double> PhiFunc = (X => X[1] - 0.5);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> UInit = X => -X[1] * (X[1] - H);

            C.InitialValues_Evaluators.Add("VelocityX#A", UInit);
            C.InitialValues_Evaluators.Add("VelocityX#B", UInit);

            #endregion

            #region boundary condition

            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", UInit);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", UInit);
            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("pressure_outlet_right");

            #endregion

            #region material
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            #endregion

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees

            C.SetDGdegree(p);

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver

            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

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

        #endregion

        #region Basic XHeat Tests

        /// <summary>
        /// Heat conduction in a Box, Singlephase
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <returns></returns>
        public static XNSE_Control BoxSinglePhase(int p = 2, int kelem = 4) {

            XNSE_Control C = new XNSE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Channel flow for BC testing";
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
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "ZeroGradient_lower");
                grd.EdgeTagNames.Add(2, "ZeroGradient_upper");

                grd.EdgeTagNames.Add(3, "ConstantTemperature_left");
                grd.EdgeTagNames.Add(4, "ConstantTemperature_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion

            #region initial condition

            Func<double[], double> PhiFunc = (X => -1);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> UInit = X => 1.0;

            C.InitialValues_Evaluators.Add("Temperature#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#A", X => Math.Sin(2 * Math.PI * X[1]));
            C.InitialValues_Evaluators.Add("VelocityY#A", X => Math.Cos(2 * Math.PI * X[0]));

            #endregion

            #region boundary condition

            C.AddBoundaryValue("ConstantTemperature_left", "Temperature#A", X => 0.0);
            C.AddBoundaryValue("ZeroGradient_lower");
            C.AddBoundaryValue("ZeroGradient_upper");
            C.AddBoundaryValue("ConstantTemperature_right", "Temperature#A", UInit);

            #endregion

            #region material
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 1.0;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 1.0;
            C.ThermalParameters.rho_A = 1.0;
            C.ThermalParameters.rho_A = 1.0;

            C.ThermalParameters.IncludeConvection = true;
            C.ThermalParameters.hVap = 0.0;
            #endregion

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees

            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            C.solveCoupledHeatEquation = true;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.LDG;
            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

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

        /// <summary>
        /// Heat conduction in a box, Dualphase, kink in Temperatureprofile at Interface
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <returns></returns>
        public static XNSE_Control BoxDualPhase(int p = 2, int kelem = 4) {

            XNSE_Control C = new XNSE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/elementalTest";
            C.ProjectDescription = "Channel flow for BC testing";
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
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "ZeroGradient_lower");
                grd.EdgeTagNames.Add(2, "ZeroGradient_upper");

                grd.EdgeTagNames.Add(3, "ConstantTemperature_left");
                grd.EdgeTagNames.Add(4, "ConstantTemperature_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                AffineTrafo ROT = AffineTrafo.Some2DRotation(Math.PI/6.0);

                var grdT = grd.Transform(ROT);
                return grdT;
            };

            #endregion

            #region initial condition

            Func<double[], double> PhiFunc = (X => X[0] - 0.6);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> UInit = X => 1.0;

            C.InitialValues_Evaluators.Add("Temperature#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("VelocityX#A", X => Math.Sin(2 * Math.PI * X[1]));
            //C.InitialValues_Evaluators.Add("VelocityY#A", X => Math.Cos(2 * Math.PI * X[0]));
            C.InitialValues_Evaluators.Add("VelocityX#A", X => -0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => -0.0);

            #endregion

            #region boundary condition

            C.AddBoundaryValue("ConstantTemperature_left", "Temperature#A", X => 0.0);
            C.AddBoundaryValue("ZeroGradient_lower");
            C.AddBoundaryValue("ZeroGradient_upper");
            C.AddBoundaryValue("ConstantTemperature_right", "Temperature#A", UInit);
            C.AddBoundaryValue("ConstantTemperature_right", "Temperature#B", UInit);

            #endregion

            #region material
            C.ThermalParameters.c_A = 10.0;
            C.ThermalParameters.c_B = 1.0;
            C.ThermalParameters.k_A = 10.0;
            C.ThermalParameters.k_B = 1.0;
            C.ThermalParameters.rho_A = 1.0;
            C.ThermalParameters.rho_A = 1.0;

            C.ThermalParameters.IncludeConvection = true;
            C.ThermalParameters.hVap = 0.0;
            C.ThermalParameters.T_sat = Double.MaxValue;
            #endregion

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees

            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            C.solveCoupledHeatEquation = true;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;
            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
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
        #endregion

        #region Basic XNSFE Tests

        /// <summary>
        /// XNSFE with evaporation, Dualphase
        /// </summary>
        /// <param name="p"></param>
        /// <param name="kelem"></param>
        /// <returns></returns>
        public static XNSE_Control BoxEvapDualPhase(int p = 3, int kelem = 1) {

            XNSE_Control C = new XNSE_Control();

            // ============================================
            #region IO

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSFE/elementalTest";
            C.ProjectDescription = "1D Evaporation";
            C.SuperSampling = 3;

            #endregion
            // ============================================

            // ============================================
            #region physics

            #region material
            C.ThermalParameters.c_A = 0.0;
            C.ThermalParameters.c_B = 0.0;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 1.0;
            C.ThermalParameters.rho_A = 1.0;
            C.ThermalParameters.rho_B = 0.1;

            C.ThermalParameters.IncludeConvection = true;
            C.ThermalParameters.hVap = 100.0;
            C.ThermalParameters.T_sat = 0.0;

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            C.PhysicalParameters.Sigma = 0.0;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;
            #endregion

            #region grid
            double angle = Math.PI * 0.0 / 180.0;
            AffineTrafo ROT = AffineTrafo.Some2DRotation(angle);
            double H = 4;
            double L = 1;
            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 4 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "freeslip_ZeroGradient_left");
                grd.EdgeTagNames.Add(2, "freeslip_ZeroGradient_right");

                grd.EdgeTagNames.Add(3, "wall_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(4, "pressure_outlet_ConstantHeatFlux_upper");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 4;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                

                var grdT = grd.Transform(ROT);
                return grdT;
            };

            #endregion           

            #region boundary condition
            double q = 25.0;
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#A", X => 0.0);
            C.AddBoundaryValue("freeslip_ZeroGradient_left");
            C.AddBoundaryValue("freeslip_ZeroGradient_right");
            C.AddBoundaryValue("pressure_outlet_ConstantHeatFlux_upper", "HeatFluxY#B", X => Math.Cos(angle) * q);
            C.AddBoundaryValue("pressure_outlet_ConstantHeatFlux_upper", "HeatFluxX#B", X => -Math.Sin(angle) * q);

            #endregion

            #region initial condition            
            double x0 = 2.5;          

            
            C.Phi0Initial = X => 1.0 / Math.Cos(angle) * (x0 + (0.1 + Math.Sin(angle)) * X);
            Func<double[], double> PhiFunc = (X => Math.Cos(angle) * (X[1] -  C.Phi0Initial(X[0])));
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("Temperature#B", X => (Math.Cos(angle) * X[1] - Math.Sin(angle) * X[0] - x0) * q/C.ThermalParameters.k_B);
            //C.InitialValues_Evaluators.Add("VelocityX#A", X => Math.Sin(2 * Math.PI * X[1]));
            //C.InitialValues_Evaluators.Add("VelocityY#A", X => Math.Cos(2 * Math.PI * X[0]));
            C.InitialValues_Evaluators.Add("VelocityX#A", X => -0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => -0.0);

            #endregion            

            #endregion
            // ============================================

            // ============================================
            #region numerics

            #region dg degrees
            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {
                Degree = p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion

            #region level set

            C.Option_LevelSetEvolution = LevelSetEvolution.SplineLS;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            #endregion

            #region solver

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;

            C.solveCoupledHeatEquation = true;
            C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.SIP;
            #endregion

            #region timestepping

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

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

        #endregion
    }
}
