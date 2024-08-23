using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.ParameterizedLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Application.XNSE_Solver.Logging;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for the testcase
    /// </summary>
    public static class LevelSetInCapillary {

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control ParameterizedLevelSetInCapillary(int p = 2, int nCells = 8, int AMRlvl = 0) {

            XNSE_Control C = new XNSE_Control();

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;
            bool steadyInterface = false;

            string _DbPath = null;

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "CapillaryParameterizedLevelSetTest";

            //C.ContinueOnIoError = false;

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
                Degree = Math.Max(2, p),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = Math.Max(2, p),
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            #endregion


            // Physical Parameters
            // ===================
            #region physics
            C.solveCoupledHeatEquation = true;
            C.PhysicalParameters.Material = false;

            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;
            //C.IncludeRecoilPressure = true;

            C.PhysicalParameters.rho_A = 600;
            C.PhysicalParameters.rho_B = 9;

            C.PhysicalParameters.mu_A = 0.013;
            C.PhysicalParameters.mu_B = 0.0001019;

            C.PhysicalParameters.Sigma = 0.0;

            C.ThermalParameters.rho_A = 600;
            C.ThermalParameters.rho_B = 9;

            C.ThermalParameters.c_A = 100;
            C.ThermalParameters.c_B = 37;

            C.ThermalParameters.k_A = 4.8;
            C.ThermalParameters.k_B = 0.0251;

            C.ThermalParameters.hVap = 10000;
            C.ThermalParameters.T_sat = 0;
            //C.ThermalParameters.T_wall = 163.75;
            double T_wall = 163.75;

            C.PhysicalParameters.betaS_A = 0.0001;
            C.PhysicalParameters.betaS_B = 0.0001;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.sliplength = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 6.0;

            #endregion

            // grid generation
            // ===============
            #region grid

            double L = 0.001;
            double xSize = L;
            double ySize = L * 5.0;


            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, 2*nCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 5*nCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "NavierSlip_Linear_ZeroGradient");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_ZeroGradient");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_ConstantTemperature");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if ((Math.Abs(X[0] - Xnodes.First()) < 1e-8) || (Math.Abs(X[0] - Xnodes.Last()) < 1e-8))
                        return 1; // walls
                    else if ((Math.Abs(X[1] - Ynodes.Last()) < 1e-8))
                        return 2; // upper border
                    else if ((Math.Abs(X[1] - Ynodes.First()) < 1e-8))
                        return 3; // bottom

                    return et;
                });

                return grd;
            };

            #endregion

            // Initial Values
            // ==============
            #region init

            C.ParameterizedLevelSetControl = new ParameterizedLevelSetControlEllipse(2.0, 2.0, 2.0); // xSemiAxis0, ySemiAxis0, yCenter0

            Func<double[], double>  PhiFunc = C.ParameterizedLevelSetControl.PhiFunc;

            Func<double[], double>  VelocityFunc = (X => (0.1 - 0.1 * X[0].Pow2()/L.Pow2()));
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add("Temperature#B", X => C.ThermalParameters.T_sat);
            C.InitialValues_Evaluators.Add("Temperature#A", X => T_wall);

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityY#A", VelocityFunc);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("NavierSlip_Linear_ZeroGradient"); // walls
            C.AddBoundaryValue("Pressure_Outlet_ZeroGradient"); // upper boundary
            C.AddBoundaryValue("Velocity_Inlet_ConstantTemperature", "Temperature#A", "X => " + T_wall, false); //bottom boundary
            C.AddBoundaryValue("Velocity_Inlet_ConstantTemperature", "VelocityX#A", "X => " + 0.0, false); //bottom boundary
            C.AddBoundaryValue("Velocity_Inlet_ConstantTemperature", "VelocityY#A", VelocityFunc); //bottom boundary
            #endregion

            // misc. solver options
            // ====================
            #region solver
            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;


            C.AdaptiveMeshRefinement = false;
            //C.AMR_startUpSweeps = 2;

            //C.LinearSolver = LinearSolverCode.automatic.GetConfig();
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.NonLinearSolver.MaxSolverIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-7;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;

            C.Option_LevelSetEvolution = LevelSetEvolution.ParameterizedLevelSet;
           
            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;


            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;

            //C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            //C.Timestepper_LevelSetHandling = (steadyInterface) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = compMode;

            double dt = 0.1; //0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1; 
            C.saveperiod = 10;

            #endregion


            return C;

        }
    }
}