using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Queries;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using static System.Math;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    static public partial class FullNSEControlExamples {





        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSEC_Control ThermodynamicEquilibrium_steadyStateTest2(int p = 2, int kelemR = 17, string _DbPath = null) {

            XNSEC_Control C = new XNSEC_Control();

            //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            bool steady = true;
            bool separated = false;

            //_DbPath = @"D:\local\local_XNSE_StudyDB";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_HeatedWall";
            _DbPath = @"C:\Databases\BoSSS_DB";
            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/ThermodynamEquilib_steady";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";

            C.ContinueOnIoError = false;

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
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            C.FieldOptions.Add("MassFraction0", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
#endregion


            // Physical Parameters
            // ===================
            #region physics
            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;

            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

 

            // Water (A: liquid, B: gaseous)
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            C.solveCoupledHeatEquation = true;
            if(separated)
                C.conductMode = Solution.XheatCommon.ConductivityInSpeciesBulk.ConductivityMode.LDG;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 0.001;
            C.ThermalParameters.k_A = 1.0;
            double kv = 0.1;
            C.ThermalParameters.k_B = kv;

            if(C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = 100.0;
                //C.ThermalParameters.hVap_B = -100.0;
            }

            double Tsat = 100.0;
            C.ThermalParameters.T_sat = Tsat;
            double pSat = 10;
            C.ThermalParameters.p_sat = pSat;


            bool includeConv = true;
            C.PhysicalParameters.IncludeConvection = includeConv;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 0.1;
            bool periodic = true;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: periodic);

                if(!steady) {
                    grd.EdgeTagNames.Add(1, "wall_ConstantHeatFlux_lower");
                    grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ZeroGradient_upper");
                } else {
                    grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature_lower");
                    grd.EdgeTagNames.Add(2, "velocity_inlet_ZeroGradient_upper");
                }
                if(!periodic)
                    grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;
                    if(!periodic) {
                        if(Math.Abs(X[0]) <= 1.0e-8 || Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 3;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double zi0 = 0.01;

            Func<double[], double> PhiFunc = (X => zi0 - X[1]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double qv = 10.0;
            C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);
            C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat + (qv / kv) * (zi0 - X[1]));

            C.prescribedMassflux_Evaluator = (X, t) => -0.1;

            if(!steady) {
                C.InitialValues_Evaluators.Add("Pressure#A", X => pSat);
                C.InitialValues_Evaluators.Add("Pressure#B", X => pSat - (0.01) * (10 - 1));

                C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.9);
            } else {
                C.InitialValues_Evaluators.Add("VelocityY#A", X => -0.1);
                C.InitialValues_Evaluators.Add("VelocityY#B", X => -1.0);
            }

            //if (separated) {
            //    C.InitialValues_Evaluators.Add("HeatFluxY#B", X => qv);
            //}

            //C.RestartInfo = new Tuple<Guid, TimestepNumber>(restartID, null);


            #endregion


            // boundary conditions
            // ===================
            #region BC

            if(!steady) {
                //C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFlux#A", (X, t) => HeatFlux);
                if(separated) {
                    C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#B", (X, t) => qv);
                } else {
                    C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#B", (X, t) => -qv);
                }

                C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper", "Pressure#A", (X, t) => pSat);
            } else {
                C.AddBoundaryValue("pressure_outlet_ConstantTemperature_lower", "Temperature#B", (X, t) => Tsat + (qv / kv) * (zi0 - X[1]));
                C.AddBoundaryValue("velocity_inlet_ZeroGradient_upper", "VelocityY#A", (X, t) => -0.1);
            }

            if(!periodic)
                C.AddBoundaryValue("slipsymmetry_ZeroGradient");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.useSolutionParamUpdate = true;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSEC_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 1;

            C.InitSignedDistance = false;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = steady ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            #endregion


            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = steady ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = steady ? AppControl._TimesteppingMode.Steady : AppControl._TimesteppingMode.Transient;
            C.dtMax = 5e-4;
            C.dtMin = 5e-4;
            C.Endtime = 10000;
            C.NoOfTimesteps = 200;
            C.saveperiod = 2;

            #endregion


            // additional parameters
            double[] param = new double[2];
            param[0] = L;           // domain height
            param[1] = L / 4.2;   // x probe

            C.AdditionalParameters = param;

            //C.LogValues = XNSE_Control.LoggingValues.EvaporationL;
            //C.LogPeriod = 2;
            //C.PostprocessingModules.Add(new EvaporationLogging() { LogPeriod = 2, mode = EvaporationLogging.Mode.LineInterface });

            return C;
        }




        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static XNSEC_Control ThermodynamicEquilibrium_steadyStateTest(int p = 1, int kelemR = 5, string _DbPath = @"C:\Databases\BoSSS_DB") {

            XNSEC_Control C = new XNSEC_Control();

            bool lowerA = true;

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/ThermodynamEquilib_steady";
            //C.ProjectDescription = "Leikonfiguration for SFB 1194";

            C.ContinueOnIoError = false;


            C.ImmediatePlotPeriod = 1;
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
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Curvature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxX", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFluxY", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("MassFraction0", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            #endregion


            // Physical Parameters
            // ===================
            #region physics
            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;

            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            // Water (A: liquid, B: gaseous)
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            C.PhysicalParameters.Sigma = 0.0;

            C.solveCoupledHeatEquation = true;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_B = C.PhysicalParameters.rho_B;
            C.ThermalParameters.c_A = 1.0;
            C.ThermalParameters.c_B = 0.001;
            C.ThermalParameters.k_A = 1.0;
            C.ThermalParameters.k_B = 1.0;
            double kv = lowerA ? C.ThermalParameters.k_A : C.ThermalParameters.k_B;

            if(C.solveCoupledHeatEquation) {
                C.ThermalParameters.hVap = 100.0;
                //C.ThermalParameters.hVap_B = -100.0;
            }

            double Tsat = 100.0;
            C.ThermalParameters.T_sat = Tsat;
            double pSat = 10;
            C.ThermalParameters.p_sat = pSat;


            bool includeConv = true;
            C.PhysicalParameters.IncludeConvection = includeConv;
            C.ThermalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double L = 0.1;
            bool periodic = true;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelemR + 1);
                double[] Ynodes = GenericBlas.Linspace(0, L, kelemR + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: periodic);

                grd.EdgeTagNames.Add(1, "wall_ConstantHeatFlux_lower");
                grd.EdgeTagNames.Add(2, "pressure_Dirichlet_ZeroGradient_upper");

                if(!periodic)
                    grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 2;
                    if(!periodic) {
                        if(Math.Abs(X[0]) <= 1.0e-8 || Math.Abs(X[0] - L) <= 1.0e-8)
                            et = 3;
                    }

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double zi0 = 0.01;

            Func<double[], double> PhiFunc = (X => lowerA ? X[1] - zi0 : zi0 - X[1]);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double qv = 10.0;

            double massflux = qv / C.ThermalParameters.hVap;
            C.prescribedMassflux_Evaluator = (X, t) => lowerA ? -massflux : massflux;

            if(lowerA) {
                C.InitialValues_Evaluators.Add("Pressure#A", X => pSat + (massflux.Pow2()) * (10 - 1));
                C.InitialValues_Evaluators.Add("Pressure#B", X => pSat);

                C.InitialValues_Evaluators.Add("VelocityY#A", X => 0.0);
                C.InitialValues_Evaluators.Add("VelocityY#B", X => massflux / C.PhysicalParameters.rho_A - massflux / C.ThermalParameters.rho_B);

                C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat);
                C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat - (qv / kv) * (zi0 - X[1]));

                C.InitialValues_Evaluators.Add("MassFraction0#A", X => 1.0);
                C.InitialValues_Evaluators.Add("MassFraction0#B", X => 1.0);

            } else {
                C.InitialValues_Evaluators.Add("Pressure#A", X => pSat);
                C.InitialValues_Evaluators.Add("Pressure#B", X => pSat - (massflux.Pow2()) * (10 - 1));

                C.InitialValues_Evaluators.Add("VelocityY#A", X => massflux / C.PhysicalParameters.rho_B - massflux / C.ThermalParameters.rho_A);
                C.InitialValues_Evaluators.Add("VelocityY#B", X => 0.0);

                C.InitialValues_Evaluators.Add("Temperature#A", X => Tsat);
                C.InitialValues_Evaluators.Add("Temperature#B", X => Tsat + (qv / kv) * (zi0 - X[1]));

                C.InitialValues_Evaluators.Add("MassFraction0#A", X => 1.0);
                C.InitialValues_Evaluators.Add("MassFraction0#B", X => 1.0);
            }


            #endregion


            // boundary conditions
            // ===================
            #region BC

            if(lowerA) {
                C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#A", (X, t) => qv);
                C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper", "Pressure#B", (X, t) => pSat);
            } else {
                C.AddBoundaryValue("wall_ConstantHeatFlux_lower", "HeatFluxY#B", (X, t) => -qv);
                C.AddBoundaryValue("pressure_Dirichlet_ZeroGradient_upper", "Pressure#A", (X, t) => pSat);
            }

            if(!periodic)
                C.AddBoundaryValue("slipsymmetry_ZeroGradient");

            #endregion


            // misc. solver options
            // ====================
            #region solver

            C.ComputeEnergyProperties = false;

            C.LSContiProjectionMethod = BoSSS.Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            //C.useSolutionParamUpdate = true;


            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = dt * 0.2;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.AdvancedDiscretizationOptions.CurvatureNeeded = true;


            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSEC_Control.RefinementStrategy.constantInterface;
            C.RefineNavierSlipBoundary = false;
            C.BaseRefinementLevel = 1;

            C.InitSignedDistance = false;

            #endregion


            // level-set
            // =========
            #region levelset

            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            C.FastMarchingPenaltyTerms = BoSSS.Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

            #endregion

            // Timestepping
            // ============
            #region time

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtMax = 1e-4;
            C.dtMin = 1e-4;
            C.Endtime = 10000;
            C.NoOfTimesteps = 200;
            C.saveperiod = 2;

            #endregion


            // additional parameters
            double[] param = new double[2];
            param[0] = L;           // domain height
            param[1] = L / 4.2;   // x probe

            C.AdditionalParameters = param;
            C.SuperSampling = 3;
            //C.LogValues = XNSE_Control.LoggingValues.EvaporationL;
            //C.LogPeriod = 2;
            //C.PostprocessingModules.Add(new EvaporationLogging() { LogPeriod = 2, mode = EvaporationLogging.Mode.LineInterface });
            //C.SkipSolveAndEvaluateResidual = true;
            return C;
        }

    }
}