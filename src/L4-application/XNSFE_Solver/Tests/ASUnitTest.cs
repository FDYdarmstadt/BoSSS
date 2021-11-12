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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using ilPSP;
using BoSSS.Solution.AdvancedSolvers.Testing;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSFE_Solver.Tests {

    /// <summary>
    /// A collection of all-up NUnit tests for the new XNSE solver (<see cref="XNSE"/>).
    /// </summary>
    [TestFixture]
    static public partial class ASUnitTest {

        /// <summary>
        /// Simple pure Heatconductivity Test with Temperature kink at the Interface
        /// </summary>
        [Test]
        public static void HeatConductivityTest(
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm) {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new HeatConductivityTest();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 3);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            XHeatSolverTest(Tst, C);
        }

        /// <summary>
        /// Simple Heatconductivity Test to ensure Energy conservation of Heat equation
        /// </summary>
        [Test]
        public static void HeatDecayTest(
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.8598)] double r,
            [Values(-130, -50, 0.0, 10, 167)] double q,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm) {
            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new HeatDecayTest(r, q);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 3);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            XHeatSolverTest(Tst, C);
        }

        /// <summary>
        /// Simple Test for Evaporation of a straight interface
        /// </summary>
        [Test]
        public static void SteadyStateEvaporationTest(
            [Values(0.0, 15.0, 45.0, 73.1264, 90.0)] double rawangle,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver) // evaporation currently only implemented with use of newton solver
            {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new SteadyStateEvaporationTest(rawangle * Math.PI / 180.0);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 2, nonlinsolver: nonlinsolver);
            XNSFESolverTest(Tst, C);
        }

        private static void XHeatSolverTest(IXHeatTest Tst, XNSFE_Control C) {
            using (var solver = new XHeat()) {

                solver.Init(C);
                solver.RunSolverMode();
                solver.OperatorAnalysis();

                //-------------------Evaluate Temperature Error ---------------------------------------- 
                var evaluator = new XHeatErrorEvaluator<XNSFE_Control>(solver);
                if (Tst.CheckT) {                    
                    double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                    double[] ErrThresh = Tst.AcceptableL2Error;
                    if (LastErrors.Length != ErrThresh.Length)
                        throw new ApplicationException();
                    for (int i = 0; i < ErrThresh.Length; i++) {
                        Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                    }
                    for (int i = 0; i < ErrThresh.Length; i++)
                        Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);
                }

                //-------------------Evaluate Energy Error ---------------------------------------- 
                if (Tst.CheckE) {
                    double LastErrors = evaluator.ComputeEnergyError(Tst.GetE(), Tst.dt);

                    double ErrThresh = Tst.AcceptableL2Error[0];   
                    
                    Console.WriteLine("relative error, '{0}': \t{1}", "total thermal Energy", LastErrors);

                    // less than one promille
                    Assert.LessOrEqual(LastErrors, 1e-3);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);
                }

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);
                
            }
        }
        private static void XNSFESolverTest(IXNSFETest Tst, XNSFE_Control C) {

            using (var solver = new XNSFE<XNSFE_Control>()) {

                solver.Init(C);
                solver.RunSolverMode();
                solver.OperatorAnalysis();

                //-------------------Evaluate Flow Error ---------------------------------------- 
                var flowevaluator = new XNSEErrorEvaluator<XNSFE_Control>(solver);
                double[] LastErrors = flowevaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length - 1) // Last Error thershold is for temperature
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length - 1; i++) {
                    Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                }
                for (int i = 0; i < LastErrors.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                //-------------------Evaluate Temperature Error ---------------------------------------- 
                var heatevaluator = new XHeatErrorEvaluator<XNSFE_Control>(solver);
                if (Tst.CheckT) {
                    LastErrors = LastErrors.Cat(heatevaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C) );

                    // Last Error belongs to temperature
                    for (int i = ErrThresh.Length - 1; i < ErrThresh.Length; i++) {
                        Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                    }
                    for (int i = ErrThresh.Length - 1; i < ErrThresh.Length; i++)
                        Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);
                }

                //-------------------Evaluate Energy Error ---------------------------------------- 
                if (Tst.CheckE) {
                    double LastError = heatevaluator.ComputeEnergyError(Tst.GetE(), Tst.dt);

                    Console.WriteLine("relative error, '{0}': \t{1}", "total thermal Energy", LastErrors);

                    // less than one promille
                    Assert.LessOrEqual(LastError, 1e-3);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);

            }
        }

        class AS_XHeat_Control : XNSFE_Control {
            public override Type GetSolverType() {
                return typeof(XHeat);
            }
        }

        static AS_XHeat_Control TstObj2CtrlObj(IXHeatTest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1) {
            AS_XHeat_Control C = new AS_XHeat_Control();
            int D = tst.SpatialDimension;


            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XHEAT/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add(VariableNames.Curvature, new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            C.GridFunc = () => tst.CreateGrid(GridResolution);

            // boundary conditions
            // ===================

            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // Physical parameters
            // ====================

            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;
            C.ThermalParameters.c_A = tst.c_A;
            C.ThermalParameters.c_B = tst.c_B;
            C.ThermalParameters.k_A = tst.k_A;
            C.ThermalParameters.k_B = tst.k_B;
            C.ThermalParameters.hVap = tst.h_vap;
            C.ThermalParameters.T_sat = tst.T_sat;
            C.ThermalParameters.IncludeConvection = tst.IncludeConvection;

            // initial values and exact solution
            // =================================

            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionTemperature.Add(spc, tst.GetT(spc));

                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, tst.GetT(spc).Convert_Xt2X(0.0));
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AgglomerationThreshold = AgglomerationTreshold;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            C.CutCellQuadratureType = CutCellQuadratureType;

            // timestepping and solver
            // =======================


            if (tst.steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            // return
            // ======

            return C;
        }

        static AS_XHeat_Control TstObj2CtrlObj(IXNSFETest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Newton) {
            AS_XHeat_Control C = new AS_XHeat_Control();
            int D = tst.SpatialDimension;


            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSFE/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========
            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = FlowSolverDegree - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add(VariableNames.Curvature, new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree ,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

            // grid
            // ====

            C.GridFunc = () => tst.CreateGrid(GridResolution);

            // boundary conditions
            // ===================

            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // Physical parameters
            // ====================

            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;
            C.ThermalParameters.c_A = tst.c_A;
            C.ThermalParameters.c_B = tst.c_B;
            C.ThermalParameters.k_A = tst.k_A;
            C.ThermalParameters.k_B = tst.k_B;
            C.ThermalParameters.hVap = tst.h_vap;
            C.ThermalParameters.T_sat = tst.T_sat;
            C.ThermalParameters.IncludeConvection = tst.IncludeConvection;


            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;
            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;
            // initial values and exact solution
            // =================================

            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
                        
            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionTemperature.Add(spc, tst.GetT(spc));
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, tst.GetT(spc).Convert_Xt2X(0.0));

                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));

                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                }

                Func<double[], double, double>[] Gravity = new Func<double[], double, double>[D];
                for (int d = 0; d < D; d++) {
                    var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AgglomerationThreshold = AgglomerationTreshold;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            C.CutCellQuadratureType = CutCellQuadratureType;

            // timestepping and solver
            // =======================


            if (tst.steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.NonLinearSolver.SolverCode = nonlinsolver;

            // return
            // ======

            return C;
        }
    }
}
