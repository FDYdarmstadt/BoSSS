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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.Tests.IsentropicVortex;
using ilPSP.Utils;
using NUnit.Framework;
using System;

namespace CNS.Tests.IBMTests {

    /// <summary>
    /// Tests based on the well-known isentropic vortex solution conceived by
    /// Lamb.
    /// </summary>
    [TestFixture]
    public class IBMIsentropicVortexTest : TestProgram<IBMControl> {


        ///// <summary>
        ///// Alternative entry point of this assembly that allows to perform
        ///// isentropic vortex tests conveniently.
        ///// </summary>
        ///// <param name="args">
        ///// Command line arguments
        ///// </param>
        //public static void Main(string[] args) {
        //    SetUp();
        //    //IBMVortexOneStepGaussAndStokesNoAgglomerationTest();
        //    IBMVortexClassicAgglomerationTest();
        //    //IBMVortexCutNextToCutAgglomerationTest();
        //    //IBMVortexClassicOptimizedHLLCAgglomerationTest();
        //    //IBMVortexLocalTimeSteppingTest();
        //}

        /// <summary>
        /// Tests the IBM CNS solver with the Rusanov flux using the
        /// example moving isentropic vortex in an ideal gas. Uses Florian's
        /// extended HMF quadrature <b>without</b> agglomeration.
        /// </summary>
        [Test]
        public static void IBMVortexOneStepGaussAndStokesNoAgglomerationTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlNoAgglomeration());
            p.RunSolverMode();

            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 4.5e-3),
                Tuple.Create("L2ErrorPressure", 5.5e-3),
                Tuple.Create("L2ErrorEntropy", 6.0e-3));
        }

        /// <summary>
        /// Control file for <see cref="IBMVortexOneStepGaussAndStokesNoAgglomerationTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlNoAgglomeration() {
            IBMControl c = ControlTemplate(dgDegree: 2, divisions: 1, levelSetPosition: -0.25);

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 10;
            c.AgglomerationThreshold = 0.0;

            return c;
        }

        /// <summary>
        /// Tests the IBM CNS solver with the Rusanov flux using the
        /// example moving isentropic vortex in an ideal gas. Uses the original
        /// HMF quadrature <b>with</b> agglomeration.
        /// </summary>
        [Test]
        public static void IBMVortexClassicAgglomerationTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlRusanovAgglomeration());
            p.RunSolverMode();

            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 4.0e-3),
                Tuple.Create("L2ErrorPressure", 5.0e-3),
                Tuple.Create("L2ErrorEntropy", 5.5e-3));
        }

        /// <summary>
        /// Control file for <see cref="IBMVortexClassicAgglomerationTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlRusanovAgglomeration() {
            IBMControl c = ControlTemplate(dgDegree: 2, divisions: 1, levelSetPosition: -0.05);

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = 5;
            c.AgglomerationThreshold = 0.1;

            return c;
        }

        /// <summary>
        /// Tests the IBM CNS solver with the HLLC flux using the
        /// example moving isentropic vortex in an ideal gas. Uses the original
        /// HMF quadrature <b>with</b> agglomeration.
        /// </summary>
        [Test]
        public static void IBMVortexClassicHLLCAgglomerationTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlHLLCAgglomeration());
            p.RunSolverMode();

            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 3.5e-3),
                Tuple.Create("L2ErrorPressure", 4.5e-3),
                Tuple.Create("L2ErrorEntropy", 5.0e-3));
        }

        /// <summary>
        /// Control file for <see cref="IBMVortexClassicHLLCAgglomerationTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlHLLCAgglomeration() {
            IBMControl c = ControlTemplate(dgDegree: 2, divisions: 1, levelSetPosition: -0.05);

            c.ConvectiveFluxType = ConvectiveFluxTypes.HLLC;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = 5;
            c.AgglomerationThreshold = 0.1;

            return c;
        }

        /// <summary>
        /// Tests the IBM CNS solver with the optimized HLLC flux using the
        /// example moving isentropic vortex in an ideal gas. Uses the original
        /// HMF quadrature <b>with</b> agglomeration.
        /// </summary>
        [Test]
        public static void IBMVortexClassicOptimizedHLLCAgglomerationTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlOptimizedHLLCAgglomeration());
            p.RunSolverMode();

            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 3.5e-3),
                Tuple.Create("L2ErrorPressure", 4.5e-3),
                Tuple.Create("L2ErrorEntropy", 5.2e-3));
        }

        /// <summary>
        /// Control file for <see cref="IBMVortexClassicOptimizedHLLCAgglomerationTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlOptimizedHLLCAgglomeration() {
            IBMControl c = ControlTemplate(dgDegree: 2, divisions: 1, levelSetPosition: -0.05);

            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = 5;
            c.AgglomerationThreshold = 0.1;

            return c;
        }

        /// <summary>
        /// Tests the IBM CNS solver with the optimized HLLC flux using the
        /// example moving isentropic vortex in an ideal gas. Uses the original
        /// HMF quadrature and local time stepping (LTS) <see cref="IBMAdamsBashforthLTS"/>
        /// </summary>
        [Test]
        public static void IBMVortexLocalTimeSteppingTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlLocalTimeStepping());
            p.RunSolverMode();

            // due to changes in HMF, the error thresholds went up in comparison to the original reference-runs
            // once the 'IBMCylinderTests.zip' database is replaced, these tolerances may finally be  restored to 
            // original, sharper thresholds 
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2*0.00297096129679995 + 1e-14),
                Tuple.Create("L2ErrorPressure", 2*0.00362102145336951 + 1e-14),
                Tuple.Create("L2ErrorEntropy", 2*0.00355348316330002 + 1e-14));
        }

        /// <summary>
        /// Control file for <see cref="IBMVortexLocalTimeSteppingTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlLocalTimeStepping() {
            IBMControl c = ControlTemplate(dgDegree: 2, divisions: 1, levelSetPosition: -0.25);

            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = 5;
            c.AgglomerationThreshold = 0.0;

            c.ExplicitOrder = 3;
            c.ExplicitScheme = ExplicitSchemes.LTS;
            c.NumberOfSubGrids = 3;

            return c;
        }

        /// <summary>
        /// Tests the edge case where two layers of cut cells are directly next
        /// to each other
        /// </summary>
        [Test]
        public static void IBMVortexCutNextToCutNoAgglomerationTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlCutNextToCutNoAgglomeration());
            p.RunSolverMode();
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2.8e-3),
                Tuple.Create("L2ErrorPressure", 4.0e-3),
                Tuple.Create("L2ErrorEntropy", 7.7e-3));
        }

        /// <summary>
        /// See <see cref="IBMVortexCutNextToCutNoAgglomerationTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlCutNextToCutNoAgglomeration() {
            return ControlCutNextToCut(agglomeration: false);
        }

        /// <summary>
        /// Tests the edge case where two layers of cut cells are directly next
        /// to each other
        /// </summary>
        [Test]
        public static void IBMVortexCutNextToCutAgglomerationTest() {
            Program<IBMControl> p = new IBMIsentropicVortexTest();
            p.Init(ControlCutNextToCutAgglomeration());
            p.RunSolverMode();

            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 3.0e-3),
                Tuple.Create("L2ErrorPressure", 4.2e-3),
                Tuple.Create("L2ErrorEntropy", 7.7e-3));
        }

        /// <summary>
        /// See <see cref="IBMVortexCutNextToCutAgglomerationTest"/>
        /// </summary>
        /// <returns></returns>
        public static IBMControl ControlCutNextToCutAgglomeration() {
            return ControlCutNextToCut(agglomeration: true);
        }

        private static IBMControl ControlTemplate(int dgDegree, int divisions, double levelSetPosition) {
            IBMControl c = new IBMControl();
            c.savetodb = false;

            double vortexSpeed = 1.0;

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.Rusanov;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                int noOfCellsPerDirection = (2 << divisions) * 10;
                GridCommons grid = Grid2D.Cartesian2DGrid(
                    GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1),
                    GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1));
                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            c.LevelSetBoundaryTag = "supersonicInlet";

            IsentropicVortexExactSolution solution = new IsentropicVortexExactSolution(c, vortexSpeed);

            c.InitialValues_Evaluators.Add(Variables.Density, X => solution.rho()(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => solution.u()(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => solution.v()(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => solution.p()(X, 0.0));

            c.LevelSetFunction = (X, t) => X[1] - levelSetPosition;

            c.AddBoundaryValue("supersonicInlet", Variables.Density, solution.rho());
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], solution.u());
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], solution.v());
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, solution.p());

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(Variables.Density, solution.rho()));
            c.Queries.Add("L2ErrorPressure", IBMQueries.L2Error(state => state.Pressure, solution.p()));
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.2;
            c.Endtime = 0.5;
            c.NoOfTimesteps = 100;

            return c;
        }

        private static IBMControl ControlCutNextToCut(bool agglomeration) {
            IBMControl c = new IBMControl();
            c.savetodb = false;

            int dgDegree = 2;
            double vortexSpeed = 1.0;

            // IBM Settings
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = 5;
            c.AgglomerationThreshold = agglomeration ? 0.3 : 0.0;

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.Rusanov;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;
            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 2);

            c.GridFunc = delegate {
                GridCommons grid = Grid2D.Cartesian2DGrid(
                    GenericBlas.Linspace(-5.0, 5.0, 21),
                    GenericBlas.Linspace(-0.5, 0.5, 3));
                grid.EdgeTagNames.Add(1, "supersonicInlet");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            c.LevelSetBoundaryTag = "supersonicInlet";

            IsentropicVortexExactSolution solution = new IsentropicVortexExactSolution(c, vortexSpeed);

            c.InitialValues_Evaluators.Add(Variables.Density, X => solution.rho()(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.xComponent, X => solution.u()(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Velocity.yComponent, X => solution.v()(X, 0.0));
            c.InitialValues_Evaluators.Add(CNSVariables.Pressure, X => solution.p()(X, 0.0));
            if (agglomeration) {
                c.LevelSetFunction = (X, t) => -((X[1] - 0.17) * (X[1] - 0.17) - 0.1);
            } else {
                c.LevelSetFunction = (X, t) => -(X[1] * X[1] - 0.2);
            }


            c.AddBoundaryValue("supersonicInlet", Variables.Density, solution.rho());
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[0], solution.u());
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Velocity[1], solution.v());
            c.AddBoundaryValue("supersonicInlet", CNSVariables.Pressure, solution.p());

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(Variables.Density, solution.rho()));
            c.Queries.Add("L2ErrorPressure", IBMQueries.L2Error(state => state.Pressure, solution.p()));
            c.Queries.Add("L2ErrorEntropy", IBMQueries.L2Error(state => state.Entropy, (X, t) => 1.0));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.2;
            c.Endtime = double.MaxValue;
            c.NoOfTimesteps = 100;

            return c;
        }
    }
}
