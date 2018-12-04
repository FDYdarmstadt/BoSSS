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
using BoSSS.Platform.LinAlg;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using ilPSP.Utils;
using NUnit.Framework;
using System;

namespace CNS.Tests.MovingIBMTests {

    public class PistonTests : TestProgram<IBMControl> {

        ///// <summary>
        ///// Alternative entry point of this assembly that allows to perform
        ///// isentropic vortex tests conveniently.
        ///// </summary>
        ///// <param name="args">
        ///// Command line arguments
        ///// </param>
        //public static void Main(string[] args) {
        //    SetUp();
        //    //SplittingIBMPiston0thOrderNoAgglomeration();
        //    //SplittingIBMPiston1stOrderNoAgglomeration();
        //    //SplittingIBMPiston0thOrderWithAgglomeration();
        //    //SplittingIBMPiston1stOrderWithAgglomeration();
        //    //MovingMeshIBMPiston0thOrderNoAgglomeration();
        //    //MovingMeshIBMPiston1stOrderNoAgglomeration();
        //    //MovingMeshIBMPiston0thOrderWithAgglomeration();
        //    //MovingMeshIBMPiston1stOrderWithAgglomeration();
        //}

        [Test]
        public static void SplittingIBMPiston0thOrderNoAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 0,
                rkDegree: 1,
                convectiveFlux: ConvectiveFluxTypes.Rusanov,
                timeSteppingStrategy: TimesteppingStrategies.LieSplitting,
                agglomerationThreshold: 0.0);

            c.dtMin = 3.01E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1e-15),
                Tuple.Create("L2ErrorXMomentum", 1e-15),
                Tuple.Create("L2ErrorYMomentum", 1e-15),
                Tuple.Create("L2ErrorPressure", 1e-15));
        }

        [Test]
        public static void SplittingIBMPiston1stOrderNoAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 1,
                rkDegree: 1,
                convectiveFlux: ConvectiveFluxTypes.Rusanov,
                timeSteppingStrategy: TimesteppingStrategies.LieSplitting,
                agglomerationThreshold: 0.0);

            c.dtMin = 1.0E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1e-13),
                Tuple.Create("L2ErrorXMomentum", 1e-13),
                Tuple.Create("L2ErrorYMomentum", 1e-13),
                Tuple.Create("L2ErrorPressure", 1e-13));
        }

        [Test]
        public static void SplittingIBMPiston0thOrderWithAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 0,
                rkDegree: 1,
                convectiveFlux: ConvectiveFluxTypes.Rusanov,
                timeSteppingStrategy: TimesteppingStrategies.LieSplitting,
                agglomerationThreshold: 0.2);

            c.dtMin = 3.01E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1e-15),
                Tuple.Create("L2ErrorXMomentum", 1e-15),
                Tuple.Create("L2ErrorYMomentum", 1e-15),
                Tuple.Create("L2ErrorPressure", 1e-15));
        }

        [Test]
        public static void SplittingIBMPiston1stOrderWithAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 1,
                rkDegree: 1,
                convectiveFlux: ConvectiveFluxTypes.Rusanov,
                timeSteppingStrategy: TimesteppingStrategies.LieSplitting,
                agglomerationThreshold: 0.2);

            c.dtMin = 1.01E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1e-14),
                Tuple.Create("L2ErrorXMomentum", 1e-14),
                Tuple.Create("L2ErrorYMomentum", 1e-14),
                Tuple.Create("L2ErrorPressure", 1e-14));
        }

        [Test]
        public static void MovingMeshIBMPiston0thOrderNoAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 0,
                rkDegree: 1,
                convectiveFlux: ConvectiveFluxTypes.MovingFrameRusanov,
                timeSteppingStrategy: TimesteppingStrategies.MovingFrameFlux,
                agglomerationThreshold: 0.0);

            c.dtMin = 3.01E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1e-14),
                Tuple.Create("L2ErrorXMomentum", 1e-14),
                Tuple.Create("L2ErrorYMomentum", 1e-14),
                Tuple.Create("L2ErrorPressure", 1e-14));
        }

        [Test]
        public static void MovingMeshIBMPiston1stOrderNoAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 1,
                rkDegree: 2,
                convectiveFlux: ConvectiveFluxTypes.MovingFrameRusanov,
                timeSteppingStrategy: TimesteppingStrategies.MovingFrameFlux,
                agglomerationThreshold: 0.0);

            c.dtMin = 1.01E-2;
            c.dtMax = c.dtMin;

            // Don't cross cell border
            c.Endtime = 0.5;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 4.5e-3),
                Tuple.Create("L2ErrorXMomentum", 4.0e-3),
                Tuple.Create("L2ErrorYMomentum", 5.9e-6),
                Tuple.Create("L2ErrorPressure", 4.5e-3));
        }

        [Test]
        public static void MovingMeshIBMPiston0thOrderWithAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 0,
                rkDegree: 1,
                convectiveFlux: ConvectiveFluxTypes.MovingFrameRusanov,
                timeSteppingStrategy: TimesteppingStrategies.MovingFrameFlux,
                agglomerationThreshold: 0.2);

            c.dtMin = 3.01E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1e-14),
                Tuple.Create("L2ErrorXMomentum", 1e-14),
                Tuple.Create("L2ErrorYMomentum", 1e-14),
                Tuple.Create("L2ErrorPressure", 1e-14));
        }

        [Test]
        public static void MovingMeshIBMPiston1stOrderWithAgglomeration() {
            IBMControl c = PistonControl(
                dgDegree: 1,
                rkDegree: 2,
                convectiveFlux: ConvectiveFluxTypes.MovingFrameRusanov,
                timeSteppingStrategy: TimesteppingStrategies.MovingFrameFlux,
                agglomerationThreshold: 0.1);

            c.dtMin = 1.0E-2;
            c.dtMax = c.dtMin;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            CheckErrorThresholds(
                solver.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 3.4e-4),
                Tuple.Create("L2ErrorXMomentum", 2.5e-4),
                Tuple.Create("L2ErrorYMomentum", 5.0e-5),
                Tuple.Create("L2ErrorPressure", 2.8e-4));
        }

        public static IBMControl PistonControl(int dgDegree, int rkDegree, ConvectiveFluxTypes convectiveFlux, TimesteppingStrategies timeSteppingStrategy, double agglomerationThreshold) {
            double pistonVelocity = 1.0;
            double initialLevelSetPosition = 0.1;
            //double initialLevelSetPosition = 0.5;

            IBMControl c = new IBMControl();
            c.DbPath = null;
            c.savetodb = false;

            c.ProjectName = "Piston";
            c.ProjectDescription = "Vertical moving at flow velocity through constant flow field";

            c.DomainType = DomainTypes.MovingImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = convectiveFlux;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.TimesteppingStrategy = timeSteppingStrategy;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = rkDegree;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(0.0, 2.0, 4);
                double[] yNodes = GenericBlas.Linspace(-1.0, 1.0, 4);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: true);
                grid.EdgeTagNames.Add(1, "adiabaticSlipWall");
                grid.EdgeTagNames.Add(2, "supersonicInlet");
                grid.DefineEdgeTags(X => 2);
                return grid;
            };

            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.LevelSetQuadratureOrder = 10;
            c.LevelSetBoundaryTag = "adiabaticSlipWall";
            c.AgglomerationThreshold = agglomerationThreshold;

            c.InitialValues_Evaluators.Add(Variables.Density, X => 1.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => pistonVelocity);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => 1.0);

            c.LevelSetFunction = delegate (double[] X, double time) {
                double newLevelSetPosition = initialLevelSetPosition + pistonVelocity * time;
                return X[0] - newLevelSetPosition;
            };
            c.LevelSetVelocity = (X, t) => new Vector(pistonVelocity, 0.0);

            c.AddBoundaryValue("adiabaticSlipWall", Variables.Velocity.xComponent, X => pistonVelocity);
            c.AddBoundaryValue("adiabaticSlipWall", Variables.Velocity.yComponent, X => 0.0);
            c.AddBoundaryValue("supersonicInlet", Variables.Density, X => 1.0);
            c.AddBoundaryValue("supersonicInlet", Variables.Velocity[0], X => pistonVelocity);
            c.AddBoundaryValue("supersonicInlet", Variables.Velocity[1], X => 0.0);
            c.AddBoundaryValue("supersonicInlet", Variables.Pressure, X => 1.0);

            c.Queries.Add("L2ErrorDensity", IBMQueries.L2Error(Variables.Density, (X, t) => 1.0));
            c.Queries.Add("L2ErrorXMomentum", IBMQueries.L2Error(Variables.Momentum.xComponent, (X, t) => 1.0));
            c.Queries.Add("L2ErrorYMomentum", IBMQueries.L2Error(Variables.Momentum.yComponent, (X, t) => 0.0));
            c.Queries.Add("L2ErrorPressure", IBMQueries.L2Error(state => state.Pressure, (X, t) => 1.0));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = 0.75;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }
    }
}
