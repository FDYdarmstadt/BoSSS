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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using CNS;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.LoadBalancing;
using CNS.MaterialProperty;
using CNS.ShockCapturing;
using CNS.Tests;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CNS_MPITests.Tests.LoadBalancing {

    public class ShockTubeLoadBalancingTests : TestProgram<CNSControl> {

        //private static CommandLineOptions commandLineOptions = new Application<CNSControl>.CommandLineOptions() {
        //    delPlt = true,
        //    ImmediatePlotPeriod = 1
        //};

        private static CommandLineOptions commandLineOptions = null;

        private static int REBALANCING_PERIOD = 5;

        private static bool twoD = false;

        public static void Main(string[] args) {
            SetUp();
            //TestRebalancingForDG0WithRK1();
            //TestRealancingForDG0WithAB1();
            //TestRebalancingForDG0WithLTS1SingleSubGrid();
            //TestRebalancingForDG0WithLTS1TwoSubGrids();
            TestRebalancingForDG2WithRK1AndAV();
            //TestRebalancingForDG2WithAB1AndAV();
            TearDown();
        }

        [Test]
        public static void TestRebalancingForDG0WithRK1() {
            int dgDegree = 0;
            ExplicitSchemes explicitScheme = ExplicitSchemes.RungeKutta;
            int explicitOrder = 1;

            var control = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            CheckRunsProduceSameResults(control);
        }

        [Test]
        public static void TestRealancingForDG0WithAB1() {
            int dgDegree = 0;
            ExplicitSchemes explicitScheme = ExplicitSchemes.AdamsBashforth;
            int explicitOrder = 1;

            var control = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            CheckRunsProduceSameResults(control);
        }

        [Test]
        public static void TestRebalancingForDG0WithLTS1SingleSubGrid() {
            int dgDegree = 0;
            ExplicitSchemes explicitScheme = ExplicitSchemes.LTS;
            int explicitOrder = 1;

            var control = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            CheckRunsProduceSameResults(control);
        }

        /// <summary>
        /// This test is currently deactivated because it fails; probably for
        /// the following reason: the reclustering delivers different results
        /// before and after load bal
        /// </summary>
        //[Test]
        public static void TestRebalancingForDG0WithLTS1TwoSubGrids() {
            int dgDegree = 0;
            ExplicitSchemes explicitScheme = ExplicitSchemes.LTS;
            int explicitOrder = 1;
            int noOfSubgrids = 2;
            double gridStretching = 1.0;

            var control = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder,
                gridStretching: gridStretching);
            control.NumberOfSubGrids = noOfSubgrids;

            // MUST be the same as rebalancing period since LTS scheme MUST
            // recluster after rebalancing (at least, it makes life much easier)
            control.ReclusteringInterval = REBALANCING_PERIOD;

            control.NoOfTimesteps = 5;
            control.dtFixed = 1.5e-3;

            CheckRunsProduceSameResults(control);
        }

        [Test]
        public static void TestRebalancingForDG2WithRK1AndAV() {
            int dgDegree = 2;
            ExplicitSchemes explicitScheme = ExplicitSchemes.RungeKutta;
            int explicitOrder = 1;

            var control = ShockTubeToro1WithAVTemplate(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            CheckRunsProduceSameResults(control);
        }

        [Test]
        public static void TestRebalancingForDG2WithAB1AndAV() {
            int dgDegree = 2;
            ExplicitSchemes explicitScheme = ExplicitSchemes.AdamsBashforth;
            int explicitOrder = 1;

            var control = ShockTubeToro1WithAVTemplate(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            CheckRunsProduceSameResults(control);
        }

        //[Test]
        public static void TestRebalancingForDG2WithLTS1AndAV() {
            throw new NotImplementedException("to do");
        }

        private static CNSControl ShockTubeToro1Template(int dgDegree, ExplicitSchemes explicitScheme, int explicitOrder, int noOfCells = 50, double gridStretching = 0.0) {
            double densityLeft = 1.0;
            double velocityLeft = 0.0;
            double pressureLeft = 1.0;
            double densityRight = 0.125;
            double velocityRight = 0.0;
            double pressureRight = 0.1;
            double discontinuityPosition = 0.5;

            CNSControl c = new CNSControl();
            c.DbPath = null;
            //c.DbPath = @"c:\bosss_db\";
            c.savetodb = false;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.ExplicitScheme = explicitScheme;
            c.ExplicitOrder = explicitOrder;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);

            c.GridFunc = delegate {
                double xMin = 0.0;
                double xMax = 1.0;
                double yMin = 0.0;
                double yMax = 1.0;

                double[] xNodes;
                double[] yNodes;
                if (gridStretching > 0.0) {
                    xNodes = Grid1D.TanhSpacing(xMin, xMax, noOfCells + 1, gridStretching, true);
                    yNodes = Grid1D.TanhSpacing(yMin, yMax, 1 + 1, gridStretching, true);
                } else {
                    xNodes = GenericBlas.Linspace(xMin, xMax, noOfCells + 1);
                    yNodes = GenericBlas.Linspace(yMin, yMax, 1 + 1);
                }

                GridCommons grid;
                if (twoD) {
                    grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                } else {
                    grid = Grid1D.LineGrid(xNodes, periodic: false);
                }

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(delegate (double[] _X) {
                    return 1;
                });
                return grid;
            };
            c.AddBoundaryCondition("AdiabaticSlipWall");

            Material material = new Material(c);
            StateVector stateLeft = StateVector.FromPrimitiveQuantities(
                material, densityLeft, new Vector3D(velocityLeft, 0.0, 0.0), pressureLeft);
            StateVector stateRight = StateVector.FromPrimitiveQuantities(
                material, densityRight, new Vector3D(velocityRight, 0.0, 0.0), pressureRight);

            c.InitialValues_Evaluators.Add(
                    Variables.Density,
                    X => stateLeft.Density + (stateRight.Density - stateLeft.Density) * (X[0] - discontinuityPosition).Heaviside());
            c.InitialValues_Evaluators.Add(
                Variables.Velocity.xComponent,
                X => stateLeft.Velocity.x + (stateRight.Velocity.x - stateLeft.Velocity.x) * (X[0] - discontinuityPosition).Heaviside());
            c.InitialValues_Evaluators.Add(
                Variables.Pressure,
                X => stateLeft.Pressure + (stateRight.Pressure - stateLeft.Pressure) * (X[0] - discontinuityPosition).Heaviside());
            if (twoD) {
                c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0);
            }

            if (!twoD) {
                var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector3D(1.0, 0.0, 0.0));
                riemannSolver.GetStarRegionValues(out double pStar, out double uStar);

                c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
                    Variables.Density,
                    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
                c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
                    Variables.Velocity.xComponent,
                    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
                c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
                    Variables.Pressure,
                    (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));
            }

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            //c.dtFixed = 1.0e-6;
            c.CFLFraction = 0.1;
            c.Endtime = 0.2;
            c.NoOfTimesteps = int.MaxValue;

            // Use METIS since ParMETIS is not installed on build server
            c.GridPartType = GridPartType.METIS;

            return c;
        }

        private static CNSControl ShockTubeToro1WithAVTemplate(int dgDegree, ExplicitSchemes explicitScheme, int explicitOrder, int noOfCells = 50) {
            Variable sensorVariable = Variables.Density;
            double sensorLimit = 1e-3;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            double endTime = 0.01;

            var c = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);
            if (twoD) {
                c.AddVariable(Variables.ArtificialViscosity, 2);
            } else {
                c.AddVariable(Variables.ArtificialViscosity, 1);
            }
            c.ActiveOperators |= Operators.ArtificialViscosity;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            c.Endtime = endTime;

            return c;
        }

        private static void CheckRunsProduceSameResults(CNSControl refControl, double differenceThreshold = 1e-16) {
            Debug.Assert(refControl.DynamicLoadBalancing_Period <= 0);
            Debug.Assert(refControl.DynamicLoadBalancing_CellCostEstimatorFactories.Count == 0);

            CNSControl loadBalControl = refControl.CloneAs();
            loadBalControl.DynamicLoadBalancing_Period = REBALANCING_PERIOD;
            loadBalControl.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 1, 10 }));
            loadBalControl.DynamicLoadBalancing_ImbalanceThreshold = 0.01;


            //// TEST ONLY SUCCEEDS IF THESE LINES ARE IN
            //loadBalControl.DynamicLoadBalancing_CellClassifier = new IndifferentCellClassifier();
            //loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.Clear();
            //loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.Add(CellCostEstimatorLibrary.AllCellsAreEqual);


            Debug.Assert(loadBalControl.DynamicLoadBalancing_Period > 0);
            Debug.Assert(loadBalControl.DynamicLoadBalancing_CellClassifier != null);
            Debug.Assert(loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.Count > 0);

            Console.WriteLine("Run WITHOUT load balancing");
            var refSolver = new ShockTubeLoadBalancingTests();
            refSolver.Init(refControl, commandLineOptions);
            refSolver.RunSolverMode();

            Console.WriteLine("Run WITH load balancing");
            var loadBalSolver = new ShockTubeLoadBalancingTests();
            loadBalSolver.Init(loadBalControl, commandLineOptions);
            loadBalSolver.RunSolverMode();

            // To be able to compare errors without using the database, we need to 
            // agree on a single grid partitioning in the end -> use ref
            Console.WriteLine("Transfering load balancing data to reference grid");
            var refPartitioning = new int[loadBalSolver.GridData.Cells.NoOfLocalUpdatedCells];
            for (int i = 0; i < refSolver.GridData.CellPartitioning.TotalLength; i++) {
                int localIndex = loadBalSolver.GridData.CellPartitioning.TransformIndexToLocal(i);
                if (localIndex >= 0 && localIndex < loadBalSolver.GridData.Cells.NoOfLocalUpdatedCells) {
                    refPartitioning[localIndex] = refSolver.GridData.CellPartitioning.FindProcess(i);
                }
            }
            loadBalSolver.MpiRedistributeAndMeshAdapt(
                int.MinValue,
                double.MinValue,
                refPartitioning,
                refSolver.GridData.CurrentGlobalIdPermutation);

            if (!twoD) {
                CompareErrors(refSolver.WorkingSet, loadBalSolver.WorkingSet, differenceThreshold);
            }
        }

        /// <summary>
        /// Note: Only to be used if <paramref name="refResults"/> and
        /// <paramref name="loadBalResults"/> have the same grid AND the same
        /// grid partitioning
        /// </summary>
        /// <param name="refResults"></param>
        /// <param name="loadBalResults"></param>
        /// <param name="differenceThreshold"></param>
        private static void CompareErrors(CNSFieldSet refResults, CNSFieldSet loadBalResults, double differenceThreshold) {
            List<Action> assertions = new List<Action>();
            {
                double densityDifference = refResults.Density.L2Error(loadBalResults.Density, overrideGridCheck: true);
                string densityMessage = String.Format(
                    "Density: {0} (Threshold is {1})",
                    densityDifference,
                    differenceThreshold);
                Console.WriteLine(densityMessage);
                assertions.Add(() => Assert.IsTrue(densityDifference < differenceThreshold, densityMessage));
            }

            for (int d = 0; d < refResults.Density.GridDat.SpatialDimension; d++) {
                double momentumDifference = refResults.Momentum[d].L2Error(loadBalResults.Momentum[d], overrideGridCheck: true);
                string momentumMessage = String.Format(
                    "Momentum[{0}]: {1} (Threshold is {2})",
                    d,
                    momentumDifference,
                    differenceThreshold);
                Console.WriteLine(momentumMessage);
                assertions.Add(() => Assert.IsTrue(momentumDifference < differenceThreshold, momentumMessage));
            }

            {
                double energyDifference = refResults.Energy.L2Error(loadBalResults.Energy, overrideGridCheck: true);
                string energyMessage = String.Format(
                    "Energy: {0} (Threshold is {1})",
                    energyDifference,
                    differenceThreshold);
                Console.WriteLine(energyMessage);
                assertions.Add(() => Assert.IsTrue(energyDifference < differenceThreshold, energyMessage));
            }

            assertions.ForEach(a => a());
        }

        /// <summary>
        /// </summary>
        [TestFixtureTearDown]
        public static void TearDown() {
            csMPI.Raw.mpiFinalize();
        }
    }
}
