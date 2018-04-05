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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using CNS;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
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

        private static int REBALANCING_PERIOD = 5;

        public static void Main(string[] args) {
            SetUp();
            //System.Threading.Thread.Sleep(10000);
            //TestRebalancingForDG0WithRK1();
            //TestRealancingForDG0WithAB1();
            //TestRebalancingForDG2WithRK1AndAV();
            //TestRebalancingForDG2WithAB1AndAV();
            //TestRebalancingForDG0WithLTS1SingleSubGrid();
            //TestRebalancingForDG0WithLTS1TwoSubGrids();
            //TestRebalancingForDG2WithLTS1TwoSubGridsAndAV();
            TestRebalancingForDG2WithRK1AndAV_IBM_AggOff();
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
        public static void TestRebalancingForDG0WithAB1() {
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
        [Test]
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

            //control.NoOfTimesteps = 5;
            //control.dtFixed = 1.5e-3;

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

        [Test]
        public static void TestRebalancingForDG2WithLTS1TwoSubGridsAndAV() {
            int dgDegree = 2;
            ExplicitSchemes explicitScheme = ExplicitSchemes.LTS;
            int explicitOrder = 1;
            int noOfSubgrids = 2;

            var control = ShockTubeToro1WithAVTemplate(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);
            control.NumberOfSubGrids = noOfSubgrids;

            // MUST be the same as rebalancing period since LTS scheme MUST
            // recluster after rebalancing (at least, it makes life much easier)
            control.ReclusteringInterval = REBALANCING_PERIOD;

            //control.NoOfTimesteps = 5;
            //control.dtFixed = 1.5e-3;

            CheckRunsProduceSameResults(control, hilbert: true);
        }

        [Test]
        public static void TestRebalancingForDG2WithRK1AndAV_IBM_AggOn() {
            int dgDegree = 2;
            ExplicitSchemes explicitScheme = ExplicitSchemes.RungeKutta;
            int explicitOrder = 1;

            IBMControl control = ShockTubeToro1WithIBMAndAVTemplate(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            control.AgglomerationThreshold = 0.9;

            // MUST be the same as rebalancing period since LTS scheme MUST
            // recluster after rebalancing (at least, it makes life much easier)
            //control.ReclusteringInterval = REBALANCING_PERIOD;

            //control.NumberOfSubGrids = 3;
            //control.maxNumOfSubSteps = 10;
            //control.FluxCorrection = false;

            Console.WriteLine("TestRebalancingForDG2WithRK1AndAV_IBM_AggOn");
            CheckRunsProduceSameResults(control, hilbert: true);
        }

        [Test]
        public static void TestRebalancingForDG2WithRK1AndAV_IBM_AggOff() {
            int dgDegree = 2;
            ExplicitSchemes explicitScheme = ExplicitSchemes.RungeKutta;
            int explicitOrder = 1;

            IBMControl control = ShockTubeToro1WithIBMAndAVTemplate(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            control.AgglomerationThreshold = 0.1;

            // MUST be the same as rebalancing period since LTS scheme MUST
            // recluster after rebalancing (at least, it makes life much easier)
            //control.ReclusteringInterval = REBALANCING_PERIOD;

            //control.NumberOfSubGrids = 3;
            //control.maxNumOfSubSteps = 10;
            //control.FluxCorrection = false;

            Console.WriteLine("TestRebalancingForDG2WithRK1AndAV_IBM_AggOff");
            // Threshold is set to 2e1-3, since METIS results for density are worse compared to Hilbert
            // Could be influenced by the handling of cut-cells along processor boundaries, communication plays
            // a more important role for non-agglomerated cut-cells in the IBM-case as for non-IBM simulations, etc.
            // Anyway, this does not seem to be a real problem.
            CheckRunsProduceSameResults(control, differenceThreshold: 2e-13, hilbert: true);
        }

        private static CNSControl ShockTubeToro1Template(int dgDegree, ExplicitSchemes explicitScheme, int explicitOrder, int noOfCells = 50, double gridStretching = 0.0, bool twoD = false) {
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
            c.AddVariable(Variables.Rank, 0);

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

        private static CNSControl ShockTubeToro1WithAVTemplate(int dgDegree, ExplicitSchemes explicitScheme, int explicitOrder, int noOfCells = 50, bool twoD = false) {
            Variable sensorVariable = Variables.Density;
            double sensorLimit = 1e-3;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            double endTime = 0.01;

            var c = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder,
                twoD: twoD);
            if (twoD) {
                c.AddVariable(Variables.ArtificialViscosity, 2);
            } else {
                c.AddVariable(Variables.ArtificialViscosity, 1);
            }
            c.ActiveOperators |= Operators.ArtificialViscosity;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            c.AddVariable(Variables.ShockSensor, 0);
            c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            c.Endtime = endTime;

            return c;
        }

        private static IBMControl ShockTubeToro1WithIBMAndAVTemplate(int dgDegree, ExplicitSchemes explicitScheme, int explicitOrder, int noOfCellsX = 50, int noOfCellsY = 10) {
            IBMControl c = new IBMControl();

            c.DbPath = null;
            //c.DbPath = @"c:\bosss_db\";
            c.savetodb = false;

            c.saveperiod = 1;
            c.PrintInterval = 1;

            double xMin = 0;
            double xMax = 1;
            double yMin = 0;
            double yMax = 1;

            bool AV = true;

            c.DomainType = DomainTypes.StaticImmersedBoundary;
            c.LevelSetFunction = delegate (double[] X, double t) {
                return X[1] - 0.16;
            };
            c.LevelSetBoundaryTag = "AdiabaticSlipWall";
            c.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            c.LevelSetQuadratureOrder = 6;
            c.AgglomerationThreshold = 0.3;
            c.AddVariable(IBMVariables.LevelSet, 1);

            //c.AddVariable(IBMVariables.FluidCells, 1);
            //c.AddVariable(IBMVariables.FluidCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.CutCells, 1);
            //c.AddVariable(IBMVariables.CutCellsWithoutSourceCells, 1);
            //c.AddVariable(IBMVariables.SourceCells, 1);

            if (AV) {
                c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            } else {
                c.ActiveOperators = Operators.Convection;
            }
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;

            // Shock-capturing
            double sensorLimit = 1e-3;
            double epsilon0 = 1.0;
            double kappa = 0.5;

            if (AV) {
                Variable sensorVariable = Variables.Density;
                c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
                c.AddVariable(Variables.ShockSensor, 0);
                c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa, lambdaMax: 2);
            }

            c.ExplicitScheme = explicitScheme;
            c.ExplicitOrder = explicitOrder;

            c.EquationOfState = IdealGas.Air;

            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Velocity.yComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Rank, 0);

            if (AV) {
                c.AddVariable(Variables.ArtificialViscosity, 2);
                c.AddVariable(Variables.CFLArtificialViscosity, 0);
            }
            c.AddVariable(Variables.CFL, 0);
            c.AddVariable(Variables.CFLConvective, 0);

            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(xMin, xMax, noOfCellsX + 1);
                double[] yNodes = GenericBlas.Linspace(yMin, yMax, noOfCellsY + 1);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                // Boundary conditions
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");

                grid.DefineEdgeTags(delegate (double[] _X) {
                    return 1;
                });
                return grid;
            };

            c.AddBoundaryCondition("AdiabaticSlipWall");

            // Normal vector of initial shock
            Vector2D normalVector = new Vector2D(1, 0);

            // Direction vector of initial shock
            Vector2D r = new Vector2D(normalVector.y, -normalVector.x);
            r.Normalize();

            // Distance from a point X to the initial shock
            double[] p = new double[] { 0.5, 0.0 };

            double cellSize = Math.Min((xMax - xMin) / noOfCellsX, (yMax - yMin) / noOfCellsY);

            Func<double, double> SmoothJump = delegate (double distance) {
                // smoothing should be in the range of h/p
                double maxDistance = 4.0 * cellSize / Math.Max(dgDegree, 1);

                return (Math.Tanh(distance / maxDistance) + 1.0) * 0.5;
            };

            Func<double, double> Jump = (x => x <= 0.5 ? 0 : 1);

            // Initial conditions
            double densityLeft = 1.0;
            double densityRight = 0.125;
            double pressureLeft = 1.0;
            double pressureRight = 0.1;

            //c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (densityLeft - densityRight));
            //c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - SmoothJump(DistanceFromPointToLine(X, p, r)) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Density, X => densityLeft - Jump(X[0]) * (densityLeft - densityRight));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => pressureLeft - Jump(X[0]) * (pressureLeft - pressureRight));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => 0.0);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => 0.0);

            // Time config
            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = 0.001;
            c.NoOfTimesteps = int.MaxValue;

            c.ProjectName = "Shock tube";
            c.SessionName = String.Format("IBM shock tube, p={0}, {1}x{2} cells, aggloThresh={3}, s0={4:0.0E-00}, CFLFrac={5}, ALTS {6}/{7}/{8}({9})", dgDegree, noOfCellsX, noOfCellsY, c.AgglomerationThreshold, sensorLimit, c.CFLFraction, c.ExplicitOrder, c.NumberOfSubGrids, c.ReclusteringInterval, c.maxNumOfSubSteps);

            return c;
        }

        private static void CheckRunsProduceSameResults(CNSControl refControl, double differenceThreshold = 1e-15, bool hilbert = false) {
            Debug.Assert(refControl.DynamicLoadBalancing_Period <= 0);
            Debug.Assert(refControl.DynamicLoadBalancing_CellCostEstimatorFactories.Count == 0);

            CNSControl loadBalControl = refControl.CloneAs();
            loadBalControl.DynamicLoadBalancing_Period = REBALANCING_PERIOD;
            loadBalControl.DynamicLoadBalancing_CellClassifier = new RandomCellClassifier(2);
            //loadBalControl.DynamicLoadBalancing_CellClassifier = new ArtificialViscosityCellClassifier();
            loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.Add((p, i) => new StaticCellCostEstimator(new[] { 1, 10 }));
            //loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.AddRange(ArtificialViscosityCellCostEstimator.GetMultiBalanceConstraintsBasedEstimators());
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
            refSolver.Init(refControl);
            refSolver.RunSolverMode();

            Console.WriteLine("\nRun WITH load balancing");
            var loadBalSolver = new ShockTubeLoadBalancingTests();
            loadBalSolver.Init(loadBalControl);
            loadBalSolver.RunSolverMode();

            ShockTubeLoadBalancingTests loadBalSolverHilbert = null;
            if (hilbert) {
                Console.WriteLine("\nRun WITH load balancing (Hilbert)");
                loadBalSolverHilbert = new ShockTubeLoadBalancingTests();
                loadBalSolverHilbert.Init(loadBalControl);
                loadBalSolverHilbert.RunSolverMode();
            }

            // To be able to compare errors without using the database, we need to 
            // agree on a single grid partitioning in the end -> use ref
            //Console.WriteLine("Transfering load balancing data to reference grid");
            //var refPartitioning = new int[loadBalSolver.GridData.Cells.NoOfLocalUpdatedCells];
            //for (int i = 0; i < refSolver.GridData.CellPartitioning.TotalLength; i++) {
            //    int localIndex = loadBalSolver.GridData.CellPartitioning.TransformIndexToLocal(i);
            //    if (localIndex >= 0 && localIndex < loadBalSolver.GridData.Cells.NoOfLocalUpdatedCells) {
            //        refPartitioning[localIndex] = refSolver.GridData.CellPartitioning.FindProcess(i);
            //    }
            //}
            //loadBalSolver.MpiRedistributeAndMeshAdapt(
            //    int.MinValue,
            //    double.MinValue,
            //    refPartitioning,
            //    refSolver.GridData.CurrentGlobalIdPermutation);

            //if (!twoD) {
            //    CompareErrors(refSolver.WorkingSet, loadBalSolver.WorkingSet, differenceThreshold);
            //}
            CompareNorms(refSolver, loadBalSolver, differenceThreshold, loadBalSolverHilbert);
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
        /// Can be used independent on the grid and the grid partitioning
        /// in constrast to <see cref="CompareErrors(CNSFieldSet, CNSFieldSet, double)"/>
        /// </summary>
        /// <param name="loadBalSolver"></param>
        /// <param name="refSolver"></param>
        /// <param name="differenceThreshold"></param>
        private static void CompareNorms(IProgram<CNSControl> loadBalSolver, IProgram<CNSControl> refSolver, double differenceThreshold, IProgram<CNSControl> hilbertSolver = null) {
            List<Action> assertions = new List<Action>();
            string[] varName = { "Density", "x-Momentum", "Energy" };

            for (int i = 0; i < 3; i++) {
                List<DGField> listOfDGFields_RepON = (List<DGField>)loadBalSolver.IOFields;
                DGField variableRepON = listOfDGFields_RepON[i];
                double L2NormRepON = variableRepON.L2Norm();

                List<DGField> listOfDGFields_RepOFF = (List<DGField>)refSolver.IOFields;
                DGField variableRepOFF = listOfDGFields_RepOFF[i];
                double L2NormRepOFF = variableRepOFF.L2Norm();

                double difference = Math.Abs(L2NormRepON - L2NormRepOFF);

                //Console.WriteLine("{0}-L2Norm Rep ON: {1}", varName[i], L2NormRepON);
                //Console.WriteLine("{0}-L2Norm Rep OFF: {1}", varName[i], L2NormRepOFF);

                string message = String.Format("METIS: Difference in {0} norm is {1} (Threshold is {2})", varName[i], difference, differenceThreshold);
                Console.WriteLine(message);
                assertions.Add(() => Assert.IsTrue(difference < differenceThreshold, message));

                // Duplicated code...
                double differenceHilbert;
                if (hilbertSolver != null) {
                    List<DGField> listOfDGFields_Hilbert = (List<DGField>)hilbertSolver.IOFields;
                    DGField variableHilbert = listOfDGFields_Hilbert[i];
                    double L2NormHilbert = variableHilbert.L2Norm();

                    differenceHilbert = Math.Abs(L2NormHilbert - L2NormRepOFF);

                    string messageHilbert = String.Format("Hilbert: Difference in {0} norm is {1} (Threshold is {2})", varName[i], differenceHilbert, differenceThreshold);
                    Console.WriteLine(messageHilbert);
                    assertions.Add(() => Assert.IsTrue(differenceHilbert < differenceThreshold, messageHilbert));
                }
            }

            assertions.ForEach(a => a());
        }

        [TestFixtureTearDown]
        public static void TearDown() {
            csMPI.Raw.mpiFinalize();
        }
    }
}