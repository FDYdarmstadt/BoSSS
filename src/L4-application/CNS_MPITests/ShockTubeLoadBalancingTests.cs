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
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CNS.Tests.LoadBalancing {

    public class ShockTubeLoadBalancingTests : TestProgram<CNSControl> {


        static CommandLineOptions bla = new Application<CNSControl>.CommandLineOptions() {
            delPlt = true,
            ImmediatePlotPeriod = 1
        };

        public static void Main(string[] args) {
            SetUp();
            TestLoadBalancingForRK1();
            csMPI.Raw.mpiFinalize();
        }

        private static void CheckRunsProduceSameResults(CNSControl refControl, CNSControl loadBalControl, double differenceThreshold = 1e-16) {
            //Debug.Assert(refControl.DynamicLoadBalancing_Period <= 0);
            //Debug.Assert(refControl.DynamicLoadBalancing_CellCostEstimatorFactories.Count == 0);

            //Debug.Assert(loadBalControl.DynamicLoadBalancing_Period > 0);
            //Debug.Assert(loadBalControl.DynamicLoadBalancing_CellClassifier != null);
            //Debug.Assert(loadBalControl.DynamicLoadBalancing_CellCostEstimatorFactories.Count > 0);

            Console.WriteLine("Run WITHOUT load balancing");
            var solver = new Program();
            solver.Init(refControl);
            //solver.Init(refControl, bla);
            solver.RunSolverMode();
            var refResults = solver.WorkingSet;

            Console.WriteLine("Run WITH load balancing");
            solver = new Program();
            solver.Init(loadBalControl);
            //solver.Init(loadBalControl, bla);
            solver.RunSolverMode();
            var loadBalResults = solver.WorkingSet;

            CompareErrors(refResults, loadBalResults, differenceThreshold);
        }

        protected static void CompareErrors(CNSFieldSet refResults, CNSFieldSet loadBalResults, double differenceThreshold) {
            List<Action> assertions = new List<Action>();

            {
                Debugger.Launch();

                refResults.Density.MPIExchange();
                loadBalResults.Density.MPIExchange();

                CellQuadratureScheme scheme = new CellQuadratureScheme();
                var rule = scheme.SaveCompile(refResults.Density.GridDat, 10);

                double[] refNorms = refResults.Density.LocalLxError((ScalarFunction)null, null, rule);
                double[] loadBalNorms = loadBalResults.Density.LocalLxError((ScalarFunction)null, null, rule);
                double maxDiff = 0.0;

                ilPSP.Environment.StdoutOnlyOnRank0 = false;

                Console.WriteLine(((GridData)refResults.Density.GridDat).MpiRank + " ref norms: " + refNorms.Length);
                Console.WriteLine(((GridData)loadBalResults.Density.GridDat).MpiRank + " load bal norms: " + loadBalNorms.Length);

                for (int i = 0; i < refResults.Density.GridDat.iLogicalCells.NoOfLocalUpdatedCells; i++) {
                    maxDiff = Math.Max(maxDiff, Math.Abs(refNorms[i] - loadBalNorms[i]));
                    //if (Math.Abs(refNorms[i] - loadBalNorms[i]) > differenceThreshold) {
                    //    throw new Exception();
                    //}
                }

                maxDiff = maxDiff.MPIMax();

                Console.WriteLine("Max diff: " + maxDiff);

                double densityDifference = refResults.Density.L2Error(loadBalResults.Density, quadratureDegree: 100);
                string densityMessage = String.Format(
                    "Density: {0} (Threshold is {1})",
                    densityDifference,
                    differenceThreshold);
                Console.WriteLine(densityMessage);
                assertions.Add(() => Assert.IsTrue(densityDifference < differenceThreshold, densityMessage));
            }

            for (int d = 0; d < refResults.Density.GridDat.SpatialDimension; d++) {
                double momentumDifference = refResults.Momentum[d].L2Error(loadBalResults.Momentum[d], quadratureDegree: 10);
                string momentumMessage = String.Format(
                    "Momentum[{0}]: {1} (Threshold is {2})",
                    d,
                    momentumDifference,
                    differenceThreshold);
                Console.WriteLine(momentumMessage);
                assertions.Add(() => Assert.IsTrue(momentumDifference < differenceThreshold, momentumMessage));
            }

            {
                double energyDifference = refResults.Energy.L2Error(loadBalResults.Energy, quadratureDegree: 10);
                string energyMessage = String.Format(
                    "Energy: {0} (Threshold is {1})",
                    energyDifference,
                    differenceThreshold);
                Console.WriteLine(energyMessage);
                assertions.Add(() => Assert.IsTrue(energyDifference < differenceThreshold, energyMessage));
            }

            assertions.ForEach(a => a());
        }

        //[Test]
        public static void TestLoadBalancingForRK1() {
            int dgDegree = 0;
            ExplicitSchemes explicitScheme = ExplicitSchemes.RungeKutta;
            int explicitOrder = 1;

            var cRef = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);

            
            cRef.NoOfTimesteps = 0;



            var cLoadBal = ShockTubeToro1Template(
                dgDegree: dgDegree,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder);
            cLoadBal.GridPartType = GridPartType.none;
            //cLoadBal.DynamicLoadBalancing_Period = 5;
            //cLoadBal.DynamicLoadBalancing_CellClassifier = new IndifferentCellClassifier();
            //cLoadBal.DynamicLoadBalancing_CellCostEstimatorFactories.Add(CellCostEstimatorLibrary.AllCellsAreEqual);
            //cLoadBal.DynamicLoadBalancing_ImbalanceThreshold = 0.01;


            cLoadBal.NoOfTimesteps = 0;



            CheckRunsProduceSameResults(cRef, cLoadBal);
        }

        //[Test]
        //public static void TestLoadBalancingForRK1WithAV() {
        //    ExplicitSchemes explicitScheme = ExplicitSchemes.RungeKutta;
        //    int explicitOrder = 1;
        //    int numOfClusters = -1;

        //    var cRef = ArtificialViscosityShockTubeTests.SetupToroTest1(
        //        explicitScheme: explicitScheme, explicitOrder: explicitOrder, numOfClusters: numOfClusters);
        //    var solver = new Program();
        //    solver.Init(cRef, null);
        //    solver.RunSolverMode();

        //    var cLoadBal = ArtificialViscosityShockTubeTests.SetupToroTest1(
        //        explicitScheme: explicitScheme, explicitOrder: explicitOrder, numOfClusters: numOfClusters);
        //    cLoadBal.DynamicLoadBalancing_Period = 5;
        //    //cLoadBal.


        //    cLoadBal.AddVariable(Variables.ArtificialViscosity, 1);
        //    int dgDegree = 99;

        //    double sensorLimit = 1e-4;
        //    double epsilon0 = 1.0;
        //    double kappa = 0.5;
        //    Variable sensorVariable = Variables.Density;
        //    cLoadBal.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
        //    cLoadBal.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(cLoadBal.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);
        //}

        private static CNSControl ShockTubeToro1Template(int dgDegree, ExplicitSchemes explicitScheme, int explicitOrder, int noOfCells = 49) {
            double densityLeft = 1.0;
            double velocityLeft = 0.0;
            double pressureLeft = 1.0;
            double densityRight = 0.125;
            double velocityRight = 0.0;
            double pressureRight = 0.1;
            double discontinuityPosition = 0.5;

            CNSControl c = new CNSControl();
            c.DbPath = null;
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
                double[] xNodes = GenericBlas.Linspace(0.0, 1.0, noOfCells + 1);
                var grid = Grid1D.LineGrid(xNodes, periodic: false);

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

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.dtFixed = 1.0e-6;
            c.Endtime = 2e-04;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }
    }
}
