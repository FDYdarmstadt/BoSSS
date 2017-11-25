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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using CNS.ShockCapturing;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace CNS.Tests.ArtificialViscosity {

    /// <summary>
    /// Tests for the shock tube problem using artificial viscosity. The results are
    /// compared to the results of the ExactRiemannSolver implementation in BoSSS. All tests
    /// are one-dimensional and use an explicit Runge-Kutta time stepping scheme. The test
    /// cases are chosen from Toro 2009, p.129, table 4.1.
    /// </summary>
    public class ArtificialViscosityShockTubeTests : TestProgram<CNSControl> {

        private static CommandLineOptions commandLineOptions = null;
        //private static CommandLineOptions commandLineOptions = new CommandLineOptions() {
        //    delPlt = true,
        //    ImmediatePlotPeriod = 1
        //};

        //public static void Main(string[] args) {
        //    SetUp();
        //    ToroTest1();
        //}

        /// <summary>
        /// Test 1 from Toro 2009, p. 119, table 4.1 using the OptimizedHLLC convective flux.
        /// The results are compared to the analytical solution that is calculated by
        /// the ExactRiemannSolver implementation in BoSSS.
        /// </summary>
        private static Dictionary<string, object> SetupToroTest1(ExplicitSchemes explicitScheme, int explicitOrder, int numOfClusters) {
            CNSControl c = GetArtificialViscosityShockTubeControlTemplate(
                convectiveFlux: ConvectiveFluxTypes.OptimizedHLLC,
                densityLeft: 1.0,
                velocityLeft: 0.0,
                pressureLeft: 1.0,
                densityRight: 0.125,
                velocityRight: 0.0,
                pressureRight: 0.1,
                discontinuityPosition: 0.5,
                explicitScheme: explicitScheme,
                explicitOrder: explicitOrder,
                numOfClusters: numOfClusters);

            c.ProjectName = "Artificial viscosity, shock tube, Toro test 1";
            c.ProjectDescription = "Toro 2009, p. 129, table 4.1, test 1";

            var solver = new Program();
            solver.Init(c, commandLineOptions);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        [Test]
        public static void ToroTest1_RK1() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.RungeKutta, explicitOrder: 1, numOfClusters: -1),
                //Tuple.Create("L2ErrorDensity", 2.133e-2),
                //Tuple.Create("L2ErrorVelocity", 1.125e-2),
                //Tuple.Create("L2ErrorPressure", 2.217e-2));
                Tuple.Create("L2ErrorDensity", 0.0213295440070947 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0112419647692771 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0221686517166847 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS1_3() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 1, numOfClusters: 3),
                Tuple.Create("L2ErrorDensity", 0.0213233520929519 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0112404986333599 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0221621831639872 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS2_3() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 2, numOfClusters: 3),
                Tuple.Create("L2ErrorDensity", 0.0213221338536003 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0112395603913243 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0221608822934621 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS3_3() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 3, numOfClusters: 3),
                Tuple.Create("L2ErrorDensity", 0.0213233152388806 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.01123851527621 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0221620445902869 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS3_4() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 3, numOfClusters: 4),
                Tuple.Create("L2ErrorDensity", 0.0213033175670039 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0112159195300983 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0221417911406063 + 1e-14));
        }

        /// <summary>
        /// Template of the control file for all shock tube tests using artificial viscosity
        /// </summary>
        /// <param name="convectiveFlux"></param>
        /// <param name="densityLeft"></param>
        /// <param name="velocityLeft"></param>
        /// <param name="pressureLeft"></param>
        /// <param name="densityRight"></param>
        /// <param name="velocityRight"></param>
        /// <param name="pressureRight"></param>
        /// <param name="discontinuityPosition"></param>
        /// <param name="explicitScheme"></param>
        /// <param name="explicitOrder"></param>
        /// <param name="numOfClusters"></param>
        /// <returns></returns>
        private static CNSControl GetArtificialViscosityShockTubeControlTemplate(ConvectiveFluxTypes convectiveFlux, double densityLeft, double velocityLeft, double pressureLeft, double densityRight, double velocityRight, double pressureRight, double discontinuityPosition, ExplicitSchemes explicitScheme, int explicitOrder, int numOfClusters) {
            CNSControl c = new CNSControl();

            c.DbPath = null;
            c.savetodb = false;

            //c.DbPath = @"c:\bosss_db\";
            //c.savetodb = true;
            //c.saveperiod = 100;
            //c.PrintInterval = 100;

            c.ActiveOperators = Operators.Convection | Operators.ArtificialViscosity;
            c.ConvectiveFluxType = convectiveFlux;

            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);
            c.ReynoldsNumber = 1.0;
            c.PrandtlNumber = 0.71;

            // Time stepping scheme
            c.ExplicitScheme = explicitScheme;
            c.ExplicitOrder = explicitOrder;
            if (explicitScheme == ExplicitSchemes.LTS) {
                c.NumberOfSubGrids = numOfClusters;
                c.ReclusteringInterval = 1;
                c.FluxCorrection = false;
            }

            int dgDegree = 3;
            int noOfCellsPerDirection = 50;

            double sensorLimit = 1e-4;
            double epsilon0 = 1.0;
            double kappa = 0.5;
            Variable sensorVariable = Variables.Density;
            c.ShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.ShockSensor, dgDegree, sensorLimit, epsilon0, kappa);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Velocity.xComponent, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Entropy, dgDegree);
            c.AddVariable(Variables.Viscosity, dgDegree);
            c.AddVariable(Variables.ArtificialViscosity, 1);
            //c.AddVariable(Variables.Sensor, dgDegree);
            c.AddVariable(Variables.LocalMachNumber, dgDegree);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(Variables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);
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
            double pStar, uStar;
            riemannSolver.GetStarRegionValues(out pStar, out uStar);

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
            c.Endtime = 5e-04;
            //c.NoOfTimesteps = 500;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }
    }
}