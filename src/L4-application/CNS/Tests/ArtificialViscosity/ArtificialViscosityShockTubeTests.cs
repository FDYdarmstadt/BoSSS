﻿/* =======================================================================
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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
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

        //public static void Main(string[] args) {
        //    SetUp();
        //    ToroTest1_ALTS3_3();
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

            var solver = new CNSProgram();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        [Test]
        public static void ToroTest1_RK1() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.RungeKutta, explicitOrder: 1, numOfClusters: -1),
                //Tuple.Create("L2ErrorDensity", 0.0213295440070947 + 1e-14),   //thresholds before changing quad order of AV
                //Tuple.Create("L2ErrorVelocity", 0.0112419647692771 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221686517166847 + 1e-14));

                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2ErrorDensity", 0.0213295144110587 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.0112424957360476 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221692599463277 + 1e-14));

                Tuple.Create("L2ErrorDensity", 0.0179908422878919 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0648935297526108 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0255974274856578 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS1_3() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 1, numOfClusters: 3),
                //Tuple.Create("L2ErrorDensity", 0.0213233520929519 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.0112404986333599 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221621831639872 + 1e-14));

                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2ErrorDensity", 0.0213233225887192 + 1e-14),     
                //Tuple.Create("L2ErrorVelocity", 0.0112410192725479 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221627888301022 + 1e-14));

                Tuple.Create("L2ErrorDensity", 0.0180310597608648 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0635360861561759 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0254089222630976 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS2_3() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 2, numOfClusters: 3),
                //Tuple.Create("L2ErrorDensity", 0.0213221338536003 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.0112395603913243 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221608822934621 + 1e-14));

                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2ErrorDensity", 0.0213221043663509 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.0112400763052089 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221614877490578 + 1e-14));

                Tuple.Create("L2ErrorDensity", 0.018111209757686 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0654910224264829 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.025763138327774 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS3_3() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 3, numOfClusters: 3),
                //Tuple.Create("L2ErrorDensity", 0.0213233152388806 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.01123851527621 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221620445902869 + 1e-14));

                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2ErrorDensity", 0.0213232857482851 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.0112390275141092 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221626504269355 + 1e-14));

                Tuple.Create("L2ErrorDensity", 0.0180871678450302 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0653868400982687 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0257302177039963 + 1e-14));
        }
        [Test]
        public static void ToroTest1_ALTS3_4() {
            CheckErrorThresholds(
                SetupToroTest1(explicitScheme: ExplicitSchemes.LTS, explicitOrder: 3, numOfClusters: 4),
                //Tuple.Create("L2ErrorDensity", 0.0213033175670039 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.0112159195300983 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221417911406063 + 1e-14));

                // Before changing Refactoring (changing LTS cell metric to harmonic sum, etc.)
                //Tuple.Create("L2ErrorDensity", 0.0213032882872332 + 1e-14),
                //Tuple.Create("L2ErrorVelocity", 0.011216449840882 + 1e-14),
                //Tuple.Create("L2ErrorPressure", 0.0221423943810654 + 1e-14));

                Tuple.Create("L2ErrorDensity", 0.0168812438223513 + 1e-14),
                Tuple.Create("L2ErrorVelocity", 0.0622376806873876 + 1e-14),
                Tuple.Create("L2ErrorPressure", 0.0247342489766164 + 1e-14));
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
            Variable sensorVariable = CompressibleVariables.Density;
            c.CNSShockSensor = new PerssonSensor(sensorVariable, sensorLimit);
            c.ArtificialViscosityLaw = new SmoothedHeavisideArtificialViscosityLaw(c.CNSShockSensor, dgDegree, sensorLimit, epsilon0, kappa);

            c.AddVariable(CompressibleVariables.Density, dgDegree);
            c.AddVariable(CompressibleVariables.Momentum.xComponent, dgDegree);
            c.AddVariable(CompressibleVariables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);
            c.AddVariable(CNSVariables.Viscosity, dgDegree);
            c.AddVariable(CNSVariables.ArtificialViscosity, 1);
            c.AddVariable(CNSVariables.ShockSensor, 0);
            c.AddVariable(CNSVariables.LocalMachNumber, dgDegree);
            if (c.ExplicitScheme.Equals(ExplicitSchemes.LTS)) {
                c.AddVariable(CNSVariables.LTSClusters, 0);
            }

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(0.0, 1.0, noOfCellsPerDirection + 1);
                var grid = Grid1D.LineGrid(xNodes, periodic: false);

                // Boundary conditions
                grid.EdgeTagNames.Add(1, "AdiabaticSlipWall");
                grid.DefineEdgeTags(delegate (Vector _X) {
                    return (byte)1;
                });
                return grid;
            };

            c.AddBoundaryValue("AdiabaticSlipWall");

            Material material = c.GetMaterial();
            StateVector stateLeft = StateVector.FromPrimitiveQuantities(
                material, densityLeft, new Vector(velocityLeft, 0.0, 0.0), pressureLeft);
            StateVector stateRight = StateVector.FromPrimitiveQuantities(
                material, densityRight, new Vector(velocityRight, 0.0, 0.0), pressureRight);

            c.InitialValues_Evaluators.Add(
                    CompressibleVariables.Density,
                    X => stateLeft.Density + (stateRight.Density - stateLeft.Density) * (X[0] - discontinuityPosition).Heaviside());
            c.InitialValues_Evaluators.Add(
                CNSVariables.Velocity.xComponent,
                X => stateLeft.Velocity.x + (stateRight.Velocity.x - stateLeft.Velocity.x) * (X[0] - discontinuityPosition).Heaviside());
            c.InitialValues_Evaluators.Add(
                CNSVariables.Pressure,
                X => stateLeft.Pressure + (stateRight.Pressure - stateLeft.Pressure) * (X[0] - discontinuityPosition).Heaviside());

            var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0, 0.0, 0.0));
            double pStar, uStar;
            riemannSolver.GetStarRegionValues(out pStar, out uStar);

            double DensityAnalytical(double[] X, double t) {
                var v = riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density;
                return v;
            }
            double VelocityAnalytical(double[] X, double t) {
                var v = riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x;
                return v;
            }
            double PressureAnalytical(double[] X, double t) {
                var v = riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure;
                return v;
            }

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(CompressibleVariables.Density, DensityAnalytical));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(CNSVariables.Velocity.xComponent, VelocityAnalytical));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(CNSVariables.Pressure, PressureAnalytical));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = 1e-02;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }
    }
}