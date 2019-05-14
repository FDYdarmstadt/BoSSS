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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace CNS.Tests.ConvectiveFlux {

    /// <summary>
    /// First five test cases from Toro 2009, p. 334, using zeroth order. All
    /// tests are simple one-dimensional shock tube tests that are integrated
    /// using explicit Euler
    /// </summary>
    public class ShockTubeTests : TestProgram<CNSControl> {

        /// <summary>
        /// All fluxes except Rusanov flux (which is significantly less
        /// accurate)
        /// </summary>
        private static ConvectiveFluxTypes[] accurateFluxes = new ConvectiveFluxTypes[] {
            ConvectiveFluxTypes.HLL,
            ConvectiveFluxTypes.HLLC,
            ConvectiveFluxTypes.OptimizedHLLC,
            ConvectiveFluxTypes.Godunov
        };

       

        //public static void Main(string[] args) {
        //    SetUp();
        //    Toro1AllButRusanovTest(ConvectiveFluxTypes.Godunov);
        //}

        /// <summary>
        /// Toro 2009, p. 334, case 1
        /// </summary>
        /// <param name="flux"></param>
        /// <returns></returns>
        private static Dictionary<string, object> Toro1Test(ConvectiveFluxTypes flux) {
            CNSControl c = GetShockTubeControlTemplate(
                convectiveFlux: flux,
                densityLeft: 1.0,
                velocityLeft: 0.75,
                pressureLeft: 1.0,
                densityRight: 0.125,
                velocityRight: 0.0,
                pressureRight: 0.1,
                discontinuityPosition: 0.3);

            c.ProjectName = "Shock tube, Toro case 1";
            c.ProjectDescription = "Toro 2009, p. 334, case 1";
            c.Endtime = 0.2;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        /// <summary>
        /// Toro 2009, p. 334, case 1, using Rusanov flux
        /// </summary>
        [Test]
        public static void Toro1RusanovTest() {
            CheckErrorThresholds(
                Toro1Test(ConvectiveFluxTypes.Rusanov),
                Tuple.Create("L2ErrorDensity", 3.6e-2),
                Tuple.Create("L2ErrorVelocity", 8.6e-2),
                Tuple.Create("L2ErrorPressure", 3.4e-2));
        }

        /// <summary>
        /// Toro 2009, p. 334, case 1, using variable flux
        /// </summary>
        [Test]
        [TestCaseSource("accurateFluxes")]
        public static void Toro1AllButRusanovTest(ConvectiveFluxTypes flux) {
            CheckErrorThresholds(
                Toro1Test(flux),
                Tuple.Create("L2ErrorDensity", 2.5e-2),
                Tuple.Create("L2ErrorVelocity", 7.9e-2),
                Tuple.Create("L2ErrorPressure", 2.0e-2));

        }

        /// <summary>
        /// Toro 2009, p. 334, case 2
        /// </summary>
        private static Dictionary<string, object> Toro2Test(ConvectiveFluxTypes flux) {
            CNSControl c = GetShockTubeControlTemplate(
                convectiveFlux: flux,
                densityLeft: 1.0,
                velocityLeft: -2.0,
                pressureLeft: 0.4,
                densityRight: 1.0,
                velocityRight: 2.0,
                pressureRight: 0.4,
                discontinuityPosition: 0.5);

            c.ProjectName = "Shock tube, Toro case 2";
            c.ProjectDescription = "Toro 2009, p. 334, case 2";
            c.Endtime = 0.15;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        /// <summary>
        /// Toro 2009, p. 334, case 2, using Rusanov flux
        /// </summary>
        [Test]
        public static void Toro2RusanovTest() {
            CheckErrorThresholds(
                Toro2Test(ConvectiveFluxTypes.Rusanov),
                Tuple.Create("L2ErrorDensity", 3.2e-2),
                Tuple.Create("L2ErrorVelocity", 7.6e-2),
                Tuple.Create("L2ErrorPressure", 1.8e-2));
        }

        /// <summary>
        /// Toro 2009, p. 334, case 2, using variable flux
        /// </summary>
        [Test]
        [TestCaseSource("accurateFluxes")]
        public static void Toro2AllButRusanovTest(ConvectiveFluxTypes flux) {
            CheckErrorThresholds(
                Toro2Test(flux),
                Tuple.Create("L2ErrorDensity", 3.2e-2),
                Tuple.Create("L2ErrorVelocity", 7.1e-2),
                Tuple.Create("L2ErrorPressure", 1.8e-2));
        }

        /// <summary>
        /// Toro 2009, p. 334, case 3
        /// </summary>
        private static Dictionary<string, object> Toro3Test(ConvectiveFluxTypes flux) {
            CNSControl c = GetShockTubeControlTemplate(
                convectiveFlux: flux,
                densityLeft: 1.0,
                velocityLeft: 0.0,
                pressureLeft: 1000.0,
                densityRight: 1.0,
                velocityRight: 0.0,
                pressureRight: 0.01,
                discontinuityPosition: 0.5);

            c.ProjectName = "Shock tube, Toro case 3";
            c.ProjectDescription = "Toro 2009, p. 334, case 3";
            c.Endtime = 0.012;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        /// <summary>
        /// Toro 2009, p. 334, case 3, using Rusanov flux
        /// </summary>
        [Test]
        public static void Toro3RusanovTest() {
            CheckErrorThresholds(
                Toro3Test(ConvectiveFluxTypes.Rusanov),
                Tuple.Create("L2ErrorDensity", 6.4e-1),
                Tuple.Create("L2ErrorVelocity", 2.0e-0),
                Tuple.Create("L2ErrorPressure", 3.8e+1));

        }

        /// <summary>
        /// Toro 2009, p. 334, case 3, using variable flux
        /// </summary>
        [Test]
        [TestCaseSource("accurateFluxes")]
        public static void Toro3AllButRusanovTest(ConvectiveFluxTypes flux) {
            CheckErrorThresholds(
                Toro3Test(flux),
                Tuple.Create("L2ErrorDensity", 6.0e-1),
                Tuple.Create("L2ErrorVelocity", 2.0e-0),
                Tuple.Create("L2ErrorPressure", 3.7e+1));
        }

        /// <summary>
        /// Toro 2009, p. 334, case 4
        /// </summary>
        private static Dictionary<string, object> Toro4Test(ConvectiveFluxTypes flux) {
            CNSControl c = GetShockTubeControlTemplate(
                convectiveFlux: flux,
                densityLeft: 5.99924,
                velocityLeft: 19.5975,
                pressureLeft: 460.894,
                densityRight: 5.99242,
                velocityRight: -6.19633,
                pressureRight: 46.0950,
                discontinuityPosition: 0.4);

            c.ProjectName = "Shock tube, Toro case 4";
            c.ProjectDescription = "Toro 2009, p. 334, case 4";
            c.Endtime = 0.035;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        /// <summary>
        /// Toro 2009, p. 334, case 4, using Rusanov flux
        /// </summary>
        [Test]
        public static void Toro4RusanovTest() {
            CheckErrorThresholds(
                Toro4Test(ConvectiveFluxTypes.Rusanov),
                Tuple.Create("L2ErrorDensity", 2.3e-0),
                Tuple.Create("L2ErrorVelocity", 1.2e-0),
                Tuple.Create("L2ErrorPressure", 1.1e+2));

        }

        /// <summary>
        /// Toro 2009, p. 334, case 4, using Rusanov flux
        /// </summary>
        [Test]
        [TestCaseSource("accurateFluxes")]
        public static void Toro4AllButRusanovTest(ConvectiveFluxTypes flux) {
            CheckErrorThresholds(
                Toro4Test(flux),
                Tuple.Create("L2ErrorDensity", 2.0e-0),
                Tuple.Create("L2ErrorVelocity", 9.0e-1),
                Tuple.Create("L2ErrorPressure", 7.0e+1));
        }

        /// <summary>
        /// Toro 2009, p. 334, case 5
        /// </summary>
        private static Dictionary<string, object> Toro5Test(ConvectiveFluxTypes flux) {
            CNSControl c = GetShockTubeControlTemplate(
                convectiveFlux: flux,
                densityLeft: 1.0,
                velocityLeft: -19.59745,
                pressureLeft: 1000.0,
                densityRight: 1.0,
                velocityRight: -19.59745,
                pressureRight: 0.01,
                discontinuityPosition: 0.8);

            c.ProjectName = "Shock tube, Toro case 5";
            c.ProjectDescription = "Toro 2009, p. 334, case 5";
            c.Endtime = 0.012;

            var solver = new Program();
            solver.Init(c);
            solver.RunSolverMode();

            return solver.QueryHandler.QueryResults;
        }

        /// <summary>
        /// Toro 2009, p. 334, case 5, using Rusanov flux
        /// </summary>
        [Test]
        public static void Toro5RusanovTest() {
            CheckErrorThresholds(
                Toro5Test(ConvectiveFluxTypes.Rusanov),
                Tuple.Create("L2ErrorDensity", 5.9e-1),
                Tuple.Create("L2ErrorVelocity", 1.9e-0),
                Tuple.Create("L2ErrorPressure", 3.5e+1));
        }

        /// <summary>
        /// Toro 2009, p. 334, case 5, using variable flux
        /// </summary>
        [Test]
        [TestCaseSource("accurateFluxes")]
        public static void Toro5AllButRusanovTest(ConvectiveFluxTypes flux) {
            CheckErrorThresholds(
                Toro5Test(flux),
                Tuple.Create("L2ErrorDensity", 5.1e-1),
                Tuple.Create("L2ErrorVelocity", 1.4e-0),
                Tuple.Create("L2ErrorPressure", 2.9e+1));
        }

        /// <summary>
        /// Template for all shock tube tests
        /// </summary>
        /// <param name="convectiveFlux"></param>
        /// <param name="densityLeft"></param>
        /// <param name="velocityLeft"></param>
        /// <param name="pressureLeft"></param>
        /// <param name="densityRight"></param>
        /// <param name="velocityRight"></param>
        /// <param name="pressureRight"></param>
        /// <param name="discontinuityPosition"></param>
        /// <returns></returns>
        private static CNSControl GetShockTubeControlTemplate(ConvectiveFluxTypes convectiveFlux, double densityLeft, double velocityLeft, double pressureLeft, double densityRight, double velocityRight, double pressureRight, double discontinuityPosition) {
            CNSControl c = new CNSControl();
            c.DbPath = null;
            c.savetodb = false;

            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = convectiveFlux;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            int dgDegree = 0;
            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);

            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Velocity.xComponent, dgDegree);

            c.GridFunc = delegate {
                double[] xNodes = GenericBlas.Linspace(0.0, 1.0, 201);
                var grid = Grid1D.LineGrid(xNodes, false);
                grid.EdgeTagNames.Add(1, "supersonicOutlet");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            // Take inner values everywhere
            c.AddBoundaryValue("supersonicOutlet");

            Material material = c.GetMaterial();
            StateVector stateLeft = StateVector.FromPrimitiveQuantities(
                material, densityLeft, new Vector(velocityLeft), pressureLeft);
            StateVector stateRight = StateVector.FromPrimitiveQuantities(
                material, densityRight, new Vector(velocityRight), pressureRight);

            c.InitialValues_Evaluators.Add(
                Variables.Density,
                X => stateLeft.Density + (stateRight.Density - stateLeft.Density) * (X[0] - discontinuityPosition).Heaviside());
            c.InitialValues_Evaluators.Add(
                CNSVariables.Velocity.xComponent,
                X => stateLeft.Velocity.x + (stateRight.Velocity.x - stateLeft.Velocity.x) * (X[0] - discontinuityPosition).Heaviside());
            c.InitialValues_Evaluators.Add(
                CNSVariables.Pressure,
                X => stateLeft.Pressure + (stateRight.Pressure - stateLeft.Pressure) * (X[0] - discontinuityPosition).Heaviside());

            var riemannSolver = new ExactRiemannSolver(stateLeft, stateRight, new Vector(1.0));
            double pStar, uStar;
            riemannSolver.GetStarRegionValues(out pStar, out uStar);

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(
                Variables.Density,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Density));
            c.Queries.Add("L2ErrorVelocity", QueryLibrary.L2Error(
                CNSVariables.Velocity.xComponent,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Velocity.x));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(
                CNSVariables.Pressure,
                (X, t) => riemannSolver.GetState(pStar, uStar, X[0] - discontinuityPosition, t).Pressure));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.5;
            c.NoOfTimesteps = int.MaxValue;

            c.PrintInterval = 50;

            return c;
        }
    }
}
