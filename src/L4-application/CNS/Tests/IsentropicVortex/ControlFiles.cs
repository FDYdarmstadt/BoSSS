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
using BoSSS.Solution.Queries;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using CNS.Solution;
using ilPSP.Utils;
using System;

namespace CNS.Tests.IsentropicVortex {

    /// <summary>
    /// Tests based on the well-known isentropic vortex exact solution
    /// </summary>
    public static class ControlFiles {
        
        /// <summary>
        /// Common settings for all tests within this set of test cases
        /// </summary>
        /// <param name="divisions"></param>
        /// <param name="dgDegree"></param>
        /// <returns></returns>
        private static VortexControl Template(int divisions, int dgDegree) {
            VortexControl c = new VortexControl();
            c.savetodb = false;

            c.ActiveOperators = Operators.Convection;
            c.TimeSteppingScheme = TimeSteppingSchemes.Explicit;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;

            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(Variables.Pressure, dgDegree);
            c.AddVariable(Variables.Entropy, dgDegree);

            c.GridFunc = delegate {
                int noOfCellsPerDirection = (2 << divisions) * 10;
                GridCommons grid = Grid2D.Cartesian2DGrid(
                    GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1),
                    GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1),
                    periodicX: true,
                    periodicY: true);
                return grid;
            };

            c.VortexSpeed = 1.0;

            c.dtMin = 0.01;
            c.dtMax = 0.01;
            c.CFLFraction = double.NaN;
            c.Endtime = 0.2;
            c.NoOfTimesteps = int.MaxValue;

            return c;
        }

        private static VortexControl IdealGasTemplate(int divisions, int dgDegree) {
            VortexControl c = Template(divisions, dgDegree);

            c.EquationOfState = IdealGas.Air;

            IsentropicVortexExactSolution solution = new IsentropicVortexExactSolution(c, c.VortexSpeed);

            c.InitialValues_Evaluators.Add(Variables.Density, X => solution.rho()(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => solution.u()(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => solution.v()(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => solution.p()(X, 0.0));

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, solution.rho()));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(Variables.Pressure, solution.p()));
            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(Variables.Entropy, (X, t) => 1.0));

            return c;
        }

        /// <summary>
        /// Test case using <see cref="IdealGas"/> in combination with the
        /// <see cref="RusanovFlux"/>.
        /// </summary>
        /// <returns></returns>
        public static VortexControl IsentropicVortexIdealGasRusanov() {
            VortexControl c = IdealGasTemplate(1, 2);
            c.ConvectiveFluxType = ConvectiveFluxTypes.Rusanov;
            return c;
        }

        /// <summary>
        /// Test case using <see cref="IdealGas"/> in combination with the
        /// <see cref="HLLCFlux"/>.
        /// </summary>
        /// <returns></returns>
        public static VortexControl IsentropicVortexIdealGasHLLC() {
            VortexControl c = IdealGasTemplate(1, 2);
            c.ConvectiveFluxType = ConvectiveFluxTypes.HLLC;
            return c;
        }

        /// <summary>
        /// Test case using <see cref="IdealGas"/> in combination with the
        /// <see cref="OptimizedHLLCFlux"/>.
        /// </summary>
        /// <returns></returns>
        public static VortexControl IsentropicVortexIdealGasOptimizedHLLC() {
            VortexControl c = IdealGasTemplate(1, 2);
            c.ConvectiveFluxType = ConvectiveFluxTypes.OptimizedHLLC;
            return c;
        }

        /// <summary>
        /// Test case using <see cref="StiffenedGas"/> in combination with the
        /// <see cref="RusanovFlux"/>.
        /// </summary>
        /// <returns></returns>
        public static VortexControl IsentropicVortexStiffenedGasRusanov() {
            VortexControl c = Template(1, 2);

            double referencePressure = 10.0;
            c.EquationOfState = new StiffenedGas(1.4, referencePressure);
            c.ConvectiveFluxType = ConvectiveFluxTypes.Rusanov;

            double gamma = c.EquationOfState.HeatCapacityRatio;
            Func<double[], double, double> x = (X, t) => X[0] - c.VortexSpeed * t;
            Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
            Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
            Func<double[], double, double> rho = (X, t) => Math.Pow(
                2.3535468936502524 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)),
                1.0 / (gamma - 1.0));
            Func<double[], double, double> uAbs = (X, t) => r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
            Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma) - referencePressure;

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => c.VortexSpeed - Math.Sin(phi(X, 0.0)) * uAbs(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => Math.Cos(phi(X, 0.0)) * uAbs(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, rho));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(Variables.Pressure, p));
            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(Variables.Entropy, (X, t) => 1.0));

            return c;
        }

        /// <summary>
        /// Test case using <see cref="CovolumeGas"/> in combination with the
        /// <see cref="RusanovFlux"/>.
        /// </summary>
        /// <returns></returns>
        public static VortexControl IsentropicVortexCovolumeGasRusanov() {
            VortexControl c = Template(1, 2);

            CovolumeGas eos = new CovolumeGas(1.4, 0.1);
            double integrationConstant = 3.6;

            c.EquationOfState = eos;
            c.ConvectiveFluxType = ConvectiveFluxTypes.Rusanov;

            Func<double[], double, StateVector> solution = (X, t) =>
                CovolumeVortexExactSolution.GetSolution(X[0], X[1], c.VortexSpeed, t, c, integrationConstant);

            c.InitialValues_Evaluators.Add(Variables.Density, X => solution(X, 0.0).Density);
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => solution(X, 0.0).Velocity[0]);
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => solution(X, 0.0).Velocity[1]);
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => solution(X, 0.0).Pressure);

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, (X, t) => solution(X, t).Density, 10));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(Variables.Pressure, (X, t) => solution(X, t).Pressure, 10));
            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(Variables.Entropy, (X, t) => 1.0, 10));

            return c;
        }
    }
}
