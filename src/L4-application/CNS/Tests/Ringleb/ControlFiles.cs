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

using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Queries;
using CNS.EquationSystem;
using System;
using static BoSSS.Solution.CompressibleFlowCommon.Boundary.ExactRinglebBoundaryState;

namespace CNS.Tests.Ringleb {

    /// <summary>
    /// Test cases based on Ringleb's exact solution for the two-dimensional
    /// Euler equations
    /// </summary>
    public static class ControlFiles {

        /// <summary>
        /// Common settings for all tests within this set of test cases
        /// </summary>
        /// <param name="dgDegree"></param>
        /// <param name="eos"></param>
        /// <returns></returns>
        private static RinglebControl GetTemplate(int dgDegree, StiffenedGas eos) {
            RinglebControl c = new RinglebControl();
            c.DbPath = "../../Tests/Ringleb/ringlebTests.zip";
            c.savetodb = false;

            c.ActiveOperators = Operators.Convection;
            c.EquationOfState = eos;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 4;

            c.MachNumber = 1 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(CNSVariables.Pressure, dgDegree);
            c.AddVariable(CNSVariables.Entropy, dgDegree);

            Func<double[], RinglebExactSolution.FlowState> solution = X => RinglebExactSolution.GetFlowState(
                X[0],
                X[1],
                eos.HeatCapacityRatio,
                eos.ReferencePressure,
                c.RinglebReferenceSpeedOfSound,
                c.RinglebReferenceTotalPressure);

            c.InitialValues_Evaluators.Add(Variables.Density, X => solution(X).Density);
            c.InitialValues_Evaluators.Add(Variables.Momentum.xComponent, X => solution(X).Momentum[0]);
            c.InitialValues_Evaluators.Add(Variables.Momentum.yComponent, X => solution(X).Momentum[1]);
            c.InitialValues_Evaluators.Add(Variables.Energy, X => solution(X).Energy);

            c.AddBoundaryValue("ringleb");

            c.Queries.Add("L2ErrorDensity", QueryLibrary.L2Error(Variables.Density, (X, t) => solution(X).Density));
            c.Queries.Add("L2ErrorPressure", QueryLibrary.L2Error(CNSVariables.Pressure, (X, t) => solution(X).Pressure));

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.5;
            c.Endtime = double.MaxValue;
            c.NoOfTimesteps = 50;

            return c;
        }

        /// <summary>
        /// Test cases using an <see cref="IdealGas"/>
        /// </summary>
        /// <returns></returns>
        public static RinglebControl RinglebIdealGasTest() {
            RinglebControl c = GetTemplate(2, new StiffenedGas(1.4, 0.0));

            c.GridGuid = new Guid("499a52ea-9a36-48c4-9c3a-a13d1414b936");
            c.ConvectiveFluxType = Convection.ConvectiveFluxTypes.Rusanov;

            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(CNSVariables.Entropy, (X, t) => 0.7142857142857142));

            return c;
        }

        /// <summary>
        /// Test cases using a <see cref="StiffenedGas"/>
        /// </summary>
        /// <returns></returns>
        public static RinglebControl RinglebStiffenedGasTest() {
            RinglebControl c = GetTemplate(2, new StiffenedGas(7.0, 10.0));

            c.GridGuid = new Guid("f5a0fea5-156d-42ce-abf2-a418469140bb");
            c.ConvectiveFluxType = Convection.ConvectiveFluxTypes.Rusanov;
            c.RinglebReferenceSpeedOfSound = 5.0;
            c.RinglebReferenceTotalPressure = 30.0;

            c.Queries.Add("L2ErrorEntropy", QueryLibrary.L2Error(CNSVariables.Entropy, (X, t) => 1.80939686135e-6));

            return c;
        }
    }
}