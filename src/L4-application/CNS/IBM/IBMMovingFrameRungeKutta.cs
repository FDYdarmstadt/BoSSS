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
using BoSSS.Solution;
using CNS.EquationSystem;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CNS.IBM {

    public class IBMMovingFrameRungeKutta : IBMRungeKutta {

        public IBMMovingFrameRungeKutta(SpatialOperator standardOperator, SpatialOperator boundaryOperator, CoordinateMapping fieldsMap, CoordinateMapping parametersMap, ImmersedSpeciesMap speciesMap, IList<TimeStepConstraint> timeStepConstraints)
            : base(standardOperator, boundaryOperator, fieldsMap, parametersMap, speciesMap, timeStepConstraints) {

            if (speciesMap.Control.DomainType != DomainTypes.MovingImmersedBoundary) {
                throw new Exception();
            }
        }


        private void RKstage(double dt, double[][] k, int s, BlockMsrMatrix MsInv, BlockMsrMatrix M0, double[] u0, double[] coefficients) {
            // Copy coordinates to temp array since SpMV (below) does not
            // support in-place computation
            double[] tempCoordinates = CurrentState.ToArray();
            M0.SpMV(1.0, u0, 0.0, tempCoordinates); // Non-agglomerated
            for (int l = 0; l < s; l++) {
                tempCoordinates.AccV(-coefficients[l] * dt, k[l]); // Non-agglomerated
            }

            speciesMap.Agglomerator.ManipulateRHS(tempCoordinates, Mapping);
            MsInv.SpMV(1.0, tempCoordinates, 0.0, CurrentState);
            speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping);
        }

        public override double Perform(double dt) {
            if (TimeStepConstraints != null) {
                dt = CalculateTimeStep();
            }

            int NoOfStages = this.Scheme.Stages;
            double[][] k = new double[NoOfStages][];
            
            Mapping.Fields.ForEach(f => f.Clear(speciesMap.SubGrid.VolumeMask.Complement()));
            UpdateEvaluatorsAndMasks();
            AgglomerateAndExtrapolateDGCoordinates();
            
            BlockMsrMatrix M0 = speciesMap.GetMassMatrixFactory(Mapping).NonAgglomeratedMassMatrix;
            double[] u0 = CurrentState.ToArray(); // Lives on non-agglomerated mesh

            // Initialize RK scheme
            k[0] = new double[Mapping.LocalLength];
            ComputeChangeRate(k[0], Time, 0.0);

            // Intermediate stages
            for (int stage = 1; stage < NoOfStages; stage++) {
                MoveLevelSetTo(Time + dt * this.Scheme.c[stage]);
                UpdateEvaluatorsAndMasks();
                AgglomerateAndExtrapolateDGCoordinates();

                BlockMsrMatrix MsInverse = speciesMap.GetMassMatrixFactory(Mapping).InverseMassMatrix;
                RKstage(dt, k, stage, MsInverse, M0, u0, this.Scheme.a.GetRow(stage));

                k[stage] = new double[Mapping.LocalLength];
                ComputeChangeRate(k[stage], Time, this.Scheme.c[stage] * dt);
            }

            // Final stage
            MoveLevelSetTo(Time + dt);
            UpdateEvaluatorsAndMasks();

            BlockMsrMatrix M1Inverse = speciesMap.GetMassMatrixFactory(Mapping).InverseMassMatrix;
            RKstage(dt, k, NoOfStages, M1Inverse, M0, u0, this.Scheme.b);

            m_Time = Time + dt;
            return dt;
        }

        protected override void ComputeChangeRate(double[] k, double AbsTime, double RelTime, double[] edgeFluxes = null) {
            using (new ilPSP.Tracing.FuncTrace()) {
                RaiseOnBeforeComputechangeRate(AbsTime, RelTime);

                Evaluator.time = AbsTime + RelTime;
                Evaluator.Evaluate(1.0, 0.0, k);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in standard terms");

                boundaryEvaluator.Value.time = AbsTime + RelTime;
                boundaryEvaluator.Value.Evaluate(1.0, 1.0, k);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in boundary terms");
            }
        }

        protected override void PerformStage(double[] y0, int s, double[][] k, double dt) {
            throw new System.Exception("Should never be called");
        }

    }
}
