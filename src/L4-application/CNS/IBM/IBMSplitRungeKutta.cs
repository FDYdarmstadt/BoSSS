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
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using CNS.EquationSystem;
using System.Diagnostics;
using System.Linq;
using ilPSP;
using System.Collections.Generic;

namespace CNS.IBM {

    public class IBMSplitRungeKutta : IBMRungeKutta {

        public IBMSplitRungeKutta(SpatialOperator standardOperator, SpatialOperator boundaryOperator, CoordinateMapping fieldsMap, CoordinateMapping parametersMap, ImmersedSpeciesMap speciesMap, IList<TimeStepConstraint> timeStepConstraints)
            : base(standardOperator, boundaryOperator, fieldsMap, parametersMap, speciesMap, timeStepConstraints) {

            if (speciesMap.Control.TimesteppingStrategy == TimesteppingStrategies.StrangSplitting) {
                throw new System.NotImplementedException();
            }
        }

        /// <summary>
        /// Required by <see cref="Perform(double)"/>
        /// </summary>
        private bool levelSetHasMoved = true;

        private CellAgglomerator.AgglomerationInfo oldAgglomerationInfo;

        /// <summary>
        /// In the first step, performs the initial agglomeration of the
        /// solution. Afterwards, the solution stays continuous for all times.
        /// After each step, the solution is extrapolated such that standard
        /// operator evaluation can be used. 
        /// </summary>
        /// <param name="dt"></param>
        public override double Perform(double dt) {
            if (speciesMap.Control.DomainType == DomainTypes.MovingImmersedBoundary) {
                // Splitting approach: Everything happens with NEW mass matrix
                MoveLevelSetTo(Time + dt);

                Mapping.Fields.ForEach(f => f.Clear(speciesMap.SubGrid.VolumeMask.Complement()));
                levelSetHasMoved = true;
            }

            if (levelSetHasMoved) {
                UpdateEvaluatorsAndMasks();

                if (!speciesMap.Agglomerator.AggInfo.Equals(oldAgglomerationInfo)) {
                    // Agglomeration pattern has changed, redo agglomeration
                    AgglomerateAndExtrapolateDGCoordinates();   // eq. (38)-(41)
                    oldAgglomerationInfo = speciesMap.Agglomerator.AggInfo;
                }

                levelSetHasMoved = false;
            }
            
            dt = base.Perform(dt);  // eq. (42)
            
            speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping); // eq. (43)

            return dt;
        }

        /// <summary>
        /// Performs a single stage of Runge-Kutta scheme (see
        /// <see cref="RungeKutta.PerformStage(double[], int, double[][], double)"/>)
        /// and extrapolates the result into agglomeration "source" cells
        /// afterwards. This ensures that the standard infrastructure for the flux
        /// evaluation can be reused.
        /// </summary>
        /// <param name="y0"></param>
        /// <param name="s"></param>
        /// <param name="k"></param>
        /// <param name="dt"></param>
        protected override void PerformStage(double[] y0, int s, double[][] k, double dt) {
            base.PerformStage(y0, s, k, dt);
            speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping);
        }

        /// <summary>
        /// Computes the change rates by sequentially evaluating the standard
        /// and immersed boundary operators. Afterwards, the change rate is
        /// agglomerated and multiplied by the inverse mass matrix (of the
        /// agglomerated basis)
        /// </summary>
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

                // Agglomerate fluxes
                speciesMap.Agglomerator.ManipulateRHS(k, Mapping);

                // Apply inverse to all cells with non-identity mass matrix
                IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
                IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, k, 0.0, k, cutAndTargetCells);
            }
        }
    }
}
