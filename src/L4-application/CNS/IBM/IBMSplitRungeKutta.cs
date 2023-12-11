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
using ilPSP.Utils;
using System;
using ilPSP.LinSolvers;

namespace CNS.IBM {

    public class IBMSplitRungeKutta : IBMRungeKutta {

        public IBMSplitRungeKutta(DifferentialOperator standardOperator, DifferentialOperator boundaryOperator, CoordinateMapping fieldsMap, CoordinateMapping parametersMap, ImmersedSpeciesMap speciesMap, IList<TimeStepConstraint> timeStepConstraints)
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
            //double dist2 = Environment.CompareTo(CurrentState);
            if (levelSetHasMoved) {
                //double dist22 = Environment.CompareTo(CurrentState);
                UpdateEvaluatorsAndMasks();
                //double dist23 = Environment.CompareTo(CurrentState);

                if (!speciesMap.Agglomerator.AggInfo.Equals(oldAgglomerationInfo)) {
                    // Agglomeration pattern has changed, redo agglomeration
                    AgglomerateAndExtrapolateDGCoordinates();   // eq. (38)-(41)
                    oldAgglomerationInfo = speciesMap.Agglomerator.AggInfo;
                }

                levelSetHasMoved = false;
            }
            //double dist3 = Environment.CompareTo(CurrentState);
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

        //int count = 0;

        /// <summary>
        /// Computes the change rates by sequentially evaluating the standard
        /// and immersed boundary operators. Afterwards, the change rate is
        /// agglomerated and multiplied by the inverse mass matrix (of the
        /// agglomerated basis)
        /// </summary>
        protected override void ComputeChangeRate(double[] k, double AbsTime, double RelTime, double[] edgeFluxes = null) {
            using (new ilPSP.Tracing.FuncTrace()) {
                RaiseOnBeforeComputechangeRate(AbsTime, RelTime);

                //var CV = new CoordinateVector(Evaluator.DomainFields);
                //double dist = Environment.CompareTo(CV);

                //(new CoordinateVector(Evaluator.DomainFields)).SaveToTextFile("inp-rk.txt");
                //(new CoordinateVector(Evaluator.Parameters.ToArray())).SaveToTextFile("para-rk.txt");

                Debug.Assert(Evaluator.DomainFields.Fields.ListEquals(boundaryEvaluator.Value.DomainFields.Fields, (a, b) => object.ReferenceEquals(b, a)));

                //var cv = new CoordinateVector(Evaluator.DomainFields);
                //Random r = new Random(666);
                //for(int ir = 0; ir < cv.Length; ir++) {
                //    cv[ir] = r.NextDouble();
                //}
                

                Evaluator.time = AbsTime + RelTime;
                Evaluator.Evaluate(1.0, 0.0, k);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in standard terms");

                //Console.WriteLine(String.Format("\nBULK: L2-Norm of change rate = {0}", k.L2Norm()));
                //k.SaveToTextFile(String.Format("k_BULK_{0}.txt", count));

                boundaryEvaluator.Value.time = AbsTime + RelTime;
                boundaryEvaluator.Value.Evaluate(1.0, 1.0, k);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in boundary terms");

                //Console.WriteLine(String.Format("BOUNDARY: L2-Norm of change rate = {0}", k.L2Norm()));
                //k.SaveToTextFile(String.Format("k_BOUND_{0}.txt", count));

                // Agglomerate fluxes
                speciesMap.Agglomerator.ManipulateRHS(k, Mapping);

                // Apply inverse to all cells with non-identity mass matrix
                IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
                IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, k, 0.0, k, cutAndTargetCells);

                //massMatrixFactory.MassMatrix.SaveToTextFileSparse("massMatrix.txt");

                //Console.WriteLine($"AGGLOMERATOR: L2-Norm of change rate = {k.L2Norm()}");
                //k.SaveToTextFile($"k_CUT_{count}.txt");
                //k.SaveToTextFile($"c:\\tmp\\cns_k_CUT_{count}.txt");

                //count++;
            }
        }
    }
}
