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

using System;
using System.Collections;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using BoSSS.Solution.Timestepping;
using System.Collections.Generic;
using ilPSP.LinSolvers;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;

namespace CNS.IBM {

    class IBMAdamsBashforth : AdamsBashforth {

        private ImmersedSpeciesMap speciesMap;

        private DifferentialOperator boundaryOperator;

        private CoordinateMapping boundaryParameterMap;

        private Lazy<IEvaluatorNonLin> boundaryEvaluator;

        private CellMask cutCells;

        private CellMask cutAndTargetCells;

        public IBMAdamsBashforth(
            DifferentialOperator standardOperator,
            DifferentialOperator boundaryOperator,
            CoordinateMapping fieldsMap,
            CoordinateMapping parametersMap,
            ISpeciesMap ibmSpeciesMap,
            IBMControl control,
            IList<TimeStepConstraint> timeStepConstraints)
            : base(standardOperator, fieldsMap, parametersMap, control.ExplicitOrder, timeStepConstraints, ibmSpeciesMap.SubGrid) {  // TO DO: I SIMPLY REMOVED PARAMETERMAP HERE; MAKE THIS MORE PRETTY

            speciesMap = ibmSpeciesMap as ImmersedSpeciesMap;
            if (this.speciesMap == null) {
                throw new ArgumentException(
                    "Only supported for species maps of type 'ImmersedSpeciesMap'",
                    "speciesMap");
            }

            this.boundaryOperator = boundaryOperator;
            this.boundaryParameterMap = parametersMap;
            agglomerationPatternHasChanged = true;

            // StarUp Phase needs also an IBM time stepper
            RungeKuttaScheme = new IBMSplitRungeKutta(
                standardOperator,
                boundaryOperator,
                fieldsMap,
                parametersMap,
                speciesMap,
                timeStepConstraints);
        }

        private void BuildEvaluatorsAndMasks() {
            CellMask fluidCells = speciesMap.SubGrid.VolumeMask;
            cutCells = speciesMap.Tracker.Regions.GetCutCellMask();
            cutAndTargetCells = cutCells.Union(speciesMap.Agglomerator.AggInfo.TargetCells);

            IBMControl control = speciesMap.Control;
            SpeciesId species = speciesMap.Tracker.GetSpeciesId(control.FluidSpeciesName);

            CellQuadratureScheme volumeScheme = speciesMap.QuadSchemeHelper.GetVolumeQuadScheme(
                species, true, fluidCells, control.LevelSetQuadratureOrder);

            // Does _not_ include agglomerated edges
            EdgeMask nonVoidEdges = speciesMap.QuadSchemeHelper.GetEdgeMask(species);
            EdgeQuadratureScheme edgeScheme = speciesMap.QuadSchemeHelper.GetEdgeQuadScheme(
                species, true, nonVoidEdges, control.LevelSetQuadratureOrder);

            this.m_Evaluator = new Lazy<IEvaluatorNonLin>(delegate () {
                this.Operator.EdgeQuadraturSchemeProvider = g => edgeScheme;
                this.Operator.VolumeQuadraturSchemeProvider = g => volumeScheme;
                var opi = this.Operator.GetEvaluatorEx(
                    Mapping,
                    //null,
                    this.boundaryParameterMap,
                    Mapping);
                //opi.ActivateSubgridBoundary(volumeScheme.Domain, subGridBoundaryTreatment: SubGridBoundaryModes.InnerEdgeLTS);
                opi.ActivateSubgridBoundary(fluidCells, subGridBoundaryTreatment: SubGridBoundaryModes.InnerEdgeLTS);
                return opi;
            });

            // Evaluator for boundary conditions at level set zero contour
            CellQuadratureScheme boundaryVolumeScheme = speciesMap.QuadSchemeHelper.GetLevelSetquadScheme(
                0, cutCells, control.LevelSetQuadratureOrder);

            this.boundaryEvaluator = new Lazy<IEvaluatorNonLin>(delegate() {
                boundaryOperator.VolumeQuadraturSchemeProvider = g => boundaryVolumeScheme; 
                boundaryOperator.EdgeQuadraturSchemeProvider = g => null; // Contains no boundary terms;
                return boundaryOperator.GetEvaluatorEx(
                    Mapping,
                    boundaryParameterMap,
                    Mapping);
            });
        }


        private static bool agglomerationPatternHasChanged = true;

        public override double Perform(double dt) {
            if (agglomerationPatternHasChanged) {
                // TO DO: Agglomerate difference between old $cutAndTargetCells and new $cutAndTargetCells only
                BuildEvaluatorsAndMasks();

                // Required whenever agglomeration pattern changes
                SpeciesId speciesId = speciesMap.Tracker.GetSpeciesId(speciesMap.Control.FluidSpeciesName);
                IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
                BlockMsrMatrix nonAgglomeratedMassMatrix = massMatrixFactory.BaseFactory.GetMassMatrix(
                    Mapping,
                    new Dictionary<SpeciesId, IEnumerable<double>>() {
                        { speciesId, Enumerable.Repeat(1.0, Mapping.NoOfVariables) } },
                    inverse: false);

                IBMUtility.SubMatrixSpMV(nonAgglomeratedMassMatrix, 1.0, CurrentState, 0.0, CurrentState, cutCells);
                speciesMap.Agglomerator.ManipulateRHS(CurrentState, Mapping);
                IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, CurrentState, 0.0, CurrentState, cutAndTargetCells);
                speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping);

                //Broadcast to RungeKutta ???
            }

            dt = base.Perform(dt);

            speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping);
            return dt;
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

                Evaluator.time = AbsTime;
                Evaluator.Evaluate(1.0, 0.0, k);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in standard terms");

                boundaryEvaluator.Value.time = AbsTime;
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
