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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace CNS.IBM {

    /// <summary>
    /// Variant of <see cref="RungeKutta"/> that accounts for the influence of
    /// an immersed boundary. That is:
    /// <list type="bullet">
    /// <item>Uses adapted quadrature for cut cells</item>
    /// <item>Evaluates an additional operator to account for BCs at the level set</item>
    /// <item>Uses cell agglomeration</item>
    /// </list>
    /// </summary>
    public abstract class IBMRungeKutta : RungeKutta {

        protected readonly ImmersedSpeciesMap speciesMap;

        protected Lazy<SpatialOperator.Evaluator> boundaryEvaluator;

        protected CellMask cutCells;

        protected CellMask cutAndTargetCells;

        protected readonly CoordinateMapping boundaryParameterMap;

        private SpatialOperator boundaryOperator;

        public IBMRungeKutta(
            SpatialOperator standardOperator,
            SpatialOperator boundaryOperator,
            CoordinateMapping fieldsMap,
            CoordinateMapping parametersMap,
            ImmersedSpeciesMap speciesMap,
            IList<TimeStepConstraint> timeStepConstraints)
            : base(RungeKutta.GetDefaultScheme(speciesMap.Control.ExplicitOrder), standardOperator, fieldsMap, parametersMap, timeStepConstraints, speciesMap.SubGrid) {

            this.speciesMap = speciesMap;
            this.boundaryOperator = boundaryOperator;
            this.boundaryParameterMap = parametersMap;
        }

        protected void UpdateEvaluatorsAndMasks() {
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

            this.m_Evaluator = new Lazy<SpatialOperator.Evaluator>(() =>
                this.Operator.GetEvaluatorEx(
                    Mapping,
                    boundaryParameterMap,
                    Mapping,
                    edgeScheme,
                    volumeScheme,
                    subGridBoundaryTreatment: SpatialOperator.SubGridBoundaryModes.InnerEdgeLTS));

            // Evaluator for boundary conditions at level set zero contour
            CellQuadratureScheme boundaryVolumeScheme = speciesMap.QuadSchemeHelper.GetLevelSetquadScheme(
                0, cutCells, control.LevelSetQuadratureOrder);

            this.boundaryEvaluator = new Lazy<SpatialOperator.Evaluator>(() =>
                boundaryOperator.GetEvaluatorEx(
                    Mapping,
                    boundaryParameterMap,
                    Mapping,
                    null, // Contains no boundary terms
                    boundaryVolumeScheme));
        }

        protected void MoveLevelSetTo(double time) {
            LevelSet levelSet = speciesMap.Tracker.LevelSets[0].As<LevelSet>();
            levelSet.Clear();
            levelSet.ProjectField(X => speciesMap.Control.LevelSetFunction(X, time));
            speciesMap.Tracker.UpdateTracker();

            cutCells = speciesMap.Tracker.Regions.GetCutCellMask();
            cutAndTargetCells = cutCells.Union(speciesMap.Agglomerator.AggInfo.TargetCells);

            // EVIL HACK SINCE UPDATE OF GRADIENT ONLY HAPPENS AFTER TIME-STEP SO FAR
            foreach (var gradientField in boundaryParameterMap.Fields) {
                int d = int.Parse(gradientField.Identification.Last().ToString());

                gradientField.Clear();
                gradientField.Derivative(
                    -1.0,
                    levelSet,
                    d,
                    cutCells);
            }
        }

        protected void AgglomerateAndExtrapolateDGCoordinates() {
            IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
            BlockMsrMatrix nonAgglomeratedMassMatrix = massMatrixFactory.NonAgglomeratedMassMatrix;

            IBMUtility.SubMatrixSpMV(nonAgglomeratedMassMatrix, 1.0, DGCoordinates, 0.0, DGCoordinates, cutCells);  // eq. (39)
            speciesMap.Agglomerator.ManipulateRHS(DGCoordinates, Mapping);  // eq. (39)
            IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, DGCoordinates, 0.0, DGCoordinates, cutAndTargetCells);   // eq. (39)
            speciesMap.Agglomerator.Extrapolate(DGCoordinates.Mapping); // eq. (41)
        }
    }
}
