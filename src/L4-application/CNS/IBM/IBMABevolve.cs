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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Timestepping;
using ilPSP.Utils;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Solution.CompressibleFlowCommon;

namespace CNS.IBM {

    class IBMABevolve : ABevolve {

        private ImmersedSpeciesMap speciesMap;

        private DifferentialOperator boundaryOperator;

        private CoordinateMapping boundaryParameterMap;

        private Lazy<IEvaluatorNonLin> boundaryEvaluator;

        private CellMask cutCells;

        private CellMask cutAndTargetCells;

        public IBMABevolve(
            DifferentialOperator standardOperator,
            DifferentialOperator boundaryOperator,
            CoordinateMapping fieldsMap,
            CoordinateMapping parametersMap,
            ISpeciesMap ibmSpeciesMap,
            int explicitOrder,
            int levelSetQuadratureOrder,
            XQuadFactoryHelper.MomentFittingVariants momentFittingVariant,
            SubGrid sgrd,
            bool adaptive = false)
            : base(standardOperator, fieldsMap, parametersMap, explicitOrder, adaptive: adaptive, sgrd: sgrd) {

            speciesMap = ibmSpeciesMap as ImmersedSpeciesMap;
            if (speciesMap == null) {
                throw new ArgumentException(
                    "Only supported for species maps of type 'ImmersedSpeciesMap'",
                    "speciesMap");
            }

            this.boundaryOperator = boundaryOperator;
            this.boundaryParameterMap = parametersMap;
            agglomerationPatternHasChanged = true;
        }

        internal void BuildEvaluatorsAndMasks() {

            CellMask fluidCells = speciesMap.SubGrid.VolumeMask.Intersect(ABSubGrid.VolumeMask);
            cutCells = speciesMap.Tracker.Regions.GetCutCellMask().Intersect(ABSubGrid.VolumeMask);
            cutAndTargetCells = cutCells.Union(speciesMap.Agglomerator.AggInfo.TargetCells).Intersect(ABSubGrid.VolumeMask);


            IBMControl control = speciesMap.Control;
            SpeciesId species = speciesMap.Tracker.GetSpeciesId(control.FluidSpeciesName);

            CellQuadratureScheme volumeScheme = speciesMap.QuadSchemeHelper.GetVolumeQuadScheme(
                species, true, fluidCells, control.LevelSetQuadratureOrder);


            // Does _not_ include agglomerated edges
            EdgeMask nonVoidEdges = speciesMap.QuadSchemeHelper.GetEdgeMask(species);
            nonVoidEdges = nonVoidEdges.Intersect(ABSubGrid.AllEdgesMask.ToGeometicalMask());
            EdgeQuadratureScheme edgeScheme = speciesMap.QuadSchemeHelper.GetEdgeQuadScheme(
                species, true, nonVoidEdges, control.LevelSetQuadratureOrder);

            this.m_Evaluator = new Lazy<IEvaluatorNonLin>(delegate () {
                this.Operator.EdgeQuadraturSchemeProvider = g => edgeScheme;
                this.Operator.VolumeQuadraturSchemeProvider = g => volumeScheme;

                var opi = this.Operator.GetEvaluatorEx(
                    Mapping,
                    boundaryParameterMap,
                    Mapping);
                opi.ActivateSubgridBoundary(ABSubGrid.VolumeMask, subGridBoundaryTreatment: SubGridBoundaryModes.InnerEdgeLTS);
                return opi;
            });

            // Evaluator for boundary conditions at level set zero contour
            CellQuadratureScheme boundaryVolumeScheme = speciesMap.QuadSchemeHelper.GetLevelSetquadScheme(
                0, cutCells, control.LevelSetQuadratureOrder);

            this.boundaryEvaluator = new Lazy<IEvaluatorNonLin>(delegate() {
                boundaryOperator.EdgeQuadraturSchemeProvider = g => null; // Contains no boundary terms --> PROBLEM??????????
                boundaryOperator.VolumeQuadraturSchemeProvider = g => boundaryVolumeScheme;
                return boundaryOperator.GetEvaluatorEx(
                    Mapping,
                    boundaryParameterMap,
                    Mapping);
            });
        }


        protected override void ComputeChangeRate(double[] k, double AbsTime, double RelTime, double[] edgeFluxes = null) {
            using (new ilPSP.Tracing.FuncTrace()) {
                RaiseOnBeforeComputechangeRate(AbsTime, RelTime);

                //(new CoordinateVector(Evaluator.DomainFields)).SaveToTextFile("inp-lts.txt");
                //(new CoordinateVector(Evaluator.Parameters.ToArray())).SaveToTextFile("para-lts.txt");

                Evaluator.time = AbsTime;
                Evaluator.Evaluate(1.0, 0.0, k, outputBndEdge: edgeFluxes);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in standard terms");

                //k.SaveToTextFile("k-lts-bulk.txt");

                boundaryEvaluator.Value.time = AbsTime;
                boundaryEvaluator.Value.Evaluate(1.0, 1.0, k);
                Debug.Assert(
                    !k.Any(f => double.IsNaN(f)),
                    "Unphysical flux in boundary terms");

                //k.SaveToTextFile("k-lts-bndy.txt");

                // Agglomerate fluxes
                speciesMap.Agglomerator.ManipulateRHS(k, Mapping);

                // Apply inverse to all cells with non-identity mass matrix
                IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
                IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, k, 0.0, k, cutAndTargetCells);
            }
        }

        internal bool agglomerationPatternHasChanged = true;

        public override double Perform(double dt) {
            if (agglomerationPatternHasChanged) {
                BuildEvaluatorsAndMasks();

                // Required whenever agglomeration pattern changes
                //SpeciesId speciesId = speciesMap.Tracker.GetSpeciesId(speciesMap.Control.FluidSpeciesName);
                //IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
                //BlockDiagonalMatrix nonAgglomeratedMassMatrix = massMatrixFactory.BaseFactory.GetMassMatrix(
                //    Mapping,
                //    new Dictionary<SpeciesId, IEnumerable<double>>() {
                //        { speciesId, Enumerable.Repeat(1.0, Mapping.NoOfVariables) } },
                //    inverse: false,
                //    VariableAgglomerationSwitch: new bool[Mapping.Fields.Count]);

                //IBMUtility.SubMatrixSpMV(nonAgglomeratedMassMatrix, 1.0, DGCoordinates, 0.0, DGCoordinates, cutCells);
                //speciesMap.Agglomerator.ManipulateRHS(DGCoordinates, Mapping);
                //IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, DGCoordinates, 0.0, DGCoordinates, cutAndTargetCells);
                //speciesMap.Agglomerator.Extrapolate(DGCoordinates.Mapping);

                agglomerationPatternHasChanged = false;
            }

            dt = base.Perform(dt);
            return dt;
        }

        protected override double[] ComputesUpdatedDGCoordinates(double[] completeChangeRate) {
            using (new ilPSP.Tracing.FuncTrace("ComputesUpdatedDGCoordinates")) {
                double[] y0 = new double[Mapping.LocalLength];
                CurrentState.CopyTo(y0, 0);
                CurrentState.axpy<double[]>(completeChangeRate, -1);

                // Speciality for IBM: Do the extrapolation
                speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping);

                double[] upDGC = CurrentState.ToArray();
                upDGC = OrderValuesBySgrd(upDGC);
                // DGCoordinates should be untouched after calling this method
                CurrentState.Clear();
                CurrentState.CopyFrom(y0, 0);

                return upDGC;
            }
        }

        /// <summary>
        /// Takes a double[] with results for the global grid and gives an
        /// array with only entries for the specific sub-grid
        /// </summary>
        /// <param name="results">Result for the complete Grid</param>
        /// <returns>Array with entries only for the sgrd-cells</returns>
        private double[] OrderValuesBySgrd(double[] results) {
            double[] ordered = new double[Mapping.LocalLength];

            for (int j = 0; j < ABSubGrid.LocalNoOfCells; j++) {
                int cell = ABSubGrid.SubgridIndex2LocalCellIndex[j];
                // cell in sgrd
                // f== each field
                // n== basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        ordered[index] = results[index];
                    }
                }
            }
            return ordered;
        }
    }
}
