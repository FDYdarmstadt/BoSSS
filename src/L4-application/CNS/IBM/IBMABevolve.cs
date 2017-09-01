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

namespace CNS.IBM {
    class IBMABevolve :ABevolve {

        private ImmersedSpeciesMap speciesMap;

        private SpatialOperator boundaryOperator;

        private CoordinateMapping boundaryParameterMap;

        private Lazy<SpatialOperator.Evaluator> boundaryEvaluator;

        private CellMask cutCells;

        private CellMask cutAndTargetCells;

        public IBMABevolve(
            SpatialOperator standardOperator,
            SpatialOperator boundaryOperator,
            CoordinateMapping fieldsMap,
            CoordinateMapping parametersMap,
            ISpeciesMap ibmSpeciesMap,
            int explicitOrder,
            int levelSetQuadratureOrder,
            XQuadFactoryHelper.MomentFittingVariants momentFittingVariant,
            SubGrid sgrd)
            : base(standardOperator, fieldsMap, null, explicitOrder, sgrd: sgrd)  {

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

            CellMask fluidCells = speciesMap.SubGrid.VolumeMask.Intersect(sgrd.VolumeMask);
            cutCells = speciesMap.Tracker._Regions.GetCutCellMask().Intersect(sgrd.VolumeMask);
            cutAndTargetCells = cutCells.Union(speciesMap.Agglomerator.AggInfo.TargetCells).Intersect(sgrd.VolumeMask);


            IBMControl control = speciesMap.Control;
            SpeciesId species = speciesMap.Tracker.GetSpeciesId(control.FluidSpeciesName);

            CellQuadratureScheme volumeScheme = speciesMap.QuadSchemeHelper.GetVolumeQuadScheme(
                species, true, fluidCells, control.LevelSetQuadratureOrder);


            // Does _not_ include agglomerated edges
            EdgeMask nonVoidEdges = speciesMap.QuadSchemeHelper.GetEdgeMask(species);
            nonVoidEdges = nonVoidEdges.Intersect(sgrd.AllEdgesMask);
            EdgeQuadratureScheme edgeScheme = speciesMap.QuadSchemeHelper.GetEdgeQuadScheme(
                species, true, nonVoidEdges, control.LevelSetQuadratureOrder);

            this.m_Evaluator = new Lazy<SpatialOperator.Evaluator>(() =>
                this.Operator.GetEvaluatorEx(
                    Mapping,
                    null, // TO DO: I SIMPLY REMOVE PARAMETERMAP HERE; MAKE THIS MORE PRETTY
                    Mapping,
                    edgeScheme,
                    volumeScheme,
                    sgrd,
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


        protected override void ComputeChangeRate(double[] k, double AbsTime, double RelTime, double[] edgeFluxes =null) {
            Evaluator.Evaluate(1.0, 0.0, k, AbsTime, outputBndEdge: edgeFluxes);
            Debug.Assert(
                !k.Any(f => double.IsNaN(f)),
                "Unphysical flux in standard terms");

            boundaryEvaluator.Value.Evaluate(1.0, 1.0, k, AbsTime);
            Debug.Assert(
                !k.Any(f => double.IsNaN(f)),
                "Unphysical flux in boundary terms");

            // Agglomerate fluxes
            speciesMap.Agglomerator.ManipulateRHS(k, Mapping);

            // Apply inverse to all cells with non-identity mass matrix
            IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
            IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, k, 0.0, k, cutAndTargetCells);
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
            double[] y0 = new double[Mapping.LocalLength];
            DGCoordinates.CopyTo(y0, 0);
            DGCoordinates.axpy<double[]>(completeChangeRate, -1);
            
            // Speciality for IBM: Do the extrapolation
            speciesMap.Agglomerator.Extrapolate(DGCoordinates.Mapping);

            double[] upDGC = DGCoordinates.ToArray();
            upDGC = orderValuesBySgrd(upDGC);
            // DGCoordinates should be untouched after calling this method
            DGCoordinates.Clear();
            DGCoordinates.CopyFrom(y0, 0);

            return upDGC;
        }

        /// <summary>
        /// Takes a double[] with results for the global grid and gives an
        /// array with only entries for the specific sub-grid
        /// </summary>
        /// <param name="results">Result for the complete Grid</param>
        /// <returns>Array with entries only for the sgrd-cells</returns>
        private double[] orderValuesBySgrd(double[] results) {
            double[] ordered = new double[Mapping.LocalLength];

            for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                int cell = sgrd.SubgridIndex2LocalCellIndex[j];
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
