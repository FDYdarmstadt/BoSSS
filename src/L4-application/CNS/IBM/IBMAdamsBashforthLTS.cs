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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.IO;

namespace CNS.IBM {

    class IBMAdamsBashforthLTS : AdamsBashforthLTS {

        private ImmersedSpeciesMap speciesMap;

        private SpatialOperator boundaryOperator;

        private CoordinateMapping boundaryParameterMap;

        private Lazy<SpatialOperator.Evaluator> boundaryEvaluator;

        private CellMask cutCells;

        private CellMask cutAndTargetCells;

        private SpatialOperator standardOperator;

        private CoordinateMapping fieldsMap;

        private IBMControl control;

        public IBMAdamsBashforthLTS(SpatialOperator standardOperator, SpatialOperator boundaryOperator, CoordinateMapping fieldsMap, CoordinateMapping boundaryParameterMap, ISpeciesMap ibmSpeciesMap, IBMControl control, IList<TimeStepConstraint> timeStepConstraints, int reclusteringInterval, bool fluxCorrection)
            : base(standardOperator, fieldsMap, boundaryParameterMap, control.ExplicitOrder, control.NumberOfSubGrids, true, timeStepConstraints, reclusteringInterval: reclusteringInterval, fluxCorrection: fluxCorrection, subGrid: ibmSpeciesMap.SubGrid) {

            this.speciesMap = ibmSpeciesMap as ImmersedSpeciesMap;
            if (this.speciesMap == null) {
                throw new ArgumentException(
                    "Only supported for species maps of type 'ImmersedSpeciesMap'",
                    "speciesMap");
            }
            this.standardOperator = standardOperator;
            this.boundaryOperator = boundaryOperator;
            this.boundaryParameterMap = boundaryParameterMap;
            this.fieldsMap = fieldsMap;
            this.control = control;

            agglomerationPatternHasChanged = true;

            cutCells = speciesMap.Tracker.Regions.GetCutCellMask();
            cutAndTargetCells = cutCells.Union(speciesMap.Agglomerator.AggInfo.TargetCells);

            // Normal LTS constructor
            NumberOfLocalTimeSteps = new List<int>(control.NumberOfSubGrids);

            clusterer = new Clusterer(this.gridData, this.TimeStepConstraints);

            CurrentClustering = clusterer.CreateClustering(control.NumberOfSubGrids, speciesMap.SubGrid);
            CurrentClustering = CalculateNumberOfLocalTS(CurrentClustering); // Might remove sub-grids when time step sizes are too similar

            ABevolver = new IBMABevolve[CurrentClustering.NumberOfClusters];

            for (int i = 0; i < ABevolver.Length; i++) {
                ABevolver[i] = new IBMABevolve(standardOperator, boundaryOperator, fieldsMap, boundaryParameterMap, speciesMap, control.ExplicitOrder, control.LevelSetQuadratureOrder, control.CutCellQuadratureType, sgrd: CurrentClustering.Clusters[i], adaptive: this.adaptive);
                ABevolver[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
            }

            GetBoundaryTopology();

#if DEBUG
            for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                Console.WriteLine("IBM AB LTS ctor: id=" + i + " -> sub-steps=" + NumberOfLocalTimeSteps[i] + " and elements=" + CurrentClustering.Clusters[i].GlobalNoOfCells);
            }
#endif

            // Start-up phase needs an IBM Runge-Kutta time stepper
            RungeKuttaScheme = new IBMSplitRungeKutta(standardOperator, boundaryOperator, fieldsMap, boundaryParameterMap, speciesMap, timeStepConstraints);
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

            this.m_Evaluator = new Lazy<SpatialOperator.Evaluator>(() =>
                this.Operator.GetEvaluatorEx(
                    Mapping,
                    null, // TO DO: I SIMPLY REMOVE PARAMETERMAP HERE; MAKE THIS MORE PRETTY
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

        protected override void ComputeChangeRate(double[] k, double AbsTime, double RelTime, double[] edgeFluxes = null) {
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

        /// <summary>
        /// Required by <see cref="Perform(double)"/>
        /// </summary>
        internal static bool agglomerationPatternHasChanged = true;

        public override double Perform(double dt) {
            if (agglomerationPatternHasChanged) {
                // TO DO: Agglomerate difference between old $cutAndTargetCells and new $cutAndTargetCells only
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

                //Broadcast to RungeKutta and ABevolve ???
                foreach (IBMABevolve evolver in ABevolver) {
                    evolver.BuildEvaluatorsAndMasks();
                    evolver.agglomerationPatternHasChanged = false;
                }
            }

            dt = base.Perform(dt);

            speciesMap.Agglomerator.Extrapolate(DGCoordinates.Mapping);
            return dt;
        }

        protected override void CreateNewABevolver() {
            // Create array of Abevolve objects based on the new clustering
            ABevolver = new IBMABevolve[CurrentClustering.NumberOfClusters];

            for (int i = 0; i < ABevolver.Length; i++) {
                ABevolver[i] = new IBMABevolve(standardOperator, boundaryOperator, fieldsMap, boundaryParameterMap, speciesMap, control.ExplicitOrder, control.LevelSetQuadratureOrder, control.CutCellQuadratureType, sgrd: CurrentClustering.Clusters[i], adaptive: this.adaptive);
                ABevolver[i].ResetTime(m_Time);
                ABevolver[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
            }
        }
    }
}
