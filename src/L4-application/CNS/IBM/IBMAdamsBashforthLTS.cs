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
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Utils;

namespace CNS.IBM {

    class IBMAdamsBashforthLTS : AdamsBashforthLTS {

        private ImmersedSpeciesMap speciesMap;

        private DifferentialOperator boundaryOperator;

        private CoordinateMapping boundaryParameterMap;

        private Lazy<IEvaluatorNonLin> boundaryEvaluator;

        private CellMask cutCells;

        private CellMask cutAndTargetCells;

        private DifferentialOperator standardOperator;

        private CoordinateMapping fieldsMap;

        private IBMControl control;

        public IBMAdamsBashforthLTS(DifferentialOperator standardOperator, DifferentialOperator boundaryOperator, CoordinateMapping fieldsMap, CoordinateMapping boundaryParameterMap, ISpeciesMap ibmSpeciesMap, IBMControl control, IList<TimeStepConstraint> timeStepConstraints, int reclusteringInterval, bool fluxCorrection, bool forceReclustering, int maxNumOfSubSteps = 0, bool logging = false, bool consoleOutput = false)
            : base(standardOperator, fieldsMap, boundaryParameterMap, control.ExplicitOrder, control.NumberOfSubGrids, true, timeStepConstraints, reclusteringInterval: reclusteringInterval, fluxCorrection: fluxCorrection, subGrid: ibmSpeciesMap.SubGrid, forceReclustering: forceReclustering, logging: logging, consoleOutput: consoleOutput) {

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

            if (consoleOutput) {
                Console.WriteLine("### This is IBM ABLTS ctor ###");
            }

            // Normal LTS constructor
            clusterer = new Clusterer(this.gridData, maxNumOfSubSteps, cellAgglomerator: speciesMap.Agglomerator);
            CurrentClustering = clusterer.CreateClustering(control.NumberOfSubGrids, this.TimeStepConstraints, speciesMap.SubGrid);
            CurrentClustering = clusterer.TuneClustering(CurrentClustering, Time, this.TimeStepConstraints); // Might remove sub-grids when time step sizes are too similar

            ABevolver = new IBMABevolve[CurrentClustering.NumberOfClusters];

            for (int i = 0; i < ABevolver.Length; i++) {
                ABevolver[i] = new IBMABevolve(standardOperator, boundaryOperator, fieldsMap, boundaryParameterMap, speciesMap, control.ExplicitOrder, control.LevelSetQuadratureOrder, control.CutCellQuadratureType, sgrd: CurrentClustering.Clusters[i], adaptive: this.adaptive);
                ABevolver[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
            }

            GetBoundaryTopology();

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

            this.m_Evaluator = new Lazy<IEvaluatorNonLin>(delegate () {
                this.Operator.EdgeQuadraturSchemeProvider = g => edgeScheme;
                this.Operator.VolumeQuadraturSchemeProvider = g => volumeScheme;
                var opi = this.Operator.GetEvaluatorEx(
                    Mapping,
                    null, // TO DO: I SIMPLY REMOVE PARAMETERMAP HERE; MAKE THIS MORE PRETTY
                          //this.boundaryParameterMap, // TO DO: I SIMPLY REMOVE PARAMETERMAP HERE; MAKE THIS MORE PRETTY
                    Mapping);
                opi.ActivateSubgridBoundary(volumeScheme.Domain, subGridBoundaryTreatment: SubGridBoundaryModes.InnerEdgeLTS);
                //opi.ActivateSubgridBoundary(fluidCells, subGridBoundaryTreatment: SubGridBoundaryModes.InnerEdgeLTS);
                return opi;
            });

            // Evaluator for boundary conditions at level set zero contour
            CellQuadratureScheme boundaryVolumeScheme = speciesMap.QuadSchemeHelper.GetLevelSetquadScheme(
                0, cutCells, control.LevelSetQuadratureOrder);

            this.boundaryEvaluator = new Lazy<IEvaluatorNonLin>(delegate() {
                boundaryOperator.EdgeQuadraturSchemeProvider = g => null; // Contains no boundary terms
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

                // Agglomerate fluxes
                speciesMap.Agglomerator.ManipulateRHS(k, Mapping);

                // Apply inverse to all cells with non-identity mass matrix
                IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
                IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, k, 0.0, k, cutAndTargetCells);
            }
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

                AgglomerateAndExtrapolateDGCoordinates();

                agglomerationPatternHasChanged = false;

                //Broadcast to RungeKutta and ABevolve ???
                foreach (IBMABevolve evolver in ABevolver) {
                    evolver.BuildEvaluatorsAndMasks();
                    evolver.agglomerationPatternHasChanged = false;
                }
            }

            dt = base.Perform(dt);

            speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping);
            return dt;
        }

        protected void AgglomerateAndExtrapolateDGCoordinates() {
            IBMMassMatrixFactory massMatrixFactory = speciesMap.GetMassMatrixFactory(Mapping);
            BlockMsrMatrix nonAgglomeratedMassMatrix = massMatrixFactory.NonAgglomeratedMassMatrix;

            IBMUtility.SubMatrixSpMV(nonAgglomeratedMassMatrix, 1.0, CurrentState, 0.0, CurrentState, cutCells);  // eq. (39)
            speciesMap.Agglomerator.ManipulateRHS(CurrentState, Mapping);  // eq. (39)
            IBMUtility.SubMatrixSpMV(massMatrixFactory.InverseMassMatrix, 1.0, CurrentState, 0.0, CurrentState, cutAndTargetCells);   // eq. (39)
            speciesMap.Agglomerator.Extrapolate(CurrentState.Mapping); // eq. (41)
        }

        protected override void CreateNewABevolver() {
            // Create array of Abevolve objects based on the new clustering
            ABevolver = new IBMABevolve[CurrentClustering.NumberOfClusters];

            for (int i = 0; i < ABevolver.Length; i++) {
                ABevolver[i] = new IBMABevolve(standardOperator, boundaryOperator, fieldsMap, boundaryParameterMap, speciesMap, control.ExplicitOrder, control.LevelSetQuadratureOrder, control.CutCellQuadratureType, sgrd: CurrentClustering.Clusters[i], adaptive: this.adaptive);
                ABevolver[i].ResetTime(m_Time, TimeInfo.TimeStepNumber);
                ABevolver[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
            }
        }
    }
}
