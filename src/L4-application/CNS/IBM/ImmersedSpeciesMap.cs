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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using CNS.MaterialProperty;
using System;
using System.IO;
using System.Linq;
using ilPSP;

namespace CNS.IBM {

    /// <summary>
    /// A species map for a single species that is embedded into domain, i.e.
    /// that encloses and/or is enclosed by a void domain. 
    /// </summary>
    public class ImmersedSpeciesMap : ISpeciesMap, IObserver<LevelSetTracker.LevelSetRegionsInfo> {

        /// <summary>
        /// The material/fluid of the represented, immersed fluid.
        /// </summary>
        private readonly Material material;

        /// <summary>
        /// The control file object
        /// </summary>
        public IBMControl Control {
            get;
            private set;
        }

        bool firstCall = true;

        private LevelSetTracker tracker;

        /// <summary>
        /// The tracker for the level set that defines the interface between
        /// the fluid and the void region
        /// </summary>
        public LevelSetTracker Tracker {
            get {
                if (firstCall) {
                    // This IS necessary... don't ask me why, though
                    tracker.UpdateTracker();
                    tracker.UpdateTracker();
                }
                firstCall = false;
                return tracker;
            }
            private set {
                this.tracker = value;
            }
        }

        /// <summary>
        /// Backing field for <see cref="QuadSchemeHelper"/>
        /// </summary>
        private XQuadSchemeHelper quadSchemeHelper;

        private CutCellMetrics lastCutCellMetrics;

        /// <summary>
        /// Quadrature scheme helper for the integration over the species with
        /// name <see cref="IBMControl.FluidSpeciesName"/>
        /// </summary>
        public XQuadSchemeHelper QuadSchemeHelper {
            get {
                if (quadSchemeHelper == null) {
                    if (Control.MomentFittingVariant == XQuadFactoryHelper.MomentFittingVariants.Classic) {
                        BoSSS.Foundation.XDG.Quadrature.HMF.LevelSetSurfaceQuadRuleFactory.UseNodesOnLevset =
                            Control.SurfaceHMF_ProjectNodesToLevelSet;
                        BoSSS.Foundation.XDG.Quadrature.HMF.LevelSetSurfaceQuadRuleFactory.RestrictNodes =
                            Control.SurfaceHMF_RestrictNodes;
                        BoSSS.Foundation.XDG.Quadrature.HMF.LevelSetSurfaceQuadRuleFactory.UseGaussNodes =
                            Control.SurfaceHMF_UseGaussNodes;

                        BoSSS.Foundation.XDG.Quadrature.HMF.LevelSetVolumeQuadRuleFactory.NodeCountSafetyFactor =
                            Control.VolumeHMF_NodeCountSafetyFactor;
                        BoSSS.Foundation.XDG.Quadrature.HMF.LevelSetVolumeQuadRuleFactory.RestrictNodes =
                            Control.VolumeHMF_RestrictNodes;
                        BoSSS.Foundation.XDG.Quadrature.HMF.LevelSetVolumeQuadRuleFactory.UseGaussNodes =
                            Control.VolumeHMF_UseGaussNodes;
                    }

                    SpeciesId[] species = new SpeciesId[] {
                        Tracker.GetSpeciesId(Control.FluidSpeciesName)
                    };

                    CutCellMetrics cutCellMetrics = new CutCellMetrics(
                        Control.MomentFittingVariant, Control.LevelSetQuadratureOrder, Tracker, species);

                    bool agglomerateNewbornAndDeceased = true;
                    var oldCCM = new CutCellMetrics[] { lastCutCellMetrics };
                    var oldAggThreshold = new double[] { Control.AgglomerationThreshold };
                    if (lastCutCellMetrics == null) {
                        agglomerateNewbornAndDeceased = false;
                        lastCutCellMetrics = cutCellMetrics;
                        oldAggThreshold = null;
                        oldCCM = null;
                    }
                    MultiphaseCellAgglomerator agglomerator = new MultiphaseCellAgglomerator(
                        cutCellMetrics,
                        Control.AgglomerationThreshold,
                        AgglomerateNewborn: agglomerateNewbornAndDeceased, AgglomerateDecased: agglomerateNewbornAndDeceased,
                        oldCcm: oldCCM,
                        oldTs__AgglomerationTreshold: oldAggThreshold
                        );


                    lastCutCellMetrics = cutCellMetrics;

                    var speciesAgglomerator = agglomerator.GetAgglomerator(
                        Tracker.GetSpeciesId(Control.FluidSpeciesName));
                    if (Control.PrintAgglomerationInfo) {

                        bool stdoutOnlyOnRank0 = ilPSP.Environment.StdoutOnlyOnRank0;
                        ilPSP.Environment.StdoutOnlyOnRank0 = false;
                        Console.WriteLine(
                            "Agglomerating {0} cells on rank {1}",
                            speciesAgglomerator.AggInfo.AgglomerationPairs.Length,
                            Tracker.GridDat.MpiRank);
                        ilPSP.Environment.StdoutOnlyOnRank0 = stdoutOnlyOnRank0;
                    }

                    if (Control.SaveAgglomerationPairs) {
                        int i = 0;
                        string fileName;
                        do {
                            fileName = String.Format(
                                "agglomerationParis_rank{0}_{1}.txt", Tracker.GridDat.MpiRank, i);
                            i++;
                        } while (File.Exists(fileName));

                        speciesAgglomerator.PlotAgglomerationPairs(fileName, includeDummyPointIfEmpty: true);
                    }

                    quadSchemeHelper = new XQuadSchemeHelper(agglomerator);
                }

                return quadSchemeHelper;
            }
        }

        /// <summary>
        /// Agglomerator for small/ill-shaped cells
        /// </summary>
        public CellAgglomerator Agglomerator {
            get {
                return QuadSchemeHelper.CellAgglomeration.GetAgglomerator(
                    Tracker.GetSpeciesId(Control.FluidSpeciesName));
            }
        }

        /// <summary>
        /// Constructs a new species map using the level set information
        /// provided by <paramref name="levelSet"/>.
        /// </summary>
        /// <param name="control"></param>
        /// <param name="levelSet"></param>
        /// <param name="material"></param>
        public ImmersedSpeciesMap(IBMControl control, LevelSet levelSet, Material material) {
            this.Control = control;
            this.material = material;

            this.tracker = new LevelSetTracker(
                (GridData) levelSet.GridDat,
                0,
                new string[] { Control.VoidSpeciesName, Control.FluidSpeciesName },
                levelSet);
            tracker.Subscribe(this);
        }

        /// <summary>
        /// Creates an instance of <see cref="IBMMassMatrixFactory"/> that is
        /// appropriate for the given <paramref name="mapping"/>
        /// </summary>
        /// <param name="mapping"></param>
        /// <returns></returns>
        public IBMMassMatrixFactory GetMassMatrixFactory(CoordinateMapping mapping) {
            if (MassMatrixFactory == null) {
                MassMatrixFactory = new IBMMassMatrixFactory(this, mapping);
            }
            return MassMatrixFactory;
        }

        private IBMMassMatrixFactory MassMatrixFactory;


        /// <summary>
        /// Adjusts h_min of the orginal mesh, due to cut-cells
        /// </summary>
        public MultidimensionalArray h_min
        {
            get
            {
                if (hMin == null) {
                    if (MassMatrixFactory == null) return null;
                    hMin = MultidimensionalArray.Create(GridData.Cells.NoOfLocalUpdatedCells);
                    MultidimensionalArray hMinUncut = GridData.Cells.h_min;
                    CellMask cutCells = this.Tracker._Regions.GetCutCellMask();
                    CellMask cutAndTargetCells = cutCells.Union(this.Agglomerator.AggInfo.TargetCells).Except(this.Agglomerator.AggInfo.SourceCells);
                    var massMatrix = this.MassMatrixFactory.MassMatrix;

                    for (int iCell = 0; iCell < hMin.Length; iCell++) {
                        // Iterates over all Chunks and looks if cell "i" is inside, i.e "i" is cut and/or target cell
                        if (cutAndTargetCells.SelectMany(c => c.Elements).Contains(iCell)) {
                            // Calculate Barycenter and smallest distance to neighbouring edge
                            
                            hMin[iCell] = CalculateHmin(IBMUtility.GetBlock( massMatrix, iCell), iCell);
                        } else {
                            hMin[iCell] = hMinUncut[iCell];
                        }
                    }

                    // Conservative assumption for agglomerated cells:
                    // Both cells have h_min of target cell, i.e
                    // the increase of h_min due to the addional source cell is only marginal.
                    foreach (var agg_pair in Agglomerator.AggInfo.AgglomerationPairs) {
                        hMin[agg_pair.jCellSource] = hMin[agg_pair.jCellTarget];
                    }
                }

                return hMin;
            }

        }

        private MultidimensionalArray hMin;

        //Assumes linear cells
        //Calculates the barycenter for each cell and then the minimal distance to each face
        private double CalculateHmin(ilPSP.Utils.IMatrix massMatrix, int cell) {
            int[] edges = GridData.Cells.Cells2Edges[cell];

            // Calculate Barycenter of Cell via MassMatrix
            double[] cellCenterLocal = new double[GridData.SpatialDimension];
            // y-coordinate
            cellCenterLocal[1] = 1.0 / Math.Sqrt(3) * massMatrix[0, 1] / massMatrix[0, 0];
            // x-coordinate
            cellCenterLocal[0] = 1.0 / Math.Sqrt(3) * massMatrix[0, 2] / massMatrix[0, 0];
            MultidimensionalArray cellCenterGlobal = MultidimensionalArray.Create(1, 2);
            GridData.TransformLocal2Global(MultidimensionalArray.CreateWrapper(cellCenterLocal, 1, 2), cellCenterGlobal, cell);

            double minimalDistance = this.Control.LevelSetFunction(new double[] { cellCenterGlobal[0, 0], cellCenterGlobal[0, 1] }, 0.0);

            // Calculate distance from barycenter with Hesse normal form
            int noOfFaces = GridData.Grid.RefElements[0].NoOfFaces;
            var edgeCenters = GridData.Grid.RefElements[0].FaceCenters;
            MultidimensionalArray edgeCentersGlobal = MultidimensionalArray.Create(edgeCenters.Lengths);
            GridData.TransformLocal2Global(edgeCenters, edgeCentersGlobal, cell);
            for (int ec = 0; ec < edges.Length; ec++) {
                int iEdge = edges[ec];
                double sign = Math.Sign(iEdge);
                iEdge = Math.Abs(iEdge) - 1;
                int iFace = GridData.Edges.FaceIndices[iEdge, 0];

                if (sign < 0.0) {
                    iFace = GridData.Edges.FaceIndices[iEdge, 1];
                }

                double xNormal = GridData.Edges.NormalsForAffine[iEdge, 0] * sign;
                double yNormal = GridData.Edges.NormalsForAffine[iEdge, 1] * sign;

                double edgePoint_x = edgeCentersGlobal[iFace, 0];
                double edgePoint_y = edgeCentersGlobal[iFace, 1];

                // Note: distance should be always positiv, because normals are pointing outward
                double distance = (edgePoint_x * xNormal + edgePoint_y * yNormal) - (cellCenterGlobal[0, 0] * xNormal + cellCenterGlobal[0, 1] * yNormal);
                System.Diagnostics.Debug.Assert(distance > 0.0, "Calculated distance should be a positive number, but was " + distance);
                if (distance < minimalDistance) {
                    minimalDistance = distance;
                }
            }

            return 2 * minimalDistance;

        }

        #region ISpeciesMap Members

        /// <summary>
        /// Information about the grid
        /// </summary>
        public GridData GridData
        {
            get
            {
                return tracker.GridDat;
            }
        }

        /// <summary>
        /// Always returns the configured equation of state (see
        /// <see cref="ImmersedSpeciesMap.ImmersedSpeciesMap"/>).
        /// </summary>
        /// <param name="levelSetValue">Ignored</param>
        /// <returns>
        /// The configured equation of state
        /// </returns>
        public Material GetMaterial(double levelSetValue) {
            return material;
        }

        /// <summary>
        /// Returns the sub-grid associated with the represented fluid.
        /// </summary>
        public SubGrid SubGrid {
            get {
                return Tracker._Regions.GetSpeciesSubGrid(Control.FluidSpeciesName);
            }
        }

        #endregion

        #region IObservable<<evelSetData> Members

        /// <summary>
        /// Discards old quadrature information
        /// </summary>
        /// <param name="value"></param>
        public void OnNext(LevelSetTracker.LevelSetRegionsInfo value) {
            this.quadSchemeHelper = null;
            this.hMin = null;
            this.MassMatrixFactory = null;
        }

        /// <summary>
        /// Not implemented
        /// </summary>
        /// <param name="error"></param>
        public void OnError(System.Exception error) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        public void OnCompleted() {
        }

        #endregion
    }
}
