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
    public class ImmersedSpeciesMap : ISpeciesMap, IObserver<LevelSetTracker.LevelSetRegions> {

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

        /// <summary>
        /// Hack required by <see cref="Tracker"/>
        /// </summary>
        private bool firstCall = true;

        /// <summary>
        /// Backing field for <see cref="Tracker"/>
        /// </summary>
        private LevelSetTracker tracker;

        /// <summary>
        /// The tracker for the level set that defines the interface between
        /// the fluid and the void region
        /// </summary>
        public LevelSetTracker Tracker {
            get {
                tracker.IncreaseHistoryLength(1);
                if (firstCall) {
                    // This IS necessary... don't ask me why, though
                    // to have one previous state available, so all the incremental updates work
                    tracker.UpdateTracker();
                    tracker.PushStacks();
                }

                firstCall = false;
                return tracker;
            }
            private set {
                this.tracker = value;
            }
        }

        /// <summary>
        /// Backing fields for <see cref="CellAgglomeration"/>
        /// </summary>
        private MultiphaseCellAgglomerator cellAgglomeration;

        /// <summary>
        /// Current cell agglomeration.
        /// </summary>
        public MultiphaseCellAgglomerator CellAgglomeration {
            get {
                //if (cellAgglomeration == null) {
                    if (Control.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Classic) {
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
                    
                    bool agglomerateNewbornAndDeceased = true;
                    var oldAggThreshold = new double[] { Control.AgglomerationThreshold };
                    if (Tracker.PopulatedHistoryLength <= 0) {
                        agglomerateNewbornAndDeceased = false;
                        oldAggThreshold = null;
                    }

                    cellAgglomeration = Tracker.GetAgglomerator(
                        new SpeciesId[] { Tracker.GetSpeciesId(Control.FluidSpeciesName) },
                        Control.LevelSetQuadratureOrder,
                        Control.AgglomerationThreshold,
                        AgglomerateNewborn: agglomerateNewbornAndDeceased,
                        AgglomerateDecased: agglomerateNewbornAndDeceased,
                        oldTs__AgglomerationTreshold: oldAggThreshold);

                    var speciesAgglomerator = cellAgglomeration.GetAgglomerator(
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
                //}

                return cellAgglomeration;
            }
        }

        /// <summary>
        /// Quadrature scheme helper for the integration over the species with
        /// name <see cref="IBMControl.FluidSpeciesName"/>
        /// </summary>
        public XQuadSchemeHelper QuadSchemeHelper {
            get {
                SpeciesId[] species = new SpeciesId[] {
                    Tracker.GetSpeciesId(Control.FluidSpeciesName)
                };

                return Tracker.GetXDGSpaceMetrics(species, Control.LevelSetQuadratureOrder, 1).XQuadSchemeHelper;
            }
        }

        /// <summary>
        /// Agglomerator for small/ill-shaped cells
        /// </summary>
        public CellAgglomerator Agglomerator {
            get {
                return this.CellAgglomeration.GetAgglomerator(
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
                (GridData)levelSet.GridDat,
                Control.CutCellQuadratureType,
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
            if (MassMatrixFactory == null || !mapping.Equals(MassMatrixFactory.Mapping)) {
                MassMatrixFactory = new IBMMassMatrixFactory(this, mapping, Control.FluidSpeciesName, Control.LevelSetQuadratureOrder);
            }
            return MassMatrixFactory;
        }

        private IBMMassMatrixFactory MassMatrixFactory;

        #region ISpeciesMap Members

        /// <summary>
        /// Information about the grid
        /// </summary>
        public GridData GridData {
            get {
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
                return Tracker.Regions.GetSpeciesSubGrid(Control.FluidSpeciesName);
            }
        }

        #endregion

        #region IObservable<<evelSetData> Members

        /// <summary>
        /// Discards old quadrature information
        /// </summary>
        /// <param name="value"></param>
        public void OnNext(LevelSetTracker.LevelSetRegions value) {
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
