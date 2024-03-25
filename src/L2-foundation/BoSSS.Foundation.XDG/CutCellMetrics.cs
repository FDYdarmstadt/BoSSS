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

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {


    public interface ICutCellMetrics {
        /// <summary>
        /// Quadrature order used to obtain the metrics
        /// </summary>
        int CutCellQuadratureOrder { get; }

        /// <summary>
        /// available species
        /// </summary>
        IEnumerable<SpeciesId> SpeciesList { get; }


        /// <summary>
        /// Volume of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        Dictionary<SpeciesId, MultidimensionalArray> CutCellVolumes { get; }


        /// <summary>
        /// Surface (Level-Set plus cut edges) of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        Dictionary<SpeciesId, MultidimensionalArray> CellSurface { get; }
    }


    /// <summary>
    /// HMF-based computation of cut-cell volumes and cut-edge areas.
    /// </summary>
    public class CutCellMetrics : ICutCellMetrics {

        /// <summary>
        /// owner object.
        /// </summary>
        public XDGSpaceMetrics XDGSpaceMetrics {
            get;
            private set;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        internal CutCellMetrics(XDGSpaceMetrics owner) {
            XDGSpaceMetrics = owner;
            ComputeNonAgglomeratedMetrics();
        }

        /// <summary>
        /// The quadrature order used for computing cell volumes and edge areas.
        /// </summary>
        public int CutCellQuadratureOrder {
            get {
                return XDGSpaceMetrics.CutCellQuadOrder;
            }
        }

        /// <summary>
        /// The kind of HMF which is be used for computing cell volumes.
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants HMFvariant {
            get {
                return XDGSpaceMetrics.CutCellQuadratureType;
            }
        }

        /// <summary>
        /// All species for which agglomeration is available.
        /// </summary>
        public IEnumerable<SpeciesId> SpeciesList {
            get {
                return XDGSpaceMetrics.SpeciesList;
            }
        }

        /// <summary>
        /// Volume of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutCellVolumes {
            get;
            private set;
        }

        /// <summary>
        /// Volume of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> InterfaceArea {
            get;
            private set;
        }

        /// <summary>
        /// Area of non-agglomerated cut edges.
        /// - key: species
        /// - index: edge index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutEdgeAreas {
            get;
            private set;
        }

        /// <summary>
        /// Surface (Level-Set plus cut edges) of non-agglomerated cut cells.
        /// - key: species
        /// - index: cell index
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CellSurface {
            get;
            private set;
        }


        /// <summary>
        /// Computes Cell-volumes and edge areas before agglomeration.
        /// </summary>
        void ComputeNonAgglomeratedMetrics() {
            using(var tr = new FuncTrace()) {
                MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD);

                var gd = XDGSpaceMetrics.GridDat;
                int JE = gd.iLogicalCells.Count;
                int J = gd.iLogicalCells.NoOfLocalUpdatedCells;
                int D = gd.SpatialDimension;
                int EE = gd.Edges.Count;
                SpeciesId[] species = this.SpeciesList.ToArray();
                int NoSpc = species.Count();
                int[,] E2C = gd.iLogicalEdges.CellIndices;

                var schH = new XQuadSchemeHelper(XDGSpaceMetrics);

                // collect all per-cell-metrics in the same MultidimArry, for MPI-exchange (only 1 exchange instead of three, saving some overhead)
                // 1st index: cell
                // 2nd index: species
                // 3rd index: 0 for interface surface per cell, 1 for cut-cell-volume, 2 for cut-cell surface
                double[] vec_cellMetrics = new double[JE * NoSpc * 3];
                MultidimensionalArray cellMetrics = MultidimensionalArray.CreateWrapper(vec_cellMetrics, JE, NoSpc, 3);

                this.CutEdgeAreas = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CutCellVolumes = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.InterfaceArea = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CellSurface = new Dictionary<SpeciesId, MultidimensionalArray>();

                //var schS = new List<CellQuadratureScheme>();
                //var rulz = new List<ICompositeQuadRule<QuadRule>>();

                tr.Info("Checkpoint1");

                // edges and volumes
                // =================
                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var cellVol = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1);
                    SpeciesId spc = species[iSpc];

                    // compute cut edge area
                    // ---------------------
                    MultidimensionalArray edgArea = MultidimensionalArray.Create(EE);
                    this.CutEdgeAreas.Add(spc, edgArea);

                    var edgeScheme = schH.GetEdgeQuadScheme(spc);
                    var edgeRule = edgeScheme.Compile(gd, this.CutCellQuadratureOrder);


                    BoSSS.Foundation.Quadrature.EdgeQuadrature.GetQuadrature(
                        new int[] { 1 }, gd,
                        edgeRule,
                        _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                        {
                            EvalResult.SetAll(1.0);
                        },
                        _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                        {
                            for (int i = 0; i < Length; i++) {
                                int iEdge = i + i0;
                                Debug.Assert(edgArea[iEdge] == 0);
                                edgArea[iEdge] = ResultsOfIntegration[i, 0];
                                Debug.Assert(!(double.IsNaN(edgArea[iEdge]) || double.IsInfinity(edgArea[iEdge])));
                            }
                        }).Execute();

                    // sum up edges for surface
                    // ------------------------

                    tr.Info("Checkpoint1.1");

                    var cellSurf = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 2);

                    for(int e = 0; e < EE; e++) {
                        double a = edgArea[e];
                        int jCell0 = E2C[e, 0];
                        int jCell2 = E2C[e, 1];
                        cellSurf[jCell0] += a;
                        if(jCell2 >= 0)
                            cellSurf[jCell2] += a;

                    }

                    tr.Info("Checkpoint1.2");

                    // compute cut cell volumes
                    // ------------------------

                    tr.Info("Species: " + XDGSpaceMetrics.Tracker.GetSpeciesName(spc));
                    var volScheme = schH.GetVolumeQuadScheme(spc);
                    tr.Info("Checkpoint1.3");
                    var volRule = volScheme.Compile(gd, this.CutCellQuadratureOrder);
                    tr.Info("Checkpoint1.4");

                    BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                        new int[] { 1 }, gd,
                        volRule,
                        _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                        {
                            EvalResult.SetAll(1.0);
                        },
                        _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                        {
                            for(int i = 0; i < Length; i++) {
                                int jCell = i + i0;
                                Debug.Assert(cellVol[jCell] == 0);
                                cellVol[jCell] = ResultsOfIntegration[i, 0];
                                //Console.WriteLine("cellVol " + cellVol[jCell] + " of cell " + jCell);
                                Debug.Assert(!(double.IsNaN(cellVol[jCell]) || double.IsInfinity(cellVol[jCell])));
                            }
                        }).Execute();
                }

                tr.Info("Checkpoint2");

                // interface surface
                // =================

                // loop over surfaces
                if (species.Length > 0) {
                    var AllSpc = XDGSpaceMetrics.TotalSpeciesList;
                    var requiredSpecies = XDGSpaceMetrics.SpeciesList;

                    // loop over all possible pairs of species
                    for (int iSpcA = 0; iSpcA < AllSpc.Count - 1; iSpcA++) {
                        var SpeciesA = AllSpc[iSpcA];
                        var SpeciesADom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesA);
                        int iLocalSpcA = requiredSpecies.IndexOf(SpeciesA);

                        // Standard XDG case, where level sets are always
                        // an interface between species
                        for (int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                            var SpeciesB = AllSpc[iSpcB];
                            int iLocalSpcB = requiredSpecies.IndexOf(SpeciesB);
                            
                            if(iLocalSpcA > -1 || iLocalSpcB > -1) {
                                var SpeciesBDom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesB);
                                var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);
                                int NoOfLs = XDGSpaceMetrics.NoOfLevelSets;
                                for (int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                                    if (schH.SpeciesAreSeparatedByLevSet(iLevSet, SpeciesA, SpeciesB)) {
                                        var LsDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
                                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                        CellQuadratureScheme SurfIntegration = schH.GetLevelSetquadScheme(iLevSet, SpeciesA, IntegrationDom);
                                        var rule = SurfIntegration.Compile(gd, this.CutCellQuadratureOrder);
                                        BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                                            new int[] { 1 }, gd,
                                            rule,
                                            _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                                                EvalResult.SetAll(1.0);
                                            },
                                            _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                for (int i = 0; i < Length; i++) {
                                                    int jCell = i + i0;
                                                    if (iLocalSpcA > -1) {
                                                        cellMetrics[jCell, iLocalSpcA, 0] += ResultsOfIntegration[i, 0];
                                                        Debug.Assert(!(double.IsNaN(cellMetrics[jCell, iLocalSpcA, 0]) || double.IsInfinity(cellMetrics[jCell, iLocalSpcA, 0])));
                                                    }
                                                    if (iLocalSpcB > -1) {
                                                        cellMetrics[jCell, iLocalSpcB, 0] += ResultsOfIntegration[i, 0];
                                                        Debug.Assert(!(double.IsNaN(cellMetrics[jCell, iLocalSpcB, 0]) || double.IsInfinity(cellMetrics[jCell, iLocalSpcB, 0])));
                                                    }
                                                }
                                            }).Execute();
                                    }
                                }
                            }
                        }
                    }
                }

                tr.Info("Checkpoint3");

                // MPI exchange & store
                // ====================

                if (species.Length > 0) {
#if DEBUG
                    int NoOfSpc = species.Length;
                    var cellMetricsB4 = cellMetrics.ExtractSubArrayShallow(new[] { 0, 0, 0 }, new[] { J - 1, NoOfSpc - 1, 1 }).CloneAs();
#endif
                    vec_cellMetrics.MPIExchange(gd);
#if DEBUG
                    var cellMetricsComp = cellMetrics.ExtractSubArrayShallow(new[] { 0, 0, 0 }, new[] { J - 1, NoOfSpc - 1, 1 }).CloneAs();
                    cellMetricsComp.Acc(-1.0, cellMetricsB4);
                    Debug.Assert(cellMetricsComp.L2Norm() == 0.0);
#endif
                }
                for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                    var spc = species[iSpc];
                    this.InterfaceArea.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 0).CloneAs());
                    this.CutCellVolumes.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1).CloneAs());
                    this.CellSurface.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 2).CloneAs());

                    this.CellSurface[spc].Acc(1.0, this.InterfaceArea[spc]);
                }
            }
        }
    }
}