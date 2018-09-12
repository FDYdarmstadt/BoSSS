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
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {


    /// <summary>
    /// HMF-based computation of cut-cell volumes and cut-edge areas.
    /// </summary>
    public class CutCellMetrics {

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
        /// Computes Cell-volumes and edge areas before agglomeration.
        /// </summary>
        void ComputeNonAgglomeratedMetrics() {
            var gd = XDGSpaceMetrics.GridDat;
            int JE = gd.Cells.Count;
            int J = gd.Cells.NoOfLocalUpdatedCells;
            int D = gd.SpatialDimension;
            int EE = gd.Edges.Count;
            SpeciesId[] species = this.SpeciesList.ToArray();
            int NoSpc = species.Count();
            int[,] E2C = gd.Edges.CellIndices;

            var schH = new XQuadSchemeHelper(XDGSpaceMetrics);

            double[] vec_cellMetrics = new double[JE * NoSpc * 2];

            // collect all per-cell-metrics in the same MultidimArry, for MPI-exchange
            // 1st index: cell
            // 2nd index: species
            // 3rd index: 0 for interface surface per cell, 1 for cut-cell-volume
            MultidimensionalArray cellMetrics = MultidimensionalArray.CreateWrapper(vec_cellMetrics, JE, NoSpc, 2);

            this.CutEdgeAreas = new Dictionary<SpeciesId, MultidimensionalArray>();
            this.CutCellVolumes = new Dictionary<SpeciesId, MultidimensionalArray>();
            this.InterfaceArea = new Dictionary<SpeciesId, MultidimensionalArray>();

            //var schS = new List<CellQuadratureScheme>();
            //var rulz = new List<ICompositeQuadRule<QuadRule>>();


            // edges and volumes
            // =================
            for (int iSpc = 0; iSpc < species.Length; iSpc++) {
                var cellVol = cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1);
                var spc = species[iSpc];

                MultidimensionalArray edgArea = MultidimensionalArray.Create(EE);
                this.CutEdgeAreas.Add(spc, edgArea);

                // compute cut edge area
                // ---------------------

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

                // compute cut cell volumes
                // ------------------------

                var volScheme = schH.GetVolumeQuadScheme(spc);
                var volRule = volScheme.Compile(gd, this.CutCellQuadratureOrder);

                //schS.Add(volScheme);
                //rulz.Add(volRule);

                BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                    new int[] { 1 }, gd,
                    volRule,
                    _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                    {
                        EvalResult.SetAll(1.0);
                    },
                    _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                    {
                        for (int i = 0; i < Length; i++) {
                            int jCell = i + i0;
                            Debug.Assert(cellVol[jCell] == 0);
                            cellVol[jCell] = ResultsOfIntegration[i, 0];
                            Debug.Assert(!(double.IsNaN(cellVol[jCell]) || double.IsInfinity(cellVol[jCell])));
                        }
                    }).Execute();
            }

            
            // interface surface
            // =================

            // loop over all possible pairs of species
            /*
            for (int iSpcA = 0; iSpcA < species.Length; iSpcA++) {
                var SpeciesA = species[iSpcA];
                var SpeciesADom = lsTrk.GetSpeciesMask(SpeciesA);
                if (SpeciesADom.NoOfItemsLocally <= 0)
                    continue;

                for (int iSpcB = iSpcA + 1; iSpcB < species.Length; iSpcB++) {
                    var SpeciesB = species[iSpcB];
                    var SpeciesBDom = lsTrk.GetSpeciesMask(SpeciesB);
                    if (SpeciesBDom.NoOfItemsLocally <= 0)
                        continue;

                    var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);
                    if (SpeciesCommonDom.NoOfItemsLocally <= 0)
                        continue;

                    // loop over level-sets
                    int NoOfLs = lsTrk.LevelSets.Count;
                    for (int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {


                        var LsDom = lsTrk.GetCutCellMask4LevSet(iLevSet);
                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                        if (IntegrationDom.NoOfItemsLocally > 0) {
                            CellQuadratureScheme SurfIntegration = schH.GetLevelSetquadScheme(iLevSet, IntegrationDom);
                            var rule = SurfIntegration.Compile(gd, this.HMForder);


                            BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(
                                new int[] { 1 }, gd,
                                rule,
                                _Evaluate: delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) //
                                {
                                    EvalResult.SetAll(1.0);
                                },
                                _SaveIntegrationResults: delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) //
                                {
                                    for (int i = 0; i < Length; i++) {
                                        int jCell = i + i0;
                                        if (cellMetrics[jCell, iSpcA, 0] != 0.0)
                                            throw new NotSupportedException("More than one zero-level-set per cell is not supported yet.");
                                        if (cellMetrics[jCell, iSpcB, 0] != 0.0)
                                            throw new NotSupportedException("More than one zero-level-set per cell is not supported yet.");

                                        cellMetrics[jCell, iSpcA, 0] = ResultsOfIntegration[i, 0];
                                        cellMetrics[jCell, iSpcB, 0] = ResultsOfIntegration[i, 0];
                                    }
                                }).Execute();
                        }
                    }
                }
            }*/

            // loop over level-sets
            if (species.Length > 0) {
                // find domain of all species: 
                CellMask SpeciesCommonDom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(species[0]);
                for (int iSpc = 1; iSpc < species.Length; iSpc++) {
                    SpeciesCommonDom = SpeciesCommonDom.Union(XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(species[iSpc]));
                }
                BitArray[] SpeciesBitMask = species.Select(spc => XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(spc).GetBitMask()).ToArray();

                int NoOfLs = XDGSpaceMetrics.NoOfLevelSets;
                int NoOfSpc = species.Length;
                for (int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                    
                    var LsDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
                    var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                    //if (IntegrationDom.NoOfItemsLocally > 0) { -> Doesn't work if Bjoerns HMF is used, eds up in an mpi dead lock
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // this level-set is actually relevant for someone in 'species'
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        CellQuadratureScheme SurfIntegration = schH.GetLevelSetquadScheme(iLevSet, IntegrationDom);
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



                                    //if (cellMetrics[jCell, iSpcA, 0] != 0.0)
                                    //    throw new NotSupportedException("More than one zero-level-set per cell is not supported yet.");
                                    //if (cellMetrics[jCell, iSpcB, 0] != 0.0)
                                    //    throw new NotSupportedException("More than one zero-level-set per cell is not supported yet.");

                                    //cellMetrics[jCell, iSpcA, 0] = ResultsOfIntegration[i, 0];
                                    //cellMetrics[jCell, iSpcB, 0] = ResultsOfIntegration[i, 0];


                                    for(int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                                        if(SpeciesBitMask[iSpc][jCell]) {
                                            cellMetrics[jCell, iSpc, 0] += ResultsOfIntegration[i, 0];
                                            Debug.Assert(!(double.IsNaN(cellMetrics[jCell, iSpc, 0]) || double.IsInfinity(cellMetrics[jCell, iSpc, 0])));
                                        }
                                
                                    }
                                }
                            }).Execute();
                    //}
                }

            }

            // MPI exchange & store
            // ====================

            if(species.Length > 0)
                vec_cellMetrics.MPIExchange(gd);

            for (int iSpc = 0; iSpc < species.Length; iSpc++) {
                var spc = species[iSpc];
                this.InterfaceArea.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 0).CloneAs());
                this.CutCellVolumes.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1).CloneAs());
            }

            //Console.WriteLine("Erinn - debug code.");
            //for(int j = 0; j < J; j++) {
            //    double totVol = gd.Cells.GetCellVolume(j);

            //    double blaVol = 0;
            //    for(int iSpc = 0; iSpc < species.Length; iSpc++) {
            //        var spc = species[iSpc];
            //        var cellVolS = this.CutCellVolumes[spc][j];
            //        blaVol += cellVolS;
            //    }

            //    Debug.Assert(Math.Abs(totVol - blaVol) / totVol < 1.0e-8);

            //}
        }
    }
}
