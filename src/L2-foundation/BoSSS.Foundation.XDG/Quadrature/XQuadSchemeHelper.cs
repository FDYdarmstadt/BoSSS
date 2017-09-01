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
using System.Globalization;
using System.IO;
using System.Linq;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Provides quadrature schemes for typical XDG-cases;
    /// </summary>
    public class XQuadSchemeHelper {

        /// <summary>
        /// Cell agglomeration object/can be null if not provided at construction time.
        /// </summary>
        public MultiphaseCellAgglomerator CellAgglomeration {
            private set;
            get;
        }

        /// <summary>
        /// 
        /// </summary>
        public CutCellMetrics NonAgglomeratedMetrics {
            get;
            private set;
        }


        /// <summary>
        /// Selected variant of the moment-fitting procedure
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant {
            get;
            private set;
        }

        /// <summary>
        /// If true, <see cref="GetEdgeGhostScheme(SpeciesId, EdgeMask)"/> is supported, otherwise not.
        /// </summary>
        /// <remarks>
        /// Currently (april2016), nobody is using this feature (<see cref="XSpatialOperator.GhostEdgesOperator"/>),
        /// therefore this is hardcoded to false.
        /// </remarks>
        public bool GhostSupport {
            get;
            private set;
        }



        /// <summary>
        /// All species for which agglomeration is available.
        /// </summary>
        public IEnumerable<SpeciesId> SpeciesList {
            get;
            private set;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="_lsTrk"></param>
        /// <param name="momentFittingVariant"></param>
        public XQuadSchemeHelper(LevelSetTracker _lsTrk, XQuadFactoryHelper.MomentFittingVariants momentFittingVariant, params SpeciesId[] __SpeciesList) {
            MPICollectiveWatchDog.Watch();
            this.lsTrk = _lsTrk;
            this.MomentFittingVariant = momentFittingVariant;
            this.CellAgglomeration = null;
            this.NonAgglomeratedMetrics = null;
            this.SpeciesList = __SpeciesList.ToList().AsReadOnly();
            this.GhostSupport = false;

            ConstructorCommon();
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="_lsTrk"></param>
        /// <param name="momentFittingVariant"></param>
        public XQuadSchemeHelper(LevelSetTracker _lsTrk, XQuadFactoryHelper.MomentFittingVariants momentFittingVariant, int order = 0, params SpeciesId[] __SpeciesList) {

            MPICollectiveWatchDog.Watch();
            this.lsTrk = _lsTrk;
            this.MomentFittingVariant = momentFittingVariant;
            this.CellAgglomeration = null;
            this.NonAgglomeratedMetrics = new CutCellMetrics(momentFittingVariant, order, lsTrk, __SpeciesList);
            this.SpeciesList = __SpeciesList.ToList().AsReadOnly();
            this.GhostSupport = false;

            ConstructorCommon();
        }



        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="momentFittingVariant"></param>
        /// <param name="order"></param>
        /// <param name="NonAggCutCellMetric">Can be null, see <see cref="CellAgglomeration"/>.</param>
        public XQuadSchemeHelper(CutCellMetrics NonAggCutCellMetric, XQuadFactoryHelper.MomentFittingVariants momentFittingVariant) {
            MPICollectiveWatchDog.Watch();
            this.lsTrk = NonAggCutCellMetric.Tracker;
            this.MomentFittingVariant = momentFittingVariant;
            this.CellAgglomeration = null;
            this.NonAgglomeratedMetrics = NonAggCutCellMetric;
            this.SpeciesList = NonAgglomeratedMetrics.SpeciesList.ToList().AsReadOnly();
            this.GhostSupport = false;

            ConstructorCommon();
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="agg">Can be null, see <see cref="CellAgglomeration"/>.</param>
        public XQuadSchemeHelper(MultiphaseCellAgglomerator agg) {
            MPICollectiveWatchDog.Watch();
            this.lsTrk = agg.Tracker;
            this.MomentFittingVariant = agg.HMFvariant;
            this.CellAgglomeration = agg;
            this.NonAgglomeratedMetrics = agg.NonAgglomeratedMetrics;
            this.SpeciesList = agg.SpeciesList.ToList().AsReadOnly();
            this.GhostSupport = false;

            ConstructorCommon();
        }

        void ConstructorCommon() {
            var Krefs = this.lsTrk.GridDat.Grid.RefElements;


            // initialize some mask's and subgrids
            // ===================================

            // since most methods of this class are non-collective,
            // an access to some Subgrid-member may cause an MPI-deadlock.
            // (this becomes even more unpredictable due to the on-demand-pattern in which most members of the Subgrid-class are implemented);

            // 
            foreach (var spId in this.SpeciesList) {
                SubGrid spSgrd = this.lsTrk._Regions.GetSpeciesSubGrid(spId);

                this.m_SpeciesSubgrid_InnerAndDomainEdges.Add(
                    spId,
                    spSgrd.InnerEdgesMask.Union(spSgrd.AllEdgesMask.Intersect(lsTrk.GridDat.BoundaryEdges)));
            }

            // all cut edges
            // =============
            m_CutCellSubgrid_InnerEdges = this.lsTrk._Regions.GetCutCellSubGrid().InnerEdgesMask;

            // cut edges sorted according to reference element and level-set
            // =============================================================
            this.m_CutEdges = new EdgeMask[Krefs.Length, this.lsTrk.LevelSets.Count()];
            this.m_HMFEdgesDomain = new EdgeMask[Krefs.Length, this.lsTrk.LevelSets.Count()];
            for (int iKref = 0; iKref < this.m_CutEdges.GetLength(0); iKref++) {

                for (int iLevSet = 0; iLevSet < this.m_CutEdges.GetLength(1); iLevSet++) {
                    EdgeMask cutEdges = this.lsTrk._Regions.GetCutCellSubgrid4LevSet(iLevSet).InnerEdgesMask.Union(
                            this.lsTrk._Regions.GetCutCellSubgrid4LevSet(iLevSet).AllEdgesMask.Intersect(lsTrk.GridDat.BoundaryEdges));

                    cutEdges = cutEdges.Intersect(this.lsTrk.GridDat.GetRefElementSubGrid(iKref).AllEdgesMask);

                    this.m_CutEdges[iKref, iLevSet] = cutEdges;

                    // (all edges of 'Kref'-elements) \cap (all edges of cells cut by 'iLevSet')  
                    this.m_HMFEdgesDomain[iKref, iLevSet] = lsTrk._Regions.GetCutCellSubgrid4LevSet(iLevSet).AllEdgesMask.Intersect(this.lsTrk.GridDat.GetRefElementSubGrid(iKref).AllEdgesMask);
                }
            }

            //if (this.GhostSupport == true) {
            //    for (int iKref = 0; iKref < Krefs.Length; iKref++) {
            //        m_Subgrid4Kref_AllEdges.Add(
            //            Krefs[iKref],
            //            this.lsTrk.GridDat.GetRefElementSubGrid(iKref).AllEdgesMask);
            //    }
            //}

            //this.GhostIntegrationDomain = new Dictionary<SpeciesId, EdgeMask>();

            
            /*
            if (this.CellAgglomeration != null) {
                if (!object.ReferenceEquals(this.CellAgglomeration.Tracker, this.lsTrk))
                    throw new ArgumentException();

                SpeciesId[] AgglomSpecies = this.CellAgglomeration.SpeciesList.ToArray();

                var Edge2Cell = this.lsTrk.GridDat.Edges.CellIndices;
                int J = this.lsTrk.GridDat.Cells.NoOfLocalUpdatedCells;
                int JE = this.lsTrk.GridDat.Cells.NoOfCells;
                int E = this.lsTrk.GridDat.Edges.Count;
                

                // eliminate "empty" edges
                for (int iSpc = 0; iSpc < AgglomSpecies.Length; iSpc++) {
                    var spId = AgglomSpecies[iSpc];
                    MultidimensionalArray EdgeArea = this.NonAgglomeratedMetrics.CutEdgeAreas[spId];
                    MultidimensionalArray Volumes = this.NonAgglomeratedMetrics.CutCellVolumes[spId];

                    
                    // compute ghost exclusions (part 1)
                    // =================================
                    List<int> GhostExclusions = new List<int>();
                    var E2C = this.lsTrk.GridDat.Edges.CellIndices;
                    foreach (int jEdge in this.GetEdgeMask(spId).ItemEnum) {
                        int jCell0 = E2C[jEdge, 0];
                        if (Volumes[jCell0] <= 1.0e-10) {
                            GhostExclusions.Add(jEdge);
                        } else {
                            int jCell1 = E2C[jEdge, 1];
                            if (jCell1 >= 0) {
                                if (Volumes[jCell1] <= 1.0e-10) {
                                    GhostExclusions.Add(jEdge);
                                }
                            }
                        }
                    }


                    // determine edges over which agglomeration is forbidden
                    // (because their measure is zero)
                    // =====================================================
                    {

                        var scheme = this.GetEdgeQuadScheme(spId).Compile(_lsTrk.GridDat, order);

                        foreach (ChunkRulePair<QuadRule> cqr in scheme) {
                            for (int iEdge = cqr.Chunk.i0; iEdge < cqr.Chunk.JE; iEdge++) {
                                var edgArea = cqr.Rule.Weights.Sum();
                                EdgeArea[iEdge] = edgArea;
                                if (edgArea <= 1.0e-12) {
                                    // edge has zero measure => should not be in the ghost scheme
                                    GhostExclusions.Add(iEdge);
                                }
                            }
                        }
                    }

                    if (GhostExclusions.Count > 0) {
                        var _newGhostIntegrationDomain = _GhostIntegrationDomain.GetBitMask().CloneAs();
                        foreach (int iEdge in GhostExclusions)
                            _newGhostIntegrationDomain[iEdge] = false;
                        _GhostIntegrationDomain = new EdgeMask(lsTrk.GridDat, _newGhostIntegrationDomain);
                    }

                    this.GhostIntegrationDomain.Add(spId, _GhostIntegrationDomain);
                }
            }
            */
        }

        /// <summary>
        /// All edges (for each reference element, for each level-set) for which the HMF must be used.
        /// 1st index: volume reference element 
        /// 2nd index: level set
        /// </summary>
        EdgeMask[,] m_HMFEdgesDomain;

        

        /// <summary>
        /// All edges which are cut by a level-set.
        /// 1st index: volume reference element 
        /// 2nd index: level set
        /// </summary>
        /// <remarks>
        /// Be aware that this migt contain also edges like <c>e</c>,
        /// which are, topologically inner edges of the cut-cell-subgrid but not cut by the level-set:
        /// <code>
        ///     o----------o-----------o
        ///     |          |           |
        ///     |     **********       |
        ///     |     *    |   *       |
        ///     o ----*----o---*-------o
        ///     |     *    |   *       |
        ///   ********     e    ************Level-Set
        ///     |          |           |
        ///     o----------o-----------o
        /// </code>         
        /// </remarks>
        EdgeMask[,] m_CutEdges;

        /// <summary>
        /// all edges which are cut by level set #<paramref name="iLevSet"/> and which belong to a cell with
        /// reference element <paramref name="Kref"/>.
        /// </summary>
        EdgeMask GetCutEdges(RefElement Kref, int iLevSet) {
            int iKref = this.lsTrk.GridDat.Grid.RefElements.IndexOf(Kref);
            return m_CutEdges[iKref, iLevSet];
        }
        

        /// <summary>
        /// all edges of cells which are cut by level set #<paramref name="iLevSet"/> and which belong to a cell with
        /// reference element <paramref name="Kref"/>.
        /// </summary>
        EdgeMask GetHMFEdgesDomain(RefElement Kref, int iLevSet) {
            int iKref = this.lsTrk.GridDat.Grid.RefElements.IndexOf(Kref, (a, b) => object.ReferenceEquals(a, b));
            return m_HMFEdgesDomain[iKref, iLevSet];
        }


        ///// <summary>
        ///// edge masks for the ghost domain, see <see cref="XSpatialOperator.GhostEdgesOperator"/>
        ///// </summary>
        //Dictionary<SpeciesId, EdgeMask> GhostIntegrationDomain;

        /// <summary>
        /// initialized by the constructor to avoid MPI-deadlocks;
        /// keys: species 'S' <br/>
        /// value: an edge-mask containing all edges that are at least partly covered by species 'S'
        /// </summary>
        Dictionary<SpeciesId, EdgeMask> m_SpeciesSubgrid_InnerAndDomainEdges = new Dictionary<SpeciesId, EdgeMask>();

        /// <summary>
        /// initialized by the constructor to avoid MPI-deadlocks;
        /// </summary>
        EdgeMask m_CutCellSubgrid_InnerEdges;

        /*
        /// <summary>
        /// initialized by the constructor to avoid MPI-deadlocks;
        /// </summary>
        Dictionary<RefElement, EdgeMask> m_Subgrid4Kref_AllEdges = new Dictionary<RefElement, EdgeMask>();
        */

        LevelSetTracker lsTrk;


        public EdgeQuadratureScheme Get_SurfaceElement_EdgeQuadScheme(SpeciesId sp) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");

            var allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].Intersect(this.m_CutCellSubgrid_InnerEdges);

            EdgeMask AggEdges = this.CellAgglomeration != null ? this.CellAgglomeration.GetAgglomerator(sp).AggInfo.AgglomerationEdges : null;
            if (AggEdges != null && AggEdges.NoOfItemsLocally > 0)
                allRelevantEdges = allRelevantEdges.Except(AggEdges);


            var edgeQrIns = new EdgeQuadratureScheme(false, allRelevantEdges);

            foreach (var Kref in lsTrk.GridDat.Grid.RefElements) {
                for (int iLevSet = 0; iLevSet < lsTrk.LevelSets.Count; iLevSet++) { // loop over level sets...
                    EdgeMask cutEdges = this.GetCutEdges(Kref, iLevSet);

                    var factory = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant).GetSurfaceElement_BoundaryRuleFactory(iLevSet, Kref);

                    edgeQrIns.AddFactory(factory, cutEdges);
                }
            }

            return edgeQrIns;
        }

        public CellQuadratureScheme Get_SurfaceElement_VolumeQuadScheme(SpeciesId sp) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");



            var spdom = lsTrk._Regions.GetSpeciesMask(sp);
            var IntegrationDom = this.lsTrk._Regions.GetCutCellMask().Intersect(spdom);

            var LevSetQrIns = new CellQuadratureScheme(false, IntegrationDom);

            foreach (var Kref in lsTrk.GridDat.Grid.RefElements) {
                for (int iLevSet = 0; iLevSet < lsTrk.LevelSets.Count; iLevSet++) { // loop over level sets...
                    var surfaceFactory = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant).GetSurfaceFactory(iLevSet, Kref);
                    LevSetQrIns = LevSetQrIns.AddFactory(surfaceFactory, this.lsTrk._Regions.GetCutCellMask4LevSet(iLevSet));
                }
            }

            return LevSetQrIns;
        }


        public EdgeQuadratureScheme GetEdgeQuadScheme(SpeciesId sp, bool UseDefaultFactories = true, EdgeMask IntegrationDomain = null, int? fixedOrder = null) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");


            // determine domain
            // ================
            var allRelevantEdges = GetEdgeMask(sp, IntegrationDomain);

            // create quadrature scheme
            // ========================
            {
                // default rules for all edges:
                EdgeQuadratureScheme edgeQrIns = new EdgeQuadratureScheme(UseDefaultFactories, allRelevantEdges);

                // overwrite with cut-cell-rules in cut-cells:
                foreach (var Kref in lsTrk.GridDat.Grid.RefElements) {
                    for (int iLevSet = 0; iLevSet < lsTrk.LevelSets.Count; iLevSet++) { // loop over level sets...
                        EdgeMask cutEdges = this.GetCutEdges(Kref, iLevSet).Intersect(allRelevantEdges);
#if DEBUG
                        CellMask difference = cutEdges.GetAdjacentCells(lsTrk.GridDat).Except(lsTrk._Regions.GetCutCellMask4LevSet(iLevSet));
                        if (difference.Count() > 0)
                            throw new ArithmeticException("Edges of the Cells" + difference.GetSummary() + " are detected as cut, but these cells are not contained in the cut Cell-Mask of the Level-Set-Tracker");
#endif


                        var jmp = IdentifyWing(iLevSet, sp);
                        var factory = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant).GetEdgeRuleFactory(iLevSet, jmp, Kref);
                        edgeQrIns.AddFactoryDomainPair(factory, cutEdges, fixedOrder);
                    }
                }

                return edgeQrIns;
            }
        }

        /// <summary>
        /// For some species <paramref name="sp"/>, 
        /// this function computes on which side/wing
        /// of level-set no. <paramref name="levSetIdx"/>)
        /// the species is located.
        /// </summary>
        /// <remarks>
        /// Nur ein Provisorium, das ganze Konzept ist noch etwas unausgereift. (Habe einen Nachmittag lang darueber nachgedacht, keinen bessere Idee gehabt, 
        /// und darum...).
        /// </remarks>
        public JumpTypes IdentifyWing(int levSetIdx, SpeciesId sp) {
            int NoOfLevSets = lsTrk.LevelSets.Count;

            if (levSetIdx < 0 || levSetIdx > NoOfLevSets)
                throw new ArgumentOutOfRangeException();

            if (NoOfLevSets == 1) {
                string[] speciesTable = (string[])lsTrk.SpeciesTable;
                Debug.Assert(speciesTable.Length == 2);

                string spN = lsTrk.GetSpeciesName(sp);

                if (spN == speciesTable[0])
                    return JumpTypes.OneMinusHeaviside;
                else if (spN == speciesTable[1])
                    return JumpTypes.Heaviside;
                else
                    throw new Exception("should not happen.");

            } else if (NoOfLevSets == 2) {
                string[,] speciesTable = (string[,])lsTrk.SpeciesTable;
                Debug.Assert(speciesTable.GetLength(0) == 2);
                Debug.Assert(speciesTable.GetLength(1) == 2);

                string spN = lsTrk.GetSpeciesName(sp);

                int cnt = 0;
                JumpTypes jmpRet = JumpTypes.Implicit;
                int[] _i = new int[2];
                for (_i[0] = 0; _i[0] < 2; _i[0]++) { // loop over signs of level-set 0 ...
                    for (_i[1] = 0; _i[1] < 2; _i[1]++) { // loop over signs of level-set 1 ...
                        if (speciesTable[_i[0], _i[1]] == spN) {
                            cnt++;

                            if (_i[levSetIdx] == 0)
                                jmpRet = JumpTypes.OneMinusHeaviside;
                            else if (_i[levSetIdx] == 1)
                                jmpRet = JumpTypes.Heaviside;
                            else
                                throw new ApplicationException();
                        }
                    }
                }

                if (cnt != 1)
                    throw new NotImplementedException("unable to identify.");


                return jmpRet;
            } else {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Quadrature scheme which is used for the penalty components in <see cref="XSpatialOperator.GhostEdgesOperator"/>, for species <paramref name="sp"/>.
        /// </summary>
        public EdgeQuadratureScheme GetEdgeGhostScheme(SpeciesId sp, EdgeMask IntegrationDomainRestriction = null) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");

            // select all edges for species sp with a positive measure
            // =======================================================
            var EdgBitMask = this.GetEdgeMask(sp, IntegrationDomainRestriction).GetBitMask().CloneAs();
            var EdgArea = this.NonAgglomeratedMetrics.CutEdgeAreas[sp];

            int L = EdgBitMask.Length;
            for (int e = 0; e < L; e++) {
                if (EdgArea[e] <= 0)
                    EdgBitMask[e] = false;
            }

            // return
            // ======
            return new EdgeQuadratureScheme(true, new EdgeMask(this.lsTrk.GridDat, EdgBitMask));


            /*

            if (!this.GhostSupport)
                throw new NotSupportedException();


            // determine domain
            // ================

            var allRelevantEdges = this.GhostIntegrationDomain[sp];
            if (IntegrationDomain != null) {
                allRelevantEdges = allRelevantEdges.Intersect(IntegrationDomain);
            }

            // special treatment for the edges of agglomerated cells
            // =====================================================
            EdgeMask fuckedEdges = null;
            {

                var AggCells = this.CellAgglomeration.GetAgglomerator(sp).AggInfo.SourceCells;
                if (AggCells != null && AggCells.NoOfItemsLocally > 0) {
                    //int EdgB4 = allRelevantEdges.NoOfItemsGlobally;

                    //allRelevantEdges = allRelevantEdges.Except(AggCellsSgrd.InnerEdgesMask); // option 1
                    //allRelevantEdges = allRelevantEdges.Except(AggCellsSgrd.AllEdgesMask);   // option 2
                    //                                                                          option 3: überhaupt nix

                    fuckedEdges = allRelevantEdges.Intersect(this.CellAgglomeration.GetAgglomerator(sp).AggInfo.SourceCellsEdges);
                }
            }

            // create Scheme
            // =============
            {
                var GhostScheme = new EdgeQuadratureScheme(true, allRelevantEdges);

                // use the HMF-rule on the edges of the agglomerated cells:
                if (fuckedEdges != null) {
                    foreach (var Kref in lsTrk.GridDat.Grid.RefElements) {
                        for (int iLevSet = 0; iLevSet < lsTrk.LevelSets.Count; iLevSet++) { // loop over level sets...
                            var cutEdges = fuckedEdges;
                            cutEdges = cutEdges.Intersect(this.m_Subgrid4Kref_AllEdges[Kref]);
                            var jmp = IdentifyWing(iLevSet, sp);
                            var factory = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant).GetEdgeRuleFactory(iLevSet, jmp, Kref);
                            GhostScheme.AddFactory(factory, cutEdges);
                        }
                    }
                }

                // return
                return GhostScheme;
            }
            */
        }

        /// <summary>
        /// All edges that have to be considered for the integration of species <paramref name="sp"/>.
        /// </summary>
        /// <param name="sp"></param>
        /// <param name="IntegrationDomainRestriction">Optional restriction to the integration domain.</param>
        public EdgeMask GetEdgeMask(SpeciesId sp, EdgeMask IntegrationDomainRestriction = null) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");

            EdgeMask allRelevantEdges;
            if (IntegrationDomainRestriction == null) {
                allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp];
            } else {
                // user provides integration domain, intersect with that
                allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].Intersect(IntegrationDomainRestriction);
            }

            // optionally, exclude edges between agglomerated cells
            EdgeMask AggEdges = (this.CellAgglomeration != null) ? (this.CellAgglomeration.GetAgglomerator(sp).AggInfo.AgglomerationEdges) : (default(EdgeMask));
            if (AggEdges != null && AggEdges.NoOfItemsLocally > 0)
                allRelevantEdges = allRelevantEdges.Except(AggEdges);
            return allRelevantEdges;
        }


        /// <summary>
        /// Fills up the volume scheme for species <paramref name="sp"/>.
        /// </summary>
        public CellQuadratureScheme GetVolumeQuadScheme(SpeciesId sp, bool UseDefaultFactories = true, CellMask IntegrationDomain = null, int? fixedOrder = null) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");

            CellMask CellMask = GetCellMask(sp, IntegrationDomain);

            /// Debugging Code
            //if (IntegrationDomain != null) {
            //    CellMask.ToTxtFile("VolDom-" + this.lsTrk.GridDat.MyRank + "of" + this.lsTrk.GridDat.Size + ".csv", false);
            //}

            // default rule for "normal" cells
            var volQrIns = (new CellQuadratureScheme(UseDefaultFactories, CellMask));

            // now: rules for the cut-cells:
            for (int iLevSet = 0; iLevSet < lsTrk.LevelSets.Count; iLevSet++) { // loop over level sets
                var cutDom = lsTrk._Regions.GetCutCellMask4LevSet(iLevSet);
                var cutCells = cutDom.Intersect(CellMask);

                var jmp = IdentifyWing(iLevSet, sp);

                for (int iKref = 0; iKref < lsTrk.GridDat.Grid.RefElements.Length; iKref++) {
                    RefElement Kref = lsTrk.GridDat.Grid.RefElements[iKref];
                    var _cutDom = cutCells.Intersect(lsTrk.GridDat.Cells.GetCells4Refelement(iKref));
                    var factory = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant).GetVolRuleFactory(iLevSet, jmp, Kref);
                    volQrIns.AddFactoryDomainPair(factory, _cutDom, fixedOrder);
                }
            }

            return volQrIns;
        }

        public CellMask GetCellMask(SpeciesId sp, CellMask IntegrationDomain) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species '" + this.lsTrk.GetSpeciesName(sp) + "' is not supported.");


            CellMask OutCellMask;

            if (IntegrationDomain == null) {
                OutCellMask = lsTrk._Regions.GetSpeciesMask(sp);
            } else {
                OutCellMask = lsTrk._Regions.GetSpeciesMask(sp).Intersect(IntegrationDomain);
            }
            return OutCellMask;
        }

        public CellQuadratureScheme GetLevelSetquadScheme(int iLevSet, CellMask IntegrationDom, int? fixedOrder = null) {

            CellQuadratureScheme LevSetQrIns = new CellQuadratureScheme(false, IntegrationDom);

            foreach (var Kref in lsTrk.GridDat.Grid.RefElements) {
                var surfaceFactory = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant).GetSurfaceFactory(iLevSet, Kref);
                LevSetQrIns.AddFactoryDomainPair(surfaceFactory, (CellMask)null, fixedOrder);
            }

            return LevSetQrIns;
        }



        

        /*
        
        static int Jmp2Idx(JumpTypes jmp) {
            switch (jmp) {
                case JumpTypes.OneMinusHeaviside:
                    return 0;
                case JumpTypes.Heaviside:
                    return 1;
                default:
                    throw new NotSupportedException();
            }
        }

        static internal double Jmp2Sign(JumpTypes jmp) {
            return (double)Math.Sign(((double)Jmp2Idx(jmp)) - 0.5);
        }
        */

        /// <summary>
        /// Writes diagnostic information about quadrature rules into csv-textfiles.
        /// </summary>
        public void RuleInfo(string spc, string volruleName, string levsetRuleName, string cellBndRule, string edgeRuleName, int order) {
            RuleInfo(lsTrk.GetSpeciesId(spc), volruleName, levsetRuleName, cellBndRule, edgeRuleName, order);
        }

        /// <summary>
        /// Writes diagnostic information about quadrature rules into csv-textfiles.
        /// </summary>
        public void RuleInfo(SpeciesId spc, string volruleName, string levsetRuleName, string cellBndRule, string edgeRuleName, int order) {

            var _Context = this.lsTrk.GridDat;
            var LsTrk = this.lsTrk;

            // test parameters
            var jmp = Foundation.XDG.Quadrature.HMF.JumpTypes.Heaviside;
            var sch = this.lsTrk.GetXQuadFactoryHelper(MomentFittingVariant);
            var spNm = this.lsTrk.GetSpeciesName(spc);
            var spId = LsTrk.GetSpeciesId(spNm);
            var DomainOfInterest = LsTrk._Regions.GetCutCellSubgrid4LevSet(0);


            if (volruleName != null) {
                int J = _Context.Cells.NoOfLocalUpdatedCells;
                double[] TotWeights = new double[J];

                foreach (var Kref in this.lsTrk.GridDat.Grid.RefElements) {
                    var volfactory = sch.GetVolRuleFactory(0, jmp, Kref);
                    var qrset = volfactory.GetQuadRuleSet(DomainOfInterest.VolumeMask, order);
                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            TotWeights[j] = chk.Rule.Weights.Sum();
                        }
                    }
                }

                DomainOfInterest.VolumeMask.SaveToTextFile(volruleName + "-" + spNm + ".csv", false,
                    delegate(double[] x, int jL, int iG) {
                        return TotWeights[jL];
                    });
            }


            if (levsetRuleName != null) {
                int J = _Context.Cells.NoOfLocalUpdatedCells;
                double[] TotWeights = new double[J];

                foreach (var Kref in this.lsTrk.GridDat.Grid.RefElements) {
                    var levsetFactory = sch.GetSurfaceFactory(0, Kref);
                    var qrset = levsetFactory.GetQuadRuleSet(DomainOfInterest.VolumeMask, order);

                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            TotWeights[j] = chk.Rule.Weights.Sum();
                        }
                    }
                }

                DomainOfInterest.VolumeMask.SaveToTextFile(levsetRuleName + "-" + spNm + ".csv", false,
                    delegate(double[] x, int jL, int iG) {
                        return TotWeights[jL];
                    });

            }



            {

                int NoEdg = this.lsTrk.GridDat.Edges.Count;
                double[] TotWeights = new double[NoEdg];

                foreach (var Kref in this.lsTrk.GridDat.Grid.RefElements) {
                    var edgeFactory = sch.GetEdgeRuleFactory(0, jmp, Kref);
                    var qrset = edgeFactory.GetQuadRuleSet(DomainOfInterest.AllEdgesMask, order);

                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            TotWeights[j] = chk.Rule.Weights.Sum();
                        }
                    }
                }

                DomainOfInterest.AllEdgesMask.SaveToTextFile(edgeRuleName + "-" + spNm + ".csv", false,
                    delegate(double[] x, int jL, int iG) {
                        return TotWeights[jL];
                    });
            } // */

            /*
            if (cellBndRule != null) {
                int J = _Context.Cells.NoOfLocalUpdatedCells;
                var FaceWeightSum = MultidimensionalArray.Create(J, this.lsTrk.GridDat.Grid.RefElements.Max(Kref => Kref.NoOfFaces));

                foreach (var Kref in this.lsTrk.GridDat.Grid.RefElements) {
                    var ruleConv = sch.GetVolumeEdgeRuleFactory(0, Kref);
                    var qrset = ruleConv.GetQuadRuleSet(DomainOfInterest.VolumeMask, order);


                    int F = Kref.NoOfFaces;


                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            var Rule = chk.Rule;
                            var NodesPerFace = Rule.NumbersOfNodesPerFace;

                            int N0 = 0;

                            for (int f = 0; f < F; f++) {
                                for (int n = N0; n < N0 + NodesPerFace[f]; n++) {
                                    FaceWeightSum[j, f] += Rule.Weights[n];
                                }
                                N0 += NodesPerFace[f];
                            }
                        }
                    }
                }

                //DomainOfInterest.VolumeMask.ToTxtFile("cellBndRule-" + jmp + ".csv", false,
                //    delegate(double[] x, int j) { return FaceWeightSum[j, 0]; },
                //    delegate(double[] x, int j) { return FaceWeightSum[j, 1]; },
                //    delegate(double[] x, int j) { return FaceWeightSum[j, 2]; },
                //    delegate(double[] x, int j) { return FaceWeightSum[j, 3]; });

                ToTxtFile(cellBndRule + ".csv", FaceWeightSum, DomainOfInterest.VolumeMask);
            }*/
        }


        /// <summary>
        /// writes the center coordinates of all cells in mask <paramref name="cm"/> to some text file
        /// </summary>
        static void ToTxtFile(GridData _Context, string fileName, MultidimensionalArray wgtSum, CellMask cm) {
            int D = _Context.SpatialDimension;
            using (var file = new StreamWriter(fileName)) {

                double[] x = new double[D];


                foreach (Chunk chunk in cm) {
                    for (int i = 0; i < chunk.Len; i++) {
                        int jCell = chunk.i0 + i;
                        var Kref = _Context.Cells.GetRefElement(jCell);
                        int NoOfFaces = Kref.NoOfFaces;

                        MultidimensionalArray globalCenters = MultidimensionalArray.Create(chunk.Len, NoOfFaces, D);
                        _Context.TransformLocal2Global(Kref.FaceCenters, jCell, 1, globalCenters, 0);

                        for (int idxF = 0; idxF < NoOfFaces; idxF++) {

                            file.Write(chunk.i0 + i);
                            for (int d = 0; d < D; d++) {
                                file.Write("\t" + globalCenters[i, idxF, d].ToString("e", NumberFormatInfo.InvariantInfo));
                            }

                            file.Write("\t" + wgtSum[jCell, idxF].ToString("e", NumberFormatInfo.InvariantInfo));

                            file.WriteLine();
                        }


                    }
                }
            }

        }


       

    }
}


