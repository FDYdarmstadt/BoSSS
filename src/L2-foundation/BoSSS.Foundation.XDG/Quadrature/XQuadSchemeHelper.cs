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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Provides quadrature schemes for typical XDG-cases;
    /// </summary>
    public class XQuadSchemeHelper {

        /// <summary>
        ///
        /// </summary>
        public CutCellMetrics NonAgglomeratedMetrics {
            get {
                return XDGSpaceMetrics.CutCellMetrics;
            }
        }

        /// <summary>
        /// Selected variant of the moment-fitting procedure
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant {
            get {
                return XDGSpaceMetrics.CutCellQuadratureType;
            }
        }

        /// <summary>
        /// If true, <see cref="GetEdgeGhostScheme(SpeciesId, EdgeMask)"/> is supported, otherwise not.
        /// </summary>
        /// <remarks>
        /// Currently (april2016), nobody is using this feature,
        /// therefore this is hard-coded to false.
        /// </remarks>
        public bool GhostSupport {
            get {
                return false;
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
        /// owner object.
        /// </summary>
        public XDGSpaceMetrics XDGSpaceMetrics {
            get;
            private set;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        internal XQuadSchemeHelper(XDGSpaceMetrics __XDGSpaceMetrics) {
            MPICollectiveWatchDog.Watch();
            this.XDGSpaceMetrics = __XDGSpaceMetrics;
            ConstructorCommon();
        }

        private void ConstructorCommon() {
            var Krefs = this.XDGSpaceMetrics.GridDat.Grid.RefElements;

            // initialize some mask's and subgrids
            // ===================================

            // since most methods of this class are non-collective,
            // an access to some Subgrid-member may cause an MPI-deadlock.
            // (this becomes even more unpredictable due to the on-demand-pattern in which most members of the Subgrid-class are implemented);

            //
            foreach (var spId in this.SpeciesList) {
                SubGrid spSgrd = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesSubGrid(spId);

                this.m_SpeciesSubgrid_InnerAndDomainEdges.Add(
                    spId,
                    spSgrd.InnerEdgesMask.Union(spSgrd.AllEdgesMask.Intersect(XDGSpaceMetrics.GridDat.BoundaryEdges)));
            }

            // all cut edges
            // =============
            m_CutCellSubgrid_InnerEdges = this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubGrid().InnerEdgesMask;

            // cut edges sorted according to reference element and level-set
            // =============================================================
            int NoOfLs = this.XDGSpaceMetrics.NoOfLevelSets;
            this.m_CutEdges = new EdgeMask[Krefs.Length, NoOfLs];
            this.m_HMFEdgesDomain = new EdgeMask[Krefs.Length, NoOfLs];
            for (int iKref = 0; iKref < this.m_CutEdges.GetLength(0); iKref++) {
                for (int iLevSet = 0; iLevSet < this.m_CutEdges.GetLength(1); iLevSet++) {
                    EdgeMask cutEdges = this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet).InnerEdgesMask.Union(
                            this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet).AllEdgesMask.Intersect(this.XDGSpaceMetrics.GridDat.BoundaryEdges));

                    cutEdges = cutEdges.Intersect(this.XDGSpaceMetrics.GridDat.GetRefElementSubGrid(iKref).AllEdgesMask);
                    cutEdges = cutEdges.ToGeometicalMask();

                    Debug.Assert(cutEdges.MaskType == MaskType.Geometrical);
                    this.m_CutEdges[iKref, iLevSet] = cutEdges;

                    // (all edges of 'Kref'-elements) \cap (all edges of cells cut by 'iLevSet')
                    this.m_HMFEdgesDomain[iKref, iLevSet] = this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet).AllEdgesMask.Intersect(this.XDGSpaceMetrics.GridDat.GetRefElementSubGrid(iKref).AllEdgesMask);
                }
            }
        }

        /// <summary>
        /// All edges (for each reference element, for each level-set) for which the HMF must be used.
        /// 1st index: volume reference element
        /// 2nd index: level set
        /// </summary>
        private EdgeMask[,] m_HMFEdgesDomain;

        /// <summary>
        /// All edges which are cut by a level-set.
        /// 1st index: volume reference element
        /// 2nd index: level set
        /// </summary>
        /// <remarks>
        /// Be aware that this might contain also edges like <c>e</c>,
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
        private EdgeMask[,] m_CutEdges;

        /// <summary>
        /// all edges which are cut by level set #<paramref name="iLevSet"/> and which belong to a cell with
        /// reference element <paramref name="Kref"/>.
        /// </summary>
        private EdgeMask GetCutEdges(RefElement Kref, int iLevSet) {
            int iKref = this.XDGSpaceMetrics.GridDat.Grid.RefElements.IndexOf(Kref);
            return m_CutEdges[iKref, iLevSet];
        }

        /// <summary>
        /// all edges of cells which are cut by level set #<paramref name="iLevSet"/> and which belong to a cell with
        /// reference element <paramref name="Kref"/>.
        /// </summary>
        private EdgeMask GetHMFEdgesDomain(RefElement Kref, int iLevSet) {
            int iKref = this.XDGSpaceMetrics.GridDat.Grid.RefElements.IndexOf(Kref, (a, b) => object.ReferenceEquals(a, b));
            return m_HMFEdgesDomain[iKref, iLevSet];
        }

        /// <summary>
        /// initialized by the constructor to avoid MPI-deadlocks;
        /// keys: species 'S' <br/>
        /// value: an edge-mask containing all edges that are at least partly covered by species 'S'
        /// </summary>
        private Dictionary<SpeciesId, EdgeMask> m_SpeciesSubgrid_InnerAndDomainEdges = new Dictionary<SpeciesId, EdgeMask>();

        /// <summary>
        /// initialized by the constructor to avoid MPI-deadlocks;
        /// </summary>
        private EdgeMask m_CutCellSubgrid_InnerEdges;

        public EdgeQuadratureScheme Get_SurfaceElement_EdgeQuadScheme(SpeciesId sp, int iLevSet) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species (id = " + sp.cntnt + ") is not supported.");
            //Default behaviour: If Species are not divided by Level Set, function should not be called
            Debug.Assert(!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp));
            //var allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].Intersect(this.m_CutCellSubgrid_InnerEdges);

            var innerCutCellEdges = this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet).InnerEdgesMask;
            var boundaryCutCellEdges = ExecutionMask.Intersect(this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubGrid().BoundaryEdgesMask, this.XDGSpaceMetrics.GridDat.BoundaryEdges);
            var allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].Intersect(ExecutionMask.Union(innerCutCellEdges, boundaryCutCellEdges));
            allRelevantEdges = allRelevantEdges.ToGeometicalMask();
            //EdgeMask AggEdges = this.CellAgglomeration != null ? this.CellAgglomeration.GetAgglomerator(sp).AggInfo.AgglomerationEdges : null;
            //if (AggEdges != null && AggEdges.NoOfItemsLocally > 0)
            //    allRelevantEdges = allRelevantEdges.Except(AggEdges);

            var edgeQrIns = new EdgeQuadratureScheme(false, allRelevantEdges);

            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                //for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                {
                    EdgeMask cutEdges = this.GetCutEdges(Kref, iLevSet);

                    var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, Kref);

                    edgeQrIns.AddFactory(factory, cutEdges);
                }
            }
            //Handle doubly cut cells
            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                    if (iLevSet != jLevSet) {
                        if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                            CellMask doublyCutCells = this.GetCutCells(iLevSet, jLevSet);
                            EdgeMask dCCEdges = doublyCutCells.AllEdges();
                            var doublyCut = dCCEdges.Intersect(allRelevantEdges);
                            if (doublyCut.Count() > 0) {
                                var jmpJ = IdentifyWingA(jLevSet, sp);
                                var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, Kref);
                                var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, jLevSet, jmpJ, Kref, backupFactory);
                                edgeQrIns.AddFactory(factory, doublyCut);
                            }
                        }
                    }
                }
            }
            return edgeQrIns;
        }

        /// <summary>
        /// Interior quadrature for the surface elements, i.e. for each cut background-cell \f$ K_j \f$ a quadrature to approximate
        /// \f[
        ///    \oint_{K_j \cap \mathfrak{I} } \ldots \mathrm{dS} .
        /// \f]
        /// </summary>
        public CellQuadratureScheme Get_SurfaceElement_VolumeQuadScheme(SpeciesId sp, int iLevSet) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species (id = " + sp.cntnt + ") is not supported.");
            //Default behaviour: If Species are not divided by Level Set, function should not be called
            Debug.Assert(!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp));

            var spdom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(sp);
            var IntegrationDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).Intersect(spdom);

            var LevSetQrIns = new CellQuadratureScheme(false, IntegrationDom);

            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                //for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                {
                    var surfaceFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);
                    LevSetQrIns = LevSetQrIns.AddFactory(surfaceFactory, XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).ToGeometicalMask());
                }
            }
            //Handle doubly cut cells
            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                    if (iLevSet != jLevSet) {
                        if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                            CellMask doublyCut = this.GetCutCells(iLevSet, jLevSet);
                            if (doublyCut.Count() > 0) {
                                var jmpJ = IdentifyWingA(jLevSet, sp);
                                var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);
                                var surfaceFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, jLevSet, jmpJ, Kref, backupFactory);
                                LevSetQrIns.AddFactory(surfaceFactory, doublyCut);
                            }
                        }
                    }
                }
            }
            return LevSetQrIns;
        }

        public CellQuadratureScheme GetContactLineQuadScheme(SpeciesId sp, int iLevSet) {
            //Find domain
            CellMask allDoublyCuts = CellMask.GetEmptyMask(XDGSpaceMetrics.GridDat, MaskType.Geometrical);

            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                    if (iLevSet != jLevSet) {
                        if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                            allDoublyCuts = allDoublyCuts.Union(GetCutCells(iLevSet, jLevSet));
                        }
                    }
                }
            }

            var spdom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(sp).ToGeometicalMask();
            var IntegrationDom = allDoublyCuts.Intersect(spdom);
            var LevSetQrIns = new CellQuadratureScheme(false, IntegrationDom);

            //Handle doubly cut cells, do it all again, this time add quadrature Factory
            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                    if (iLevSet != jLevSet) {
                        if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                            CellMask doublyCut = this.GetCutCells(iLevSet, jLevSet);
                            if (doublyCut.Count() > 0) {
                                var jmpJ = IdentifyWingA(jLevSet, sp);
                                var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, Kref);
                                var surfaceFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetIntersectionRuleFactory(iLevSet, jLevSet, Kref, backupFactory);
                                LevSetQrIns.AddFactory(surfaceFactory, doublyCut);
                            }
                        }
                    }
                }
            }
            return LevSetQrIns;
        }

        /// <summary>
        /// Boundary quadrature for the surface elements, i.e. for each cut background-cell \f$ K_j \f$ a quadrature to approximate
        /// \f[
        ///    \int_{\partial K_j \cap \mathfrak{I} } \ldots \mathrm{dS} .
        /// \f]
        /// </summary>
        public EdgeQuadratureScheme GetEdgeQuadScheme(SpeciesId sp, bool UseDefaultFactories = true, EdgeMask IntegrationDomain = null, int? fixedOrder = null) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species (id = " + sp.cntnt + ") is not supported.");

            // determine domain
            // ================
            var allRelevantEdges = GetEdgeMask(sp, IntegrationDomain);

            // create quadrature scheme
            // ========================
            {
                // default rules for all edges:
                EdgeQuadratureScheme edgeQrIns = new EdgeQuadratureScheme(UseDefaultFactories, allRelevantEdges);

                // overwrite with cut-cell-rules in cut-cells:
                foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                    for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                        if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {
                            EdgeMask cutEdges = this.GetCutEdges(Kref, iLevSet).Intersect(allRelevantEdges);
#if DEBUG
                            CellMask difference = cutEdges.GetAdjacentCells().Except(XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet));
                            if (difference.Count() > 0)
                                throw new ArithmeticException("Edges of the Cells" + difference.GetSummary() + " are detected as cut, but these cells are not contained in the cut Cell-Mask of the Level-Set-Tracker");
#endif
                            var jmp = IdentifyWingA(iLevSet, sp);
                            var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(iLevSet, jmp, Kref);
                            edgeQrIns.AddFactoryDomainPair(factory, cutEdges, fixedOrder);
                        }
                    }
                }

                // overwrite with double-cut-cell-rules in double-cut-cells:
                foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                    for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                        if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {
                            EdgeMask cutEdges = this.GetCutEdges(Kref, iLevSet).Intersect(allRelevantEdges);
                            var jmp = IdentifyWingA(iLevSet, sp);
                            //handle rules for cells/edges where two levelsets are present
                            for (int jLevSet = iLevSet + 1; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                                if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                                    EdgeMask doublyCut = cutEdges.Intersect(GetCutEdges(Kref, jLevSet));
                                    if (doublyCut.Count() > 0) {
                                        var jmpJ = IdentifyWingA(jLevSet, sp);
                                        var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(iLevSet, jmp, Kref);
                                        var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(iLevSet, jmp, jLevSet, jmpJ, Kref, backupFactory);
                                        edgeQrIns.AddFactoryDomainPair(twoLSFactory, doublyCut, fixedOrder);
                                    }
                                }
                            }
                        }
                    }
                }

                return edgeQrIns;
            }
        }

        private CellMask GetCutCells(int iLevSet, int jLevSet) {
            CellMask iCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
            CellMask jCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(jLevSet);
            return iCells.Intersect(jCells).ToGeometicalMask();
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
            int NoOfLevSets = this.XDGSpaceMetrics.NoOfLevelSets;

            if (levSetIdx < 0 || levSetIdx > NoOfLevSets)
                throw new ArgumentOutOfRangeException();

            if (NoOfLevSets == 1) {
                string[] speciesTable = (string[])(this.XDGSpaceMetrics.LevelSetRegions.SpeciesTable);
                Debug.Assert(speciesTable.Length == 2);

                string spN = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(sp);

                if (spN == speciesTable[0])
                    return JumpTypes.OneMinusHeaviside;
                else if (spN == speciesTable[1])
                    return JumpTypes.Heaviside;
                else
                    throw new Exception("should not happen.");
            } else if (NoOfLevSets == 2) {
                string[,] speciesTable = (string[,])(this.XDGSpaceMetrics.LevelSetRegions.SpeciesTable);
                Debug.Assert(speciesTable.GetLength(0) == 2);
                Debug.Assert(speciesTable.GetLength(1) == 2);

                string spN = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(sp);

                int[] LevSetSigns;
                bool foundCell;
                {
                    // we need the signs of other level sets to identify the wing
                    // so we search for some cut cell of Level-Set #levSetIdx
                    // which actually contains species 'sp'

                    int J = this.XDGSpaceMetrics.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                    int[] LenToNextchange = this.XDGSpaceMetrics.LevelSetRegions.m_LenToNextChange;
                    LevSetSigns = new int[NoOfLevSets];
                    foundCell = false;

                    for (int j = 0; j < J; j += LenToNextchange[j]) {
                        ushort code = this.XDGSpaceMetrics.LevelSetRegions.m_LevSetRegions[j];

                        int dist = LevelSetTracker.DecodeLevelSetDist(code, levSetIdx);
                        if (dist != 0)
                            continue;
                        // cut by Level-Set

                        bool present = this.XDGSpaceMetrics.LevelSetRegions.IsSpeciesPresentInCell(sp, j);
                        if (!present)
                            continue;
                        // contains species 'sp'

                        var Signs = this.XDGSpaceMetrics.LevelSetRegions.GetCellSignCode(j);
                        for (int iLs = 0; iLs < NoOfLevSets; iLs++) {
                            if (iLs != levSetIdx) {
                                var s = Signs.GetSign(iLs);

                                if (s == LevelsetSign.Both) {
                                    continue;
                                } else if (s == LevelsetSign.Negative) {
                                    foundCell = true;
                                    LevSetSigns[iLs] = 0;
                                } else if (s == LevelsetSign.Positive) {
                                    foundCell = true;
                                    LevSetSigns[iLs] = 1;
                                } else {
                                    throw new NotImplementedException();
                                }
                            }
                        }
                    }
                }
                int cnt = 0;
                JumpTypes jmpRet = JumpTypes.Implicit;

                /*
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
                */

                if (foundCell == false)
                    return JumpTypes.Heaviside; // no cell on this proc, so anyway pretty irrelevant

                for (int i = 0; i < 2; i++) {
                    LevSetSigns[levSetIdx] = i;

                    if (speciesTable[LevSetSigns[0], LevSetSigns[1]] == spN) {
                        cnt++;
                        if (LevSetSigns[levSetIdx] == 0)
                            jmpRet = JumpTypes.OneMinusHeaviside;
                        else if (LevSetSigns[levSetIdx] == 1)
                            jmpRet = JumpTypes.Heaviside;
                        else
                            throw new ApplicationException();
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
        /// Returns the Jumptype of the species.
        /// Handle with care! Does this really work?
        /// </summary>
        /// <param name="levSetIdx"></param>
        /// <param name="sp"></param>
        /// <returns></returns>
        public JumpTypes IdentifyWingA(int levSetIdx, SpeciesId sp) {
            if (this.XDGSpaceMetrics.NoOfLevelSets != 2) {
                return IdentifyWing(levSetIdx, sp);
            }
            //JumpTypes oho = IdentifyWing(levSetIdx, sp);

            string[,] speciesTable = (string[,])(this.XDGSpaceMetrics.LevelSetRegions.SpeciesTable);
            string spName = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(sp);

            //Find sp indices, only 2LS
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    if (speciesTable[i, j] == spName) {
                        int[] indice = new[] { i, j };
                        if (indice[levSetIdx] == 0) {
                            //if (JumpTypes.OneMinusHeaviside != oho) {
                            //    throw new Exception("oh noes");
                            //}
                            return JumpTypes.OneMinusHeaviside;
                        } else {
                            //if (JumpTypes.Heaviside != oho) {
                            //    throw new Exception("oh noes");
                            //}
                            return JumpTypes.Heaviside;
                        }
                    }
                }
            }
            throw new Exception("Species not found");
        }

        /// <summary>
        /// Does levSet separate species A and B?
        /// </summary>
        public bool SpeciesAreSeparatedByLevSet(int levSet, SpeciesId A, SpeciesId B) {
            switch (this.XDGSpaceMetrics.NoOfLevelSets) {
                case 1:
                return AreSeparatedByLevSet1LS(levSet, A, B);

                case 2:
                return AreSeparatedByLevSet2LS(levSet, A, B);

                default:
                throw new NotSupportedException();
            }
        }

        private bool AreSeparatedByLevSet1LS(int levSet, SpeciesId A, SpeciesId B) {
            string nameA = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(A);
            string nameB = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(B);
            string[] speciesTable = (string[])this.XDGSpaceMetrics.LevelSetRegions.SpeciesTable;
            if (speciesTable[0] == nameA && speciesTable[1] == nameB) {
                return true;
            } else if (speciesTable[1] == nameA && speciesTable[0] == nameB) {
                return true;
            } else {
                return false;
            }
        }

        private bool AreSeparatedByLevSet2LS(int levSet, SpeciesId A, SpeciesId B) {
            string nameA = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(A);
            string nameB = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(B);

            string[,] speciesTable2LS = (string[,])(this.XDGSpaceMetrics.LevelSetRegions.SpeciesTable);

            int j = levSet;

            if (speciesTable2LS[(0 + j) % 2, 0] == nameA && speciesTable2LS[1, (0 + j) % 2] == nameB
                || speciesTable2LS[(0 + j) % 2, 0] == nameB && speciesTable2LS[1, (0 + j) % 2] == nameA) {
                return true;
            } else if (speciesTable2LS[0, (1 + j) % 2] == nameA && speciesTable2LS[(1 + j) % 2, 1] == nameB
                  || speciesTable2LS[0, (1 + j) % 2] == nameB && speciesTable2LS[(1 + j) % 2, 1] == nameA) {
                return true;
            }
            return false;
        }

        /// <summary>
        /// Quadrature scheme which is used for the penalty components on ghost edges, for species <paramref name="sp"/>.
        /// </summary>
        public EdgeQuadratureScheme GetEdgeGhostScheme(SpeciesId sp, EdgeMask IntegrationDomainRestriction = null) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species (" + sp.cntnt + ") is not supported.");

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
            return new EdgeQuadratureScheme(true, new EdgeMask(this.XDGSpaceMetrics.GridDat, EdgBitMask));

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
                throw new ArgumentException("Given species ( id = " + sp.cntnt + ") is not supported.");

            EdgeMask allRelevantEdges;
            if (IntegrationDomainRestriction == null) {
                allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].ToGeometicalMask();
            } else {
                // user provides integration domain, intersect with that

                if (IntegrationDomainRestriction.MaskType == MaskType.Geometrical)
                    allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].ToGeometicalMask().Intersect(IntegrationDomainRestriction);
                else
                    allRelevantEdges = this.m_SpeciesSubgrid_InnerAndDomainEdges[sp].Intersect(IntegrationDomainRestriction).ToGeometicalMask();
            }

            //// optionally, exclude edges between agglomerated cells
            //EdgeMask AggEdges = (this.CellAgglomeration != null) ? (this.CellAgglomeration.GetAgglomerator(sp).AggInfo.AgglomerationEdges) : (default(EdgeMask));
            //if (AggEdges != null && AggEdges.NoOfItemsLocally > 0)
            //    allRelevantEdges = allRelevantEdges.Except(AggEdges);
            return allRelevantEdges;
        }

        /// <summary>
        /// Fills up the volume scheme for species <paramref name="sp"/>.
        /// </summary>
        public CellQuadratureScheme GetVolumeQuadScheme(SpeciesId sp, bool UseDefaultFactories = true, CellMask IntegrationDomain = null, int? fixedOrder = null) {
                ilPSP.MPICollectiveWatchDog.Watch();
                if (!this.SpeciesList.Contains(sp))
                    throw new ArgumentException("Given species (id = " + sp.cntnt + ") is not supported.");

                if (IntegrationDomain != null)
                    if (IntegrationDomain.MaskType != MaskType.Logical)
                        throw new ArgumentException();

                CellMask CellMask = GetCellMask(sp, IntegrationDomain);

                // Debugging Code
                //if (IntegrationDomain != null) {
                //    CellMask.ToTxtFile("VolDom-" + this.lsTrk.GridDat.MyRank + "of" + this.lsTrk.GridDat.Size + ".csv", false);
                //}

                // default rule for "normal" cells
                var volQrIns = (new CellQuadratureScheme(UseDefaultFactories, CellMask));

                // now: rules for the cut-cells:
                for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets
                    if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {
                        var cutDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).ToGeometicalMask();
                        var cutCells = cutDom.Intersect(CellMask);

                        var jmp = IdentifyWingA(iLevSet, sp);

                        for (int iKref = 0; iKref < XDGSpaceMetrics.GridDat.Grid.RefElements.Length; iKref++) {
                            RefElement Kref = XDGSpaceMetrics.GridDat.Grid.RefElements[iKref];
                            var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetVolRuleFactory(iLevSet, jmp, Kref);
                            var _cutDom = cutCells.Intersect(XDGSpaceMetrics.GridDat.Cells.GetCells4Refelement(iKref));
                            volQrIns.AddFactoryDomainPair(factory, _cutDom, fixedOrder);
                        }
                    }
                }

                //now: rules for the doubly cut-cells
                for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets
                    if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {
                        var cutDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).ToGeometicalMask();
                        var cutCells = cutDom.Intersect(CellMask);

                        var jmp = IdentifyWingA(iLevSet, sp);

                        for (int iKref = 0; iKref < XDGSpaceMetrics.GridDat.Grid.RefElements.Length; iKref++) {
                            RefElement Kref = XDGSpaceMetrics.GridDat.Grid.RefElements[iKref];
                            var _cutDom = cutCells.Intersect(XDGSpaceMetrics.GridDat.Cells.GetCells4Refelement(iKref));

                            //handle rules for cells/edges where two levelsets are present
                            for (int jLevSet = iLevSet + 1; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                                if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                                    CellMask doublyCut = GetCutCells(iLevSet, jLevSet).Intersect(_cutDom);
                                    if (doublyCut.Count() > 0) {
                                        var jmpJ = IdentifyWingA(jLevSet, sp);
                                        var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetVolRuleFactory(iLevSet, jmp, Kref);
                                        var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetVolRuleFactory(iLevSet, jmp, jLevSet, jmpJ, Kref, backupFactory);
                                        volQrIns.AddFactoryDomainPair(twoLSFactory, doublyCut, fixedOrder);
                                    }
                                }
                            }
                        }
                    }
                }
                return volQrIns;
        }

        private CellMask GetCellMask(SpeciesId sp, CellMask IntegrationDomain) {
            if (!this.SpeciesList.Contains(sp))
                throw new ArgumentException("Given species (sp = " + sp.cntnt + ") is not supported.");
            Debug.Assert(IntegrationDomain == null || IntegrationDomain.MaskType == MaskType.Logical);

            CellMask OutCellMask;

            if (IntegrationDomain == null) {
                OutCellMask = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(sp);
                Debug.Assert(OutCellMask.MaskType == MaskType.Logical);
            } else {
                OutCellMask = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(sp).Intersect(IntegrationDomain);
            }

            OutCellMask = OutCellMask.ToGeometicalMask();

            return OutCellMask;
        }

        /// <summary>
        /// Quadrature scheme for the integration over the level-set, i.e. for each cut background-cell \f$ K_j \f$ a quadrature to approximate
        /// \f[
        ///    \oint_{K_j \cap \mathfrak{I} } \ldots \mathrm{dS} .
        /// \f]
        /// </summary>
        /// <param name="iLevSet"></param>
        /// <param name="IntegrationDom"></param>
        /// <param name="fixedOrder"></param>
        /// <returns></returns>
        public CellQuadratureScheme GetLevelSetquadScheme(int iLevSet, CellMask IntegrationDom, int? fixedOrder = null) {
            if (IntegrationDom.MaskType == MaskType.Logical)
                IntegrationDom = IntegrationDom.ToGeometicalMask();

            CellQuadratureScheme LevSetQrIns = new CellQuadratureScheme(false, IntegrationDom);

            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                var surfaceFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);
                LevSetQrIns.AddFactoryDomainPair(surfaceFactory, (CellMask)null, fixedOrder);
            }

            return LevSetQrIns;
        }

        public CellQuadratureScheme GetLevelSetquadScheme(int iLevSet, SpeciesId spA, CellMask IntegrationDom, int? fixedOrder = null) {
            if (IntegrationDom.MaskType == MaskType.Logical)
                IntegrationDom = IntegrationDom.ToGeometicalMask();

            CellQuadratureScheme LevSetQrIns = new CellQuadratureScheme(false, IntegrationDom);
            foreach (var Kref in XDGSpaceMetrics.GridDat.Grid.RefElements) {
                var surfaceFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);

                //handle rules for cells/edges where two levelsets are present
                for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                    if (jLevSet != iLevSet) {
                        CellMask doublyCut = GetCutCells(iLevSet, jLevSet);
                        if (doublyCut.Count() > 0) {
                            IntegrationDom = IntegrationDom.Except(doublyCut);

                            //Ist das so in Ordnung? Scheint keine Invariante zu verletzen.
                            var jmpA = IdentifyWingA(jLevSet, spA);
                            //Debug.Assert(jmpA == IdentifyWing(jLevSet, spB));
                            IQuadRuleFactory<QuadRule> backupFactory;
                            if (iLevSet == 1)
                                backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(jLevSet, jmpA, Kref);
                            else
                                backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);
                            var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, jLevSet, jmpA, Kref, backupFactory);
                            LevSetQrIns.AddFactoryDomainPair(twoLSFactory, doublyCut, fixedOrder);
                        }
                    }
                }

                LevSetQrIns.AddFactoryDomainPair(surfaceFactory, IntegrationDom, fixedOrder);
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

        ///// <summary>
        ///// Writes diagnostic information about quadrature rules into csv-textfiles.
        ///// </summary>
        //public void RuleInfo(string spc, string volruleName, string levsetRuleName, string cellBndRule, string edgeRuleName, int order) {
        //    RuleInfo(lsTrk.GetSpeciesId(spc), volruleName, levsetRuleName, cellBndRule, edgeRuleName, order);
        //}

        /// <summary>
        /// Writes diagnostic information about quadrature rules into csv-textfiles.
        /// </summary>
        public void RuleInfo(SpeciesId spc, string volruleName, string levsetRuleName, string cellBndRule, string edgeRuleName, int order, int iLevSet) {
            var _Context = this.XDGSpaceMetrics.GridDat;

            // test parameters
            var jmp = Foundation.XDG.Quadrature.HMF.JumpTypes.Heaviside;
            var sch = this.XDGSpaceMetrics.XQuadFactoryHelper;
            var spNm = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(spc);
            //var spId = LsTrk.GetSpeciesId(spNm);
            var DomainOfInterest = XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet);

            if (volruleName != null) {
                int J = _Context.Cells.NoOfLocalUpdatedCells;
                double[] TotWeights = new double[J];

                foreach (var Kref in this.XDGSpaceMetrics.GridDat.Grid.RefElements) {
                    var volfactory = sch.GetVolRuleFactory(0, jmp, Kref);
                    var qrset = volfactory.GetQuadRuleSet(DomainOfInterest.VolumeMask, order);
                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            TotWeights[j] = chk.Rule.Weights.Sum();
                        }
                    }
                }

                DomainOfInterest.VolumeMask.SaveToTextFile(volruleName + "-" + spNm + ".csv", false,
                    delegate (double[] x, int jL, int iG) {
                        return TotWeights[jL];
                    });
            }

            if (levsetRuleName != null) {
                int J = _Context.Cells.NoOfLocalUpdatedCells;
                double[] TotWeights = new double[J];

                foreach (var Kref in this.XDGSpaceMetrics.GridDat.Grid.RefElements) {
                    var levsetFactory = sch.GetSurfaceFactory(0, Kref);
                    var qrset = levsetFactory.GetQuadRuleSet(DomainOfInterest.VolumeMask, order);

                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            TotWeights[j] = chk.Rule.Weights.Sum();
                        }
                    }
                }

                DomainOfInterest.VolumeMask.SaveToTextFile(levsetRuleName + "-" + spNm + ".csv", false,
                    delegate (double[] x, int jL, int iG) {
                        return TotWeights[jL];
                    });
            }

            {
                int NoEdg = this.XDGSpaceMetrics.GridDat.Edges.Count;
                double[] TotWeights = new double[NoEdg];

                foreach (var Kref in this.XDGSpaceMetrics.GridDat.Grid.RefElements) {
                    var edgeFactory = sch.GetEdgeRuleFactory(0, jmp, Kref);
                    var qrset = edgeFactory.GetQuadRuleSet(DomainOfInterest.AllEdgesMask, order);

                    foreach (var chk in qrset) {
                        for (int j = chk.Chunk.i0; j < chk.Chunk.JE; j++) {
                            TotWeights[j] = chk.Rule.Weights.Sum();
                        }
                    }
                }

                DomainOfInterest.AllEdgesMask.SaveToTextFile(edgeRuleName + "-" + spNm + ".csv", false,
                    delegate (double[] x, int jL, int iG) {
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
        private static void ToTxtFile(GridData _Context, string fileName, MultidimensionalArray wgtSum, CellMask cm) {
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