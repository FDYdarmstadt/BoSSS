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
using BoSSS.Foundation.XDG.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Provides quadrature schemes 
    /// This is a driver-class for the creation of various quadrature schemes.
    /// All of these have different, involved, constructor calls.
    /// Via this utility class, a somewhat unified access to the creation of these rules is provided.
    /// (see also <see cref="CutCellQuadratureMethod"/>).
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
        public CutCellQuadratureMethod CutCellQuadratureMethod {
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

        LevelSetTracker Tracker => this.XDGSpaceMetrics.Tracker;

        IGridData gdat => this.Tracker.GridDat;

        /// <summary>
        /// ctor.
        /// </summary>
        internal XQuadSchemeHelper(XDGSpaceMetrics __XDGSpaceMetrics) {
            MPICollectiveWatchDog.Watch();
            this.XDGSpaceMetrics = __XDGSpaceMetrics;
            ConstructorCommon();
        }

        /// <summary>
        /// creates all objects which require MPI-collective calls at the beginning, 
        /// when the call happens on all ranks simultaneously.
        /// </summary>
        private void ConstructorCommon() {
            var Krefs = this.gdat.iGeomCells.RefElements;
            var KrefsEdges = this.gdat.iGeomEdges.EdgeRefElements;

            // initialize some mask's and subgrids
            // ===================================

            // since most methods of this class are non-collective,
            // an access to some Subgrid-member may cause an MPI-deadlock.
            // (this becomes even more unpredictable due to the on-demand-pattern in which most members of the Subgrid-class are implemented);

            //
            foreach(var spId in this.SpeciesList) {
                SubGrid spSgrd = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesSubGrid(spId);

                this.m_SpeciesSubgrid_InnerAndDomainEdges.Add(
                    spId,
                    spSgrd.InnerEdgesMask.Union(spSgrd.AllEdgesMask.Intersect(gdat.GetBoundaryEdgeMask())));
            }



            // cut edges sorted according to reference element and level-set
            // =============================================================
            int NoOfLs = this.XDGSpaceMetrics.NoOfLevelSets;

            this.m_CutEdges = new EdgeMask[KrefsEdges.Length, NoOfLs];
            for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                var cutCellSubGrid = this.XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet);
                var innerCut = cutCellSubGrid.InnerEdgesMask;
                var boundaryCut = cutCellSubGrid.AllEdgesMask.Intersect(gdat.GetBoundaryEdgeMask());

                EdgeMask cutEdges = innerCut.Union(boundaryCut).ToGeometicalMask();

                for(int iKref = 0; iKref < KrefsEdges.Length; iKref++) {
                    cutEdges = cutEdges.Intersect(gdat.GetEdges4RefElement(iKref));

                    Debug.Assert(cutEdges.MaskType == MaskType.Geometrical);
                    this.m_CutEdges[iKref, iLevSet] = cutEdges;
                }
            }

            // intersection of species boundary and level-set
            // ==============================================
            SurfaceElementEdges = new Dictionary<(SpeciesId SpeciesA, SpeciesId SpeciesB, int iLevelSet), EdgeMask>();
            LevelSetSeparationLayer = new Dictionary<(SpeciesId SpeciesA, SpeciesId SpeciesB, int iLevelSet), CellMask>();
            foreach(var SpeciesA in this.SpeciesList) { // for species A, we only take the species from our list
                var SpeciesADom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesA);

                foreach(var SpeciesB in this.XDGSpaceMetrics.Tracker.SpeciesIdS) { // for species B, we allow all possible species
                    if(SpeciesB ==  SpeciesA) 
                        continue;
                    for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                        if(SurfaceElementEdges.ContainsKey((SpeciesA, SpeciesB, iLevSet)))
                            continue;

                        var SpeciesBDom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(SpeciesB);
                        var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);

                        CellMask levelSetdomain = this.XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);

                        CellMask LevelSetSeparationLayer_etc = levelSetdomain.Intersect(SpeciesCommonDom);
                        LevelSetSeparationLayer.Add((SpeciesA, SpeciesB, iLevSet), LevelSetSeparationLayer_etc);
                        LevelSetSeparationLayer.Add((SpeciesB, SpeciesA, iLevSet), LevelSetSeparationLayer_etc);

                        EdgeMask edgesInSeparationLayer = InnerAndDomainBoundary(LevelSetSeparationLayer_etc);
                        SurfaceElementEdges.Add((SpeciesA, SpeciesB, iLevSet), edgesInSeparationLayer);
                        SurfaceElementEdges.Add((SpeciesB, SpeciesA, iLevSet), edgesInSeparationLayer);
                        //CellBoundaryLayer_Bitmask.Add((spId, iLevSet), cellBoundaryLayer.GetBitMaskWithExternal());
                    }
                }
            }
        }

        /// <summary>
        /// Note: this is allocates in the constructor, since it requires MPI-communication
        /// </summary>
        Dictionary<(SpeciesId SpeciesA, SpeciesId SpeciesB, int iLevelSet), EdgeMask> SurfaceElementEdges;
        
        /// <summary>
        /// Note: we pre-compute this, since we need it for <see cref="SurfaceElementEdges"/> anyway
        /// </summary>
        Dictionary<(SpeciesId SpeciesA, SpeciesId SpeciesB, int iLevelSet), CellMask> LevelSetSeparationLayer;




        ///// <summary>
        ///// intersection of species boundary and level-set;
        ///// To include also external/ghost cells, this must be initialized in the constructor
        ///// </summary>
        //Dictionary<(SpeciesId sp, int iLevSet), BitArray> CellBoundaryLayer_Bitmask;

        /// <summary>
        /// All edges which are cut by a level-set.
        /// - 1st index: edge reference element, i.e., index into <see cref="IGeometricalEdgeData.EdgeRefElements"/>
        /// - 2nd index: level set
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
        /// reference element <paramref name="KrefEdge"/>.
        /// </summary>
        private EdgeMask GetCutEdges(RefElement KrefEdge, int iLevSet) {
            int iKref = gdat.iGeomEdges.EdgeRefElements.IndexOf(KrefEdge);
            //int iKref = this.gdat.iGeomCells.RefElements.IndexOf(Kref);
            return m_CutEdges[iKref, iLevSet];
        }

        /// <summary>
        /// initialized by the constructor to avoid MPI-deadlocks;
        /// - keys: species 'S' 
        /// - value: an edge-mask containing all edges that are at least partly covered by species 'S'
        /// </summary>
        private Dictionary<SpeciesId, EdgeMask> m_SpeciesSubgrid_InnerAndDomainEdges = new Dictionary<SpeciesId, EdgeMask>();

        /*
        EdgeMask __InnerAndDomainBoundary(CellMask cellRestriction, SpeciesId sp, int iLevSet) {
            //var subGrd = new SubGrid(cells);
            //var innerEdges = subGrd.InnerEdgesMask;
            //var bndy = cells.GridData.GetBoundaryEdgeMask().Intersect(subGrd.BoundaryEdgesMask);
            //return bndy.Union(innerEdges);

            if(cellRestriction != null && cellRestriction.MaskType != MaskType.Logical)
                throw new ArgumentException("expecting a logical cell mask", nameof(cellRestriction));
            
            var bCellRstm = cellRestriction?.GetBitMask();
            var bCellMask = this.CellBoundaryLayer_Bitmask[(sp, iLevSet)];
            
            int[,] E2C = this.gdat.iLogicalEdges.CellIndices;
            int NoOfEdges = this.gdat.iLogicalEdges.Count;
            var bEdgeMask = new BitArray(NoOfEdges);
            int Jup = this.gdat.iLogicalCells.NoOfLocalUpdatedCells;



            for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                int jCell0 = E2C[iEdge, 0];
                int jCell1 = E2C[iEdge, 1];

                if(bCellRstm != null) {
                    if(jCell1 >= 0)
                        if(bCellRstm[jCell0] == false && bCellRstm[jCell1] == false)
                            continue;

                    if(jCell1 < 0 || jCell1 >= Jup)
                        if(bCellRstm[jCell0] == false)
                            continue;
                }


                if(jCell1 >= 0) {
                    // inner egde: both cells must be part of the bit-mask
                    if(bCellMask[jCell0] && bCellMask[jCell1])
                        bEdgeMask[iEdge] = true;
                } else {
                    // boundary edge: only in-cell must be part of the mask
                    if(bCellMask[jCell0])
                        bEdgeMask[iEdge] = true;
                }


            }

            return new EdgeMask(this.gdat, bEdgeMask, mt: MaskType.Logical);
        }
        */

        EdgeMask InnerAndDomainBoundary(CellMask cells) {
            //var subGrd = new SubGrid(cells);
            //var innerEdges = subGrd.InnerEdgesMask;
            //var bndy = cells.GridData.GetBoundaryEdgeMask().Intersect(subGrd.BoundaryEdgesMask);
            //return bndy.Union(innerEdges);

            if(cells.MaskType != MaskType.Logical)
                throw new ArgumentException("expecting a logical cell mask", nameof(cells));



            
            int[,] E2C = cells.GridData.iLogicalEdges.CellIndices;
            var bCellMask = cells.GetBitMaskWithExternal();
            int NoOfEdges = cells.GridData.iLogicalEdges.Count;
            var bEdgeMask = new BitArray(NoOfEdges);



            for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                if(E2C[iEdge, 1] >= 0) {
                    // inner egde: both cells must be part of the bit-mask
                    if(bCellMask[E2C[iEdge, 0]] && bCellMask[E2C[iEdge, 1]])
                        bEdgeMask[iEdge] = true;
                } else {
                    // boundary edge: only in-cell must be part of the mask
                    if(bCellMask[E2C[iEdge, 0]])
                        bEdgeMask[iEdge] = true;
                }


            }

            return new EdgeMask(cells.GridData, bEdgeMask, mt: MaskType.Logical);
        }


        /// <summary>
        /// Edge quadrature for the surface elements, i.e., a $D-2$ dimensional integral.
        /// For each cut background-cell $ K_j $, the surface element is $ K_j \cap \mathfrak{I} $ and its boundaries is  $ \partial K_j \cap \mathfrak{I} $; 
        /// Therefore, one need to integrate over the line integrals:
        /// \[
        ///    \int_{\Gamma \cap \mathfrak{I}_i \cap \partial \mathfrak{A} \cap \partial \mathfrak{B} } \ldots \mathrm{dl} ,
        /// \]
        /// where: 
        /// - $ \Gamma $ is the union of all mesh edges, i.e., $ \Gamma = \bigcup_j \partial K_j $,
        /// - $ i $ is the level set index <paramref name="iLevSet"/>, i.e., $ \mathfrak{I}_i $ is the respective interface, 
        /// - $ \partial \mathfrak{A} $ and $ \cap \partial \mathfrak{B} $ are the boundaries of species <paramref name="spA"/> and <paramref name="spB"/>, respectively.
        /// 
        /// These rules are used, e.g., by the (edge parts of) surface element operator <see cref="XDifferentialOperatorMk2.SurfaceElementOperator_Ls0"/>
        /// </summary>
        public EdgeQuadratureScheme Get_SurfaceElement_EdgeQuadScheme(SpeciesId spA, SpeciesId spB, int iLevSet) {
            var tracker = this.XDGSpaceMetrics.Tracker;
            IGridData gdat = tracker.GridDat;
            var levSetRegions = this.XDGSpaceMetrics.LevelSetRegions;


            //CellMask cellRestriction = null;

            if(!this.SpeciesList.Contains(spA))
                throw new ArgumentException($"Given species (id = {tracker.GetSpeciesName(spA)}) is not supported.");

            EdgeMask allRelevantEdges = SurfaceElementEdges[(spA, spB, iLevSet)];
            /*{
                var allRelevantEdges2 = __InnerAndDomainBoundary(cellRestriction, sp, iLevSet);


                CellMask cellBoundaryLayer = levSetRegions.GetSpeciesMask(sp).Intersect(levSetRegions.GetCutCellMask4LevSet(iLevSet));
                if(cellRestriction != null)
                    cellBoundaryLayer = cellBoundaryLayer.Intersect(cellRestriction);

                allRelevantEdges = InnerAndDomainBoundary(cellBoundaryLayer);


                if(!allRelevantEdges2.Equals(allRelevantEdges)) {
                    allRelevantEdges.SaveToTextFile("edgesCorr.csv");
                    allRelevantEdges2.SaveToTextFile("edgesWrong.csv");
                    throw new ApplicationException("ich hasse dich");
                }
            }*/
            allRelevantEdges = allRelevantEdges.ToGeometicalMask();

            var edgeQrIns = new EdgeQuadratureScheme(
                scaling: new SurfaceElementEdgeIntegrationMetric(this.XDGSpaceMetrics.LevelSetData[iLevSet]),
                UseDefaultFactories: false, 
                domain: allRelevantEdges);

            // Regular case: cells cut by one level-set
            // ========================================

            foreach (var KrefEdge in gdat.iGeomEdges.EdgeRefElements) {
                //for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                {
                    EdgeMask cutEdges = this.GetCutEdges(KrefEdge, iLevSet);
                    var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, KrefEdge);
                    edgeQrIns.AddFactory(factory, cutEdges);
                }
            }


            // Handle doubly cut cells
            // =======================
            for(int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; jLevSet++) {
                if(iLevSet != jLevSet) {
                    //if(!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) 
                    {
                        for(int iKrefEdge = 0; iKrefEdge < gdat.iGeomEdges.EdgeRefElements.Length; iKrefEdge++) {
                            var KrefEdge = gdat.iGeomEdges.EdgeRefElements[iKrefEdge];

                            var cutDom1 = m_CutEdges[iKrefEdge, iLevSet];
                            var cutDom2 = m_CutEdges[iKrefEdge, jLevSet];
                            var doublyCut = cutDom1.Intersect(cutDom2);

                            if(doublyCut.NoOfItemsLocally > 0) {
                                var jmpJ = IdentifyWingA(jLevSet, spA);
                                foreach(var Kref in gdat.iGeomCells.RefElements) {

                                    var _doublyCut = doublyCut.Intersect(gdat.GetEdges4RefElement(Kref, KrefEdge));
                                    if(_doublyCut.NoOfItemsLocally <= 0)
                                        continue;

                                    var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, KrefEdge);
                                    var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, jLevSet, jmpJ, Kref, backupFactory);
                                    edgeQrIns.AddFactory(factory, _doublyCut);
                                }
                            }
                        }
                    }
                }
            }


            // Special case: coinciding level-sets
            // ===================================


            if(this.XDGSpaceMetrics.LevelSetRegions.LevSetCoincidingFaces != null || this.XDGSpaceMetrics.LevelSetRegions.LevSetCoincidingCoFaces != null) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Special case handling:
                // Some  level-set is (more-or-less) exactly on a cell edge;
                // therefore, we need some stable handling of such cases;
                // in the respective cells, the level-set surface quadrature will be overwritten
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                int D = gdat.SpatialDimension;

                foreach(var edgeKref in gdat.iGeomEdges.EdgeRefElements) {
                    int iKrefEdge = Array.IndexOf(gdat.iGeomEdges.EdgeRefElements, edgeKref);


                    EdgeMask IntegrationDom = Quadrature.LevelSetOnEdge.SurfaceElementBoundaryIntegration.ComputeMask(XDGSpaceMetrics.LevelSetRegions, iLevSet, iKrefEdge);
                    Debug.Assert(IntegrationDom.MaskType == MaskType.Geometrical);
                    var fact = new Quadrature.LevelSetOnEdge.SurfaceElementBoundaryIntegration(
                        edgeKref,
                        XDGSpaceMetrics.LevelSetData.ToArray(),
                        XDGSpaceMetrics.Tracker.GetLevelSetSignCodes(spA),
                        iLevSet,
                        this.XDGSpaceMetrics.LevelSetRegions.LevSetCoincidingFaces,
                        this.XDGSpaceMetrics.LevelSetRegions.LevSetCoincidingCoFaces
                        );
                    edgeQrIns.AddFactoryDomainPair(fact, IntegrationDom);


                }
            }


            

            //edgeQrIns.AddFactory(new QRoverride(), QRoverride.GetEM(gdat));


            return edgeQrIns;
        }



        /// <summary>
        /// Interior quadrature for the surface elements, i.e. for each cut background-cell $ K_j $ a quadrature to approximate
        /// \[
        ///    \int_{ K_j \cap \mathfrak{I}_i  \cap \partial \mathfrak{A} \cap \partial \mathfrak{B}  } \ldots \mathrm{dS} ,
        /// \]
        /// where: 
        /// - $ i $ is the level set index <paramref name="iLevSet"/>, i.e., $ \mathfrak{I}_i $ is the respective interface, 
        /// - $ \partial \mathfrak{A} $ and $ \cap \partial \mathfrak{B} $ are the boundaries of species <paramref name="spA"/> and <paramref name="spB"/>, respectively.   
        /// 
        /// These rules are used, e.g., by the (volume parts of) surface element operator <see cref="XDifferentialOperatorMk2.SurfaceElementOperator_Ls0"/>
        /// </summary>
        public CellQuadratureScheme Get_SurfaceElement_VolumeQuadScheme(SpeciesId spA, SpeciesId spB, int iLevSet) {
            if (!this.SpeciesList.Contains(spA))
                throw new ArgumentException("Given species (id = " + spA.cntnt + ") is not supported.");
            //Default behaviour: If Species are not divided by Level Set, function should not be called
            //Debug.Assert(!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp));

            //var spdom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(sp);
            //var IntegrationDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).Intersect(spdom);
            var IntegrationDom = this.LevelSetSeparationLayer[(spA, spB, iLevSet)];

            return GetLevelSetQuadScheme(iLevSet, spA, spB, IntegrationDom);
        }


        /// <summary>
        /// Intersection of two level-sets <paramref name="iLevSet0"/> and <paramref name="iLevSet1"/>
        /// \[
        ///    \int_{K_j \cap \mathfrak{I}_0 \cap \mathfrak{I}_1 \cap \partial \mathfrak{s} } \ldots \mathrm{dl} ,
        /// \]
        /// where $ \partial \mathfrak{s} $ is the boundary of species <paramref name="sp"/>
        /// </summary>
        public CellQuadratureScheme GetContactLineQuadScheme(SpeciesId sp, int iLevSet0, int iLevSet1) {
            if(iLevSet0 == iLevSet1)
                throw new ArgumentException("Both level-set indices are the same; self-intersection is not allowed.");
            if(iLevSet0 < 0 || iLevSet0 >= XDGSpaceMetrics.NoOfLevelSets)
                throw new ArgumentOutOfRangeException($"iLevSet0 index out-of-range: is {iLevSet0}, number of level sets: {XDGSpaceMetrics.NoOfLevelSets}");
            if(iLevSet1 < 0 || iLevSet1 >= XDGSpaceMetrics.NoOfLevelSets)
                throw new ArgumentOutOfRangeException($"iLevSet1 index out-of-range: is {iLevSet0}, number of level sets: {XDGSpaceMetrics.NoOfLevelSets}");

            //Find domain
            CellMask allDoublyCuts = CellMask.GetEmptyMask(gdat, MaskType.Geometrical);

            foreach(var Kref in gdat.iGeomCells.RefElements) {
                //for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                //    if (iLevSet != jLevSet) {
                if(!SpeciesAreSeparatedByLevSet(iLevSet1, sp, sp)) {
                    allDoublyCuts = allDoublyCuts.Union(GetDoubleCutCells(iLevSet0, iLevSet1));
                }
            }

            var spdom = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(sp).ToGeometicalMask();
            var IntegrationDom = allDoublyCuts.Intersect(spdom);
            var LevSetQrIns = new CellQuadratureScheme(
                new LevelSetIntersectionIntegrationMetric(this.XDGSpaceMetrics.LevelSetData[iLevSet0], this.XDGSpaceMetrics.LevelSetData[iLevSet1]),
                false, IntegrationDom);

            foreach(var Kref in gdat.iGeomCells.RefElements) {
                //for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                //    if (iLevSet != jLevSet) {
                if(!SpeciesAreSeparatedByLevSet(iLevSet1, sp, sp)) {
                    CellMask doublyCut = this.GetDoubleCutCells(iLevSet0, iLevSet1);
                    if(doublyCut.Count() > 0) {
                        var jmpJ = IdentifyWingA(iLevSet1, sp);
                        var KrefEdge = gdat.iGeomEdges.EdgeRefElements.Single();
                        var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet0, KrefEdge);
                        
                        var intersectionFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetIntersectionRuleFactory(iLevSet0, iLevSet1, Kref, backupFactory);
                        LevSetQrIns.AddFactory(intersectionFactory, doublyCut);
                    }
                }
            }

            //var rule1 = LevSetQrIns.Compile(gdat, 4);
            //var rule2 = (new QRoverride2()).GetQuadRuleSet(LevSetQrIns.Domain, 4)
            //LevSetQrIns.AddFactory(new QRoverride2(), IntegrationDom);

            // special case handling
            // =====================

            if(XDGSpaceMetrics.LevelSetRegions.LevSetCoincidingFaces != null) {
                int iKref = 0;
                foreach(var Kref in gdat.iGeomCells.RefElements) {
                    var doublyCutAndCoinciding = Quadrature.LevelSetOnEdge.InterfacesIntersection.ComputeCellMask(XDGSpaceMetrics.LevelSetRegions, iLevSet0, iLevSet1, iKref, IntegrationDom);


                    var intersectionFactory = new Quadrature.LevelSetOnEdge.InterfacesIntersection(
                        Kref, XDGSpaceMetrics.LevelSetRegions,
                        iLevSet0, iLevSet1,
                        this.XDGSpaceMetrics.XQuadFactoryHelper._GetSurfaceElement_BoundaryRuleFactory(iLevSet0, Kref),
                        this.XDGSpaceMetrics.XQuadFactoryHelper._GetSurfaceElement_BoundaryRuleFactory(iLevSet1, Kref),
                        this.XDGSpaceMetrics.CutCellMetrics.CutCellVolumes[sp]
                        );


                    LevSetQrIns.AddFactory(intersectionFactory, doublyCutAndCoinciding);

                    iKref++;
                }
            }


            //var rule = LevSetQrIns.Compile(gdat, XDGSpaceMetrics.CutCellQuadOrder);


            return LevSetQrIns;
        }

        /*
        class QRoverride2 : IQuadRuleFactory<QuadRule> {

            

            public RefElement RefElement => Square.Instance;

            public int[] GetCachedRuleOrders() {
                return [];
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {

                Console.WriteLine("go fuck yourself2");

                var ret = new List<ChunkRulePair<QuadRule>>();
                void AddQr(int iEdge, double x, double y, double w) {
                    var r = QuadRule.CreateBlank(Square.Instance, 1, 2, false);
                    r.Nodes[0, 0] = x;
                    r.Nodes[0, 1] = y;
                    r.Nodes.LockForever();
                    r.Weights[0] = w;
                    r.OrderOfPrecision = order + 12;
                    ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge), r));
                }

                foreach(int iCell in mask.ItemEnum) {
                    switch(iCell) {
                        case 6: {
                            AddQr(iCell, 1, 0.25, 1);
                            break;
                        }
                        default: {
                            AddQr(iCell, 0, 0, 0);
                            break;
                        }
                    }

                }
                return ret;
            }
        }
        */




        /*
        /// <summary>
        /// All cells, in which level-set <paramref name="iLevSet"/> intersects with some other level-set.
        /// </summary>
        public CellMask GetDoubleCutCells(int iLevSet) {
            CellMask allDoublyCuts = CellMask.GetEmptyMask(gdat, MaskType.Logical);

            //foreach (var Kref in gdat.Grid.RefElements) {
            for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                if (iLevSet != jLevSet) {
                    allDoublyCuts = allDoublyCuts.Union(GetDoubleCutCells(iLevSet, jLevSet));

                }
            }
            //}

            return allDoublyCuts;
        }
        */

        /// <summary>
        /// Quadrature for edges, i.e. for each cut background-cell $` K_j `$ and each species $` \mathfrak{s} `$ a quadrature to approximate
        /// ```math
        ///    \int_{\partial K_j \cap \mathfrak{s} } \ldots \mathrm{dS} .
        /// ```
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
                foreach (var KrefEdge in gdat.iGeomEdges.EdgeRefElements) {
                    for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                        if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {
                            
                            EdgeMask cutEdges = this.GetCutEdges(KrefEdge, iLevSet).Intersect(allRelevantEdges);
#if DEBUG
                            CellMask difference = cutEdges.GetAdjacentCells().Except(XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet));
                            if (difference.Count() > 0)
                                throw new ArithmeticException("Edges of the Cells" + difference.GetSummary() + " are detected as cut, but these cells are not contained in the cut Cell-Mask of the Level-Set-Tracker");
#endif
                            var jmp = IdentifyWingA(iLevSet, sp);
                            var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(iLevSet, jmp, KrefEdge);
                            edgeQrIns.AddFactoryDomainPair(factory, cutEdges, fixedOrder);
                        }
                    }
                }

                // overwrite with double-cut-cell-rules in double-cut-cells:
                //foreach(var Kref in gdat.iGeomCells.RefElements) {
                for(int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets...
                    if(!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {



                        var jmp = IdentifyWingA(iLevSet, sp);
                        //handle rules for cells/edges where two levelsets are present
                        for(int jLevSet = iLevSet + 1; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                            if(!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {


                                foreach(var Kref in gdat.iGeomCells.RefElements) {
                                    foreach(var KrefEdge in gdat.iGeomEdges.EdgeRefElements) {
                                        EdgeMask cutEdges = this.GetCutEdges(KrefEdge, iLevSet).Intersect(allRelevantEdges);
                                        EdgeMask doublyCut = cutEdges.Intersect(GetCutEdges(KrefEdge, jLevSet));
                                        
                                        if(doublyCut.NoOfItemsLocally > 0) {
                                            var jmpJ = IdentifyWingA(jLevSet, sp);

                                            var _doublyCut = doublyCut.Intersect(gdat.GetEdges4RefElement(Kref, KrefEdge));
                                            if(_doublyCut.NoOfItemsLocally <= 0)
                                                continue;


                                            var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(iLevSet, jmp, KrefEdge);
                                            var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(iLevSet, jmp, jLevSet, jmpJ, Kref, backupFactory);
                                            edgeQrIns.AddFactoryDomainPair(twoLSFactory, _doublyCut, fixedOrder);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                //}
                

                return edgeQrIns;
            }
        }


        private CellMask GetDoubleCutCells(int iLevSet, int jLevSet) {
            CellMask iCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
            CellMask jCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(jLevSet);
            return iCells.Intersect(jCells).ToGeometicalMask();
        }
        //private CellMask GetDoubleCutCells(int iLevSet, int jLevSet) {
        //    CellMask iCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);
        //    CellMask jCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(jLevSet);
        //    return iCells.Intersect(jCells);
        //}

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

                    int J = this.gdat.iLogicalCells.NoOfLocalUpdatedCells;
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
        /// <see cref="XDifferentialOperatorMk2.GhostEdgesOperator"/>
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
            return new EdgeQuadratureScheme(true, new EdgeMask(this.gdat, EdgBitMask));

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

               

                // default rule for "normal" cells
                var volQrIns = (new CellQuadratureScheme(UseDefaultFactories, CellMask));

                // now: rules for the cut-cells:
                for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets
                    if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {                    // it seems that with this usage it checks if both sides of this level set are different species (which can happen like {A,A})
                        var cutDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).ToGeometicalMask();
                        var cutCells = cutDom.Intersect(CellMask);

                        var jmp = IdentifyWingA(iLevSet, sp);

                        for (int iKref = 0; iKref < gdat.iGeomCells.RefElements.Length; iKref++) {
                            RefElement Kref = gdat.iGeomCells.RefElements[iKref];
                            var factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetVolRuleFactory(iLevSet, jmp, Kref);
                            var _cutDom = cutCells.Intersect(gdat.iGeomCells.GetCells4Refelement(Kref));
                            volQrIns.AddFactoryDomainPair(factory, _cutDom, fixedOrder);
                        }
                    }
                }

                //now: rules for the doubly cut-cells
                for (int iLevSet = 0; iLevSet < XDGSpaceMetrics.NoOfLevelSets; iLevSet++) { // loop over level sets
                    if (!SpeciesAreSeparatedByLevSet(iLevSet, sp, sp)) {                    // it seems that with this usage it checks if both sides of this level set are different species (which can happen like {A,A})
					var cutDom = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet).ToGeometicalMask();
                        var cutCells = cutDom.Intersect(CellMask);

                        var jmp = IdentifyWingA(iLevSet, sp);

                        for (int iKref = 0; iKref < gdat.iGeomCells.RefElements.Length; iKref++) {
                            RefElement Kref = gdat.iGeomCells.RefElements[iKref];
                            var _cutDom = cutCells.Intersect(gdat.iGeomCells.GetCells4Refelement(Kref));

                            //handle rules for cells/edges where two level sets are present
                            for (int jLevSet = iLevSet + 1; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                                if (!SpeciesAreSeparatedByLevSet(jLevSet, sp, sp)) {
                                    CellMask doublyCut = GetDoubleCutCells(iLevSet, jLevSet).Intersect(_cutDom);

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


        /*
        public CellQuadratureScheme bal(int iLevSet) {
            var CutCells = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask4LevSet(iLevSet);

            foreach(var Kref in gdat.iGeomCells.RefElements) {
                this.XDGSpaceMetrics.XQuadFactoryHelper._GetSurfaceElement_BoundaryRuleFactory()
            }
        }
        */






        /// <summary>
        /// Quadrature scheme for the integration over the level-set <paramref name="iLevSet"/>, i.e. for each cut background-cell \f$ K_j \f$ a quadrature to approximate
        /// \f[
        ///    \oint_{K_j \cap \mathfrak{I} \cap \partial \mathfrak{A} } \ldots \mathrm{dS} ,
        /// \f]
        /// where \f$ \mathfrak{A} \f$ is the first species in <see cref="SpeciesList"/>.
        /// </summary>
        public CellQuadratureScheme GetLevelSetQuadScheme(int iLevSet, CellMask IntegrationDom, int? fixedOrder = null) {

            var PossibleSpecies = Tracker.GetSpeciesSeparatedByLevSet(iLevSet).Select(Tracker.GetSpeciesId);

            return GetLevelSetQuadScheme(iLevSet, PossibleSpecies.First(), PossibleSpecies.ElementAt(1), IntegrationDom, fixedOrder);
        }

        /// <summary>
        /// Quadrature scheme for the integration over the level-set <paramref name="iLevSet"/>, i.e. for each cut background-cell \f$ K_j \f$ a quadrature to approximate
        /// \f[
        ///    \oint_{K_j \cap \mathfrak{I} \cap \partial \mathfrak{A} \cap \partial \mathfrak{b} } \ldots \mathrm{dS} ,
        /// \f]
        /// where \f$ \mathfrak{A} \f$ and \f$ \mathfrak{B} \f$  are the domains of species <paramref name="spA"/> and <paramref name="spB"/>, respectively.
        /// </summary>
        public CellQuadratureScheme GetLevelSetQuadScheme(int iLevSet, SpeciesId spA, SpeciesId spB, CellMask IntegrationDom, int? fixedOrder = null) {
            if (IntegrationDom.MaskType == MaskType.Logical)
                IntegrationDom = IntegrationDom.ToGeometicalMask();
            if(gdat.iGeomCells.RefElements.Length > 1)
                throw new NotImplementedException("more than one reference element is currently not supported");

            


            CellQuadratureScheme LevSetQrIns = new CellQuadratureScheme(
                scaling: new LevelSetIntegrationMetric(this.XDGSpaceMetrics.LevelSetData[iLevSet]),
                UseDefaultFactories: false, domain: IntegrationDom);

            // "ordinary" cut cells for level set `iLevSet`
            // ============================================
            foreach(var Kref in gdat.iGeomCells.RefElements) {
                var surfaceFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);

                // exclude doubly-cut cells;
                // (should not be necessary, because for a quadrature scheme, the rules which come later would overwrite the earier ones)
                CellMask modIntegrationDom = IntegrationDom;
                for(int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) {
                    if(jLevSet != iLevSet) {
                        CellMask doublyCut = GetDoubleCutCells(iLevSet, jLevSet);
                        if(doublyCut.NoOfItemsLocally > 0) {
                            modIntegrationDom = modIntegrationDom.Except(doublyCut);
                        }
                    }
                }

                LevSetQrIns.AddFactoryDomainPair(surfaceFactory, modIntegrationDom, fixedOrder);
            }

            // special case: double-cut cut cells
            // ===================================
            foreach(var Kref in gdat.iGeomCells.RefElements) {
                //handle rules for cells/edges where two levelsets are present

                CellMask modIntegrationDom = IntegrationDom;
                for (int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; ++jLevSet) { // loop over all other 'trimming' level-sets...
                    if (jLevSet != iLevSet) {
                        CellMask doublyCut = GetDoubleCutCells(iLevSet, jLevSet);
                        if (doublyCut.NoOfItemsLocally > 0) {
                            modIntegrationDom = modIntegrationDom.Except(doublyCut);

                            //var jmpA = IdentifyWingA(jLevSet, spA);
                            var jmpA = GetTrimmingLevelSetSign(jLevSet);

                            var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, jLevSet, jmpA, Kref, null);
                            LevSetQrIns.AddFactoryDomainPair(twoLSFactory, doublyCut, fixedOrder);

                            /*

                            if(gdat.iGeomEdges.EdgeRefElements.Length > 1)
                                    throw new NotImplementedException("Not implemented for more than one edge ref element");

                            //Debug.Assert(jmpA == IdentifyWing(jLevSet, spB));
                            if (iLevSet == 1) {

                                var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(jLevSet, jmpA, gdat.iGeomEdges.EdgeRefElements[0]);

                                //var _doublyCut = doublyCut.Intersect(gdat.GetEdges4RefElement(Kref, KrefEdge));
                                var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, jLevSet, jmpA, Kref, backupFactory);
                                LevSetQrIns.AddFactoryDomainPair(twoLSFactory, doublyCut, fixedOrder);

                            } else {
                                var backupFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, Kref);
                                var twoLSFactory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetSurfaceFactory(iLevSet, jLevSet, jmpA, Kref, backupFactory);
                                LevSetQrIns.AddFactoryDomainPair(twoLSFactory, doublyCut, fixedOrder);
                            }
                            */
                        }
                    }
                }
            }

            // special case: level set `iLevSet` (might) coincide with some edges
            // ===================================================================

            if(this.XDGSpaceMetrics.LevelSetRegions.LevSetCoincidingFaces != null) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Special case handling:
                // Some  level-set is (more-or-less) exactly on a cell edge;
                // therefore, we need some stable handling of such cases;
                // in the respective cells, the level-set surface quadrature will be overwritten
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                foreach(var Kref in gdat.iGeomCells.RefElements) {
                    int iKref = Array.IndexOf(gdat.iGeomCells.RefElements, Kref);
                    
                    var coincidingCells = Quadrature.LevelSetOnEdge.InterfaceIntegration.ComputeCellMask(XDGSpaceMetrics.LevelSetRegions, iLevSet, iKref);
                    var mask = coincidingCells.Intersect(IntegrationDom);

                    var fact = new Quadrature.LevelSetOnEdge.InterfaceIntegration(Kref, this.XDGSpaceMetrics.LevelSetData[iLevSet]);
                    LevSetQrIns.AddFactoryDomainPair(fact, mask, fixedOrder);


                    // coinciding level-sets which are cut by other level-sets...
                    for(int jLevSet = 0; jLevSet < XDGSpaceMetrics.NoOfLevelSets; jLevSet++) { // loop over all other level-sets...
                        if(jLevSet != iLevSet) {
                            CellMask doublyCut = GetDoubleCutCells(iLevSet, jLevSet);
                            var doublyCut_onEdge = doublyCut.Intersect(mask);

                            if(doublyCut_onEdge.NoOfItemsLocally > 0) {
                                //var jmpA = IdentifyWingA(jLevSet, spA);
                                var jmpA = GetTrimmingLevelSetSign(jLevSet);

                                var KrefEdge = gdat.iGeomEdges.EdgeRefElements.Single();
                                var trimming_factory = this.XDGSpaceMetrics.XQuadFactoryHelper.GetEdgeRuleFactory(jLevSet, jmpA, KrefEdge);

                                var dblCutEdges = this.GetCutEdges(KrefEdge, jLevSet);

                                var fact_dbl_cut = new Quadrature.LevelSetOnEdge.InterfaceIntegration4DoubleCut(
                                    Kref, 
                                    this.XDGSpaceMetrics.LevelSetData[iLevSet], 
                                    trimming_factory, dblCutEdges,
                                    XDGSpaceMetrics.LevelSetData.ToArray(),
                                    XDGSpaceMetrics.Tracker.GetLevelSetSignCodes(spA));

                                LevSetQrIns.AddFactoryDomainPair(fact_dbl_cut, doublyCut_onEdge, fixedOrder);
                            }
                        }
                    }
                }
            }
            return LevSetQrIns;

            JumpTypes bool2JumpType(bool b) => b ? JumpTypes.Heaviside : JumpTypes.OneMinusHeaviside;

            JumpTypes GetTrimmingLevelSetSign(int iTrimmingLevSet) {
                var spnA = Tracker.GetSpeciesName(spA);
                var spnB = Tracker.GetSpeciesName(spB);
                var SignsA = Tracker.GetLevelSetSignCodes(spnA).Select(sign_code => sign_code.GetSign(iTrimmingLevSet));
                if(SignsA.Count() == 1)
                    return bool2JumpType(SignsA.Single());
                var SignsB = Tracker.GetLevelSetSignCodes(spnB).Select(sign_code => sign_code.GetSign(iTrimmingLevSet));
                if(SignsB.Count() == 1)
                    return bool2JumpType(SignsB.Single());
                throw new ApplicationException($"unable to determine required jump type/sign of trimming level-set {iTrimmingLevSet} for integration over Species {spnA}-{spnB} intersection and level-set {iLevSet}.");
            }
        }

        /// <summary>
        /// Writes diagnostic information about quadrature rules into csv-textfiles.
        /// </summary>
        public void RuleInfo(SpeciesId spc, string volruleName, string levsetRuleName, string cellBndRule, string edgeRuleName, int order, int iLevSet) {


            // test parameters
            var jmp = Foundation.XDG.Quadrature.JumpTypes.Heaviside;
            var sch = this.XDGSpaceMetrics.XQuadFactoryHelper;
            var spNm = this.XDGSpaceMetrics.LevelSetRegions.GetSpeciesName(spc);
            //var spId = LsTrk.GetSpeciesId(spNm);
            var DomainOfInterest = XDGSpaceMetrics.LevelSetRegions.GetCutCellSubgrid4LevSet(iLevSet);

            if (volruleName != null) {
                int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
                double[] TotWeights = new double[J];

                foreach (var Kref in this.gdat.iGeomCells.RefElements) {
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
                int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
                double[] TotWeights = new double[J];

                foreach (var Kref in this.gdat.iGeomCells.RefElements) {
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
                int NoEdg = this.gdat.iGeomEdges.Count;
                double[] TotWeights = new double[NoEdg];

                foreach (var Kref in this.gdat.iGeomCells.RefElements) {
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

            
        }
        
    }
}