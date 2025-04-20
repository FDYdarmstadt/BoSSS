using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature {

    /// <summary>
    /// Generation of quadrature rules for level-sets which coincide with cell faces,
    /// centrally triggered by
    /// <see cref="LevelSetTracker.LevelSetRegions.LevSetCoincidingFaces"/>
    /// </summary>
    class LevelSetOnEdgeRuleFactory : IQuadRuleFactory<QuadRule> {


        static public CellMask ComputeMask(LevelSetTracker tracker, int iLevelSet, int iKref) {
            var CoincidFaces = tracker.Regions.LevSetCoincidingFaces;
            var gdat = tracker.GridDat;
            
            if(CoincidFaces == null) {
                return CellMask.GetEmptyMask(gdat, MaskType.Geometrical);
            }

            int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
            BitArray bitMask = new BitArray(J);

            for(int j = 0; j < J; j++) {
                if(gdat.iGeomCells.GetRefElementIndex(j) != iKref)
                    continue;

                if(CoincidFaces[j] == null) 
                    continue;

                foreach(var t in CoincidFaces[j]) {
                    if(t.iLevSet == iLevelSet)
                        bitMask[j] = true;
                }

            }


            return new CellMask(gdat, bitMask, MaskType.Geometrical);
        }



        public LevelSetOnEdgeRuleFactory(RefElement _RefElement, LevelSetTracker.LevelSetData levelSetData) {
            RefElement = _RefElement;
            m_LevelSetData = levelSetData;
            m_LevelSetOnEdgeRule = new LevelSetOnEdgeRule(m_LevelSetData);
        }

        public RefElement RefElement {
            get;
            private set;
        }

        readonly LevelSetTracker.LevelSetData m_LevelSetData;

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        readonly LevelSetOnEdgeRule m_LevelSetOnEdgeRule;

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if(mask.MaskType != MaskType.Geometrical) {
                throw new ArgumentException($"wrong type of mask ({mask.MaskType})");
            }

            var ret = new List<ChunkRulePair<QuadRule>>();

            foreach(int jCell in mask.ItemEnum) {
                ret.Add(m_LevelSetOnEdgeRule.SurfaceQuadRule(order, jCell));
            }

            return ret;
        }
    }


    /// <summary>
    /// Generation of quadrature rules for level-sets which coincide with cell faces,
    /// centrally triggered by
    /// <see cref="LevelSetTracker.LevelSetRegions.LevSetCoincidingFaces"/>
    /// </summary>
    class LevelSetBoundaryOnEdgeRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {


        static IEnumerable<int> GeometricalCellFace2Edges(IGridData gdat, int jCellGeom, int iFace) {

            int jCellLog;
            if(gdat.iGeomCells.GeomCell2LogicalCell != null)
                jCellLog = gdat.iGeomCells.GeomCell2LogicalCell[jCellGeom];
            else
                jCellLog = jCellGeom;


            int[] LogEdges = gdat.iLogicalCells.Cells2Edges[jCellLog];

            var ret = new List<int>();
            foreach(int __iLogEdge in LogEdges) {
                int iLogEdge = Math.Abs(__iLogEdge) - 1;

                foreach(int iGeomEdge in gdat.GetGeometricEdgeIndices(iLogEdge)) {
                    int match = -1;
                    if(gdat.iGeomEdges.CellIndices[iGeomEdge, 0] == jCellGeom)
                        match = 0;
                    else if(gdat.iGeomEdges.CellIndices[iGeomEdge, 1] == jCellGeom)
                        match = 1;

                    if(match >= 0) {
                        if(gdat.iGeomEdges.FaceIndices[iGeomEdge, match] == iFace) {
                            ret.Add(iGeomEdge);
                        }

                    }
                }
            }

            return ret;




        }

        static public EdgeMask ComputeMask(LevelSetTracker tracker, int iLevelSet, int iKref) {
            var CoincidFaces = tracker.Regions.LevSetCoincidingFaces;
            IGridData gdat = tracker.GridDat;

            if(CoincidFaces == null) {
                return EdgeMask.GetEmptyMask(gdat, MaskType.Geometrical);
            }

            int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
            int E = gdat.iGeomEdges.Count; ;
            BitArray bitMask = new BitArray(E);

            for(int j = 0; j < J; j++) {
                if(CoincidFaces[j] == null)
                    continue;


                foreach(var t in CoincidFaces[j]) {
                    if(t.iLevSet == iLevelSet) {

                        foreach(var e in GeometricalCellFace2Edges(gdat, j, t.iFace)) {
                            if(gdat.iGeomEdges.GetRefElementIndex(e) == iKref) {

                                bitMask[e] = true;
                            }
                        }

                    }
                }

            }


            return new EdgeMask(gdat, bitMask, MaskType.Geometrical);
        }


        public LevelSetBoundaryOnEdgeRuleFactory(RefElement _RefElement, LevelSetTracker.LevelSetData levelSetData) {
            RefElement = _RefElement;
            m_LevelSetData = levelSetData;
            m_LevelSetOnEdgeRule = new LevelSetOnEdgeRule(m_LevelSetData);
        }

        public RefElement RefElement {
            get;
            private set;
        }

        readonly LevelSetTracker.LevelSetData m_LevelSetData;

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        readonly LevelSetOnEdgeRule m_LevelSetOnEdgeRule;

        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if(mask.MaskType != MaskType.Geometrical) {
                throw new ArgumentException($"wrong type of mask ({mask.MaskType})");
            }

            var ret = new List<ChunkRulePair<CellBoundaryQuadRule>>();

            foreach(int jCell in mask.ItemEnum) {
                ret.Add(m_LevelSetOnEdgeRule.BoundaryLineRule(order, jCell));
                Debug.Assert(ret.Last().Rule.NumbersOfNodesPerFace != null, "Not a valid boundary rule");
            }

            return ret;
        }
    }



    class SurfaceElementBounaryFactoryForCoincidingEdges : IQuadRuleFactory<QuadRule> {


        /// <summary>
        /// conversion (<paramref name="jCellGeom"/>,<paramref name="iFace"/>) to grid edges
        /// </summary>
        static IEnumerable<int> GeometricalCellFace2Edges(IGridData gdat, int jCellGeom, int iFace) {

            int jCellLog;
            if(gdat.iGeomCells.GeomCell2LogicalCell != null)
                jCellLog = gdat.iGeomCells.GeomCell2LogicalCell[jCellGeom];
            else
                jCellLog = jCellGeom;


            int[] LogEdges = gdat.iLogicalCells.Cells2Edges[jCellLog];

            var ret = new List<int>();
            foreach(int __iLogEdge in LogEdges) {
                int iLogEdge = Math.Abs(__iLogEdge) - 1;

                foreach(int iGeomEdge in gdat.GetGeometricEdgeIndices(iLogEdge)) {
                    int match = -1;
                    if(gdat.iGeomEdges.CellIndices[iGeomEdge, 0] == jCellGeom)
                        match = 0;
                    else if(gdat.iGeomEdges.CellIndices[iGeomEdge, 1] == jCellGeom)
                        match = 1;

                    if(match >= 0) {
                        if(gdat.iGeomEdges.FaceIndices[iGeomEdge, match] == iFace) {
                            ret.Add(iGeomEdge);
                        }

                    }
                }
            }

            return ret;




        }

        /// <summary>
        /// conversion (<paramref name="jCellGeom"/>,<paramref name="iFace"/>) to grid edges
        /// </summary>
        static IEnumerable<int> GeometricalCellCoFace2Edges(IGridData gdat, int jCellGeom, int iCoFace) {

            int jCellLog;
            if(gdat.iGeomCells.GeomCell2LogicalCell != null)
                jCellLog = gdat.iGeomCells.GeomCell2LogicalCell[jCellGeom];
            else
                jCellLog = jCellGeom;


            var Kref = gdat.iGeomCells.GetRefElement(jCellGeom);
            int face0 = Kref.CoFaceToFaceIndices[iCoFace, 0];
            int face1 = Kref.CoFaceToFaceIndices[iCoFace, 0];



            return GeometricalCellFace2Edges(gdat, jCellGeom, face0).SetUnion(GeometricalCellFace2Edges(gdat, jCellGeom, face0));
        }




        /// <summary>
        /// on which edges this class can actually provide working quadrature rules?
        /// </summary>
        static public EdgeMask ComputeMask(LevelSetTracker tracker, int iLevelSet, int iKrefEdge) {
            var CoincidFaces = tracker.Regions.LevSetCoincidingFaces;
            var CoincidCoFcs = tracker.Regions.LevSetCoincidingCoFaces;
            IGridData gdat = tracker.GridDat;

            if(CoincidFaces == null && CoincidCoFcs == null) {
                return EdgeMask.GetEmptyMask(gdat, MaskType.Geometrical);
            }

            int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
            int E = gdat.iGeomEdges.Count; 
            BitArray bitMask = new BitArray(E);

            for(int j = 0; j < J; j++) {
                if(CoincidFaces[j] != null) {
                    foreach(var t in CoincidFaces[j]) {
                        if(t.iLevSet == iLevelSet) {

                            foreach(var e in GeometricalCellFace2Edges(gdat, j, t.iFace)) {
                                if(gdat.iGeomEdges.GetRefElementIndex(e) == iKrefEdge) {

                                    bitMask[e] = true;
                                }
                            }

                        }
                    }
                }

                if(CoincidCoFcs[j] != null) {
                    foreach(var t in CoincidCoFcs[j]) {
                        if(t.iLevSet == iLevelSet) {
                            foreach(var e in GeometricalCellCoFace2Edges(gdat, j, t.iCoFace)) {
                                if(gdat.iGeomEdges.GetRefElementIndex(e) == iKrefEdge) {

                                    bitMask[e] = true;
                                }
                            }
                        }
                    }
                }
            }


            return new EdgeMask(gdat, bitMask, MaskType.Geometrical);
        }



        public SurfaceElementBounaryFactoryForCoincidingEdges(RefElement _RefElement, LevelSetTracker.LevelSetData[] __levelSetDataS, LevelSetSignCode[] __allSignCodes, int __iLevSet) {
            RefElement = _RefElement;
            if(Array.IndexOf(__levelSetDataS[0].GridDat.iGeomEdges.EdgeRefElements, this.RefElement) < 0) {
                throw new ArgumentException($"{_RefElement} is not an edge reference element of the grid.");
            }
            m_LevelSetDataS = __levelSetDataS;
            m_allSignCodes = __allSignCodes;
            m_iLevSet = __iLevSet;
        }

        public RefElement RefElement {
            get;
            private set;
        }
#pragma warning disable 169,649
        readonly LevelSetTracker.LevelSetData[] m_LevelSetDataS;
        readonly LevelSetSignCode[] m_allSignCodes;
        readonly int m_iLevSet;

        /// <summary>
        /// <see cref="LevelSetTracker.LevelSetRegions.LevSetCoincidingFaces"/>
        /// </summary>
        readonly (int iLevSet, int iFace)[][] m_LevSetCoincidingFaces;

        /// <summary>
        /// <see cref="LevelSetTracker.LevelSetRegions.LevSetCoincidingCoFaces"/>
        /// </summary>
        readonly (int iLevSet, int iFace)[][] LevSetCoincidingFaces;


        /// <summary>
        /// <see cref="LevSetCoincidingCoFaces"/>>
        /// </summary>
        internal (int iLevSet, int iCoFace)[][] LevSetCoincidingCoFaces;


        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if(!(mask is EdgeMask edgeMask)) {
                throw new ArgumentException("expecting an edge mask");
            }
            if(!(mask.MaskType != MaskType.Geometrical)) {
                throw new ArgumentException("expecting an geometrical mask");
            }
            var gdat = mask.GridData;
            int iKref = Array.IndexOf(gdat.iGeomEdges.EdgeRefElements, this.RefElement);

            var fulCoFaceRule = RefElement.FaceRefElement.GetQuadratureRule(order);
            var emptyFaceRule = QuadRule.CreateBlank(RefElement, 1, Math.Max(RefElement.FaceRefElement.SpatialDimension, 1), true);
            emptyFaceRule.OrderOfPrecision = order;
            emptyFaceRule.Nodes.LockForever();

            int[,] Edge2Cell = gdat.iGeomEdges.CellIndices;
            int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
            byte[,] Edge2Face = gdat.iGeomEdges.FaceIndices;

            int[][,] CoFaceToFace = gdat.iGeomCells.RefElements.Select(Kref => Kref.CoFaceToFaceIndices).ToArray();


            var compRule = new ChunkRulePair<QuadRule>[mask.NoOfItemsLocally];
            int cnt = -1;
            foreach(int iEdge in mask.ItemEnum) {
                cnt++;
                if(gdat.iGeomEdges.GetRefElementIndex(iEdge) != iKref)
                    throw new ArgumentException("mask violates the element");

                bool bFound = false;
                for(int inOt = 0; inOt < 2; inOt++) {
                    if(bFound)
                        continue;
                    int jCell = Edge2Cell[iEdge, inOt];
                    if(jCell < 0)
                        continue;
                    int iCellFace = Edge2Face[jCell, inOt];



                    if(this.LevSetCoincidingFaces != null && this.LevSetCoincidingFaces[jCell] != null) {
                        if(this.LevSetCoincidingFaces[jCell].Any(tt => (tt.iLevSet == this.m_iLevSet && tt.iFace == iCellFace))) {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // entire face coincides wit the level-set
                            // in this case, we don't want any co-co-dim quadrature 
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++
                            bFound = true;
                            compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge), emptyFaceRule);
                            break;
                        }
                    }
                    
                    //var gdat.iGeomCells.GetRefElement(jCell);

                    if(this.LevSetCoincidingCoFaces != null && this.LevSetCoincidingCoFaces[jCell] != null) {
                        foreach(var tt in this.LevSetCoincidingCoFaces[jCell]) {
                            if(tt.iLevSet != this.m_iLevSet)
                                // different level-set
                                continue;

                            int[,] cf2f = CoFaceToFace[gdat.iGeomCells.GetRefElementIndex(jCell)];

                            if(iCellFace != cf2f[iCellFace, 0] && iCellFace != cf2f[iCellFace, 1])
                                // some face that we are not interested in
                                continue;



                        }

                    }

                }


            }


            return compRule;
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }
    }



    /// <summary>
    /// Generation of quadrature rules for level-sets which coincide with cell faces,
    /// centrally triggered by
    /// <see cref="LevelSetTracker.LevelSetRegions.LevSetCoincidingFaces"/>
    /// </summary>
    class LevelSetOnEdgeRule {
        GridData grddat;

        LevelSetTracker.LevelSetData levelSetData;

        int D;

        public LevelSetOnEdgeRule(LevelSetTracker.LevelSetData levelSetData) {
            this.levelSetData = levelSetData;
            grddat = levelSetData.GridDat;
            D = grddat.SpatialDimension;
        }

        public bool IsSpecialCell(int j) {
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData.Region.LevSetCoincidingFaces;
            if (CoIncFaces == null)
                return false;
            if (CoIncFaces[j] == null)
                return false;

            foreach (var t in CoIncFaces[j]) {
                int levelSetIndex = levelSetData.LevelSetIndex;
                if (t.iLevSet == levelSetIndex)
                    return true; // jetzt geht der Spass los!
            }
            return false;
        }

        int GetSpecialFaceIndex(int j) {
            (int iLevSet, int iFace)[][] CoIncFaces = levelSetData.Region.LevSetCoincidingFaces;
            foreach (var t in CoIncFaces[j]) {
                int levelSetIndex = levelSetData.LevelSetIndex;
                if (t.iLevSet == levelSetIndex)
                    return t.iFace; // jetzt geht der Spass los!
            }
            throw new Exception("Face does not have registered special Face");
        }

        /// <summary>
        /// **This is a special case: level-set coincides with a cell-face; 
        /// cell is either full or empty:**
        /// 
        /// Call this function when the <see cref="IsSpecialCell"/> returns true.
        /// 
        /// </summary>
        public (ChunkRulePair<QuadRule> surfaceRule, ChunkRulePair<QuadRule> volumeRule) ComboQuadRule(int intOrder, int jCell, int speciesSign = 1) {
            ChunkRulePair<QuadRule> surfaceRule = SurfaceQuadRule(intOrder, jCell);
            ChunkRulePair<QuadRule> volumeRule = VolumeQuadRule(intOrder, jCell, speciesSign);
            return (surfaceRule, volumeRule);
        }

        /// <summary>
        /// **This is a special case: level-set coincides with a cell-face; 
        /// cell is either full or empty:**
        /// 
        /// Call this function when the <see cref="IsSpecialCell"/> returns true.
        /// 
        /// Note:
        /// Only one of the two cells on the special edge should have a surface rule, the other cell should be empty
        /// Otherwise, the level-set components are integrates twice (or never)!
        ///
        /// The convention which we use is:
        ///   - the cell with the **lower global index** has the full rule
        ///   - the other cell is empty
        /// We use the global index here, so that the result is "stable" even if we are at an MPI boundary.
        /// </summary>
        /// <param name="intOrder"></param>
        /// <param name="jCell"></param>
        /// <returns></returns>
        public ChunkRulePair<QuadRule> SurfaceQuadRule(int intOrder, int jCell) {
            int SpecialFace = GetSpecialFaceIndex(jCell);
            RefElement refElement = grddat.Cells.GetRefElement(jCell);

            // determine whether this cell is 
            // inside or outside 
            // w.r.t. the edge that corresponds with face 
            int iEdge = grddat.GetEdgesForFace(jCell, SpecialFace, out int InOrOut, out int[] FurtherEdges);
            if (FurtherEdges != null && FurtherEdges.Length > 0) {
                // throw new NotSupportedException("Hanging node on a edge which coincides with the level set - this should be avoided.");
                // Console.WriteLine("Hanging node on a edge which coincides with the level set - this should be avoided.");
            }
            int J = grddat.CellPartitioning.LocalLength;

            // surface rule -- classification, using con-formality of edges
            // always use conformal cell, if both conformal choose lower global index
            bool EmptyOrFool; // true: full rule; false: empty rule
            {
                int OtherCell = grddat.Edges.CellIndices[iEdge, InOrOut == 0 ? 1 : 0];
                bool jConform = InOrOut == 0 ? grddat.Edges.IsEdgeConformalWithCell1(iEdge) : grddat.Edges.IsEdgeConformalWithCell2(iEdge);
                bool OtherConform = InOrOut == 0 ? grddat.Edges.IsEdgeConformalWithCell2(iEdge) : grddat.Edges.IsEdgeConformalWithCell1(iEdge);

                // if both cell faces have no hanging nodes use globally lower index
                if (jConform & OtherConform) {
                    long jCellGlob = jCell + grddat.CellPartitioning.i0;
                    long OtherCellGlob = OtherCell < J ? OtherCell + grddat.CellPartitioning.i0 : grddat.Parallel.GlobalIndicesExternalCells[OtherCell - J];
                    EmptyOrFool = jCellGlob < OtherCellGlob;
                } else if (jConform & !OtherConform) {
                    EmptyOrFool = true; // this cells face consists of one edge only, empty rule
                } else if (!jConform & OtherConform) {
                    EmptyOrFool = false; // this cells face consists of more than one edge only, use full rule, so that each partial edge contains the rule of full degree!
                } else {
                    throw new NotSupportedException(String.Format("Error in cell {0}, {1}: Only one cell should have the hanging node.", jCell, OtherCell));
                }
            }

            ChunkRulePair<QuadRule> surfaceRule;
            if (EmptyOrFool) {

                var FaceRule = refElement.FaceRefElement.GetQuadratureRule(intOrder);
                int K = FaceRule.NoOfNodes;
                NodeSet VolumeNodes = new NodeSet(refElement, K, D, true);
                refElement.TransformFaceCoordinates(SpecialFace, FaceRule.Nodes, VolumeNodes);
                VolumeNodes.LockForever();

                double gTrF = refElement.FaceTrafoGramianSqrt[SpecialFace];
                //var metrics = this.levelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(VolumeNodes, jCell, 1);

                QuadRule qr_l = new QuadRule() {
                    OrderOfPrecision = FaceRule.OrderOfPrecision,
                    Weights = MultidimensionalArray.Create(K),
                    Nodes = VolumeNodes
                };

                for (int k = 0; k < K; k++) {
                    qr_l.Weights[k] = FaceRule.Weights[k] * gTrF;// / metrics[0, k];
                }

                surfaceRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);

            } else {
                // use empty rule
                QuadRule empty = new QuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(1),
                    Nodes = refElement.Center
                };

                surfaceRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), empty);
            }
            return surfaceRule;
        }


        static IEnumerable<int> ContainingFaces(RefElement R, Vector V, double eps = 1.0e-10) {
            var ret = new List<int>();

            for(int iFace = 0; iFace < R.NoOfFaces; iFace++) {

                var vFace = R.GetInverseEmbeddedFaceTrafo(iFace).Transform(V);
                double dist = vFace[vFace.Dim - 1];
                var vFacePrj = new Vector(vFace.Take(vFace.Dim - 1).ToArray());
                //Console.WriteLine("f" + iFace + ": " + vFace + " dist = " + dist + " within? " + );
                bool onFace = (R.FaceRefElement.IsWithin(vFacePrj, eps) && (dist.Abs() < eps));
                if(onFace)
                    ret.Add(iFace);

            }

            return ret.ToArray();
        }


        /// <summary>
        /// **This is a special case: level-set coincides with a cell-face; 
        /// cell is either full or empty:**
        /// 
        /// Call this function when the <see cref="IsSpecialCell"/> returns true.
        /// 
        /// Note:
        /// Only one of the two cells on the special edge should have a surface rule, the other cell should be empty
        /// Otherwise, the level-set components are integrates twice (or never)!
        ///
        /// The convention which we use is:
        ///   - the cell with the **lower global index** has the full rule
        ///   - the other cell is empty
        /// We use the global index here, so that the result is "stable" even if we are at an MPI boundary.
        /// </summary>
        /// 
        /// <param name="intOrder"></param>
        /// <param name="jCell"></param>
        /// <returns></returns>
        public ChunkRulePair<CellBoundaryQuadRule> BoundaryLineRule(int intOrder, int jCell) {
            int SpecialFace = GetSpecialFaceIndex(jCell);
            RefElement refElement = grddat.Cells.GetRefElement(jCell);

            // determine whether this cell is 
            // inside or outside 
            // w.r.t. the edge that corresponds with face 
            int iEdge = grddat.GetEdgesForFace(jCell, SpecialFace, out int InOrOut, out int[] FurtherEdges);
            if(FurtherEdges != null && FurtherEdges.Length > 0) {
                // throw new NotSupportedException("Hanging node on a edge which coincides with the level set - this should be avoided.");
                // Console.WriteLine("Hanging node on a edge which coincides with the level set - this should be avoided.");
            }
            int J = grddat.CellPartitioning.LocalLength;

            // surface rule -- classification, using con-formality of edges
            // always use conformal cell, if both conformal choose lower global index
            bool EmptyOrFool; // true: full rule; false: empty rule
            {
                int OtherCell = grddat.Edges.CellIndices[iEdge, InOrOut == 0 ? 1 : 0];
                bool jConform = InOrOut == 0 ? grddat.Edges.IsEdgeConformalWithCell1(iEdge) : grddat.Edges.IsEdgeConformalWithCell2(iEdge);
                bool OtherConform = InOrOut == 0 ? grddat.Edges.IsEdgeConformalWithCell2(iEdge) : grddat.Edges.IsEdgeConformalWithCell1(iEdge);

                // if both cell faces have no hanging nodes use globally lower index
                if(jConform & OtherConform) {
                    long jCellGlob = jCell + grddat.CellPartitioning.i0;
                    long OtherCellGlob = OtherCell < J ? OtherCell + grddat.CellPartitioning.i0 : grddat.Parallel.GlobalIndicesExternalCells[OtherCell - J];
                    EmptyOrFool = jCellGlob < OtherCellGlob;
                } else if(jConform & !OtherConform) {
                    EmptyOrFool = true; // this cells face consists of one edge only, empty rule
                } else if(!jConform & OtherConform) {
                    EmptyOrFool = false; // this cells face consists of more than one edge only, use full rule, so that each partial edge contains the rule of full degree!
                } else {
                    throw new NotSupportedException(String.Format("Error in cell {0}, {1}: Only one cell should have the hanging node.", jCell, OtherCell));
                }
            }

            ChunkRulePair<CellBoundaryQuadRule> boundaryQuadRule;
            if(EmptyOrFool) {
                // ++++++++++++++
                // non-empty rule
                // ++++++++++++++

                //
                // the `SpecialFace` coincides with the level-set;
                // a boundary rule on this face would give us a correct BoundaryLineRule;
                // The trick part is that we need to associate the quadrature nodes with the neighbor faces of `SpecialFace`
                // Otherwise, then the cell-boundary-rule to edge-rule conversion will produce double integration points
                //


                var FaceRule = refElement.FaceRefElement.GetBoundaryQuadRule(intOrder);
                int K = FaceRule.NoOfNodes;

                
                MultidimensionalArray VolumeNodes = MultidimensionalArray.Create(K, D);
                refElement.TransformFaceCoordinates(SpecialFace, FaceRule.Nodes, VolumeNodes);

                double gTrF = refElement.FaceTrafoGramianSqrt[SpecialFace];
                
                int[] NumbersOfNodesPerFace = new int[refElement.NoOfFaces];
                List<Vector>[] pointsOnFace = refElement.NoOfFaces.ForLoop(i => new List<Vector>());
                List<double>[] weightOnFace = refElement.NoOfFaces.ForLoop(i => new List<double>());
                for(int k = 0; k < K; k++) {
                    var pt = VolumeNodes.GetRowPt(k);
                    double w = FaceRule.Weights[k] * gTrF;
                    var allFaces = ContainingFaces(refElement, pt).ToList();
                    allFaces.Remove(SpecialFace);
                    var iFaceOther = allFaces[0];
                    NumbersOfNodesPerFace[iFaceOther]++;
                    pointsOnFace[iFaceOther].Add(pt);
                    weightOnFace[iFaceOther].Add(w);
                }


                CellBoundaryQuadRule qr_l = new CellBoundaryQuadRule() {
                    OrderOfPrecision = FaceRule.OrderOfPrecision,
                    Weights = MultidimensionalArray.Create(K),
                    Nodes = new NodeSet(refElement, K, D, false),
                    NumbersOfNodesPerFace = NumbersOfNodesPerFace
                };
                
                int _k = 0;
                for(int iFace = 0; iFace < refElement.NoOfFaces; iFace++) {
                    for(int __k = 0; __k < NumbersOfNodesPerFace[iFace]; __k++) {
                        qr_l.Nodes.SetRowPt(_k, pointsOnFace[iFace][__k]);
                        qr_l.Weights[_k] = weightOnFace[iFace][__k];
                        _k++;
                    }
                }

                VolumeNodes.LockForever();

                boundaryQuadRule = new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);
            } else {
                // ++++++++++++++
                // use empty rule
                // ++++++++++++++
                CellBoundaryQuadRule empty = new CellBoundaryQuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(1),
                    Nodes = refElement.Center,
                    NumbersOfNodesPerFace = new int[refElement.NoOfFaces]
                };

                boundaryQuadRule = new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), empty);
            }
            return boundaryQuadRule;
        }




        /// <summary>
        /// **This is a special case: level-set coincides with a cell-face; 
        /// cell is either full or empty:**
        /// 
        /// Call this function when the <see cref="IsSpecialCell"/> returns true.
        /// 
        /// Note:
        /// Only one of the two cells on the special edge should have a surface rule, the other cell should be empty
        /// Otherwise, the level-set components are integrates twice (or never)!
        ///
        /// The convention which we use is:
        ///   - the cell with the **lower global index** has the full rule
        ///   - the other cell is empty
        /// We use the global index here, so that the result is "stable" even if we are at an MPI boundary.
        /// </summary>
        public ChunkRulePair<QuadRule> VolumeQuadRule(int intOrder, int jCell, int speciesSign) {
            RefElement refElement = grddat.Cells.GetRefElement(jCell);

            //
            // Volume rule
            // -----------

            ChunkRulePair<QuadRule> volumeRule;
            // we assume that the cell is either completely empty or completely occupied
            bool EmptyOrFoolVol; // true: full rule; false: empty rule
            {
                var LsVals = levelSetData.GetLevSetValues(refElement.Center, jCell, 1);
                EmptyOrFoolVol = LsVals[0, 0] * speciesSign > 0;
            }

            if (EmptyOrFoolVol) {

                var VolRule = refElement.GetQuadratureRule(intOrder);

                volumeRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), VolRule);

            } else {
                // use empty rule
                QuadRule empty = new QuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(1),
                    Nodes = refElement.Center
                };

                volumeRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), empty);
            }
            return volumeRule;
        }
    }
}
