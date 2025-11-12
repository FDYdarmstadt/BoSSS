using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;


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
        /// conversion (<paramref name="jCellGeom"/>,<paramref name="iCoFace"/>) to grid edges
        /// </summary>
        static IEnumerable<int> GeometricalCellCoFace2Edges(IGridData gdat, int jCellGeom, int iCoFace) {

            int jCellLog;
            if(gdat.iGeomCells.GeomCell2LogicalCell != null)
                jCellLog = gdat.iGeomCells.GeomCell2LogicalCell[jCellGeom];
            else
                jCellLog = jCellGeom;


            var Kref = gdat.iGeomCells.GetRefElement(jCellGeom);
            int face0 = Kref.CoFaceToFaceIndices[iCoFace, 0];
            int face1 = Kref.CoFaceToFaceIndices[iCoFace, 1];



            return GeometricalCellFace2Edges(gdat, jCellGeom, face0).SetUnion(GeometricalCellFace2Edges(gdat, jCellGeom, face1));
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
                if(CoincidFaces != null && CoincidFaces[j] != null) {
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

                if(CoincidCoFcs != null && CoincidCoFcs[j] != null) {
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



        public SurfaceElementBounaryFactoryForCoincidingEdges(RefElement _RefElement, LevelSetTracker.LevelSetData[] __levelSetDataS, LevelSetSignCode[] __allSignCodes, int __iLevSet, (int iLevSet, int iFace)[][] __LevSetCoincidingFaces, (int iLevSet, int iFace)[][] __LevSetCoincidingCoFaces) {
            RefElement = _RefElement;
            if(Array.IndexOf(__levelSetDataS[0].GridDat.iGeomEdges.EdgeRefElements, this.RefElement) < 0) {
                throw new ArgumentException($"{_RefElement} is not an edge reference element of the grid.");
            }
            m_LevelSetDataS = __levelSetDataS;
            m_allSignCodes = __allSignCodes;
            m_iLevSet = __iLevSet;
            LevSetCoincidingFaces = __LevSetCoincidingFaces;
            LevSetCoincidingCoFaces = __LevSetCoincidingCoFaces;
        }

        public RefElement RefElement {
            get;
            private set;
        }

        readonly LevelSetTracker.LevelSetData[] m_LevelSetDataS;
        readonly LevelSetSignCode[] m_allSignCodes;
        readonly int m_iLevSet;

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
            if(mask.MaskType != MaskType.Geometrical) {
                throw new ArgumentException("expecting an geometrical mask");
            }
            var gdat = mask.GridData;
            int iKref = Array.IndexOf(gdat.iGeomEdges.EdgeRefElements, this.RefElement);

            var fulCoFaceRule = RefElement.FaceRefElement.GetQuadratureRule(order);
            var emptyFaceRule = QuadRule.CreateBlank(RefElement, 1, Math.Max(RefElement.FaceRefElement.SpatialDimension, 1), true);
            emptyFaceRule.OrderOfPrecision = order;
            emptyFaceRule.Nodes.LockForever();

            var edgeTest = RefElement.GetQuadratureRule(order);


            int[,] Edge2Cell = gdat.iGeomEdges.CellIndices;
            int J = gdat.iGeomCells.NoOfLocalUpdatedCells;
            byte[,] Edge2Face = gdat.iGeomEdges.FaceIndices;

            int[][,] CoFaceToFace = gdat.iGeomCells.RefElements.Select(Kref => Kref.CoFaceToFaceIndices).ToArray();

            int NoOfLevSets = m_LevelSetDataS.Length;


            var compRule = new ChunkRulePair<QuadRule>[mask.NoOfItemsLocally];
            int cnt = -1;
            foreach(int iEdge in mask.ItemEnum) {
                cnt++;
                if(gdat.iGeomEdges.GetRefElementIndex(iEdge) != iKref)
                    throw new ArgumentException("mask violates the element");

                bool bFound = false;
                for(int inOt = 0; inOt < 2; inOt++) { // loop over both cells bound to edge `iEdge`
                    if(bFound)
                        continue;
                    int jCell = Edge2Cell[iEdge, inOt];
                    if(jCell < 0)
                        continue;
                    if(jCell >= J)
                        continue;
                    if(!gdat.iGeomEdges.IsEdgeConformal(iEdge, inOt))
                        continue;

                    int iCellFace = Edge2Face[iEdge, inOt];

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

                            int iCoFace = tt.iCoFace;

                            int[,] cf2f = CoFaceToFace[gdat.iGeomCells.GetRefElementIndex(jCell)];

                            int inOtCoface = -1;
                            if(iCellFace == cf2f[iCoFace, 0]) {
                                inOtCoface = 0;
                            } else if(iCellFace == cf2f[iCoFace, 1]) {
                                inOtCoface = 1;
                            } else {
                                // some face that we are not interested in
                                continue;
                            }


                            // we have a edge that is NOT coinciding with the level-set,
                            // but an entire co-edge coincides with the level-set.
                            // Now, the question is, whether the quad-rule should be full or empty 

                            var KrefVol = gdat.iGeomCells.GetRefElement(jCell);
                            Debug.Assert(KrefVol.FaceRefElement == this.RefElement);

                            int TrafoIdx = gdat.iGeomEdges.Edge2CellTrafoIndex[iEdge, inOt];
                            var Edge2CellTrafo = gdat.iGeomEdges.Edge2CellTrafos[TrafoIdx];
                            var cellTest = new NodeSet(KrefVol, Edge2CellTrafo.Transform(edgeTest.Nodes), false);
                            cellTest.LockForever();

#if DEBUG
                            Vector FaceCenter = gdat.TransformLocal2Global(KrefVol.GetFaceCenter(iCellFace).GetRowPt(0), jCell);
                            Vector EdgeCenter = gdat.iGeomEdges.GetCenter(iEdge);
                            Debug.Assert(EdgeCenter.Dist(FaceCenter) < gdat.iGeomCells.h_min[jCell] * 1.0e-8, $"Some mismatch in face versus edge center ({FaceCenter} vs. {EdgeCenter})");
#endif



                            var LevelSetValues = MultidimensionalArray.Create(NoOfLevSets, cellTest.NoOfNodes);
                            for(int iLevSet = 0; iLevSet < NoOfLevSets; iLevSet++) {
                                LevelSetValues.ExtractSubArrayShallow(iLevSet, -1).Set(
                                    this.m_LevelSetDataS[iLevSet].GetLevSetValues(cellTest, jCell, 1).ExtractSubArrayShallow(0, -1)
                                    );
                            }

                            //var codes = new LevelSetSignCode[NoOfLevSets];
                            bool completelyEmty = true, completelyFull = true;
                            for(int k = 0; k < NoOfLevSets; k++) {
                                var code = LevelSetSignCode.ComputeLevelSetBytecode(LevelSetValues.GetColumn(k));
                                bool nonvoid = Array.IndexOf(this.m_allSignCodes, code) >= 0;

                                completelyFull &= nonvoid;
                                completelyEmty &= (!nonvoid);
                            }

                            if(completelyFull == completelyEmty)
                                throw new ArithmeticException("cannot decide");

                            if(completelyEmty) {
                                bFound = true;
                                compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge), emptyFaceRule);
                                break;
                            } else if(completelyFull) {
                                bFound = true;

                                int iFaceOfEdge = KrefVol.CoFaceToFaceFaceIndex[iCoFace, inOtCoface];
                                //#if DEBUG
                                //                                Vector A = this.RefElement.FaceCenters.GetRowPt(iFaceOfEdge);

                                //#endif



                                //var __Trafo = this.RefElement.GetFaceTrafo(iFace);
                                double scale = this.RefElement.FaceTrafoGramianSqrt[iFaceOfEdge];
                                //var fullRule = new NodeSet(this.RefElement, __Trafo.Transform(fulCoFaceRule.Nodes), false);
                                //fullRule.LockForever();

                                var volNodes = KrefVol.GetCoFaceTrafo(iCoFace).Transform(fulCoFaceRule.Nodes);
                                var fullRule2 = new NodeSet(this.RefElement, Edge2CellTrafo.TransformInverse(volNodes), false);
                                fullRule2.LockForever();

                                // if the `volNodes` are really in the plane defined by edge `iEdge`, the `Edge2CellTrafo` (which reduces the spatial dimension by 1) is still an identity
                                Debug.Assert(volNodes.L2Dist(Edge2CellTrafo.Transform(fullRule2)) < 1.0e-8, "Nodes probably not on the right edge.");



                                var __weights = fulCoFaceRule.Weights.CloneAs();
                                __weights.Scale(scale);


                                compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge),
                                    new QuadRule() {
                                        Nodes = fullRule2,
                                        Weights = __weights,
                                        OrderOfPrecision = fulCoFaceRule.OrderOfPrecision
                                    }
                                    );
                                break;
                            } else {
                                throw new ArithmeticException("should never reach this point");
                            }
                        }

                    }

                }

                if(bFound == false) {
                    throw new ArithmeticException("Unable to create quadrature role for edge " + iEdge + ".");
                }
            }


            return compRule;
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }
    }


}
