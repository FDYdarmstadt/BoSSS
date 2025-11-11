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
using System.Text;


namespace BoSSS.Foundation.XDG.Quadrature.LevelSetOnEdge {

    /*
    public static class Misc {


       


        /// <summary>
        /// on which edges this class can actually provide working quadrature rules?
        /// </summary>
        static public EdgeMask ComputeEdgeMask(LevelSetTracker tracker, int iLevelSet, int iKrefEdge) {
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


    }
    */


    internal class InterfaceIntegration : IQuadRuleFactory<QuadRule> {


        /// <summary>
        /// Returns all cells which have faces that conincide with the level-set
        /// </summary>
        static public CellMask ComputeCellMask(LevelSetTracker tracker, int iLevelSet, int iKref) {
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



        public InterfaceIntegration(RefElement _RefElement, LevelSetTracker.LevelSetData levelSetData) {
            RefElement = _RefElement;
            m_LevelSetData = levelSetData;
            grddat = levelSetData.GridDat;
        }

        public RefElement RefElement {
            get;
            private set;
        }

        readonly LevelSetTracker.LevelSetData m_LevelSetData;
        readonly IGridData grddat;

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }



        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if(mask.MaskType != MaskType.Geometrical) {
                throw new ArgumentException($"wrong type of mask ({mask.MaskType})");
            }

            var ret = new List<ChunkRulePair<QuadRule>>();

            foreach(int jCell in mask.ItemEnum) {
                ret.Add(SurfaceQuadRule(order, jCell));
            }

            return ret;
        }

        int GetSpecialFaceIndex(int j) {
            (int iLevSet, int iFace)[][] CoIncFaces = m_LevelSetData.Region.LevSetCoincidingFaces;
            foreach(var t in CoIncFaces[j]) {
                int levelSetIndex = m_LevelSetData.LevelSetIndex;
                if(t.iLevSet == levelSetIndex)
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
            RefElement refElement = grddat.iGeomCells.GetRefElement(jCell);
            int D = refElement.SpatialDimension;

            bool EmptyOrFool; // true: full rule; false: empty rule
            EmptyOrFool = IsNotEmpty(jCell);

            ChunkRulePair<QuadRule> surfaceRule;
            if(EmptyOrFool) {
                // +++++++++++++++++++++++++++++++++++++++++++++++
                // Level-Set integration should happen in `jCell`
                // +++++++++++++++++++++++++++++++++++++++++++++++

                var FaceRule = refElement.FaceRefElement.GetQuadratureRule(intOrder);
                int K = FaceRule.NoOfNodes;
                NodeSet VolumeNodes = new NodeSet(refElement, K, D, true);
                refElement.TransformFaceCoordinates(SpecialFace, FaceRule.Nodes, VolumeNodes);
                VolumeNodes.LockForever();

                double gTrF = refElement.FaceTrafoGramianSqrt[SpecialFace];
                
                QuadRule qr_l = new QuadRule() {
                    OrderOfPrecision = FaceRule.OrderOfPrecision,
                    Weights = MultidimensionalArray.Create(K),
                    Nodes = VolumeNodes
                };

                for(int k = 0; k < K; k++) {
                    qr_l.Weights[k] = FaceRule.Weights[k] * gTrF;
                }

                surfaceRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);

            } else {
                // ++++++++++++++++++++
                // use empty rule
                // ++++++++++++++++++++
                QuadRule empty = new QuadRule() {
                    OrderOfPrecision = intOrder,
                    Weights = MultidimensionalArray.Create(1),
                    Nodes = refElement.Center
                };

                surfaceRule = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), empty);
            }
            return surfaceRule;
        }

        private bool IsNotEmpty(int jCell) {
            int SpecialFace = GetSpecialFaceIndex(jCell);


            // determine whether this cell is inside or outside 
            // w.r.t. the edge that corresponds with face 
            int iEdge = grddat.GetEdgesForFace(jCell, SpecialFace, out int InOrOut, out int[] FurtherEdges);
            if(FurtherEdges != null && FurtherEdges.Length > 0) {
                // throw new NotSupportedException("Hanging node on a edge which coincides with the level set - this should be avoided.");
                // Console.WriteLine("Hanging node on a edge which coincides with the level set - this should be avoided.");
            }
            int J = grddat.CellPartitioning.LocalLength;

            // surface rule -- classification, using con-formality of edges
            // always use conformal cell, if both conformal choose lower global index



            bool EmptyOrFool;
            {
                int OtherCell = grddat.iGeomEdges.CellIndices[iEdge, InOrOut == 0 ? 1 : 0];
                bool ThisjConform = InOrOut == 0 ? grddat.iGeomEdges.IsEdgeConformalWithCell1(iEdge) : grddat.iGeomEdges.IsEdgeConformalWithCell2(iEdge);
                bool OtherConform = InOrOut == 0 ? grddat.iGeomEdges.IsEdgeConformalWithCell2(iEdge) : grddat.iGeomEdges.IsEdgeConformalWithCell1(iEdge);

                // if both cell faces have no hanging nodes use globally lower index
                if(ThisjConform && OtherConform) {
                    long jCellGlob = jCell + grddat.CellPartitioning.i0;
                    long OtherCellGlob = OtherCell < J ? OtherCell + grddat.CellPartitioning.i0 : grddat.iParallel.GlobalIndicesExternalCells[OtherCell - J];
                    EmptyOrFool = jCellGlob < OtherCellGlob;
                } else if(ThisjConform && !OtherConform) {
                    EmptyOrFool = true; // this cells face consists of one edge only, empty rule
                } else if(!ThisjConform && OtherConform) {
                    EmptyOrFool = false; // this cells face consists of more than one edge only, use full rule, so that each partial edge contains the rule of full degree!
                } else {
                    throw new NotSupportedException($"Error in cell pair {jCell}, {OtherCell}: Only one cell should have the hanging node.");
                }
            }

            return EmptyOrFool;
        }
    }





}
