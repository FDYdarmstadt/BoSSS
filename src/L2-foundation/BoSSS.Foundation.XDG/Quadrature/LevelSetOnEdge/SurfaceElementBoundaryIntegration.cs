using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature.LevelSetOnEdge {
    
    
    internal class SurfaceElementBoundaryIntegration : IQuadRuleFactory<QuadRule> {


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
        static public EdgeMask ComputeMask(LevelSetTracker.LevelSetRegions regions, int iLevelSet, int iKrefEdge) {
            var CoincidFaces = regions.LevSetCoincidingFaces;
            var CoincidCoFcs = regions.LevSetCoincidingCoFaces;
            IGridData gdat = regions.GridDat;

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



        public SurfaceElementBoundaryIntegration(RefElement _RefElement, LevelSetTracker.LevelSetData[] __levelSetDataS, LevelSetSignCode[] __allSignCodes, int __iLevSet, (int iLevSet, int iFace)[][] __LevSetCoincidingFaces, (int iLevSet, int iFace)[][] __LevSetCoincidingCoFaces) {
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
            var emptyFaceRule = QuadRule.CreateBlank(RefElement, 1, Math.Max(RefElement.SpatialDimension, 1), true);
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

                if(!gdat.iGeomEdges.IsEdgeConformal(iEdge, 0)) {
                    // in-cell is non-conformal => hope for the out-cell

                    int jCell2 = Edge2Cell[iEdge, 1];
                    if(jCell2 < 0 || jCell2 >= J)
                        throw new NotSupportedException("unable to obtain quadrature rule, since the inner ");
                }

                    
                    
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
                            // entire face coincides with the level-set
                            // in this case, we don't want any co-co-dim quadrature
                            // 
                            //           coinciding face
                            //    =====*=================*===== Level-Set
                            //         |                 |
                            //         |                 |
                            //         |    Cell         |
                            //         |                 |
                            //         -------------------
                            //
                            //    the edge elements (*) should be assigned to the 
                            //    two VERTICAL edges
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++
                            bFound = true;
                            compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge), emptyFaceRule);
                            break;
                        }
                    }

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
                            //
                            //       
                            //
                            //                   Level-Set       
                            //       Species B   /               
                            //           -------*   Species A             
                            //           |     /|                
                            //           |    / |     
                            //           |   /  |     Whether an edge is assigned an empty rule or non-empty rule is determined on basis of the species:           
                            //           |  /   |      Note that this factory is always created for some particular species;
                            //           | /    |      If the edge is within the domain of the respective species, the quadrature rule will be non-empty.
                            //   cell 2  |/ cl 3|      If the edge is on the opposite side, an empty quadrature rule will be assigned.
                            //      ---- *------|                
                            //   cell 1 /|cell 4|      This means: for species A, the coupling is: cell 1 <-> cell 4 <-> cell 3
                            //         / |      |                  for species A, the coupling is: cell 1 <-> cell 2 <-> cell 3
                            //        /  |-------
                            //
                            //   
                            //

                            

                            var KrefVol = gdat.iGeomCells.GetRefElement(jCell);
                            Debug.Assert(KrefVol.FaceRefElement == this.RefElement);
#if DEBUG
                            Vector FaceCenter = gdat.TransformLocal2Global(KrefVol.GetFaceCenter(iCellFace).GetRowPt(0), jCell);
                            Vector EdgeCenter = gdat.iGeomEdges.GetCenter(iEdge);
                            Debug.Assert(EdgeCenter.Dist(FaceCenter) < gdat.iGeomCells.h_min[jCell] * 1.0e-8, $"Some mismatch in face versus edge center ({FaceCenter} vs. {EdgeCenter})");
#endif
                            int TrafoIdx = gdat.iGeomEdges.Edge2CellTrafoIndex[iEdge, inOt];
                            var Edge2CellTrafo = gdat.iGeomEdges.Edge2CellTrafos[TrafoIdx];


                            NodeSet GetCellNodes(MultidimensionalArray edgeNodes) {
                                var cellTest = new NodeSet(KrefVol, Edge2CellTrafo.Transform(edgeNodes), false);
                                cellTest.LockForever();

                                return cellTest;
                            }
                          

                            void CheckLevelSets(out bool _completelyEmty, out bool _completelyFull, NodeSet cellTest, double LsEps = 0) {
                                var LevelSetValues = MultidimensionalArray.Create(NoOfLevSets, cellTest.NoOfNodes);
                                for(int iLevSet = 0; iLevSet < NoOfLevSets; iLevSet++) {
                                    LevelSetValues.ExtractSubArrayShallow(iLevSet, -1).Set(
                                        this.m_LevelSetDataS[iLevSet].GetLevSetValues(cellTest, jCell, 1).ExtractSubArrayShallow(0, -1)
                                        );
                                }

                                _completelyEmty = true;
                                _completelyFull = true;
                                for(int k = 0; k < cellTest.NoOfNodes; k++) {
                                    var code = LevelSetSignCode.ComputeLevelSetBytecode(LevelSetValues.GetColumn(k));
                                    bool nonvoid = Array.IndexOf(this.m_allSignCodes, code) >= 0;

                                    if(LsEps != 0) {
                                        var LsModVals = LevelSetValues.GetColumn(k);
                                        LsModVals[m_iLevSet] += LsEps;
                                        var code3 = LevelSetSignCode.ComputeLevelSetBytecode(LsModVals);
                                        nonvoid |= Array.IndexOf(this.m_allSignCodes, code3) >= 0;
                                        
                                        LsModVals[m_iLevSet] -= 2*LsEps;
                                        var code4 = LevelSetSignCode.ComputeLevelSetBytecode(LsModVals);
                                        nonvoid |= Array.IndexOf(this.m_allSignCodes, code4) >= 0;
                                    }

                                    _completelyFull &= nonvoid;
                                    _completelyEmty &= (!nonvoid);
                                }
                            }


                            CheckLevelSets(out bool completelyFull, out bool completelyEmty, GetCellNodes(edgeTest.Nodes));


                            if(completelyFull && completelyEmty) {
                                // ++++++++++++++++++++++++++
                                // something really weird
                                // ++++++++++++++++++++++++++
                                //EdgeMask em = new EdgeMask(gdat, Chunk.GetSingleElementChunk(iEdge), MaskType.Logical);
                                //em.SaveToTextFile("fubar2.csv", WriteHeader: false);
                                throw new ArithmeticException("cannot decide");
                            } else if(!completelyEmty && !completelyFull) {
                                // ++++++++++++++++++++++++++
                                // double-cut 
                                // should only happen in 3D, so the Co-Face element is a Line
                                // this line must be split up into empty and non-empty parts
                                // ++++++++++++++++++++++++++

                                if(this.RefElement.FaceRefElement != Line.Instance)
                                    throw new ApplicationException();


                                var LineInVol = KrefVol.GetCoFaceTrafo(iCoFace).Transform(this.RefElement.FaceRefElement.Vertices);
                                //var global = gdat.TransformLocal2Global(LineInVol, jCell);
                                //global.SaveToTextFile("line.csv");

                                if(this.m_LevelSetDataS.Length != 2)
                                    throw new NotSupportedException("should not happen with one level-set, more than 2 level-sets currently not supported.");
                                ILevelSet otherLevelSet;
                                if(this.m_iLevSet == 0)
                                    otherLevelSet = m_LevelSetDataS[1].LevelSet;
                                else
                                    otherLevelSet = m_LevelSetDataS[0].LevelSet;

                                Vector lineStart = LineInVol.GetRowPt(0);
                                Vector lineEnd = LineInVol.GetRowPt(1);

                                var FullCoFace = new HMF.LineSegment(3, KrefVol, lineStart, lineEnd);
                                FullCoFace.ProjectBasisPolynomials((otherLevelSet as LevelSet).Basis);
                                var roots = FullCoFace.GetRoots(otherLevelSet, jCell, gdat.iGeomCells.GetRefElementIndex(jCell));

                                var parts = FullCoFace.Split(roots);
                                
                                bool IsPartNonEmpty(LineSegment ls) {
                                    var center = 0.5 * (ls.Start + ls.End);
                                    CheckLevelSets(out _, out bool nonVoid, new NodeSet(KrefVol, center, false), 1.0);
                                    return nonVoid;
                                }

                                var nonEmpty_parts = parts.Where(IsPartNonEmpty).ToArray();
                                if(nonEmpty_parts.Length <= 0) {
                                    compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge), emptyFaceRule);
                                    bFound = true;
                                    Debug.Assert(compRule[cnt].Rule.SpatialDim == this.RefElement.SpatialDimension, "produced a quadrature rule with wrong spatial dimension");
                                    break;
                                } 

                                MultidimensionalArray volNodes = MultidimensionalArray.Create(nonEmpty_parts.Length*fulCoFaceRule.NoOfNodes, KrefVol.SpatialDimension);
                                MultidimensionalArray __weights = MultidimensionalArray.Create(volNodes.GetLength(0));

                                //gdat.TransformLocal2Global(LineInVol, jCell).SaveToTextFile("lineInVol.csv");
                                int k = 0;
                                foreach(var part in nonEmpty_parts) {
                                    double metric = Math.Abs(part.EndCoord - part.StartCoord)/2.0;

                                    for(int i = 0; i < fulCoFaceRule.NoOfNodes; i++) {
                                        double alpha = fulCoFaceRule.Nodes[i, 0];
                                        double beta = (alpha + 1) * (part.EndCoord - part.StartCoord) / 2.0 + part.StartCoord;

                                        var volNode = (beta + 1) * (lineEnd - lineStart) * 0.5 + lineStart;

                                        volNodes.SetRowPt(k, volNode);
                                        __weights[k] = fulCoFaceRule.Weights[k] * metric;

                                        k++;
                                    }
                                }

                                //gdat.TransformLocal2Global(volNodes, jCell).SaveToTextFile("line.csv");

                                var fullRule3 = new NodeSet(this.RefElement, Edge2CellTrafo.TransformInverse(volNodes), false);
                                fullRule3.LockForever();
                                compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge),
                                    new QuadRule() {
                                        Nodes = fullRule3,
                                        Weights = __weights,
                                        OrderOfPrecision = fulCoFaceRule.OrderOfPrecision
                                    });
                                Debug.Assert(compRule[cnt].Rule.SpatialDim == this.RefElement.SpatialDimension, "produced a quadrature rule with wrong spatial dimension");
                                bFound = true;
                                break;
                            } else if(completelyEmty) {
                                // ++++++++++++++++++++++++++
                                // empty - assign empty rule
                                // ++++++++++++++++++++++++++


                                compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge), emptyFaceRule);
                                Debug.Assert(compRule[cnt].Rule.SpatialDim == this.RefElement.SpatialDimension, "produced a quadrature rule with wrong spatial dimension");
                                bFound = true;
                                break;
                            } else {
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // not empty - assign boundary integral of surface element to this edge
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                                int iFaceOfEdge = KrefVol.CoFaceToFaceFaceIndex[iCoFace, inOtCoface];
                                double scale = this.RefElement.FaceTrafoGramianSqrt[iFaceOfEdge];

                                var volNodes = KrefVol.GetCoFaceTrafo(iCoFace).Transform(fulCoFaceRule.Nodes);
                                var fullRule2 = new NodeSet(this.RefElement, Edge2CellTrafo.TransformInverse(volNodes), false);
                                fullRule2.LockForever();

                                // if the `volNodes` are really in the plane defined by edge `iEdge`,
                                // the `Edge2CellTrafo` (which reduces the spatial dimension by 1) is still an identity
                                //Debug.Assert(volNodes.L2Dist(Edge2CellTrafo.Transform(fullRule2)) < 1.0e-8, "Nodes probably not on the right edge.");
                                if(volNodes.L2Dist(Edge2CellTrafo.Transform(fullRule2)) >= 1.0e-8) {
                                    throw new ArithmeticException("Nodes probably not on the right edge.");
                                }

                                var __weights = fulCoFaceRule.Weights.CloneAs();
                                __weights.Scale(scale);


                                compRule[cnt] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(iEdge),
                                    new QuadRule() {
                                        Nodes = fullRule2,
                                        Weights = __weights,
                                        OrderOfPrecision = fulCoFaceRule.OrderOfPrecision
                                    }
                                    );
                                Debug.Assert(compRule[cnt].Rule.SpatialDim == this.RefElement.SpatialDimension, "produced a quadrature rule with wrong spatial dimension");
                                bFound = true;
                                break;
                            } 
                            //else {
                            //    throw new ArithmeticException("should never reach this point");
                            //}
                        }

                    }

                }

                if(bFound == false) {
                    //(new EdgeMask(gdat, Chunk.GetSingleElementChunk(iEdge), MaskType.Geometrical)).SaveToTextFile("fubar.csv", false);
                    throw new ArithmeticException("Unable to create quadrature rule for edge " + iEdge + ".");
                }

                Debug.Assert(compRule[cnt].Rule.SpatialDim == this.RefElement.SpatialDimension, "produced a quadrature rule with wrong spatial dimension");
            }



            return compRule;
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }
    }

}
