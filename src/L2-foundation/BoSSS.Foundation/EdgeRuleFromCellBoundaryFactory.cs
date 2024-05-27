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
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using MPI.Wrappers;
using System.Threading.Tasks;
using ilPSP.Utils;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// converts a cell boundary quadrature rule into an edge quadrature rule;
    /// Therefore, the individual faces of a cell (<see cref="CellBoundaryQuadRule.NumbersOfNodesPerFace"/>) are split up into the respective edges.
    /// 
    /// Note: the conversion of cell boundary rules into edge rules implicitly assumes that the 
    /// level-set-fields is at least C^0 continuous, i.e. the level-set field is continuous across the cell boundaries.
    /// </summary>
    public class EdgeRuleFromCellBoundaryFactory : IQuadRuleFactory<QuadRule> {

        /// <summary>
        /// constructor
        /// </summary>
        public EdgeRuleFromCellBoundaryFactory(IGridData g, RefElement edgeRefElement, IEnumerable<IQuadRuleFactory<CellBoundaryQuadRule>> cellBndQF, EdgeMask maxDomain) {
            m_cellBndQF = new Dictionary<RefElement, IQuadRuleFactory<CellBoundaryQuadRule>>();
            foreach(var f in cellBndQF) {
                m_cellBndQF.Add(f.RefElement, f);
            }

            if (g.iGeomEdges.EdgeRefElements.IndexOf(edgeRefElement) < 0)
                throw new ArgumentException("edgeRefElement is not an edge ref element of the grid.");
            
            RefElement = edgeRefElement;
            grd = g;
            if (maxDomain.MaskType != MaskType.Geometrical)
                throw new ArgumentException("expecting a geometrical mask");
            m_maxDomain = maxDomain;

            if(!maxDomain.IsSubMaskOf(g.iGeomEdges.GetEdges4RefElement(edgeRefElement)))
                throw new ArgumentException("maxDomain is not a submask of the edges of the edgeRefElement.");
        }

        /// <summary>
        /// the edge simplex
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        int m_quadorder;

        public void CreateRulesAndMPIExchgange(int __quadorder) {
            //
            m_quadorder = __quadorder;
            CreateInternal(__quadorder);
        }

        IDictionary<RefElement, IQuadRuleFactory<CellBoundaryQuadRule>> m_cellBndQF;

        readonly IGridData grd;

        readonly EdgeMask m_maxDomain;

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new[] { m_quadorder };
        }


        void CreateInternal(int order) {
            // init & checks
            // =============

            var mask = m_maxDomain;
            Console.Error.WriteLine("   ----------------   EdgeRuleFromCellBoundaryFactory.GetQuadRuleSet: order = " + order);


            if (!(mask is EdgeMask))
                throw new ArgumentException("Expecting an edge mask.");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

#if DEBUG
            var maskBitMask = mask.GetBitMask();
#endif

            var Edg2Cel = this.grd.iGeomEdges.CellIndices;
            var Edg2Fac = this.grd.iGeomEdges.FaceIndices;
            int J = this.grd.iGeomCells.NoOfLocalUpdatedCells;
            QuadRule DefaultRule = this.RefElement.GetQuadratureRule(order);
            var iGeomEdges = this.grd.iGeomEdges;
            var iGeomCells = this.grd.iGeomCells;
            var iLogiCells = this.grd.iLogicalCells;

            int myIKrfeEdge = this.grd.iGeomEdges.EdgeRefElements.IndexOf(this.RefElement, (a, b) => object.ReferenceEquals(a, b));
            if (myIKrfeEdge < 0)
                throw new ApplicationException("fatal error");


            int[] EdgeIndices = mask.ItemEnum.ToArray();
            int NoEdg = EdgeIndices.Length;

            // return value
            ChunkRulePair<QuadRule>[] Ret = new ChunkRulePair<QuadRule>[NoEdg];

            // find cells
            // ==========

            int NoOfRefElements = iGeomCells.RefElements.Length;
            BitArray[] CellBitMask = NoOfRefElements.ForLoop(iKref => new BitArray(J));

            var QuadRulesForOtherProc = new Dictionary<int, List<(CellBoundaryQuadRule rule, long jGlob, int iPart, Vector Center)>>();

            int[] Cells = new int[NoEdg]; // mapping: Edge Index --> Cell Index (both geometrical)
            int[] Faces = new int[NoEdg]; // mapping: Edge Index --> Face 
            BitArray MaxDomainMask = m_maxDomain.GetBitMask();
            {
                for (int i = 0; i < NoEdg; i++) {

                    int iEdge = EdgeIndices[i];
                    if (iGeomEdges.GetRefElementIndex(iEdge) != myIKrfeEdge)
                        throw new ArgumentException("illegal edge mask");


                    int jCell0 = Edg2Cel[iEdge, 0];
                    int jCell1 = Edg2Cel[iEdge, 1];
                    int iKref0 = iGeomCells.GetRefElementIndex(jCell0);
                    int iKref1 = iGeomCells.GetRefElementIndex(jCell1);


                    bool conf0 = grd.iGeomEdges.IsEdgeConformalWithCell1(iEdge);
                    bool conf1 = grd.iGeomEdges.IsEdgeConformalWithCell2(iEdge);

                    if (!(conf0 || conf1)) {
                        throw new NotSupportedException("For an edge that is not conformal with at least one cell, no edge rule can be created from a cell boundary rule.");
                    }

                    CellBitMask[iKref0][jCell0] = true;
                    if(jCell1 < J)
                        CellBitMask[iKref1][jCell1] = true;

                    if (conf0) {
                        // use the quadrature rule from the IN-cell
                        Cells[i] = jCell0;
                        Faces[i] = Edg2Fac[iEdge, 0];
                    } else if(conf1 && jCell1 < J) {
                        // use the quadrature rule from the OUT-cell
                        Cells[i] = jCell1;
                        Faces[i] = Edg2Fac[iEdge, 1];
                    } else {
                        // the out-cell is conformal, but since it is an external cell, we must obtain it from the MPI communication
                        Debug.Assert(conf1 == true);
                        Debug.Assert(jCell1 >= J);
                        Cells[i] = jCell1;
                        Faces[i] = Edg2Fac[iEdge, 1]; 
                    }

                


                    if (conf0 == true && conf1 == false && jCell1 >= J) {
                        // +++++++++++++++++++++++++++++++++++++++++++++
                        // cell must be communicated to other processor 
                        // +++++++++++++++++++++++++++++++++++++++++++++

                        int jCell1log = grd.GetLogicalCellIndex(jCell1);
                        long jCell2glob = grd.iParallel.GlobalIndicesExternalCells[jCell1log - J];
                        int otherProcRank = grd.CellPartitioning.FindProcess(jCell2glob);

                        // index of the geometrical cell within the logical cell
                        int iPart = grd.GetGeometricalPartIndex(jCell1);

                        // center of the cell, just to check whether we really exchanged the right cell 
                        Vector CenterCheck = grd.iGeomCells.GetCenter(jCell1);


                        if(!QuadRulesForOtherProc.ContainsKey(otherProcRank)) {
                            QuadRulesForOtherProc.Add(otherProcRank, new List<(CellBoundaryQuadRule rule, long jGlob, int iPart, Vector Center)>());
                        }

                        QuadRulesForOtherProc[otherProcRank].Add((null, jCell2glob, iPart, CenterCheck));
                    }

                }
                
            }

            // get cell boundary rule
            // ======================

            CellMask[] CellMasks = CellBitMask.Select( cbm => new CellMask(this.grd, cbm, MaskType.Geometrical)).ToArray();

            // Cell Boundary Quadrature Rules for each reference element
            // - 1st index: index of the reference element
            // - 2nd index: index of the chunk
            IChunkRulePair<CellBoundaryQuadRule>[][] cellBndRuleS = NoOfRefElements.ForLoop(iKref =>
                 this.m_cellBndQF[iGeomCells.RefElements[iKref]].GetQuadRuleSet(CellMasks[iKref], order).ToArray());


            int[] jCellGeom2Chunk = new int[iGeomCells.Count];
            jCellGeom2Chunk.SetAll(-1);
            for (int iKref = 0; iKref < NoOfRefElements; iKref++) {
                var cellBndRule = cellBndRuleS[iKref];
                int L = cellBndRule.Length;

                for(int iChunk = 0; iChunk < L; iChunk++) {
                    int j0 = cellBndRule[iChunk].Chunk.i0;
                    int JE = cellBndRule[iChunk].Chunk.JE;
                    for (int jCell = j0; jCell < JE; jCell++) {
                        jCellGeom2Chunk[jCell] = iChunk;
                    }
                }
            }


            // MPI exchange
            // (required for non-conformal boundaries)
            // ========================================

            // Step 1: fill in the quadrature rules
            foreach (var kv in QuadRulesForOtherProc) {
                int iTargetRank = kv.Key;
                var list = kv.Value;
                for (int iItem = 0; iItem < list.Count; iItem++) {
                    var ttt = list[iItem];

                    int jCell = grd.CellPartitioning.TransformIndexToLocal(ttt.jGlob);
                    int jCellGeom = grd.GetGeometricCellIndex(jCell, ttt.iPart);
                    int iKref = iGeomCells.GetRefElementIndex(jCellGeom);

                    ttt.rule = cellBndRuleS[iKref][jCellGeom].Rule;

                }
            }

            // Step 2: MPI exchange
            var QuadRulesFromOtherProc = SerialisationMessenger.ExchangeData(QuadRulesForOtherProc, csMPI.Raw._COMM.WORLD);


            // Step 3: re-arange the received data
            int Jgeom = iGeomCells.NoOfLocalUpdatedCells;
            CellBoundaryQuadRule[] RulesForExternalCells = new CellBoundaryQuadRule[iGeomCells.Count - Jgeom];

            foreach (var kv in QuadRulesFromOtherProc) {
                int iSourceRank = kv.Key;
                var list = kv.Value;
                for (int iItem = 0; iItem < list.Count; iItem++) {
                    var ttt = list[iItem];
                    if(!grd.CellPartitioning.IsInLocalRange(ttt.jGlob))
                        throw new ApplicationException("fatal error; should only receive quadrature rules for external cells.");

                    int jLoc = grd.iParallel.Global2LocalIdx[ttt.jGlob];
                    int jCellGeom = grd.GetGeometricCellIndex(jLoc, ttt.iPart);


                    RulesForExternalCells[jCellGeom - Jgeom] = ttt.rule;
                }
            }



            // build rule
            // ==========
            {
                m_EdgesToChunk = new int[grd.iGeomEdges.Count];
                m_EdgesToChunk.SetAll(-1111);
                m_QuadRule = new ChunkRulePair<QuadRule>[NoEdg];


                for (int i = 0; i < NoEdg; i++) { // loop over edges in mask
                    //if (MaxDomainMask[Cells[i]] == false)
                    //    Debugger.Break();

                    int jCell = Cells[i];
                    int iFace = Faces[i];

                    CellBoundaryQuadRule CellBndR;
                    if(jCell < Jgeom) {
                        // use locally computed rule

                        CellBndR = cellBndRuleS[iGeomCells.GetRefElementIndex(jCell)][jCellGeom2Chunk[jCell]].Rule;
                    } else {
                        // use externally computed rule

                        CellBndR = RulesForExternalCells[jCell - Jgeom];
                    }

                    //if (Cells[i] >= 0) {

                    QuadRule qrEdge = null;
                    qrEdge = this.CombineQr(null, CellBndR, iFace);
                    

                    qrEdge.Nodes.LockForever();
                    qrEdge.Weights.LockForever();
                    m_QuadRule[i] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(EdgeIndices[i]), qrEdge);
                }
            }

            // return
            // ======

#if DEBUG
            for (int i = 0; i < m_QuadRule.Length; i++) {
                Chunk c = m_QuadRule[i].Chunk;
                for (int e = c.i0; e < c.JE; e++) {
                    Debug.Assert(maskBitMask[e] == true, "Internal error: rule created for edge " + e + ", which was never requested.");
                    maskBitMask[e] = false;
                }
            }

            for (int i = 0; i < maskBitMask.Length; i++) {
                Debug.Assert(maskBitMask[i] == false, "Internal error: rule for edge " + i + " was not created.");
            }
#endif
        }


        ChunkRulePair<QuadRule>[] m_QuadRule;

        /// <summary>
        /// mapping from edge index into <see cref="m_QuadRule"/>
        /// </summary>
        int[] m_EdgesToChunk;

        /// <summary>
        /// produces an edge quadrature rule
        /// </summary>
        /// <param name="mask">an edge mask</param>
        /// <param name="order">desired order</param>
        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {

            // init & checks
            // =============

            Console.Error.WriteLine("   ----------------   EdgeRuleFromCellBoundaryFactory.GetQuadRuleSet: order = " + order);


            if (!(mask is EdgeMask))
                throw new ArgumentException("Expecting an edge mask.");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

            if(!mask.IsSubMaskOf(m_maxDomain))
                throw new ArgumentException("Cannot return quadrature rule for mask larger than initially specified");
#if DEBUG
            var checkBitMask = new BitArray(grd.iGeomEdges.Count);
#endif


            var Ret = new List<ChunkRulePair<QuadRule>>();


            foreach (var chunk in mask) {


                for (int i = chunk.i0; i < chunk.JE; i++) {
                    int iQuadRule = m_EdgesToChunk[chunk.i0];
                    var pair = m_QuadRule[iQuadRule];
                    if (pair.Chunk.i0 == i && pair.Chunk.Len == (chunk.JE - i)) {
                        // the pair is an exact fit for the rest of the chunk
                        Ret.Add(pair);
                        i += pair.Chunk.Len - 1;
                    } else {
                        var newChunk = new Chunk() {
                            i0 = i,
                            Len = Math.Min(pair.Chunk.JE, chunk.JE) - i
                        };
                        Ret.Add(new ChunkRulePair<QuadRule>(newChunk, pair.Rule));
                        i += newChunk.Len - 1;
                    }
#if DEBUG
                    for (int i2 = Ret.Last().Chunk.i0; i2 < Ret.Last().Chunk.JE; i2++) {
                        Debug.Assert(checkBitMask[i2] == false, "rule for respective edge is already set.");
                        checkBitMask[i2] = true;
                    }
#endif
                }

            }


#if DEBUG
            var maskBitMask = mask.GetBitMask();
            for(int i = 0; i < grd.iLogicalEdges.Count; i++) {
                Debug.Assert(checkBitMask[i] == maskBitMask[i], "not all edges are covered by the rule.");
            }
#endif

            return Ret;
        }

        private QuadRule CombineQr(QuadRule qrEdge, CellBoundaryQuadRule givenRule, int iFace) {
            int D = grd.SpatialDimension;
            var volSplx = givenRule.Nodes.RefElement;
            int coD = D - 1;
            Debug.Assert(this.RefElement.SpatialDimension == coD);

            // extract edge rule
            // -----------------

            int i0 = 0;
            for (int i = 0; i < iFace; i++)
                i0 += givenRule.NumbersOfNodesPerFace[i];
            int iE = i0 + givenRule.NumbersOfNodesPerFace[iFace] - 1;

            if (iE < i0) {
                // rule is empty (measure is zero).

                if (qrEdge == null) {
                    QuadRule ret = new QuadRule();
                    ret.OrderOfPrecision = int.MaxValue - 1;
                    ret.Nodes = new NodeSet(this.RefElement, 1, Math.Max(1, D - 1), false);
                    ret.Weights = MultidimensionalArray.Create(1);  // this is an empty rule, since the weight is zero!
                    // (rules with zero nodes may cause problems at various places.)
                    return ret;
                } else {
                    qrEdge.Nodes.Scale(0.5);
                    qrEdge.Weights.Scale(0.5);
                    return qrEdge;
                }
            }

            MultidimensionalArray NodesVol = givenRule.Nodes.ExtractSubArrayShallow(new int[] { i0, 0 }, new int[] { iE, D - 1 });
            MultidimensionalArray Weigts = givenRule.Weights.ExtractSubArrayShallow(new int[] { i0 }, new int[] { iE }).CloneAs();
            NodeSet Nodes = new NodeSet(this.RefElement, iE - i0 + 1, coD, qrEdge == null);

            volSplx.GetInverseFaceTrafo(iFace).Transform(NodesVol, Nodes);
            Nodes.LockForever();

            //Debug.Assert((Weigts.Sum() - grd.Grid.GridSimplex.EdgeSimplex.Volume).Abs() < 1.0e-6, "i've forgotten the gramian");

            // combine 
            // -------
            if (qrEdge == null) {
                // no rule defined yet - just set the one we have got
                // ++++++++++++++++++++++++++++++++++++++++++++++++++
                qrEdge = new QuadRule();
                qrEdge.Weights = Weigts;
                qrEdge.Nodes = Nodes;
                qrEdge.OrderOfPrecision = givenRule.OrderOfPrecision;

            } else {
                // take the mean of already defined and new rule
                // +++++++++++++++++++++++++++++++++++++++++++++

                int L1 = qrEdge.Nodes.GetLength(0);
                int L2 = Nodes.GetLength(0);
                Debug.Assert(coD == qrEdge.Nodes.GetLength(1));

                NodeSet newNodes = new NodeSet(this.RefElement, L1 + L2, coD, true);
                newNodes.SetSubArray(qrEdge.Nodes, new int[] { 0, 0 }, new int[] { L1 - 1, coD - 1 });
                newNodes.SetSubArray(Nodes, new int[] { L1, 0 }, new int[] { L1 + L2 - 1, coD - 1 });
                newNodes.LockForever();

                MultidimensionalArray newWeights = MultidimensionalArray.Create(L1 + L2);
                newWeights.AccSubArray(0.5, qrEdge.Weights, new int[] { 0 }, new int[] { L1 - 1 });
                newWeights.AccSubArray(0.5, Weigts, new int[] { L1 }, new int[] { L1 + L2 - 1 });

                double oldSum = qrEdge.Weights.Sum();
                double newSum = Weigts.Sum();
                WeightInbalance += Math.Abs(oldSum - newSum);


                qrEdge.Nodes = newNodes;
                qrEdge.Weights = newWeights;
                qrEdge.OrderOfPrecision = Math.Min(qrEdge.OrderOfPrecision, givenRule.OrderOfPrecision);

            }

            // return
            // ------
            return qrEdge;

        }

        /*
        // Basically the same, but drop nodes, that are not located on the edge!
        private QuadRule CombineQrNonConformal(QuadRule qrEdge, CellBoundaryQuadRule givenRule, int iFace, int iEdge, int jCell) {
            int D = grd.SpatialDimension;
            var volSplx = m_cellBndQF.RefElement;
            int coD = D - 1;
            Debug.Assert(this.RefElement.SpatialDimension == coD);


            //if (grd.MpiRank == 1 && this.m_cellBndQF.GetType().ToString().ToLowerInvariant().Contains("pointqrf")) {
            //    Debugger.Launch();
            //}


            // extract edge rule
            // -----------------

            int i0 = 0, iE = 0;
            for (int i = 0; i < iFace; i++)
                i0 += givenRule.NumbersOfNodesPerFace[i];
            iE = i0 + givenRule.NumbersOfNodesPerFace[iFace] - 1;

            if (iE < i0) {
                // rule is empty (measure is zero).

                if (qrEdge == null) {
                    QuadRule ret = new QuadRule();
                    ret.OrderOfPrecision = int.MaxValue - 1;
                    ret.Nodes = new NodeSet(this.RefElement, 1, Math.Max(1, D - 1), false);
                    ret.Weights = MultidimensionalArray.Create(1);  // this is an empty rule, since the weight is zero!
                    // (rules with zero nodes may cause problems at various places.)
                    return ret;
                } else {
                    qrEdge.Nodes.Scale(0.5);
                    qrEdge.Weights.Scale(0.5);
                    return qrEdge;
                }
            }

            // dbg_launch();

            MultidimensionalArray NodesVol = givenRule.Nodes.ExtractSubArrayShallow(new int[] { i0, 0 }, new int[] { iE, D - 1 });
            MultidimensionalArray Weigts = givenRule.Weights.ExtractSubArrayShallow(new int[] { i0 }, new int[] { iE }).CloneAs();
            NodeSet Nodes = new NodeSet(this.RefElement, iE - i0 + 1, coD, false);

            // as before, nodes from cell transformed to face
            volSplx.GetInverseFaceTrafo(iFace).Transform(NodesVol, Nodes);
            Nodes.LockForever();

            // get edge endpoints 
            var EdgeNodes = this.RefElement.Vertices;
            NodeSet FaceNodes = new NodeSet(this.RefElement, this.RefElement.Vertices.NoOfNodes, coD, false);
            // Transform to cell
            int trf;
            if (grd.Edges.CellIndices[iEdge, 0] == jCell) {
                trf = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
            } else {
                trf = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];
            }
            var CellNodes = EdgeNodes.GetVolumeNodeSet(grd, trf, false);
            // Transform to face
            volSplx.GetInverseFaceTrafo(iFace).Transform(CellNodes, FaceNodes);

            // build transformation from face to edge
            var Trafo = AffineTrafo.FromPoints(FaceNodes.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { this.RefElement.SpatialDimension, this.RefElement.SpatialDimension -1 }), EdgeNodes.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { this.RefElement.SpatialDimension, this.RefElement.SpatialDimension - 1 }));
            double scale = Trafo.Matrix.Determinant();

            // construct the NodeSet in edge local coordinates
            NodeSet NodesEdge = new NodeSet(this.RefElement, Trafo.Transform(Nodes), false);

            double[] NodesOnEdge = new double[0];
            double[] WeightsOnEdge = new double[0];
            int NoOfNodes = 0;
            for (int i = 0; i < iE - i0 + 1; i++) {
                // remove nodes, that are not located on this edge
                double[] point = NodesEdge.ExtractSubArrayShallow(i, -1).To1DArray();
                
                double weight = Weigts[i] * scale;
                if (this.RefElement.IsWithin(point)) {
                    NodesOnEdge = NodesOnEdge.Concat(point).ToArray();
                    WeightsOnEdge = WeightsOnEdge.Concat(weight).ToArray();
                    NoOfNodes++;
                }
            }

            // create empty if none of the points lies on this edge
            if(NoOfNodes == 0) {
                NodesOnEdge = NodesOnEdge.Concat(this.RefElement.Center.ExtractSubArrayShallow(0,-1).To1DArray()).ToArray();
                WeightsOnEdge = WeightsOnEdge.Concat(0.0).ToArray();
                NoOfNodes++;
            }

            // override the Node set with the "true" nodes located on the edge and in edge local coordinates
            Nodes = new NodeSet(this.RefElement, MultidimensionalArray.CreateWrapper(NodesOnEdge, new int[] { NoOfNodes, this.RefElement.SpatialDimension }), qrEdge == null);
            Weigts = MultidimensionalArray.CreateWrapper(WeightsOnEdge, new int[] { NoOfNodes });
            Nodes.LockForever();

            //Debug.Assert((Weigts.Sum() - grd.Grid.GridSimplex.EdgeSimplex.Volume).Abs() < 1.0e-6, "i've forgotten the gramian");

            // combine 
            // -------
            if (qrEdge == null) {
                // no rule defined yet - just set the one we have got
                // ++++++++++++++++++++++++++++++++++++++++++++++++++
                qrEdge = new QuadRule();
                qrEdge.Weights = Weigts;
                qrEdge.Nodes = Nodes;
                qrEdge.OrderOfPrecision = givenRule.OrderOfPrecision;

            } else {
                // take the mean of already defined and new rule
                // +++++++++++++++++++++++++++++++++++++++++++++

                int L1 = qrEdge.Nodes.GetLength(0);
                int L2 = Nodes.GetLength(0);
                Debug.Assert(coD == qrEdge.Nodes.GetLength(1));

                NodeSet newNodes = new NodeSet(this.RefElement, L1 + L2, coD, true);
                newNodes.SetSubArray(qrEdge.Nodes, new int[] { 0, 0 }, new int[] { L1 - 1, coD - 1 });
                newNodes.SetSubArray(Nodes, new int[] { L1, 0 }, new int[] { L1 + L2 - 1, coD - 1 });
                newNodes.LockForever();

                MultidimensionalArray newWeights = MultidimensionalArray.Create(L1 + L2);
                newWeights.AccSubArray(0.5, qrEdge.Weights, new int[] { 0 }, new int[] { L1 - 1 });
                newWeights.AccSubArray(0.5, Weigts, new int[] { L1 }, new int[] { L1 + L2 - 1 });

                double oldSum = qrEdge.Weights.Sum();
                double newSum = Weigts.Sum();
                WeightInbalance += Math.Abs(oldSum - newSum);


                qrEdge.Nodes = newNodes;
                qrEdge.Weights = newWeights;
                qrEdge.OrderOfPrecision = Math.Min(qrEdge.OrderOfPrecision, givenRule.OrderOfPrecision);

            }

            // return
            // ------
            return qrEdge;

        }
        */

        /// <summary>
        /// total difference between quadrature weights on both sides of an edge
        /// </summary>
        public double WeightInbalance {
            get;
            private set;
        }

       
    }
}
