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

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// converts a cell boundary quadrature rule into an edge quadrature rule
    /// </summary>
    public class EdgeRuleFromCellBoundaryFactory : IQuadRuleFactory<QuadRule> {

        /// <summary>
        /// constructor
        /// </summary>
        public EdgeRuleFromCellBoundaryFactory(Grid.Classic.GridData g, IQuadRuleFactory<CellBoundaryQuadRule> cellBndQF, CellMask maxDomain) {
            m_cellBndQF = cellBndQF;
            grd = g;
            m_maxDomain = maxDomain;
        }

        IQuadRuleFactory<CellBoundaryQuadRule> m_cellBndQF;

        Grid.Classic.GridData grd;

        CellMask m_maxDomain;

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return this.m_cellBndQF.GetCachedRuleOrders();
        }

        
        /// <summary>
        /// produces an edge quadrature rule
        /// </summary>
        /// <param name="mask">an edge mask</param>
        /// <param name="order">desired order</param>
        /// <returns></returns>
        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {

            // init & checks
            // =============

            if (!(mask is EdgeMask))
                throw new ArgumentException("Expecting an edge mask.");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

#if DEBUG
            var maskBitMask = mask.GetBitMask();

#endif

            var Edg2Cel = this.grd.iGeomEdges.CellIndices;
            var Edg2Fac = this.grd.iGeomEdges.FaceIndices;
            int J = this.grd.Cells.NoOfLocalUpdatedCells;
            QuadRule DefaultRule = this.RefElement.GetQuadratureRule(order); ;

            int myIKrfeEdge = this.grd.Edges.EdgeRefElements.IndexOf(this.RefElement, (a, b) => object.ReferenceEquals(a, b));
            if (myIKrfeEdge < 0)
                throw new ApplicationException("fatal error");


            int[] EdgeIndices = mask.ItemEnum.ToArray();
            int NoEdg = EdgeIndices.Length;

            // return value
            ChunkRulePair<QuadRule>[] Ret = new ChunkRulePair<QuadRule>[NoEdg];

            // find cells
            // ==========

            BitArray CellBitMask = new BitArray(J);
            int[] Cells = new int[NoEdg]; // mapping: Edge Index --> Cell Index (both geometrical)
            int[] Faces = new int[NoEdg]; // mapping: Edge Index --> Face 
            BitArray MaxDomainMask = m_maxDomain.GetBitMask();
            {
                for (int i = 0; i < NoEdg; i++) {

                    int iEdge = EdgeIndices[i];
                    if (this.grd.Edges.GetRefElementIndex(iEdge) != myIKrfeEdge)
                        throw new ArgumentException("illegal edge mask");

                    if (!(grd.Edges.IsEdgeConformalWithCell1(iEdge) || grd.Edges.IsEdgeConformalWithCell2(iEdge))) {
                        throw new NotSupportedException("For an edge that is not conformal with at least one cell, no edge rule can be created from a cell boundary rule.");
                    }

                    int jCell0 = Edg2Cel[iEdge, 0];
                    int jCell1 = Edg2Cel[iEdge, 1];
                    bool conf0 = grd.Edges.IsEdgeConformalWithCell1(iEdge);
                    bool conf1 = grd.Edges.IsEdgeConformalWithCell2(iEdge);

                    // this gives no errors for surface elements in 3D
                    bool Allow0 = MaxDomainMask[jCell0];
                    bool Allow1 = (jCell1 >= 0 && jCell1 < J) ? MaxDomainMask[jCell1] : false;

                    // //this is required for MPI parallel calculations
                    //bool Allow0 = true;// AllowedCells[jCell0];
                    //bool Allow1 = (jCell1 >= 0 && jCell1 < J);// ? AllowedCells[jCell1] : false;

                    if (!Allow0 && !Allow1) {
                        // fallback onto default rule, if allowed

                        //if (this.m_DefaultRuleFallbackAllowed) {
                        //    Cells[i] = -1; // by a negative index, we mark that we take the default rule
                        //    Faces[i] = -1;
                        //    Ret[i] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(EdgeIndices[i]), DefaultRule);


                        //} else {
                        throw new ArgumentException("unable to find a cell from which the edge rule can be taken.");
                        //}
                    } else {
                        Debug.Assert(Allow0 || Allow1);

                        if (conf0 && Allow0) {
                            // cell 0 is allowed and conformal:
                            // take this, it won't get better

                            CellBitMask[jCell0] = true;
                            Faces[i] = Edg2Fac[iEdge, 0] + 333; // arbitrary shift, as we need the sign to distinguish a case later on, otherwise Face = 0 is undefined and may lead to wrong results
                            Cells[i] = jCell0;
                        } else if (conf1 && Allow1) {
                            // cell 1 is allowed and conformal:
                            // take this, it won't get better

                            CellBitMask[jCell1] = true;
                            Faces[i] = Edg2Fac[iEdge, 1] + 333;
                            Cells[i] = jCell1;

                        } else if (Allow0) {
                            // cell 0 is allowed, but NOT conformal

                            CellBitMask[jCell0] = true;
                            Faces[i] = -(Edg2Fac[iEdge, 0] + 333); // by a negative index, we mark a non-conformal edge
                            Cells[i] = jCell0;
                        } else if (Allow1) {
                            // cell 1 is allowed, but NOT conformal

                            CellBitMask[jCell1] = true;
                            Faces[i] = -(Edg2Fac[iEdge, 1] + 333); // by a negative index, we mark a non-conformal edge
                            Cells[i] = jCell1;
                        }

                    }
                }
            }


            // get cell boundary rule
            // ======================
            var CellMask = new CellMask(this.grd, CellBitMask, MaskType.Geometrical);


            IChunkRulePair<CellBoundaryQuadRule>[] cellBndRule = this.m_cellBndQF.GetQuadRuleSet(CellMask, order).ToArray();
            int[] jCell2PairIdx = new int[J];
            for (int i = 0; i < cellBndRule.Length; i++) {
                var chk = cellBndRule[i].Chunk; // cell chunk
                for (int jCell = chk.JE - 1; jCell >= chk.i0; jCell--) {
                    jCell2PairIdx[jCell] = i + 555;
                }
            }

            int[] iChunk = new int[NoEdg]; // which chunk for which edge?
            for (int i = 0; i < NoEdg; i++) {
                if (Cells[i] >= 0) {
                    iChunk[i] = jCell2PairIdx[Cells[i]] - 555;
                } else {
                    iChunk[i] = int.MinValue;
                }
            }
            // build rule
            // ==========
            {
                for (int i = 0; i < NoEdg; i++) { // loop over edges
                    //if (MaxDomainMask[Cells[i]] == false)
                    //    Debugger.Break();

                    //if (Cells[i] >= 0) {
                    var CellBndR = cellBndRule[iChunk[i]].Rule;
                    QuadRule qrEdge = null;

                    if(Faces[i] >= 0)
                        qrEdge = this.CombineQr(null, CellBndR, Faces[i] - 333);
                    else {                       
                        qrEdge = this.CombineQrNonConformal(null, CellBndR, -Faces[i] - 333, EdgeIndices[i], Cells[i]);                        
                    }

                    qrEdge.Nodes.LockForever();
                    qrEdge.Weights.LockForever();
                    Ret[i] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(EdgeIndices[i]), qrEdge);
                    //} else {
                    //    Debug.Assert(Ret[i] != null);
                    //}
                }
            }

            // return
            // ======

#if DEBUG
            for (int i = 0; i < Ret.Length; i++) {
                Chunk c = Ret[i].Chunk;
                for (int e = c.i0; e < c.JE; e++) {
                    Debug.Assert(maskBitMask[e] == true);
                }
            }
#endif

            return Ret;

            
        }

        private QuadRule CombineQr(QuadRule qrEdge, CellBoundaryQuadRule givenRule, int iFace) {
            int D = grd.SpatialDimension;
            var volSplx = m_cellBndQF.RefElement;
            int coD = D - 1;
            Debug.Assert(this.RefElement.SpatialDimension == coD);

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

        // Basically the same, but drop nodes, that are not located on the edge!
        private QuadRule CombineQrNonConformal(QuadRule qrEdge, CellBoundaryQuadRule givenRule, int iFace, int iEdge, int jCell) {
            int D = grd.SpatialDimension;
            var volSplx = m_cellBndQF.RefElement;
            int coD = D - 1;
            Debug.Assert(this.RefElement.SpatialDimension == coD);

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

        /// <summary>
        /// total difference between quadrature weights on both sides of an edge
        /// </summary>
        public double WeightInbalance {
            get;
            private set;
        }

        /// <summary>
        /// the edge simplex
        /// </summary>
        public RefElement RefElement {
            get {
                return m_cellBndQF.RefElement.FaceRefElement;
            }
        }
    }
}
