using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation.Grid;
using System.Collections;
using MPI.Wrappers;
using System.Diagnostics;
using BoSSS.Platform;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// converts a cell boundary quadrature rule into an edge quadrature rule
    /// </summary>
    public class DoubleEdgeRuleFromCellBoundaryFactory : IQuadRuleFactory<DoubleEdgeQuadRule> {

        /// <summary>
        /// constructor
        /// </summary>
        public DoubleEdgeRuleFromCellBoundaryFactory(GridData g, IQuadRuleFactory<CellBoundaryQuadRule> cellBndQF) {
            m_cellBndQF = cellBndQF;
            grd = g;
        }
        
        IQuadRuleFactory<CellBoundaryQuadRule> m_cellBndQF;

        GridData grd;

        /// <summary>
        /// produces an edge quadrature rule
        /// </summary>
        /// <param name="mask">an edge mask</param>
        /// <param name="order">desired order</param>
        /// <returns></returns>
        public IEnumerable<IChunkRulePair<DoubleEdgeQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (!(mask is EdgeMask))
                throw new ArgumentException();

            EdgeMask edgMask = mask as EdgeMask;
            var ret = new Dictionary<Chunk, DoubleEdgeQuadRule>(mask.NoOfItemsLocally);
#if DEBUG
            var EdgesOfInterest = edgMask.GetBitMask();
            var EdgeTouched = (BitArray)EdgesOfInterest.Clone();
#endif

            WeightInbalance = 0.0;

            // find all cells that are 'touched' by the edge mask
            // --------------------------------------------------
            int J = grd.NoOfLocalUpdatedCells;
            BitArray TouchedCells = new BitArray(J, false);

            var Edg2Cell = grd.Edges;

            int chunkCnt = 0;
            foreach (var chnk in mask) {
                int EE = chnk.JE;
                for (int e = chnk.i0; e < EE; e++) {
                    int j1 = Edg2Cell[e, 0], j2 = Edg2Cell[e, 1];
                    TouchedCells[j1] = true;
                    if (j2 >= 0)
                        TouchedCells[j2] = true;


                    Chunk singleEdgeChunk;
                    singleEdgeChunk.i0 = e;
                    singleEdgeChunk.Len = 1;

                    ret.Add(singleEdgeChunk, null);
                }
                chunkCnt++;
            }


            CellMask celMask = (new CellMask(grd, TouchedCells));

            // create cell boundary rule!
            IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> cellBndRule = m_cellBndQF.GetQuadRuleSet(celMask, order);

            // do MPI communication (get rules for external cells)
            {
                int size, rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (size > 1)
                    throw new NotSupportedException("currently no MPI support");
            }

            // assign the cell boundary rule to edges
            // --------------------------------------
            var volSplx = m_cellBndQF.Simplex;
            int NoOfFaces = volSplx.NoOfEdges;
            var Cells2Edge = grd.LocalCellIndexToEdges;

            foreach (var kv in cellBndRule) { // loop over cell chunks (in the cell boundary rule)...
                Chunk chk = kv.Chunk;
                CellBoundaryQuadRule qr = kv.Rule;

                int JE = chk.JE;
                for (int j = chk.i0; j < JE; j++) { // loop over all cells in chunk...
                    Debug.Assert(qr.NumbersOfNodesPerEdge.Length == NoOfFaces);

                    for (int e = 0; e < NoOfFaces; e++) { // loop over faces of cell...

                        //if (qr.NumbersOfNodesPerEdge[e] <= 0)
                        //    // no contribution from this edge
                        //    continue;

                        

                        int iEdge = Math.Abs(Cells2Edge[j, e]) - 1;

                        Chunk singleEdgeChunk = Chunk.GetSingleElementChunk(iEdge);

                        DoubleEdgeQuadRule qrEdge;
                        if (ret.TryGetValue(singleEdgeChunk, out qrEdge)) {
                            // we are interested in this edge!
                            
#if DEBUG
                            Debug.Assert(EdgesOfInterest[iEdge] == true);
                            EdgeTouched[iEdge] = false;
                            
                            var vtx = this.Simplex.Vertices;
                            MultidimensionalArray _vtx = MultidimensionalArray.Create(vtx.GetLength(0), vtx.GetLength(1));
                            _vtx.SetA2d(vtx);
                            var VolSimplex = m_cellBndQF.Simplex;

                            var RefCoord = MultidimensionalArray.Create(vtx.GetLength(0), vtx.GetLength(1) + 1);
                            var PhysCoord = MultidimensionalArray.Create(1, vtx.GetLength(0), vtx.GetLength(1) + 1);

                            VolSimplex.EdgeToVolumeCoordinates(e, _vtx, RefCoord);
                            grd.TransformLocal2Global(RefCoord, PhysCoord, j, 1, 0);
                            
#endif

                            qrEdge = CombineQr(qrEdge, qr, e, j);

                            Debug.Assert(qrEdge != null);
                            ret[singleEdgeChunk] = qrEdge;
                        } else {
                            // nop: the edge is not in the 'edgMask'!
                            continue;
                        }
                    }

                }
            }

#if DEBUG
            for (int i = EdgeTouched.Length - 1; i >= 0; i--)
                Debug.Assert(EdgeTouched[i] == false);
#endif

            return ret.Select(p => new ChunkRulePair<DoubleEdgeQuadRule>(p.Key, p.Value));
        }

        private DoubleEdgeQuadRule CombineQr(DoubleEdgeQuadRule qrEdge, CellBoundaryQuadRule givenRule, int iFace, int jCell) {
            int D = grd.SpatialDimension;
            var volSplx = m_cellBndQF.Simplex;
            int coD = grd.Grid.GridSimplex.EdgeSimplex.SpatialDimension;

            // extract edge rule
            // -----------------

            int i0 = 0, iE = 0;
            for (int i = 0; i < iFace; i++)
                i0 += givenRule.NumbersOfNodesPerEdge[i];
            iE = i0 + givenRule.NumbersOfNodesPerEdge[iFace] - 1;

            if (iE < i0) {
                // rule is empty (measure is zero).

                if (qrEdge == null) {
                    DoubleEdgeQuadRule ret = new DoubleEdgeQuadRule();
                    ret.OrderOfPrecision = int.MaxValue - 1;
                    ret.Nodes = MultidimensionalArray.Create(1, Math.Max(1,D-1));
                    ret.Weights = MultidimensionalArray.Create(1);  // this is an empty rule, since the weight is zero!
                    ret.Median = 1;
                    // (rules with zero nodes may cause problems at various places.)
                    return ret;
                } else {
                    Debug.Assert(qrEdge.Median == qrEdge.NoOfNodes);
                    return qrEdge;
                }
            }

            MultidimensionalArray NodesVol = givenRule.Nodes.ExtractSubArrayShallow(new int[] { i0, 0 }, new int[] { iE, D - 1 });
            MultidimensionalArray Weigts = givenRule.Weights.ExtractSubArrayShallow(new int[] { i0 }, new int[] { iE }).CloneAs();
            MultidimensionalArray Nodes = MultidimensionalArray.Create(iE - i0 + 1, coD);

            volSplx.VolumeToEdgeCoordinates(iFace, NodesVol, Nodes);

            //Debug.Assert((Weigts.Sum() - grd.Grid.GridSimplex.EdgeSimplex.Volume).Abs() < 1.0e-6, "i've forgotten the gramian");

            // combine 
            // -------
            if (qrEdge == null) {
                // no rule defined yet - just set the one we have got
                // ++++++++++++++++++++++++++++++++++++++++++++++++++
                qrEdge = new DoubleEdgeQuadRule();
                qrEdge.Weights = Weigts;
                qrEdge.Nodes = Nodes;
                qrEdge.Median = qrEdge.NoOfNodes;
                qrEdge.OrderOfPrecision = givenRule.OrderOfPrecision;

            } else { 
                // take the mean of already defined and new rule
                // +++++++++++++++++++++++++++++++++++++++++++++

                int L1 = qrEdge.Nodes.GetLength(0);
                int L2 = Nodes.GetLength(0);
                Debug.Assert(coD == qrEdge.Nodes.GetLength(1));

                MultidimensionalArray newNodes = MultidimensionalArray.Create(L1 + L2, coD);
                newNodes.Set(qrEdge.Nodes, new int[] { 0, 0 }, new int[] { L1 - 1, coD - 1 });
                newNodes.Set(Nodes, new int[] { L1, 0 }, new int[] { L1 + L2 - 1, coD - 1 });
                
                MultidimensionalArray newWeights = MultidimensionalArray.Create(L1 + L2);
                newWeights.Acc(1.0, qrEdge.Weights, new int[] { 0 }, new int[] { L1 - 1 });
                newWeights.Acc(1.0, Weigts, new int[] { L1 }, new int[] { L1 + L2 - 1 });

                double oldSum = qrEdge.Weights.Sum();
                double newSum = Weigts.Sum();
                WeightInbalance += Math.Abs(oldSum - newSum);


                qrEdge.Nodes = newNodes;
                qrEdge.Weights = newWeights;
                qrEdge.Median = L1;
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
        public Simplex Simplex {
            get {
                return grd.Grid.GridSimplex.EdgeSimplex;
            }
        }
    }
}
