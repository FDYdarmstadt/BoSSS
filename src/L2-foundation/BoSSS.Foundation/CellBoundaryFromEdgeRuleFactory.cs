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
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// Wrapper that creates quadrature rules in the volume coordinate
    /// system from rules given by a <see cref="IQuadRuleFactory{T}"/>s that
    /// create rules on the edges of a <see cref="RefElement"/>.
    /// </summary>
    public class CellBoundaryFromEdgeRuleFactory<T> : IQuadRuleFactory<T> where T : CellBoundaryQuadRule, new() {

        /// <summary>
        /// Omnipresent context
        /// </summary>
        private IGridData context;

        /// <summary>
        /// The factory to be wrapped
        /// </summary>
        private IQuadRuleFactory<QuadRule> edgeRuleFactory;



        /// <summary>
        /// Constructs a wrapper for <paramref name="edgeRuleFactory"/>.
        /// </summary>
        /// <param name="context">Omnipresent context</param>
        /// <param name="edgeRuleFactory">
        /// The factory to be wrapped
        /// </param>
        /// <param name="s">
        /// some reference element of the grid, see <see cref="GridCommons.RefElements"/>
        /// </param>
        public CellBoundaryFromEdgeRuleFactory(IGridData context, RefElement s, IQuadRuleFactory<QuadRule> edgeRuleFactory) {
            this.context = context;
            this.edgeRuleFactory = edgeRuleFactory;

            if (!context.iGeomEdges.EdgeRefElements.Contains(edgeRuleFactory.RefElement))
                throw new ArgumentException("Edge rule factory required");
            if (!context.iGeomCells.RefElements.Contains(s))
                throw new ArgumentException("unknown simplex.");

            this.RefElement = s;
        }

        #region IQuadRuleFactory<T> Members

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return this.edgeRuleFactory.GetCachedRuleOrders();
        }

        /// <summary>
        /// Uses the given edge rule factory to create quadrature rules for
        /// the boundaries of all cells in <paramref name="mask"/>.
        /// </summary>
        /// <param name="mask">
        /// A <b><see cref="CellMask"/></b> containing all cells to be
        /// integrated over.
        /// </param>
        /// <param name="order">
        /// <see cref="IQuadRuleFactory{T}.GetQuadRuleSet"/>
        /// </param>
        /// <returns>
        /// Quadrature rules for the boundary of all cells in
        /// <paramref name="mask"/> in the volume coordinate system.
        /// </returns>
        public IEnumerable<IChunkRulePair<T>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (!(mask is CellMask))
                throw new ArgumentException("Expecting a cell/volume mask.");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");
            if (!object.ReferenceEquals(this.context, mask.GridData))
                throw new ArgumentException();

            // helper vars
            int D = this.context.SpatialDimension;
            var CellToEdge = this.context.iLogicalCells.Cells2Edges;
            var EdgeToCell = this.context.iLogicalEdges.CellIndices;
            var FaceIndices = this.context.iGeomEdges.FaceIndices;
            var TrafoIdx = this.context.iGeomEdges.Edge2CellTrafoIndex;
            var EdgeToCellTrafos = this.context.iGeomEdges.Edge2CellTrafos;
            int NoOfFaces = RefElement.NoOfFaces;
            var Scalings = context.iGeomEdges.Edge2CellTrafos_SqrtGramian;
            int J = this.context.iLogicalCells.NoOfLocalUpdatedCells;

            // output data structure, temporary
            Dictionary<int, T> cellBndRuleMap = new Dictionary<int, T>();

            // Get edge quad rules
            BitArray edgeBitMask = new BitArray(context.iGeomEdges.Count);
            foreach (Chunk chunk in mask) {
                foreach (int cell in chunk.Elements) {
                    var LocalCellIndexToEdges_cell = CellToEdge[cell];
                    for (int e = LocalCellIndexToEdges_cell.Length - 1; e >= 0; e--) {
                        int edge = Math.Abs(LocalCellIndexToEdges_cell[e]) - 1;
                        edgeBitMask[edge] = true;
                    }
                    cellBndRuleMap.Add(cell, new T() {
                        NumbersOfNodesPerFace = new int[NoOfFaces]
                    });
                }
            }
            IEnumerable<IChunkRulePair<QuadRule>> edgeRuleMap = edgeRuleFactory.GetQuadRuleSet(
                    new EdgeMask(context, edgeBitMask, MaskType.Geometrical), order);

            
            // build cell boundary rule
            var CellBitMask = mask.GetBitMask();
            foreach (var edgeRule in edgeRuleMap) {
                QuadRule rule = edgeRule.Rule;
                Debug.Assert(rule.Nodes.GetLength(1) == Math.Max(D - 1, 1));
                int iEdge0 = edgeRule.Chunk.i0;
                int iEdgeE = edgeRule.Chunk.JE;

                for (int iEdge = iEdge0; iEdge < iEdgeE; iEdge++) {
                    for (int kk = 0; kk < 2; kk++) { // loop over in and out
                        int jCell = EdgeToCell[iEdge, kk];
                        if (jCell < 0 || jCell >= J)
                            continue;
                        if (!CellBitMask[jCell])
                            continue;
                        
                        int iTrafo = TrafoIdx[iEdge, kk];
                        var trafo = EdgeToCellTrafos[iTrafo];
                        int iFace = FaceIndices[iEdge, kk];
                        double scl = Scalings[iTrafo];

                        NodeSet VolumeNodes = new NodeSet(this.RefElement, rule.NoOfNodes, D);
                        trafo.Transform(rule.Nodes, VolumeNodes);
                        VolumeNodes.LockForever();

                        var cellBndRule = cellBndRuleMap[jCell];
                        AddToCellBoundaryRule(cellBndRule, VolumeNodes, rule.Weights, iFace, scl);
                    }
                }
            }


            // return 
            // ======
            {
                int[] Cells = cellBndRuleMap.Keys.ToArray();
                Array.Sort(Cells);

                ChunkRulePair<T>[] R = new ChunkRulePair<T>[Cells.Length];
                for (int i = 0; i < R.Length; i++) {
                    int jCell = Cells[i];
                    cellBndRuleMap[jCell].Nodes.LockForever();
                    R[i] = new ChunkRulePair<T>(Chunk.GetSingleElementChunk(jCell), cellBndRuleMap[jCell]);
                }


                return R;
            }
        }

        /// <summary>
        /// <see cref="IQuadRuleFactory{T}.RefElement"/>
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        #endregion

        static void AddToCellBoundaryRule(T cbqr, NodeSet nodes, MultidimensionalArray weights, int iFace, double scl) {

            Debug.Assert(nodes.GetLength(0) == weights.GetLength(0));
            if (nodes.GetLength(0) <= 0)
                return;

            if (cbqr.Nodes == null) {
                cbqr.Nodes = nodes;
                Debug.Assert(cbqr.Weights == null);
                cbqr.Weights = weights.CloneAs();
                cbqr.Weights.Scale(scl);
                cbqr.NumbersOfNodesPerFace[iFace] = nodes.GetLength(0);
            } else {
                int D = nodes.GetLength(1);
                Debug.Assert(D == cbqr.Nodes.GetLength(1));

                int L1 = cbqr.Nodes.GetLength(0);
                int L2 = nodes.GetLength(0);

                NodeSet newNodes = new NodeSet(cbqr.Nodes.RefElement, L1 + L2, D);
                MultidimensionalArray newWeights = MultidimensionalArray.Create(L1 + L2);

                int K = 0; // total number of nodes on all faces;
                for (int k = 0; k < iFace; k++) {
                    K += cbqr.NumbersOfNodesPerFace[k];
                }

                if (K > 0) {
                    newNodes.ExtractSubArrayShallow(
                        new int[] { 0, 0 },
                        new int[] { K - 1, D - 1 }).Set(
                            cbqr.Nodes.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, D - 1 }));
                    newWeights.ExtractSubArrayShallow(
                        new int[] { 0 },
                        new int[] { K - 1 }).Set(
                        cbqr.Weights.ExtractSubArrayShallow(new int[] { 0 }, new int[] { K - 1 }));
                }

                newNodes.ExtractSubArrayShallow(new int[] { K, 0 }, new int[] { K + L2 - 1, D - 1 }).Set(nodes);
                newWeights.ExtractSubArrayShallow(new int[] { K }, new int[] { K + L2 - 1 }).Acc(scl, weights);

                if (K < L1) {
                    newNodes.ExtractSubArrayShallow(
                        new int[] { K + L2, 0 },
                        new int[] { L1 + L2 - 1, D - 1 }).Set(
                            cbqr.Nodes.ExtractSubArrayShallow(new int[] { K, 0 }, new int[] { L1 - 1, D - 1 }));
                    newWeights.ExtractSubArrayShallow(
                        new int[] { K + L2 },
                        new int[] { L1 + L2 - 1 }).Set(
                        cbqr.Weights.ExtractSubArrayShallow(new int[] { K }, new int[] { L1 - 1 }));
                }

                cbqr.Nodes = newNodes;
                cbqr.Nodes.LockForever();
                cbqr.Weights = newWeights;
                cbqr.NumbersOfNodesPerFace[iFace] += L2;
            }
        }


    }
}
