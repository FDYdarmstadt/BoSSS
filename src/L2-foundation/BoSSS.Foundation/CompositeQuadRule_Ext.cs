using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation {

    /// <summary>
    /// extension methods
    /// </summary>
    public static class ICompositeQuadRule_Ext {

        /// <summary>
        /// Saves the sum of weights of each edge rule in the given
        /// <paramref name="compositeRule"/> together with the coordinates of
        /// the corresponding edge center into a text file.
        /// </summary>
        public static void SumOfWeightsToTextFileEdge(this ICompositeQuadRule<QuadRule> compositeRule, IGridData g, string filename) {
            int E = g.iLogicalEdges.Count;
            var bMask = new System.Collections.BitArray(E);
            double[] wSum = new double[E];
            foreach(IChunkRulePair<QuadRule> crp in compositeRule) {
                for(int iEdge = crp.Chunk.i0; iEdge < crp.Chunk.JE; iEdge++) {
                    bMask[iEdge] = true;
                    wSum[iEdge] = crp.Rule.Weights.Sum();
                }
            }

            var mask = new EdgeMask(g, bMask);
            mask.SaveToTextFile(filename, false, (X, i, ii) => wSum[i]);
        }

        /// <summary>
        /// Saves the sum of weights of each volume rule in the given
        /// <paramref name="compositeRule"/> together with the coordinates of
        /// the corresponding cell center into a text file.
        /// </summary>
        public static void SumOfWeightsToTextFileVolume(this ICompositeQuadRule<QuadRule> compositeRule, IGridData g, string filename) {
            int J = g.iLogicalCells.NoOfLocalUpdatedCells;
            var bMask = new System.Collections.BitArray(J);
            double[] wSum = new double[J];
            foreach(IChunkRulePair<QuadRule> crp in compositeRule) {
                for(int iCell = crp.Chunk.i0; iCell < crp.Chunk.JE; iCell++) {
                    bMask[iCell] = true;
                    wSum[iCell] = crp.Rule.Weights.Sum();
                }
            }

            var mask = new CellMask(g, bMask);
            mask.SaveToTextFile(filename, false, (X, i, ii) => wSum[i]);
        }

        /// <summary>
        /// Saves the location and weight associated with each node in
        /// <paramref name="compositeRule"/> into a text file
        /// </summary>
        public static void ToTextFileCell(this ICompositeQuadRule<QuadRule> compositeRule, IGridData gridData, string filename) {
            int D = gridData.SpatialDimension;
            string[] dimensions = new string[] { "x", "y", "z" };

            using(var file = new StreamWriter(filename)) {
                file.WriteLine(String.Format(
                    "Cell\t{0}\tWeight",
                    dimensions.Take(D).Aggregate((s, t) => s + "\t" + t)));

                foreach(IChunkRulePair<QuadRule> pair in compositeRule) {
                    MultidimensionalArray globalNodes = gridData.GlobalNodes.GetValue_Cell(pair.Rule.Nodes, pair.Chunk.i0, pair.Chunk.Len);
                    foreach(var cell in pair.Chunk.Elements.AsSmartEnumerable()) {
                        for(int n = 0; n < pair.Rule.NoOfNodes; n++) {
                            file.Write(cell.Value);

                            for(int d = 0; d < D; d++) {
                                file.Write("\t{0}", globalNodes[cell.Index, n, d].ToString("E", NumberFormatInfo.InvariantInfo));
                            }

                            file.WriteLine("\t{0}", pair.Rule.Weights[n].ToString("E", NumberFormatInfo.InvariantInfo));
                        }
                    }

                    file.Flush();
                }
            }
        }

        /// <summary>
        /// Saves the location and weight associated with each node in
        /// <paramref name="compositeRule"/> into a text file
        /// </summary>
        public static void ToTextFileEdge(this ICompositeQuadRule<QuadRule> compositeRule, IGridData gridData, string filename) {
            int D = gridData.SpatialDimension;
            string[] dimensions = new string[] { "x", "y", "z" };

            using(var file = new StreamWriter(filename)) {
                file.WriteLine(String.Format(
                    "Cell\t{0}\tWeight",
                    dimensions.Take(D).Aggregate((s, t) => s + "\t" + t)));

                foreach(IChunkRulePair<QuadRule> pair in compositeRule) {
                    MultidimensionalArray globalNodes = gridData.GlobalNodes.GetValue_EdgeSV(pair.Rule.Nodes, pair.Chunk.i0, pair.Chunk.Len);
                    foreach(var cell in pair.Chunk.Elements.AsSmartEnumerable()) {
                        for(int n = 0; n < pair.Rule.NoOfNodes; n++) {
                            file.Write(cell.Value);

                            for(int d = 0; d < D; d++) {
                                file.Write("\t{0}", globalNodes[cell.Index, n, d].ToString("E", NumberFormatInfo.InvariantInfo));
                            }

                            file.WriteLine("\t{0}", pair.Rule.Weights[n].ToString("E", NumberFormatInfo.InvariantInfo));
                        }
                    }

                    file.Flush();
                }
            }
        }

        /// <summary>
        /// Write the edge quadrature rules into vtp files (importable to paraview)
        /// </summary>
        /// <param name="chunRulePairList">the lsit of edge quadrature rules</param>
        /// <param name="gd">grid data</param>
        /// <param name="filename">filename header</param>
        public static void ToVtpFilesEdge(this ICompositeQuadRule<QuadRule> chunRulePairList, IGridData gd, string filename) {
            foreach(var chunkRulePair in chunRulePairList) {
                foreach(int edge in chunkRulePair.Chunk.Elements) {
                    var loopEdgeRule = chunkRulePair.Rule;
                    int iTrafo = gd.iGeomEdges.Edge2CellTrafoIndex[edge, 0];
                    int localEdge = gd.iGeomEdges.FaceIndices[edge, 0];
                    int jCell = gd.iGeomEdges.CellIndices[edge, 0];

                    // Create cell-based quadrature rule (to express them in cell coordinates)
                    RefElement KrefCell = gd.iGeomCells.GetRefElement(jCell);
                    var cellNodes = loopEdgeRule.Nodes.GetVolumeNodeSet(gd, iTrafo, false);
                    var cellRule = QuadRule.CreateEmpty(KrefCell, cellNodes.NoOfNodes, KrefCell.SpatialDimension);
                    cellRule.Nodes = cellNodes;
                    cellRule.Weights = loopEdgeRule.Weights;

                    // Transform cell-based local coordinates to global coordinates
                    cellRule.TransformLocal2Global(gd, jCell);
                    cellRule.OutputQuadratureRuleAsVtpXML(filename + "For" + "j" + jCell + "e" + edge + ".vtp");
                }
            }
        }

        /// <summary>
        /// Write the cell quadrature rules (volume or surface) into vtp files (importable to paraview)
        /// </summary>
        /// <param name="chunRulePairList">the lsit of edge quadrature rules</param>
        /// <param name="gd">grid data</param>
        /// <param name="filename">filename header</param>
        public static void ToVtpFilesCell(this ICompositeQuadRule<QuadRule> chunRulePairList, IGridData gd, string filename) {
            foreach(var chunkRulePair in chunRulePairList) {
                foreach(int jCell in chunkRulePair.Chunk.Elements) {
                    var loopCellRule = chunkRulePair.Rule;

                    // Transform cell-based local coordinates to global coordinates
                    loopCellRule.TransformLocal2Global(gd, jCell);
                    loopCellRule.OutputQuadratureRuleAsVtpXML(filename + "For" + "j" + jCell + ".vtp");
                }
            }
        }
    }
}
