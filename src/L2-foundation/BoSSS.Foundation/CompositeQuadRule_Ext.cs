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
using System.Xml;
using System.Security.Claims;
using BoSSS.Foundation.Caching;
using System.Collections;

namespace BoSSS.Foundation {

    /// <summary>
    /// extension methods
    /// </summary>
    public static class ICompositeQuadRule_Ext {

        /// <summary>
        /// Creates a <see cref="CellMask"/> containing all cells covered by this composite quadrature rule.
        /// </summary>
        /// <param name="compositeRule">The composite quadrature rule to convert</param>
        /// <returns>A cell mask matching the quadrature rule's coverage</returns>
        public static CellMask GetCellMask(this ICompositeQuadRule<QuadRule> compositeRule) {
            // Get upper bound for logical cells
            var gridData = compositeRule.GridData;
            int J = gridData.iGeomCells.NoOfLocalUpdatedCells;
            var RefElements = gridData.iGeomCells.RefElements;

            // Create bitmask
            BitArray cellMask = new BitArray(J, false);
            
            // Mark cells covered by quadrature rule
            foreach (var chunkPair in compositeRule) {
                Chunk chunk = chunkPair.Chunk;
                if(Array.IndexOf(RefElements, chunkPair.Rule.RefElement) < 0)
                    throw new ArgumentException($"unknown reference element for cell: {chunkPair.Rule.RefElement}");

                for (int j = chunk.i0; j < chunk.JE; j++) {
                    if (j >= J)
                        throw new ArgumentException($"Quadrature rule contains invalid cell index {j} (max allowed: {J-1})");
                    cellMask[j] = true;
                }
            }
            
            return new CellMask(gridData, cellMask, MaskType.Geometrical);
        }

        /// <summary>
        /// Restricts a composite quadrature rule to an execution mask.
        /// </summary>
        /// <param name="compositeRule">Original quadrature rule.</param>
        /// <param name="executionMask">Mask to restrict to.</param>
        /// <returns>New composite rule containing only parts of the original rule within the mask.</returns>
        /// <exception cref="ArgumentException">Thrown if the mask contains items not in the original rule.</exception>
        public static ICompositeQuadRule<QuadRule> RestrictToMask(
            this ICompositeQuadRule<QuadRule> compositeRule, 
            ExecutionMask executionMask) {
            
            // Get mask as BitArray
            BitArray maskBits = executionMask.GetBitMask();
            int upperBound = maskBits.Length;

            // Check coverage of original rule
            BitArray originalCoverage = new BitArray(upperBound, false);
            foreach (var chunkPair in compositeRule) {
                for (int i = chunkPair.Chunk.i0; i < chunkPair.Chunk.JE; i++) {
                    if (i >= upperBound) 
                        throw new ArgumentException("Composite rule contains indices beyond the mask's range.");
                    originalCoverage[i] = true;
                }
            }

            // Validate mask is fully covered
            for (int i = 0; i < maskBits.Length; i++) {
                if (maskBits[i] && !originalCoverage[i])
                    throw new ArgumentException($"Mask contains item {i} not present in the composite rule.");
            }

            // Filter and split chunks
            List<IChunkRulePair<QuadRule>> filteredChunks = new List<IChunkRulePair<QuadRule>>();
            foreach (var chunkPair in compositeRule) {
                Chunk chunk = chunkPair.Chunk;
                QuadRule rule = chunkPair.Rule;

                List<Chunk> validSubChunks = new List<Chunk>();
                int currentStart = -1;
                
                for (int i = chunk.i0; i < chunk.JE; i++) {
                    if (maskBits[i]) {
                        currentStart = (currentStart == -1) ? i : currentStart;
                    } else if (currentStart != -1) {
                        validSubChunks.Add(new Chunk { i0 = currentStart, Len = i - currentStart });
                        currentStart = -1;
                    }
                }
                if (currentStart != -1) 
                    validSubChunks.Add(new Chunk { i0 = currentStart, Len = chunk.JE - currentStart });

                // Add filtered chunks with original rule
                foreach (var subChunk in validSubChunks) {
                    filteredChunks.Add(new ChunkRulePair<QuadRule>(subChunk, rule));
                }
            }

            return new CompositeQuadRule<QuadRule>(compositeRule.IntegrationMetric, compositeRule.GridData) { 
                chunkRulePairs = filteredChunks
            };
        }

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
        public static void SaveToTextFileCell(this ICompositeQuadRule<QuadRule> compositeRule, IGridData gridData, string filename, bool writeHeader = true) {
            _SaveToTextFileCell<QuadRule>(compositeRule, gridData, filename, writeHeader);
        }

        /// <summary>
        /// Saves the location and weight associated with each node in
        /// <paramref name="compositeRule"/> into a text file
        /// </summary>
        public static void SaveToTextFileCellBoundary(this ICompositeQuadRule<CellBoundaryQuadRule> compositeRule, IGridData gridData, string filename, bool writeHeader = true) {
            _SaveToTextFileCell<CellBoundaryQuadRule>(compositeRule, gridData, filename, writeHeader);
        }

        static void _SaveToTextFileCell<Q>(this ICompositeQuadRule<Q> compositeRule, IGridData gridData, string filename, bool writeHeader = true) where Q : QuadRule {
            int D = gridData.SpatialDimension;

            using(var file = new StreamWriter(filename)) {
                if(writeHeader) {
                    string[] dimensions = new string[] { "x", "y", "z" };
                    file.WriteLine(String.Format(
                        "Cell\t{0}\tWeight",
                        dimensions.Take(D).Aggregate((s, t) => s + "\t" + t)));
                }

                foreach(IChunkRulePair<QuadRule> pair in compositeRule) {
                    MultidimensionalArray globalNodes = gridData.GlobalNodes.GetValue_Cell(pair.Rule.Nodes, pair.Chunk.i0, pair.Chunk.Len);
                    foreach(var cell in pair.Chunk.Elements.AsSmartEnumerable()) {

                        var cellNodes = globalNodes.ExtractSubArrayShallow(cell.Index, -1, -1);
                        /*
                        Vector center = gridData.iGeomCells.GetCenter(cell.Value);

                        int[] node2Face;
                        if(pair.Rule is CellBoundaryQuadRule cbr) {
                            node2Face = new int[cbr.NoOfNodes];
                            int n = 0;
                            int iFace = 0;
                            foreach(int NNF in cbr.NumbersOfNodesPerFace) {
                                for(int nf = 0; nf < NNF; nf++) {
                                    node2Face[n] = iFace;
                                    n++;
                                }
                                iFace++;
                            }

                        } else {
                            node2Face = null;
                        }

                        var Kref = gridData.iGeomCells.GetRefElement(cell.Value);
                        int NoOfFaces = Kref.NoOfFaces;
                        Vector[] offCenter = new Vector[NoOfFaces];
                        for(int iFce = 0; iFce < NoOfFaces; iFce++) {
                            offCenter[iFce] = gridData.TransformLocal2Global(Kref.GetFaceNormal(iFce), cell.Value);
                        }
                        */

                        for(int n = 0; n < pair.Rule.NoOfNodes; n++) {
                            file.Write(cell.Value);

                            var node = cellNodes.GetRowPt(n);

                            //var ToCen = node2Face != null ? offCenter[node2Face[n]] - node : center;
                            //var nodeMod = node + ToCen * 0.1;
                            var nodeMod = node;

                            for(int d = 0; d < D; d++) {
                                file.Write("\t{0}", nodeMod[d].ToString("E", NumberFormatInfo.InvariantInfo));
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
        public static void SaveToTextFileEdge(this ICompositeQuadRule<QuadRule> compositeRule, IGridData gridData, string filename, bool writeHeader = true) {
            int D = gridData.SpatialDimension;

            using(var file = new StreamWriter(filename)) {
                if(writeHeader) {
                    string[] dimensions = new string[] { "x", "y", "z" };
                    file.WriteLine(String.Format(
                        "edge\t{0}\tWeight",
                        dimensions.Take(D).Aggregate((s, t) => s + "\t" + t)));
                }

                foreach(IChunkRulePair<QuadRule> pair in compositeRule) {
                    
                    MultidimensionalArray globalNodes = gridData.GlobalNodes.GetValue_EdgeSV(pair.Rule.Nodes, pair.Chunk.i0, pair.Chunk.Len);
                    foreach(var edge in pair.Chunk.Elements.AsSmartEnumerable()) {

                        //Vector center = gridData.iGeomEdges.GetCenter(edge.Value);
                        var edgeNodes = globalNodes.ExtractSubArrayShallow(edge.Index, -1, -1);

                        for(int n = 0; n < pair.Rule.NoOfNodes; n++) {
                            file.Write(edge.Value);

                            Vector node = edgeNodes.GetRowPt(n);
                            //var ToCen = center - node;
                            //var nodeMod = node + ToCen * 0.1;
                            var nodeMod = node;

                            for(int d = 0; d < D; d++) {
                                file.Write("\t{0}", (nodeMod[d]).ToString("E", NumberFormatInfo.InvariantInfo));
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
                    //int iTrafo = gd.iGeomEdges.Edge2CellTrafoIndex[edge, 0];
                    //int localEdge = gd.iGeomEdges.FaceIndices[edge, 0];
                    int jCell = gd.iGeomEdges.CellIndices[edge, 0];

                    
                    // Transform cell-based local coordinates to global coordinates
                    var globNodes = loopEdgeRule.Nodes.TransformLocal2Global(gd, edge);
                    
                    globNodes.OutputQuadratureRuleAsVtpXML(loopEdgeRule.Weights, filename + "For" + "j" + jCell + "e" + edge + ".vtp");
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
                    var globNodes = loopCellRule.Nodes.TransformLocal2Global(gd, jCell);
                    OutputQuadratureRuleAsVtpXML(globNodes, loopCellRule.Weights, filename + "For" + "j" + jCell + ".vtp");
                }
            }
        }

        /// <summary>
		/// Writes a xml file for visualization (Use .vtp extension for Paraview)
		/// </summary>
		static public void OutputQuadratureRuleAsVtpXML(this QuadRule qr, string filePath) {
            OutputQuadratureRuleAsVtpXML(qr.Nodes, qr.Weights, filePath);
        }

        /// <summary>
        /// Writes a xml file for visualization (Use .vtp extension for Paraview)
        /// </summary>
        static public void OutputQuadratureRuleAsVtpXML(this MultidimensionalArray nodes, MultidimensionalArray Weights, string filePath) {
            if(nodes.Dimension != 2)
                throw new ArgumentException("nodes is expected to be a 2D array (node index, coordinate index)");
            int SptialDim = nodes.GetLength(1);
            if(SptialDim != 2 && SptialDim != 3) {
                Console.Error.WriteLine("XML output is supported only for 2D and 3D schemes.");
            }
            if(Weights.Dimension != 1)
                throw new ArgumentException("Weights is expected to be a 1D array");
            if(Weights.GetLength(0) != nodes.GetLength(0))
                throw new ArgumentException("mismatch between number of Weights and number of nodes");


            try {
                using(XmlWriter writer = XmlWriter.Create(filePath, new XmlWriterSettings { Indent = true })) {
                    writer.WriteStartDocument();
                    writer.WriteStartElement("VTKFile");
                    writer.WriteAttributeString("type", "PolyData");
                    writer.WriteAttributeString("version", "0.1");
                    writer.WriteAttributeString("byte_order", "LittleEndian");

                    writer.WriteStartElement("PolyData");
                    writer.WriteStartElement("Piece");
                    writer.WriteAttributeString("NumberOfPoints", Weights.Length.ToString());
                    writer.WriteAttributeString("NumberOfVerts", Weights.Length.ToString());
                    writer.WriteAttributeString("NumberOfLines", "0");
                    writer.WriteAttributeString("NumberOfStrips", "0");
                    writer.WriteAttributeString("NumberOfPolys", "0");

                    // Points
                    writer.WriteStartElement("Points");
                    writer.WriteStartElement("DataArray");
                    writer.WriteAttributeString("type", "Float32");
                    writer.WriteAttributeString("Name", "Points");
                    writer.WriteAttributeString("NumberOfComponents", "3");
                    writer.WriteAttributeString("format", "ascii");

                    for(int i = 0; i < Weights.Length; i++) {
                        writer.WriteString($"{nodes[i, 0]} {nodes[i, 1]} {(SptialDim == 3 ? nodes[i, 2] : 0.0)}\n");
                    }

                    writer.WriteEndElement(); // DataArray
                    writer.WriteEndElement(); // Points

                    // Verts
                    writer.WriteStartElement("Verts");
                    writer.WriteStartElement("DataArray");
                    writer.WriteAttributeString("type", "Int32");
                    writer.WriteAttributeString("Name", "connectivity");
                    writer.WriteAttributeString("format", "ascii");

                    for(int i = 0; i < Weights.Length; i++) {
                        writer.WriteString($"{i}\n");
                    }

                    writer.WriteEndElement(); // DataArray

                    writer.WriteStartElement("DataArray");
                    writer.WriteAttributeString("type", "Int32");
                    writer.WriteAttributeString("Name", "offsets");
                    writer.WriteAttributeString("format", "ascii");

                    for(int i = 1; i <= Weights.Length; i++) {
                        writer.WriteString($"{i}\n");
                    }

                    writer.WriteEndElement(); // DataArray
                    writer.WriteEndElement(); // Verts

                    // PointData
                    writer.WriteStartElement("PointData");
                    writer.WriteAttributeString("Scalars", "w");

                    writer.WriteStartElement("DataArray");
                    writer.WriteAttributeString("type", "Float32");
                    writer.WriteAttributeString("Name", "w");
                    writer.WriteAttributeString("NumberOfComponents", "1");
                    writer.WriteAttributeString("format", "ascii");

                    for(int i = 0; i < Weights.Length; i++) {
                        writer.WriteString($"{Weights[i]}\n");
                    }

                    writer.WriteEndElement(); // DataArray
                    writer.WriteEndElement(); // PointData

                    writer.WriteEndElement(); // Piece
                    writer.WriteEndElement(); // PolyData
                    writer.WriteEndElement(); // VTKFile

                    writer.WriteEndDocument();
                }
            } catch(Exception ex) {
                Console.WriteLine("Error opening file: " + ex.Message);
            }
        }
    }
}
