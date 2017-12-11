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
using System.Collections.Generic;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Tracing;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// A subdivision strategy that adaptively subdivides a node into sub-nodes
    /// if it is cut the interface.
    /// </summary>
    public class AdaptiveSubdivisionStrategy : ISubdivisionStrategy {
        

        /// <summary>
        /// The maximum number of subdivisions of the simplex
        /// </summary>
        private readonly int maxDivisions;

        /// <summary>
        /// Vertex set containing the vertices of the root simplex.
        /// </summary>
        private readonly NestedVertexSet baseVertexSet;

        /// <summary>
        /// The subdivision tree that handles the actual subdivisions of the
        /// simplex.
        /// </summary>
        private readonly SimplexSubdivisionTree subdivisionTree;


        /// <summary>
        /// level-set evaluation
        /// </summary>
        LevelSetTracker.LevelSetData LevelSetData;

        
        /// <summary>
        /// The omnipresent context.
        /// </summary>
        private GridData gridData {
            get {
                return LevelSetData.GridDat;
            }
        }


        /// <summary>
        /// Initializes the subdivision of the given simplex. Initially, the
        /// simplex remains undivided.
        /// </summary>
        /// <param name="refElement">
        /// The simplex to subdivide.
        /// </param>
        /// <param name="tracker">
        /// The level set tracker. Allows for the check if a simplex ist cut.
        /// </param>
        /// <param name="levSetIndex">Index of the level set</param>
        /// <param name="maxDivisions">
        /// The maximum number of subdivisions of the simplex
        /// </param>
        public AdaptiveSubdivisionStrategy(RefElement refElement, LevelSetTracker.LevelSetData levelSetData, int maxDivisions) {
            this.RefElement = refElement;
            this.maxDivisions = maxDivisions;
            this.baseVertexSet = new NestedVertexSet(refElement.SpatialDimension);
            this.LevelSetData = LevelSetData;

            int verticesPerCell = refElement.Vertices.GetLength(0);
            int[] simplexVertices = new int[verticesPerCell];
            for (int i = 0; i < verticesPerCell; i++) {
                double[] vertex = refElement.Vertices.GetRow(i);
                simplexVertices[i] = baseVertexSet.RegisterVertex(vertex);
            }
            this.subdivisionTree = new SimplexSubdivisionTree(
                refElement, baseVertexSet, simplexVertices);

            this.subdivisionTree.SetSavePoint();
        }

        /// <summary>
        /// Estimates the distance below which we consider an element (i.e.,
        /// cell or edge) cut. This is important in the case of curved
        /// interfaces.
        /// </summary>
        /// <param name="element">
        /// The element (edge or cell) in question.
        /// </param>
        /// <param name="mask">
        /// The mask containing the given element <paramref name="element"/>.
        /// </param>
        /// <returns>
        /// The minimal distance below which consider element
        /// <paramref name="element"/> cut.
        /// </returns>
        /// <remarks>
        /// Uses different estimates depending on the element type. The case of
        /// edges is straightforward (since we don't have enough information to
        /// do something sophisticated). For cells, we use a modified version of
        /// a formular given by MinGibou2007:
        /// \f$ 
        /// \max |\Phi(v)| \leq 0.5 * \mathrm{Lip}(\Phi) h_\mathrm{max}.
        /// \f$ 
        /// In particular, we assume that the Lipschitz constant is close to 1
        /// which is true if the level set is signed distance. A more
        /// sophisticated method would evaluate the gradients at the vertices
        /// in order to estimate the Lipschitz constant, but this would be less
        /// efficient.
        /// </remarks>
        private double EstimateMinDistance(int element, ExecutionMask mask) {
            double minDistance = double.NaN;

            if (mask is CellMask) {
                // Assume norm of gradient \approx 1
                // -> Lipschitz constant is in the order of 1
                double hMaxEdge = -1.0;
                int[] cells2Edges = gridData.Cells.Cells2Edges[element];
                for (int e = 0; e < cells2Edges.Length; e++) {
                    int edgeIndex = Math.Abs(cells2Edges[e]) - 1;
                    hMaxEdge = Math.Max(gridData.Edges.h_max_Edge[edgeIndex], hMaxEdge);
                }
                minDistance = 0.5 * hMaxEdge;
            } else if (mask is EdgeMask) {
                //minDistance = 0.5 * gridData.Edges.h_max_Edge[element];
                minDistance = 1.0 * gridData.Edges.h_max_Edge[element];
            } else {
                throw new NotImplementedException("Unknown mask type");
            }

            return minDistance;
        }

        /// <summary>
        /// Determines the vertices associated with the given element (either
        /// a cell or an edge). In case of a cell, they're directly stored in
        /// the given set <paramref name="set"/>. In case of an edge, the
        /// stored vertices have to be transformed from the edge coordinate
        /// system to the volume coordinate system.
        /// </summary>
        /// <param name="element">
        /// The index of either a cell or an edge
        /// </param>
        /// <param name="set">
        /// The set containing the vertices
        /// </param>
        /// <param name="mask">
        /// The mask containing <paramref name="element"/>.
        /// </param>
        /// <param name="cell">
        /// On exit: The cell the given element <paramref name="element"/>
        /// belongs to.
        /// </param>
        /// <returns>
        /// The vertices associated with <paramref name="element"/>.
        /// </returns>
        private NodeSet GetVertices(int element, NestedVertexSet set, ExecutionMask mask, out int cell) {
            //NodeSet vertices = new NodeSet(this.RefElement, set.LocalVertices);
            //cell = element;

            //if (vertices == null) {
            //    cell = -1;
            //    return null;
            //}

            //if (mask is EdgeMask) {
            //    // This might be dangerous if level set is discontinuous
            //    cell = gridData.Edges.CellIndices[element, 0];

            //    NodeSet volumeVertices = new NodeSet(this.RefElement, vertices.GetLength(0), gridData.SpatialDimension);

            //    int localEdge;
            //    if (gridData.Edges.CellIndices[element, 0] == cell) {
            //        localEdge = gridData.Edges.FaceIndices[element, 0];
            //    } else {
            //        localEdge = gridData.Edges.FaceIndices[element, 1];
            //    }

            //    gridData.Grid.RefElements[0].TransformFaceCoordinates(localEdge, vertices, volumeVertices);
            //    volumeVertices.LockForever();

            //    vertices = volumeVertices;
            //}

            //return vertices;

            MultidimensionalArray vertices = set.LocalVertices;
            cell = element;

            if(vertices == null) {
                cell = -1;
                return null;
            }
            
            if(mask is EdgeMask) {
                // This might be dangerous if level set is discontinuous
                cell = gridData.Edges.CellIndices[element, 0];

                MultidimensionalArray volumeVertices = MultidimensionalArray.Create(
                    vertices.GetLength(0), gridData.SpatialDimension);

                int localEdge;
                if(gridData.Edges.CellIndices[element, 0] == cell) {
                    localEdge = gridData.Edges.FaceIndices[element, 0];
                } else {
                    localEdge = gridData.Edges.FaceIndices[element, 1];
                }

                gridData.Grid.RefElements[0].TransformFaceCoordinates(localEdge, vertices, volumeVertices);
                vertices = volumeVertices;
            }

            return new NodeSet(this.gridData.Cells.GetRefElement(cell), vertices);
        }

        #region ISubdivisionStrategy Members

        /// <summary>
        /// <see cref="ISubdivisionStrategy.RefElement"/>.
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        /// <summary>
        /// Constructs the adapted subdivisions for all cells in
        /// <paramref name="mask"/>. The size of each chunk will be exactly one
        /// since, in general, we cannot expect subsequent cells to have the
        /// same subdivisions.
        /// </summary>
        /// <param name="mask">
        /// <see cref="ISubdivisionStrategy.GetSubdivisionNodes"/>
        /// </param>
        /// <returns>
        /// The subdivisions for all elements in <paramref name="mask"/>,
        /// </returns>
        public IEnumerable<KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>> GetSubdivisionNodes(ExecutionMask mask) {
            using(new FuncTrace()) {
                foreach(int iElement in mask.ItemEnum) { // loop over all cells/edges in mask

                    this.subdivisionTree.ResetToSavePoint();

                    double minDistance = EstimateMinDistance(iElement, mask);

                    NestedVertexSet currentSet = baseVertexSet;
                    for(int i = 0; i < maxDivisions; i++) {
                        if(currentSet.NumberOfLocalVertices == 0) {
                            // No new vertices were added during last subdivision
                            break;
                        }

                        int cell;
                        NodeSet vertices = GetVertices(iElement, currentSet, mask, out cell);

                        MultidimensionalArray levelSetValues = LevelSetData.GetLevSetValues(vertices, cell, 1);

                        subdivisionTree.ReadDistances(levelSetValues.ExtractSubArrayShallow(0, -1));

                        currentSet = new NestedVertexSet(currentSet);
                        subdivisionTree.Subdivide(
                            currentSet,
                            true,
                            minDistance / Math.Pow(2.0, i));
                    }

                    // Finally, read level set values in leaves (needed by IsCut)
                    {
                        int cell;
                        NodeSet vertices = GetVertices(iElement, currentSet, mask, out cell);
                        if(vertices != null) {
                            MultidimensionalArray levelSetValues = LevelSetData.GetLevSetValues(vertices, cell, 1);

                            subdivisionTree.ReadDistances(levelSetValues.ExtractSubArrayShallow(0, -1));
                        }
                    }

                    yield return new KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>(
                        new Chunk() {
                            i0 = iElement,
                            Len = 1
                        },
                        subdivisionTree.Leaves);
                }
            }

        }

        #endregion
    }
}
