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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// Represents a subdivision-tree for a given <see cref="RefElement"/>. The
    /// tree is a $n-tree where $n is the number of subdivisions for the
    /// simplex (see <see cref="RefElement.GetSubdivisionTree"/>).
    /// </summary>
    public class SimplexSubdivisionTree {

        /// <summary>
        /// The simplex to be subdivided
        /// </summary>
        private RefElement refElement;

        /// <summary>
        /// The root node of the subdivision tree
        /// </summary>
        private Node treeRoot;

        /// <summary>
        /// The current depth of the tree
        /// </summary>
        private int depth;

        /// <summary>
        /// The depth of the tree at the last savepoint (see
        /// <see cref="SetSavePoint"/>).
        /// </summary>
        private int savedDepth;

        /// <summary>
        /// Constructs a new tree with depth zero.
        /// </summary>
        /// <param name="refElement">
        /// The simplex to be subdivided
        /// </param>
        /// <param name="vertexSet">
        /// A container for the vertices of the subdivisions simplices. Must
        /// already contain the vertices of the root simplex.
        /// </param>
        /// <param name="rootVertexIndices">
        /// Indices of the vertices of the root simplex in
        /// <paramref name="vertexSet"/>.
        /// </param>
        public SimplexSubdivisionTree(RefElement refElement, NestedVertexSet vertexSet, int[] rootVertexIndices) {
            this.refElement = refElement;
            this.treeRoot = new Node(this, 0, vertexSet, rootVertexIndices, AffineTrafo.Identity(refElement.SpatialDimension));
            this.depth = 0;
        }

        /// <summary>
        /// Takes a snapshot of the current tree. If
        /// <see cref="ResetToSavePoint"/> is called, this state will be
        /// restored. This can be useful if a minimal number of subdivisions
        /// of a simplex should be present in all cells.
        /// </summary>
        public void SetSavePoint() {
            savedDepth = depth;
        }

        /// <summary>
        /// Restores the state after the last call to
        /// <see cref="SetSavePoint"/>.
        /// </summary>
        public void ResetToSavePoint() {
            depth = savedDepth;
            treeRoot.Reset(savedDepth);
        }

        /// <summary>
        /// Distributes the given distance information to the leaves of this
        /// tree (see <see cref="Leaves"/>). A "distance", in this context,
        /// denotes the signed distance to the interface to which this tree
        /// should adapt.
        /// </summary>
        /// <param name="distances">
        /// A one-dimensional array containing the distances of the nodes of
        /// the current leaves of the tree (see <see cref="Leaves"/>). The
        /// ordering of the distances corresponds to the ordering defined by
        /// the <see cref="NestedVertexSet"/> associated with the leaves.
        /// </param>
        public void ReadDistances(MultidimensionalArray distances) {
            treeRoot.ReadDistances(depth, distances);
        }

        /// <summary>
        /// Subdivides all leaves of the current tree (see
        /// <see cref="Leaves"/>) that can be considered "cut" by the
        /// interface. A node is considered cut if the distances (provided via
        /// <see cref="ReadDistances"/>) have different signs at the vertices
        /// of a leave _or_ if the distance of at least one vertex is
        /// smaller than <paramref name="minDistance"/>.
        /// </summary>
        /// <param name="refinedVertexSet">
        /// The vertex set that should hold the vertices of the newly created
        /// nodes associated with the new leaves of the tree.
        /// </param>
        /// <param name="checkIsCut">
        /// If true, leaves will only subdivided of it is considered cut. If
        /// false, all leaves will be subdivided indifferently.
        /// </param>
        /// <param name="minDistance">
        /// The minimum distance a vertx of a leave may have before the leave
        /// is considered cut.
        /// </param>
        /// <remarks>
        /// The parameter <paramref name="minDistance"/> can be used to ensure
        /// that all cut cells are captured in the case of a non-linear
        /// interface. Here, the simple check if all vertex-distances have the
        /// same sign is not sufficient since a curved interface can pass
        /// through a cell without enclosing a vertex.
        /// </remarks>
        public void Subdivide(NestedVertexSet refinedVertexSet, bool checkIsCut, double minDistance) {
            if (checkIsCut) {
                bool divided = treeRoot.SubdivideCutLeaves(depth, refinedVertexSet, minDistance);
                if (divided) {
                    depth++;
                }
            } else {
                treeRoot.SubdivideLeaves(depth, refinedVertexSet);
                depth++;
            }
        }

        /// <summary>
        /// All leaves of this tree.
        /// </summary>
        public List<SubdivisionNode> Leaves {
            get {
                List<SubdivisionNode> result = new List<SubdivisionNode>();
                CollectLeavesRecursively(treeRoot, result);
                return result;
            }
        }

        /// <summary>
        /// Recursively collects all leaves of this tree.
        /// </summary>
        /// <param name="node">
        /// The root of the current branch of the recursion.
        /// </param>
        /// <param name="list">
        /// On exit: Contains all leaves of this tree.
        /// </param>
        private void CollectLeavesRecursively(Node node, List<SubdivisionNode> list) {
            if (node.Children.Count == 0) {
                list.Add(node);
            } else {
                foreach (Node childNode in node.Children) {
                    CollectLeavesRecursively(childNode, list);
                }
            }
        }

        /// <summary>
        /// A node of the subdivisions tree.
        /// </summary>
        private class Node : SubdivisionNode {

            /// <summary>
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </summary>
            private int level;

            /// <summary>
            /// The distances of the vertices of this node to the interface.
            /// The i-th entry belongs to the vertex with global index
            /// <see cref="globalVertexIndices"/>[i] in<see cref="vertexSet"/>.
            /// </summary>
            private double[] distances;

            /// <summary>
            /// The owner of this node.
            /// </summary>
            private SimplexSubdivisionTree owner;

            /// <summary>
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </summary>
            private int[] globalVertexIndices;

            /// <summary>
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </summary>
            private NestedVertexSet vertexSet;

            /// <summary>
            /// Constructs a root node of a tree.
            /// </summary>
            /// <param name="owner">
            /// The creator of this object.
            /// </param>
            /// <param name="level">
            /// The level of this node in the subdivision tree, i.e. the number
            /// of ancestors.
            /// </param>
            /// <param name="vertexSet">
            /// The set containing all vertices of <b>all</b> nodes on this
            /// level of the subdivision tree (i.e., of this node and its
            /// siblings).
            /// </param>
            /// <param name="vertexIndices">
            /// The indices of the vertices of this node in
            /// <see cref="vertexSet"/>.
            /// </param>
            /// <param name="transformationFromRoot">
            /// The affine transformation that transforms a vertex of the root
            /// simplex to a vertex of this node.
            /// </param>
            public Node(SimplexSubdivisionTree owner, int level, NestedVertexSet vertexSet, int[] vertexIndices, AffineTrafo transformationFromRoot)
                : this(owner, null, level, vertexSet, vertexIndices, transformationFromRoot) {
            }

            /// <summary>
            /// Creates a new child of <paramref name="parrent"/>.
            /// </summary>
            /// <param name="owner">
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </param>
            /// <param name="parrent">
            /// The parrent of this node.
            /// </param>
            /// <param name="level">
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </param>
            /// <param name="vertexSet">
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </param>
            /// <param name="vertexIndices">
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </param>
            /// <param name="transformationFromRoot">
            /// <see cref="Node.Node(SimplexSubdivisionTree, Node, int, NestedVertexSet, int[], AffineTrafo)"/>
            /// </param>
            private Node(SimplexSubdivisionTree owner, Node parrent, int level, NestedVertexSet vertexSet, int[] vertexIndices, AffineTrafo transformationFromRoot)
                : base(transformationFromRoot) {
                Debug.Assert(vertexIndices.Length == owner.refElement.NoOfVertices, "Wrong number of vertices");
                Debug.Assert(vertexIndices.All((i) => i >= 0), "All vertex indices must be positive");
                Debug.Assert(vertexIndices.Distinct().Count() == vertexIndices.Length, "Vertex indices must be unique");

                this.owner = owner;
                this.level = level;
                this.vertexSet = vertexSet;
                this.globalVertexIndices = vertexIndices;

                Children = new List<Node>();
                distances = new double[globalVertexIndices.Length];
                ArrayTools.SetAll(distances, double.NaN);

                if (parrent != null) {
                    for (int i = 0; i < globalVertexIndices.Length; i++) {
                        int parrentIndex = -1;
                        for (int j = 0; j < parrent.globalVertexIndices.Length; j++) {
                            if (parrent.globalVertexIndices[j] == globalVertexIndices[i]) {
                                parrentIndex = j;
                                break;
                            }
                        }

                        if (parrentIndex >= 0) {
                            distances[i] = parrent.distances[parrentIndex];
                        }
                    }
                }
            }

            /// <summary>
            /// All children of this nodes in the tree.
            /// </summary>
            public IList<Node> Children {
                get;
                private set;
            }

            /// <summary>
            /// Tells the node to subdivide itself if it is a leave of the
            /// tree. If it is not a tree, the command is passed down to all
            /// children.
            /// </summary>
            /// <param name="leavesLevel">
            /// The current depth of tree.
            /// </param>
            /// <param name="refinedVertexSet">
            /// The vertex set that, on exit, holds the newly added vertices
            /// of the new nodes.
            /// </param>
            public void SubdivideLeaves(int leavesLevel, NestedVertexSet refinedVertexSet) {
                SubdivideLeaves(leavesLevel, refinedVertexSet, false, double.NaN);
            }

            /// <summary>
            /// Tells the node to subdivide itself if it is a leave of the
            /// tree. If it is not a tree, the command is passed down to all
            /// children. In contrast to
            /// <see cref="SubdivideLeaves(int, NestedVertexSet)"/>, nodes
            /// are only cut if they are (potentially) cut be the interface.
            /// </summary>
            /// <param name="leavesLevel">
            /// <see cref="SubdivideLeaves(int, NestedVertexSet)"/>
            /// </param>
            /// <param name="refinedVertexSet">
            /// <see cref="SubdivideLeaves(int, NestedVertexSet)"/>
            /// </param>
            /// <param name="minDistance">
            /// The minimum distance a node may have from the interface before
            /// the node is considered cut (even though the signs of the
            /// distances in all vertices are the same). In the case of curved
            /// interfaces, this is necessary to fix certain corner cases.
            /// </param>
            /// <returns>
            /// True, if at least one node has been subdivided. Otherwise,
            /// false is returned.
            /// </returns>
            public bool SubdivideCutLeaves(int leavesLevel, NestedVertexSet refinedVertexSet, double minDistance) {
                return SubdivideLeaves(leavesLevel, refinedVertexSet, true, minDistance);
            }

            /// <summary>
            /// Passes the given distances to the leaves of tree where the
            /// information is stored and used by <see cref="IsCut"/>.
            /// </summary>
            /// <param name="leavesLevel">
            /// The current depth of the tree.
            /// </param>
            /// <param name="newDistances">
            /// A one-dimensional array containing the distances of the nodes of
            /// the current leaves of the tree. The ordering of the distances
            /// corresponds to the ordering defined by the
            /// <see cref="NestedVertexSet"/> (see
            /// <see cref="Node.Node(SimplexSubdivisionTree, int, NestedVertexSet, int[], AffineTrafo)"/>)
            /// associated with the leaves.
            /// </param>
            public void ReadDistances(int leavesLevel, MultidimensionalArray newDistances) {
                if (level < leavesLevel) {
                    if (Children.Count > 0) {
                        foreach (Node child in Children) {
                            child.ReadDistances(leavesLevel, newDistances);
                        }
                    }
                } else if (level == leavesLevel) {
                    for (int i = 0; i < globalVertexIndices.Length; i++) {
                        int globalIndex = globalVertexIndices[i];
                        if (vertexSet.VertexIsContainedInLocalSet(globalIndex)) {
                            distances[i] = newDistances[vertexSet.GetLocalFromGlobalVertexIndex(globalIndex)];
                        }
                    }

                    Debug.Assert(!distances.Any((d) => double.IsNaN(d)), "Some level set values remain undefined");
                } else {
                    throw new Exception("This should not have happened.");
                }
            }

            /// <summary>
            /// Clears all information stored in this node (i.e., the distances
            /// to the interface). Additionally:
            /// <list type="bullet">
            ///     <item>
            ///     For nodes with level smaller than
            ///     <paramref name="cutOffDepth"/>: Resets all children
            ///     </item>
            ///     <item>
            ///     For nodes with level equal to <paramref name="cutOffDepth"/>:
            ///     Deletes all children
            ///     </item>
            /// </list>
            /// </summary>
            /// <param name="cutOffDepth">
            /// The minimal depth the tree must have after the reset. In other
            /// words: After the reset, will have this depth.
            /// </param>
            public void Reset(int cutOffDepth) {
                distances = new double[globalVertexIndices.Length];
                ArrayTools.SetAll(distances, double.NaN);

                if (this.level == cutOffDepth) {
                    Children = new List<Node>();
                } else if (this.level < cutOffDepth) {
                    foreach (Node child in Children) {
                        child.Reset(cutOffDepth);
                    }
                } else {
                    throw new ArgumentException("Given cut-off depth is bigger than the level of this node", "cutOffDepth");
                }
            }

            /// <summary>
            /// <see cref="IsPotentiallyCut"/>(0.0)
            /// </summary>
            /// <remarks>
            /// Assuming minDistance = 0.0 is potentially harmful but so far,
            /// no negative (or negligible) effects have been observed. A safer
            /// way to do this (that would even relief us from the burden of
            /// having to evaluate the level set in the leaves), would be to
            /// simply consider all leaves as cut. However, numerical tests
            /// suggest this is more expensive (since more leaves are
            /// considered cut) and does not really improve the accuracy.
            /// </remarks>
            public override bool IsCut {
                get {
                    //return (level == owner.depth);
                    return IsPotentiallyCut(0.0);
                }
            }

            /// <summary>
            /// Checks whether this node is cut by the interface. That is, it
            /// checks if the signs of the distances in all vertices are the
            /// same _and_ if all nodes have an absolute distance greater than
            /// <paramref name="minDistance"/>.
            /// </summary>
            /// <param name="minDistance">
            /// The minimal absolute distance a node may have before it is
            /// considered cut.
            /// </param>
            /// <returns>
            /// True, if the distances of the signs of the distances of the
            /// vertices vary _or_ at least one vertex has an absolute distance
            /// smaller than <paramref name="minDistance"/>. Otherwise, false
            /// is returned.
            /// </returns>
            private bool IsPotentiallyCut(double minDistance) {
                int sign = Math.Sign(distances.FirstOrDefault());
                return distances.Any((value) => Math.Sign(value) != sign)
                    || distances.Any(value => Math.Abs(value) < minDistance);
            }

            /// <summary>
            /// Subdivides the leaves of thee tree.
            /// </summary>
            /// <param name="leavesLevel">
            /// <see cref="SubdivideCutLeaves"/>.
            /// </param>
            /// <param name="refinedVertexSet">
            /// <see cref="SubdivideCutLeaves"/>.
            /// </param>
            /// <param name="checkIsCut">
            /// If true, leaves will only be subdivided if <see cref="IsCut"/>
            /// returns true.
            /// </param>
            /// <param name="minDistance">
            /// See <see cref="IsCut"/>.
            /// </param>
            /// <returns>
            /// True, if at least one node has been subdivided. Otherwise,
            /// false is returned.
            /// </returns>
            private bool SubdivideLeaves(int leavesLevel, NestedVertexSet refinedVertexSet, bool checkIsCut, double minDistance) {
                bool result = false;
                if (level < leavesLevel) {
                    if (Children.Count > 0) {
                        foreach (Node child in Children) {
                            result |= child.SubdivideLeaves(leavesLevel, refinedVertexSet, checkIsCut, minDistance);
                        }
                    }
                } else if (level == leavesLevel) {
                    Debug.Assert(Children.Count == 0, "Can only subdivide leaves");

                    if (!checkIsCut || (checkIsCut && IsPotentiallyCut(minDistance))) {
                        AffineTrafo[] transformations = owner.refElement.GetSubdivision();

                        int N = globalVertexIndices.Length;
                        int D = owner.refElement.SpatialDimension;
                        foreach (AffineTrafo elementaryTransformation in transformations) {
                            AffineTrafo combinedTransformation = Transformation * elementaryTransformation;

                            int[] transfomedVertexIndices = new int[N];
                            for (int i = 0; i < N; i++) {
                                double[] vertex = owner.refElement.Vertices.GetRow(i);
                                vertex = combinedTransformation.Transform(vertex);

                                transfomedVertexIndices[i] = refinedVertexSet.RegisterVertex(vertex);
                            }
                            Children.Add(new Node(
                                owner, this, level + 1, refinedVertexSet, transfomedVertexIndices, combinedTransformation));
                            result = true;
                        }
                    } else {
                        return false;
                    }
                } else {
                    throw new ArgumentException("Given leaves level is higher than the depth of the tree", "leavesLevel");
                }

                return result;
            }
        }
    }
}
