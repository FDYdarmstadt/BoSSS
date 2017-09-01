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

using System.Collections.Generic;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.LinAlg;
using System;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// Defines the interface for algorithms that determine how to subdivide a
    /// simplex into sub-simplices.
    /// </summary>
    public interface ISubdivisionStrategy {

        /// <summary>
        /// The simplex to be subdivided.
        /// </summary>
        RefElement RefElement {
            get;
        }

        /// <summary>
        /// Constructs the nodes that divide the subdivisions for all elements
        /// of <paramref name="mask"/>. Depending on the strategy, the
        /// subdivision may or may not vary from element to element.
        /// </summary>
        /// <param name="mask">
        /// A mask containing all elements (cells or edges) for which the nodes
        /// should be computed.
        /// </param>
        /// <returns>
        /// The subdivision nodes for each element in <paramref name="mask"/>,
        /// where elements with the same set of subdivision nodes are grouped
        /// into a chunk. The chunks are subsets are subsets of the chunk of
        /// the given mask.
        /// </returns>
        IEnumerable<KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>> GetSubdivisionNodes(ExecutionMask mask);
    }

    /// <summary>
    /// A subdivision node that maps the nodes of a <see cref="RefElement"/>
    /// (see <see cref="ISubdivisionStrategy.RefElement"/>) to some sub-simplex.
    /// The set of all subdivision nodes for one element (i.e., either cell or
    /// edge) must cover the whole simplex _and_ the sub-simplices must be
    /// disjoint.
    /// </summary>
    public class SubdivisionNode : IEquatable<SubdivisionNode> {

        /// <summary>
        /// Cache for <see cref="NullSubdivisionNode"/>.
        /// </summary>
        private static SubdivisionNode[] nullSubdivisionsNodes = new SubdivisionNode[] {
            new SubdivisionNode(AffineTrafo.Identity(2)), // 1D
            new SubdivisionNode(AffineTrafo.Identity(3)), // 2D
            new SubdivisionNode(AffineTrafo.Identity(4)), // 3D
        };

        /// <summary>
        /// Constructs an uncut subdivision node defined by the given
        /// transformation.
        /// </summary>
        /// <param name="transformation">
        /// The transformation that transforms the vertices of the original
        /// simplex to the vertices of the sub-simplex represented by this
        /// node.
        /// </param>
        public SubdivisionNode(AffineTrafo transformation)
            : this(transformation, false) {
        }

        /// <summary>
        /// Constructs a subdivision node defined by the given transformation.
        /// </summary>
        /// <param name="transformation">
        /// The transformation that transforms the vertices of the original
        /// simplex to the vertices of the sub-simplex represented by this
        /// node.
        /// </param>
        /// <param name="isCut">
        /// Defines whether this node is considered is cut by the interface.
        /// </param>
        public SubdivisionNode(AffineTrafo transformation, bool isCut) {
            Transformation = transformation;
            IsCut = isCut;
        }

        /// <summary>
        /// The transformation that transforms the vertices of the original
        /// simplex to the vertices of the sub-simplex represented by this
        /// node.
        /// </summary>
        public AffineTrafo Transformation {
            get;
            private set;
        }

        /// <summary>
        /// Returns true of the given subdivision node is (potentially) cut by
        /// the interface
        /// </summary>
        /// <returns>
        /// True if the nodes is potentially cut, false otherwise
        /// </returns>
        public virtual bool IsCut {
            get;
            private set;
        }

        /// <summary>
        /// The (<paramref name="spatialDimension"/> + 1)-dimensional identity
        /// transformation representing a node which does not subdivide the
        /// original simplex.
        /// </summary>
        /// <param name="spatialDimension">
        /// The spatial dimension of the simplex to (not) divide.
        /// </param>
        /// <returns>
        /// A subdivision node with an identity transformation.
        /// </returns>
        public static SubdivisionNode NullSubdivisionNode(int spatialDimension) {
            return nullSubdivisionsNodes[spatialDimension - 1];
        }

        /// <summary>
        /// See <see cref="Equals(SubdivisionNode)"/>
        /// </summary>
        /// <param name="obj">See <see cref="object.Equals(object)"/></param>
        /// <returns>See <see cref="Equals(SubdivisionNode)"/></returns>
        public override bool Equals(object obj) {
            if (obj is SubdivisionNode) {
                return Equals((SubdivisionNode)obj);
            } else {
                return false;
            }
        }

        /// <summary>
        /// Creates a hash code based on the hash codes of <see cref="IsCut"/>
        /// and <see cref="Transformation"/>.
        /// </summary>
        /// <returns>A hash code</returns>
        public override int GetHashCode() {
            return IsCut.GetHashCode() ^ Transformation.GetHashCode();
        }

        #region IEquatable<SubdivisionNode> Members

        /// <summary>
        /// Equality check.
        /// </summary>
        /// <param name="other">
        /// The object to be checked against.
        /// </param>
        /// <returns>
        /// True, if both objects have the same <see cref="Transformation"/>
        /// and the same value for <see cref="IsCut"/>. False, otherwise.
        /// </returns>
        public bool Equals(SubdivisionNode other) {
            if (other == null) {
                return false;
            }

            if (ReferenceEquals(this, other)) {
                return true;
            }

            if (IsCut != other.IsCut) {
                return false;
            }

            return other.Transformation.Equals(other.Transformation);
        }

        #endregion
    }

}
