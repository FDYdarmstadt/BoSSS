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
using System.Linq;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// a container for quadrature rules
    /// </summary>
    public class QuadRule : ICloneable, IEquatable<QuadRule> {


        /// <summary>
        /// creates a empty (i.e. all entries set to 0.0) <see cref="QuadRule"/> object 
        /// with <paramref name="noOfNodes"/> quadrature nodes.
        /// </summary>
        /// <param name="noOfNodes"></param>
        /// <param name="D">spatial dimension</param>
        /// <param name="Kref">reference element for which this quadrature rule is valid</param>
        /// <returns>an empty (i.e. all weights are 0.0) quadrature rule</returns>
        public static QuadRule CreateEmpty(RefElement Kref, int noOfNodes, int D) {
            QuadRule ret = new QuadRule();
            ret.Nodes = new NodeSet(Kref, noOfNodes, D);
            ret.Weights = MultidimensionalArray.Create(noOfNodes);
            ret.OrderOfPrecision = 0;
            return ret;
        }

        /// <summary>
        /// the reference element on which this rule is valid
        /// </summary>
        public RefElement RefElement {
            get {
                return this.Nodes.RefElement;
            }
        }
        
        /// <summary>
        /// number of nodes 
        /// </summary>
        public int NoOfNodes {
            get {
                Debug.Assert(Weights.Length == this.Nodes.GetLength(0));
                return Weights.Length;
            }
        }

        /// <summary>
        /// Nodes (i.e. Points in the domain at which the integrand is evaluated);<br/>
        /// 1st index: node index; <br/>
        /// 2nd index: spatial coordinate index,
        /// 0 for 1D and 0,1 for 2D and 0,1,2 for 3D;
        /// </summary>
        public NodeSet Nodes;
        
        /// <summary>
        /// quadrature weights; 
        /// 
        /// index corresponds with 1st index of <see cref="Nodes"/>;
        /// </summary>
        public MultidimensionalArray Weights;

        /// <summary>
        /// if applicable, the maximum degree of a polynomial which is integrated
        /// exactly by this quadrature rule.
        /// </summary>
        public int OrderOfPrecision;

        /// <summary>
        /// spatial dimension of quad rule.
        /// </summary>
        public int SpatialDim {
            get {
                if (OrderOfPrecision == int.MaxValue) {
                    return 0;
                } else {
                    return Nodes.GetLength(1);
                }
            }
        }

        /// <summary>
        /// See <see cref="Equals(QuadRule)"/>
        /// </summary>
        /// <param name="obj">See <see cref="object.Equals(object)"/></param>
        /// <returns>See <see cref="Equals(QuadRule)"/></returns>
        public override bool Equals(object obj) {
            if (obj is QuadRule) {
                return Equals((QuadRule)obj);
            } else {
                return false;
            }
        }

        /// <summary>
        /// Creates a hash based on the hashes of <see cref="Weights"/> and
        /// <see cref="Nodes"/>
        /// </summary>
        /// <returns>
        /// A hash code for this object
        /// </returns>
        public override int GetHashCode() {
            unchecked {
                int hash = 17;
                hash = hash * 23 + Weights.GetHashCode();
                hash = hash * 23 * Nodes.GetHashCode();
                return hash;
            }
        }

        #region ICloneable Members

        /// <summary>
        /// Creates a deep copy of this object.
        /// </summary>
        /// <returns>
        /// An independent clone of this object.
        /// </returns>
        public virtual object Clone() {
            return new QuadRule() {
                Nodes = this.Nodes.CloneAs(),
                OrderOfPrecision = this.OrderOfPrecision,
                Weights = this.Weights.CloneAs()
            };
        }

        #endregion

        #region IEquatable<QuadRule> Members

        /// <summary>
        /// Checks if two quad rules can be considered equal. This is true, if
        /// they have the same <see cref="Weights"/> and the same
        /// <see cref="Nodes"/>
        /// </summary>
        /// <param name="other">
        /// The quad rule for which equality should be checked
        /// </param>
        /// <returns>
        /// True, if <see cref="Weights"/> and <see cref="Nodes"/> are equal.
        /// Otherwise, false is returned.
        /// </returns>
        public bool Equals(QuadRule other) {

            // possibly the fastest test for equality
            if (object.ReferenceEquals(this, other))
                return true;

            if (other == null) {
                return false;
            }

            if (!this.RefElement.Equals(other.RefElement))
                return false;

            if (NoOfNodes != other.NoOfNodes) {
                return false;
            }

            // try to detect inequality fast:
            if (this.Weights[0] != other.Weights[0])
                return false;

            // try to detect inequality fast:
            if (this.Nodes[0, 0] != other.Nodes[0, 0])
                return false;

            if (Weights.Equals(other.Weights) && Nodes.Equals(other.Nodes)) {
                return true;
            }

            return false;

        }

        #endregion
    }
}
