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
using System.Linq;
using System.Text;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Quadrature {
    
    /// <summary>
    /// Variant of <see cref="QuadRule"/> that additionally stores the
    /// number of nodes per edge of the edge rules that have been used to
    /// construct this rule.
    /// </summary>
    public class CellBoundaryQuadRule : QuadRule {

        /// <summary>
        /// Creates an empty rule.
        /// </summary>
        public static CellBoundaryQuadRule CreateEmpty(RefElement Kref, int noOfNodes, int D, int noOfEdges) {
            CellBoundaryQuadRule rule = new CellBoundaryQuadRule();
            rule.Nodes = new NodeSet(Kref, noOfNodes, D);
            rule.Weights = MultidimensionalArray.Create(noOfNodes);
            rule.OrderOfPrecision = 0;
            rule.NumbersOfNodesPerFace = new int[noOfEdges];
            return rule;
        }

        /// <summary>
        /// For each face of the reference element: The number of nodes
        /// of the edge rules that has been used to construct this volume rule.
        /// </summary>
        public int[] NumbersOfNodesPerFace {
            get;
            set;
        }

        /// <summary>
        /// Creates a deep copy of this object
        /// </summary>
        /// <returns>A deep copy of this object</returns>
        public override object Clone() {
            return new CellBoundaryQuadRule() {
                Nodes = this.Nodes.CloneAs(),
                Weights = this.Weights.CloneAs(),
                NumbersOfNodesPerFace = this.NumbersOfNodesPerFace.CloneAs(),
                OrderOfPrecision = this.OrderOfPrecision,
            };
        }
    }


}
