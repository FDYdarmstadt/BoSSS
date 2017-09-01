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
    /// Variant of <see cref="QuadRule"/> for edge integrals, 
    /// where different nodes/weights can be used for both sides
    /// </summary>
    public class DoubleEdgeQuadRule : QuadRule {


        /// <summary>
        /// Creates an empty rule.
        /// </summary>
        new public static DoubleEdgeQuadRule CreateEmpty(RefElement Kref, int noOfNodes, int D) {
            DoubleEdgeQuadRule rule = new DoubleEdgeQuadRule();
            rule.Nodes = new NodeSet(Kref, noOfNodes, D);
            rule.Weights = MultidimensionalArray.Create(noOfNodes);
            rule.OrderOfPrecision = 0;
            return rule;
        }

        /// <summary>
        /// nodes/weights from index 0 to <see cref="Median"/>-1 regard the in-cell,
        /// indices from <see cref="Median"/> to <see cref="QuadRule.NoOfNodes"/>-1 regard the out cell.
        /// </summary>
        public int Median;


        /// <summary>
        /// Creates a deep copy of this object
        /// </summary>
        /// <returns>A deep copy of this object</returns>
        public override object Clone() {
            return new DoubleEdgeQuadRule() {
                Nodes = this.Nodes.CloneAs(),
                Weights = this.Weights.CloneAs(),
                OrderOfPrecision = this.OrderOfPrecision,
            };
        }
    }
}
