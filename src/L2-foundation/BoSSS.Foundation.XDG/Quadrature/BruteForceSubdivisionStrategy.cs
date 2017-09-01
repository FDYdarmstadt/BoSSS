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
using System.Linq;
using BoSSS.Foundation.Grid;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// Subdivision strategy that blindly subdivides a simplex into
    /// sub-simplices up to a given level. As a result, the subdivisions are
    /// the same for all elements (i.e., cells or edges).
    /// </summary>
    public class BruteForceSubdivisionStrategy : ISubdivisionStrategy {

        /// <summary>
        /// The number of subdivisions of the original simplex to perform.
        /// </summary>
        private readonly int numberOfSubdivisions;

        /// <summary>
        /// Initializes the strategy.
        /// </summary>
        /// <param name="refElement">
        /// The simplex to be subdivided.
        /// </param>
        /// <param name="numberOfSubdivisions">
        /// The number of subdivisions, see
        /// <see cref="Grid.RefElement.SubdivisionTree"/>.
        /// </param>
        public BruteForceSubdivisionStrategy(RefElement refElement, int numberOfSubdivisions) {
            this.RefElement = refElement;
            this.numberOfSubdivisions = numberOfSubdivisions;
        }

        #region ISubdivisionStrategy Members

        /// <summary>
        /// The simplex to be subdivided.
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        /// <summary>
        /// Uses <see cref="Grid.RefElement.SubdivisionTree"/> to create the brute
        /// force subdivision which is the same for all cells. As a result, the
        /// chunks of <paramref name="mask"/> are not altered.
        /// </summary>
        /// <param name="mask">
        /// <see cref="ISubdivisionStrategy.GetSubdivisionNodes"/>
        /// </param>
        /// <returns>
        /// <see cref="ISubdivisionStrategy.GetSubdivisionNodes"/>
        /// </returns>
        /// <remarks>
        /// Currently, the nodes returned by
        /// <see cref="Grid.RefElement.SubdivisionTree"/> are unnecessarily
        /// boxed by this method. This should be changed by changing the return
        /// type of <see cref="Grid.RefElement.SubdivisionTree"/> at some
        /// point.
        /// </remarks>
        public IEnumerable<KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>> GetSubdivisionNodes(ExecutionMask mask) {
            using (new FuncTrace()) {
                RefElement.SubdivisionTree[] leaves = RefElement.GetSubdivisionTree(numberOfSubdivisions).GetLeaves();

                foreach (Chunk chunk in mask) {
                    yield return new KeyValuePair<Chunk, IEnumerable<SubdivisionNode>>(
                        chunk,
                        leaves.Select((leave) => new SubdivisionNode(leave.TrafoFromRoot, true)));
                }
            }
        }

        #endregion
    }
}
