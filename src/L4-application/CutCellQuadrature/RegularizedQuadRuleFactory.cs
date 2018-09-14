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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace CutCellQuadrature {

    class RegularizedQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        private IQuadRuleFactory<QuadRule> baseFactory;

        private LevelSetTracker tracker;

        private int levSetIndex;

        private RegularizationPolynomoial polynomial;

        private double width;

        public RegularizedQuadRuleFactory(IQuadRuleFactory<QuadRule> baseFactory, LevelSetTracker tracker, RegularizationPolynomoial polynomial, double width)
            : this(baseFactory, tracker, 0, polynomial, width) { }

        public RegularizedQuadRuleFactory(IQuadRuleFactory<QuadRule> baseFactory, LevelSetTracker tracker, int levSetIndex, RegularizationPolynomoial polynomial, double width) {
            this.baseFactory = baseFactory;
            this.polynomial = polynomial;
            this.width = width;
            this.tracker = tracker;
            if (tracker.LevelSets.Count <= levSetIndex)
                throw new ArgumentOutOfRangeException("Please specify a valid index for the level set.");
            this.levSetIndex=levSetIndex;
        }

        #region IQuadRuleFactory Members

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        public RefElement RefElement {
            get {
                return baseFactory.RefElement;
            }
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");
            var baseSet = baseFactory.GetQuadRuleSet(mask, order);

            var result = new List<ChunkRulePair<QuadRule>>();
            foreach (var chunkRulePair in baseSet) {
                MultidimensionalArray levelSetValues = tracker.DataHistories[levSetIndex].Current.GetLevSetValues(
                    chunkRulePair.Rule.Nodes, chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len);
                MultidimensionalArray gradientValues = tracker.DataHistories[levSetIndex].Current.GetLevelSetGradients(
                     chunkRulePair.Rule.Nodes, chunkRulePair.Chunk.i0, chunkRulePair.Chunk.Len);
                
                foreach (int cell in chunkRulePair.Chunk.Elements) {
                    QuadRule rule = chunkRulePair.Rule.CloneAs();

                    for (int i = 0; i < rule.NoOfNodes; i++) {
                        rule.Weights[i] *= polynomial.Evaluate(levelSetValues[cell - chunkRulePair.Chunk.i0, i], width);

                        // Scale by the norm of the gradient in order to be on
                        // the safe side with non-signed-distance functions
                        double norm = 0.0;
                        for (int j = 0; j < RefElement.SpatialDimension; j++) {
                            int index = cell - chunkRulePair.Chunk.i0;
                            norm += gradientValues[index, i, j] * gradientValues[index, i, j];
                        }
                        rule.Weights[i] *= Math.Sqrt(norm);
                    }

                    result.Add(new ChunkRulePair<QuadRule>(
                        Chunk.GetSingleElementChunk(cell), rule));
                }
            }

            return result;
        }

        #endregion
    }
}
