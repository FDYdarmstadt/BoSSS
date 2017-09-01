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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// Helper class that wraps a given quadrature rule into a
    /// <see cref="IQuadRuleFactory{QuadRule}"/>. Mainly exists to support
    /// legacy code. New code should directly provide a quad rule factory.
    /// </summary>
    public class FixedRuleFactory<T> : IQuadRuleFactory<T> where T : QuadRule {

        /// <summary>
        /// Constructs the factory.
        /// </summary>
        public FixedRuleFactory(T qr) {
            RefElement = qr.RefElement;
            m_qr = qr;
        }

        T m_qr;

        #region IQuadRuleFactory<QuadRule> Members

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        /// <summary>
        /// Not implemented since it isn't required.
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        /// <summary>
        /// Returns the quad rule supplied to the constructor for every
        /// single chunk in <paramref name="mask"/>.
        /// </summary>
        /// <param name="mask">
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </param>
        /// <param name="order">
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </param>
        /// <returns>
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </returns>
        public IEnumerable<IChunkRulePair<T>> GetQuadRuleSet(ExecutionMask mask, int order) {
            var result = new List<IChunkRulePair<T>>(mask.Count());
            var quadRule = m_qr;
            foreach (Chunk chunk in mask) {
                result.Add(new ChunkRulePair<T>(chunk, quadRule));
            }
            return result;
        }

        #endregion
    }
}
