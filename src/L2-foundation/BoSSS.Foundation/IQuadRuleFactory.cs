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

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// Common interface for all providers of quadrature rules.
    /// </summary>
    /// <typeparam name="TQuadRule">
    /// The type of the quadrature rules to be created by the factory.
    /// </typeparam>
    /// <remarks>
    /// The interface has been designed such that it is covariant in its type
    /// parameter (hence the "out" keyword). As a result, assignments of the
    /// type
    /// <code>
    /// IQuadRuleFactory&lt;CellBndQuadRulegt; factory = someFactory;
    /// IQuadRuleFactory&lt;QuadRulegt variableUsingBasetype = factory;
    /// </code>
    /// are perfectly valid. To achieve this, the signature of
    /// <see cref="GetQuadRuleSet"/> had to been changed to a type that is
    /// also covariant in its type parameter (which is not true, for instance,
    /// in case of <see cref="IDictionary{S,T}"/>.
    /// </remarks>
    public interface IQuadRuleFactory<out TQuadRule> where TQuadRule : QuadRule {

        /// <summary>
        /// The simplex to be integrated over.
        /// </summary>
        Grid.RefElements.RefElement RefElement {
            get;
        }

        /// <summary>
        /// Creates the quadrature rules for all elements (i.e., cell or edges)
        /// in <paramref name="mask"/>. The used quadrature rules have an order
        /// of <paramref name="order"/>. The chunks used as keys of the
        /// resulting dictionary are subsets of the chunks of
        /// <paramref name="mask"/>. Depending on the implementation, the
        /// quadrature rule may or may not vary from cell to cell which implies
        /// that the chunks may or may not be equal to the chunks of
        /// <paramref name="mask"/>.
        /// </summary>
        /// <param name="mask">
        /// Mask containing the elements for which the quadrature rules are
        /// needed.
        /// </param>
        /// <param name="order">
        /// The desired order of the integration rules.
        /// </param>
        /// <returns>
        /// A mapping between chunks of elements and the integration rule
        /// proposed by the specific implementation.
        /// </returns>
        IEnumerable<IChunkRulePair<TQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order);


        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        int[] GetCachedRuleOrders();
    }

    /// <summary>
    /// A tuple consisting of a chunk (a subset of a domain of integration) and
    /// an associated quadrature rule.
    /// </summary>
    /// <typeparam name="TQuadRule">
    /// The type of the quadrature rule
    /// </typeparam>
    /// <remarks>
    /// Introduced since <see cref="KeyValuePair{S,T}"/> is not covariant in
    /// the type parameter T.
    /// </remarks>
    public interface IChunkRulePair<out TQuadRule> where TQuadRule : QuadRule {

        /// <summary>
        /// A chunk of integration elements.
        /// </summary>
        Chunk Chunk {
            get;
        }

        /// <summary>
        /// A quadrature rule valid for all elements of <see name="Chunk"/>.
        /// </summary>
        TQuadRule Rule {
            get;
        }
    }

    /// <summary>
    /// Represents an immutable tuple of a chunk and an associated quad rule.
    /// </summary>
    /// <typeparam name="TQuadRule">
    /// The type of quadrature rule to be represented
    /// </typeparam>
    public class ChunkRulePair<TQuadRule> : IChunkRulePair<TQuadRule> where TQuadRule : QuadRule {

        /// <summary>
        /// Just stores the given values
        /// </summary>
        public ChunkRulePair(Chunk chunk, TQuadRule rule) {
            this._Chunk = chunk;
            this.Rule = rule;
        }

        /// <summary>
        /// Half-Shallow Clone
        /// </summary>
        public ChunkRulePair(IChunkRulePair<TQuadRule> org) {
            this._Chunk = org.Chunk;
            this.Rule = org.Rule;
        }

        #region IChunkRulePair<TQuadRule> Members

        /// <summary>
        /// <see cref="IChunkRulePair{TQuadRule}.Chunk"/>
        /// </summary>
        public Chunk Chunk {
            get { return _Chunk; }
        }

        internal Chunk _Chunk;

        /// <summary>
        /// <see cref="IChunkRulePair{TQuadRule}.Rule"/>
        /// </summary>
        public TQuadRule Rule {
            get;
            internal set;
        }

        #endregion
    }
}
