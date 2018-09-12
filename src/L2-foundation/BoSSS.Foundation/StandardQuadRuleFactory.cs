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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using System;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// Standard quad rule factory that constructs quadrature rules that are
    /// suitable for continuous integrands.
    /// </summary>
    public class StandardQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        /// <summary>
        /// Constructs a new factory for the given simplex.
        /// </summary>
        /// <param name="simplex">
        /// The simplex to be integrated over.
        /// </param>
        public StandardQuadRuleFactory(RefElement simplex) {
            RefElement = simplex;
        }

        #region IQuadRuleFactory<QuadRule> Members

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        /// <summary>
        /// Uses <see cref="Grid.RefElement.GetQuadratureRule"/> to create a quad rule
        /// (i.e., the quad rule is the same for all elements of
        /// <paramref name="mask"/>)
        /// </summary>
        /// <param name="mask">
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </param>
        /// <param name="order">
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </param>
        /// <returns>
        /// <see cref="Grid.RefElement.GetQuadratureRule"/>
        /// </returns>
        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

            QuadRule rule = RefElement.GetQuadratureRule(order);
            Debug.Assert(rule.Nodes.IsLocked, "Error: non-locked quad rule from reference element.");
            Debug.Assert(object.ReferenceEquals(rule.Nodes.RefElement, RefElement), "Error: quad rule from reference element has not the right reference elment assigned.");

            return mask.Select(chunk => new ChunkRulePair<QuadRule>(chunk, rule));
        }

        /// <summary>
        /// <see cref="IQuadRuleFactory{QuadRule}.RefElement"/>
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        #endregion
    }


    /// <summary>
    /// Standard quad rule factory that constructs quadrature rules that are
    /// suitable for continuous integrands.
    /// </summary>
    public class StandardCellBoundaryQuadRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {

        /// <summary>
        /// Constructs a new factory for the given simplex.
        /// </summary>
        /// <param name="simplex">
        /// The simplex to be integrated over.
        /// </param>
        public StandardCellBoundaryQuadRuleFactory(RefElement simplex) {
            RefElement = simplex;
        }

        #region IQuadRuleFactory<QuadRule> Members

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        /// <summary>
        /// Uses <see cref="Grid.RefElement.GetQuadratureRule"/> to create a quad rule
        /// (i.e., the quad rule is the same for all elements of
        /// <paramref name="mask"/>)
        /// </summary>
        /// <param name="mask">
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </param>
        /// <param name="order">
        /// <see cref="IQuadRuleFactory{QuadRule}.GetQuadRuleSet"/>
        /// </param>
        /// <returns>
        /// <see cref="Grid.RefElement.GetQuadratureRule"/>
        /// </returns>
        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            CellBoundaryQuadRule rule = RefElement.GetBoundaryQuadRule(order);
            Debug.Assert(rule.Nodes.IsLocked, "Error: non-locked quad rule from reference element.");
            Debug.Assert(object.ReferenceEquals(rule.Nodes.RefElement, RefElement), "Error: quad rule from reference element has not the right reference element assigned.");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");


            return mask.Select(chunk => new ChunkRulePair<CellBoundaryQuadRule>(chunk, rule));
        }

        /// <summary>
        /// <see cref="IQuadRuleFactory{QuadRule}.RefElement"/>
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        #endregion
    }
}
