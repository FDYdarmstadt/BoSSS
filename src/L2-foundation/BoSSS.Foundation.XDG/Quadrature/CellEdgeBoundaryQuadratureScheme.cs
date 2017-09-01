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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// A quadrature scheme for <see cref="CellEdgeBoundaryQuadRule"/>
    /// </summary>
    public class CellEdgeBoundaryQuadratureScheme : QuadratureScheme<CellEdgeBoundaryQuadRule, CellMask> {

        /// <summary>
        /// Default constructor.
        /// </summary>
        /// <param name="useDefaultFactories"></param>
        /// <param name="domain"></param>
        public CellEdgeBoundaryQuadratureScheme(bool useDefaultFactories, CellMask domain = null)
            : base(useDefaultFactories, domain) {
        }

        /// <summary>
        /// Alternative constructor that saves a single call to
        /// <see cref="QuadratureScheme{S,T}.AddFactoryDomainPair"/>
        /// </summary>
        /// <param name="useDefaultFactories"></param>
        /// <param name="factory"></param>
        /// <param name="domain"></param>
        public CellEdgeBoundaryQuadratureScheme(
            bool useDefaultFactories,
            IQuadRuleFactory<CellEdgeBoundaryQuadRule> factory,
            CellMask domain = null)
            : base(useDefaultFactories, domain) {
            AddFactoryDomainPair(factory, domain);
        }

        /// <summary>
        /// Returns <see cref="CellMask.GetFullMask"/>
        /// </summary>
        /// <param name="gridData"></param>
        /// <returns></returns>
        protected override CellMask GetDefaultDomain(IGridData gridData) {
            return CellMask.GetFullMask(gridData);
        }

        /// <summary>
        /// Returns an instance of
        /// <see cref="CellBoundaryFromEdgeRuleFactory{CellEdgeBoundaryQuadRule}"/>
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="elem"></param>
        /// <returns></returns>
        protected override IQuadRuleFactory<CellEdgeBoundaryQuadRule> GetDefaultRuleFactory(
            IGridData gridData, RefElement elem) {
            return new CellBoundaryFromEdgeRuleFactory<CellEdgeBoundaryQuadRule>(
                 gridData, elem, new StandardQuadRuleFactory(elem.FaceRefElement));
        }
    }
}
