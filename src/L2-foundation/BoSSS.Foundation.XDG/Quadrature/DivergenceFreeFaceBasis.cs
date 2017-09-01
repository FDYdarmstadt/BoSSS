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

using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    internal class DivergenceFreeFaceBasis : DivergenceFreeBasis {

        private int localEdge;

        public DivergenceFreeFaceBasis(GridData gridData, RefElement VolKref, int degree, int localEdge)
            : base(gridData, VolKref.FaceRefElement, degree) {
            this.localEdge = localEdge;
        }

        public MultidimensionalArray Evaluate(NodeSet nodes) {
            //Debug.Assert(object.ReferenceEquals(nodes.RefElement, base.RefElement));
            NodeSet edgeNodes = new NodeSet(base.RefElement, nodes.RefElement.GetInverseFaceTrafo(localEdge).Transform(nodes));

            //MultidimensionalArray monomials = Polynomial.GetMonomials(
            //    edgeNodes, m_GridDat.Grid.RefElements[0].FaceRefElement.SpatialDimension, Degree);



            MultidimensionalArray ret = MultidimensionalArray.Create(
                edgeNodes.GetLength(0), base.Count);
            //for(int m = 0; m < base.Count; m++) {
            //    base[m].Evaluate(
            //        ret.ExtractSubArrayShallow(-1, m),
            //        edgeNodes,
            //        monomials);
            //}
            base.Evaluate(edgeNodes, ret);

            return ret;
        }
    }
}
