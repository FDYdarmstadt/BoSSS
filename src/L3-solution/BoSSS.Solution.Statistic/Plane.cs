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
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Statistic.QuadRules {

    /// <summary>
    /// 
    /// </summary>
    public class Plane : Surface, ISurfaceEvaluation {
        
        /// <summary>
        /// Creates a plane through z=0
        /// </summary>
        /// <param name="cont">The context object</param>
        /// <param name="b">Basis for the field representation</param>
        /// <param name="deltax">Parameter for the approximation of the signum function</param>
        public Plane(GridData cont, Basis b, double deltax)
            : base(cont, b, deltax) {
            m_Field.ProjectField(LevSetInit);
            NormalVec(deltax);
        }
        /// <summary>
        /// Representation as a scalar function 
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="outp"></param>
        public override void LevSetInit(MultidimensionalArray inp, MultidimensionalArray outp) {
            for (int i = 0; i < inp.GetLength(0); i++) {

                double x = inp[i, 0];
                double y = inp[i, 1];
                double z = inp[i, 2];

                outp[i] = z;

            }

        }
        /// <summary>
        /// Method that creates (<paramref name="nw"/>)^2 quadrature nodes and weights on a plate of measure [-1,1]x[-1,1]
        /// through z=0
        /// </summary>
        /// <param name="testnodes">Nodes on the surface</param>
        /// <param name="quadwghts">Corresponding quadrature weights</param>
        /// <param name="nw">Square root of the number of nodes, or weights respectively</param>
        public override void CreateNodesAndWeights(out double[,] testnodes, out double[] quadwghts, int nw) {

            testnodes = new double[nw * nw, 3];
            quadwghts = new double[nw * nw];
            int c = 0;
            for (int i = 0; i < nw; i++) {
                for (int k = 0; k < nw; k++) {
                    testnodes[c + i, 0] = 0.0;
                    testnodes[c + i, 1] = -1.0 + (double)i * 2.0 / (double)nw;
                    testnodes[c + i, 2] = -1.0 + (double)k * 2.0 / (double)nw;
                    quadwghts[c + i] = (double)4.0 / (double)(nw * nw);
                    c += nw;
                }
            }

        }

    }
}
