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
    public class Torus : Surface, ISurfaceEvaluation {
        /// <summary>
        /// Major radius of the torus
        /// </summary>
        double R;

        /// <summary>
        /// Minor radius of the torus
        /// </summary>
        double r;

        /// <summary>
        /// Creates a torus object
        /// </summary>
        /// <param name="cont">The context object</param>
        /// <param name="b">Basis for the field representation</param>
        /// <param name="deltax">Parameter for the approximation of the signum function</param>
        /// <param name="majorRadius">Major radius of the torus</param>
        /// <param name="minorRadius">Minor radius of the torus</param>
        public Torus(GridData cont, Basis b, double deltax, double majorRadius, double minorRadius)
            : base(cont, b, deltax) {

            if ((r < 0.0) || (R < 0.0))
                throw new ArgumentOutOfRangeException("Both of the radii have to be assigned a positive values");

            R = majorRadius;
            r = minorRadius;
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
                outp[i] = (Math.Sqrt((Math.Sqrt(x * x + y * y) - R) * (Math.Sqrt(x * x + y * y) - R) + z * z) - r);
            }

        }
        /// <summary>
        /// Method that creates (<paramref name="nw"/>)^2 quadrature nodes and weights on the surface of the torus
        /// </summary>
        /// <param name="testnodes">Nodes on the surface</param>
        /// <param name="quadwghts">Corresponding quadrature weights</param>
        /// <param name="nw">Square root of the number of nodes, or weights respectively</param>
        public override void CreateNodesAndWeights(out double[,] testnodes, out double[] quadwghts, int nw) {

            testnodes = new double[nw * nw, 3];
            quadwghts = new double[nw * nw];
            int counts = 0;
            for (int i = 0; i < nw; i++) {

                for (int k = 0; k < nw; k++) {

                    testnodes[counts + k, 0] = (R + r * Math.Sin(k * 2.0 * Math.PI / (double)nw)) * Math.Cos(2.0 * (double)i * Math.PI / (double)nw);
                    testnodes[counts + k, 1] = (R + r * Math.Sin(k * 2.0 * Math.PI / (double)nw)) * Math.Sin(2.0 * (double)i * Math.PI / (double)nw);
                    testnodes[counts + k, 2] = r * Math.Cos(k * 2.0 * Math.PI / (double)nw);
                    quadwghts[k + counts] = Math.Abs((2.0 * (r * r) / (double)nw) * Math.PI * (Math.Sin(((double)k + 1.0) * 2.0 * Math.PI / (double)nw) -
                        Math.Sin((double)k * 2.0 * Math.PI / (double)nw))
                         + 4.0 * Math.PI * r * R * Math.PI / (double)(nw * nw));


                }
                counts += nw;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public double majorRadius {

            get {
                return R;
            }
            set {
                if (value < 0.0)
                    throw new ArgumentOutOfRangeException("The radius of the sphere has to be assigned a positive value");
                R = value;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public double minorRadius {

            get {
                return r;
            }
            set {
                if (value < 0.0)
                    throw new ArgumentOutOfRangeException("The radius of the sphere has to be assigned a positive value");
                r = value;
            }
        }



    }
}
