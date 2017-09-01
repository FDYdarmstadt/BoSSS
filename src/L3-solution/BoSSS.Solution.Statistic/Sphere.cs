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
    public class Sphere : Surface, ISurfaceEvaluation {
        /// <summary>
        /// Radius of the sphere with center (0,0,0)
        /// </summary>
        double r;
        /// <summary>
        /// Creates a sphere with options for evaluation on the surface
        /// </summary>
        /// <param name="cont">The context object</param>
        /// <param name="b">Basis for the field representation</param>
        /// <param name="deltax">Parameter for the approximation of the signum function</param>
        /// <param name="radius">Radius of the sphere- note that the center is (0,0,0) </param>

        public Sphere(GridData cont, Basis b, double deltax, double radius)
            : base(cont, b, deltax) {
            if (r < 0.0)
                throw new ArgumentOutOfRangeException("The radius of the sphere has to be assigned a positive value");
            r = radius;
            m_Field.ProjectField(LevSetInit);
            NormalVec(deltax);
        }

        /// <summary>
        /// 
        /// </summary>
        public double Radius {

            get {
                return r;
            }
            set {
                if (value < 0.0)
                    throw new ArgumentOutOfRangeException("The radius of the sphere has to be assigned a positive value");
                r = value;
            }
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
                outp[i] = r - Math.Sqrt(x * x + y * y + z * z);
            }

        }

        /// <summary>
        /// Method that creates (<paramref name="nw"/>)^2 quadrature nodes and weights on the surface of the sphere
        /// </summary>
        /// <param name="testnodes">Nodes on the surface</param>
        /// <param name="quadwghts">Corresponding quadrature weights</param>
        /// <param name="nw">Square root of the number of nodes, or weights respectively</param>
        public override void CreateNodesAndWeights(out double[,] testnodes, out double[] quadwghts, int nw) {


            testnodes = new double[nw * nw, 3];
            int counts = 0;
            quadwghts = new double[nw * nw];
            for (int i = 0; i < nw; i++) {

                for (int k = 0; k < nw; k++) {
                    double kd = System.Convert.ToDouble(k);
                    double id = System.Convert.ToDouble(i);
                    double nwd = System.Convert.ToDouble(nw);
                    testnodes[counts + k, 0] = r * Math.Sin(kd * Math.PI / nwd) * Math.Cos(2.0 * id * Math.PI / nwd);
                    testnodes[counts + k, 1] = r * Math.Sin(kd * Math.PI / nwd) * Math.Sin(2.0 * id * Math.PI / nwd);
                    testnodes[counts + k, 2] = r * Math.Cos(kd * Math.PI / nwd);
                    /*   if (i == nw - 1) {
                           quadwghts[counts + k] = (2.0 * r * r / nwd) * Math.PI * Math.Abs(Math.Abs(Math.Cos(kd * Math.PI / nwd)) - Math.Abs(Math.Cos((kd - 0.5) * Math.PI / nwd)));
                       }
                       else if (i == nw - 2) {
                           quadwghts[counts + k] = (2.0 * r * r / nwd) * Math.PI * Math.Abs(Math.Abs(Math.Cos((kd + 0.5) * Math.PI / nwd)) - Math.Abs(Math.Cos(kd * Math.PI / nwd)));
                       }
                       else {*/
                    quadwghts[counts + k] = (2.0 * r * r / nwd) * Math.PI * r * Math.Sin(kd * Math.PI / nwd) * (Math.Abs(Math.Abs(Math.Cos((kd + 1.0) * Math.PI / nwd)) - Math.Abs(Math.Cos(kd * Math.PI / nwd))));

                    //}

                }
                counts += nw;
            }
        }

        /// <summary>
        /// Method that creates (<paramref name="nw"/>)^2 quadrature nodes and weights on the surface of a sphere shifting its center....
        /// this method probably dosn#t work out. Come up with something else .......
        /// </summary>
      
        public  void CreateNodesAndWeights(out double[,] testnodes, out double[] quadwghts, int nw, double[] shiftedCenter) {


            testnodes = new double[nw * nw, 3];
            int counts = 0;
            quadwghts = new double[nw * nw];
            for (int i = 0; i < nw; i++) {

                for (int k = 0; k < nw; k++) {
                    double kd = System.Convert.ToDouble(k);
                    double id = System.Convert.ToDouble(i);
                    double nwd = System.Convert.ToDouble(nw);
                    testnodes[counts + k, 0] = r * Math.Sin(kd * Math.PI / nwd) * Math.Cos(2.0 * id * Math.PI / nwd)-shiftedCenter[0];
                    testnodes[counts + k, 1] = r * Math.Sin(kd * Math.PI / nwd) * Math.Sin(2.0 * id * Math.PI / nwd)-shiftedCenter[1];
                    testnodes[counts + k, 2] = r * Math.Cos(kd * Math.PI / nwd)-shiftedCenter[2];
                 
                    quadwghts[counts + k] = (2.0 * r * r / nwd) * Math.PI * Math.Abs(Math.Abs(Math.Cos((kd + 1.0) * Math.PI / nwd)) - Math.Abs(Math.Cos(kd * Math.PI / nwd)));

                    

                }
                counts += nw;
            }
        }



    }
}
