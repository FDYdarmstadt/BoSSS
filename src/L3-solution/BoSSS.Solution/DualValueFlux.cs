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
using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// An abstract baseclass that helps with the implementation of 
    /// <see cref="IDualValueFlux"/>;
    /// Instead of implementing the vectorized functions in the interface, which is a bit complicated but offers better performance,
    /// the user only has to override scalar versions;
    /// </summary>
    public abstract class DualValueFlux : IDualValueFlux {
        #region IDualValueFlux Members

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// </summary>
        public void BorderEdgeFlux(double time, int jEdge, MultidimensionalArray x, MultidimensionalArray normal, bool flipNormal, byte[] EdgeTags, int EdgeTagsOffset, MultidimensionalArray[] Uin, int Offset, int Lenght, MultidimensionalArray Output) {

            // preparation
            // -----------
            int D = normal.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = Uin[0].GetLength(1);
            int L = m_Uin.Length;
            double normSign = flipNormal ? -1 : 1;

            // loop
            // ----

            for (int e = 0; e < Lenght; e++) {
                byte EdgeTag = EdgeTags[e + EdgeTagsOffset];

                for (int n = 0; n < NoOfNodes; n++) {

                    for (int d = 0; d < D; d++) {
                        m_x[d] = x[e + Offset, n, d];
                        m_normal[d] = normal[e + Offset, n, d] * normSign;
                    }
                    for (int l = 0; l < L; l++) m_Uin[l] = Uin[l][e + Offset, n];

                    Output[e + Offset, n] += BorderEdgeFlux(time, m_x, m_normal, EdgeTag, m_Uin, e + jEdge);
                }
            }
        }


        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// </summary>
        public void InnerEdgeFlux(double time, int jEdge, MultidimensionalArray x, MultidimensionalArray normal, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, int Offset, int Lenght, MultidimensionalArray Output_InCell, MultidimensionalArray Output_OutCell) {

            // preparation
            // -----------
            int D = normal.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = Uin[0].GetLength(1);
            int L = m_Uin.Length;

            // loop
            // ----

            for (int e = 0; e < Lenght; e++) {
                for (int n = 0; n < NoOfNodes; n++) {

                    for (int d = 0; d < D; d++) {
                        m_normal[d] = normal[e + Offset, n, d];
                        m_x[d] = x[e + Offset, n, d];
                    }
                    for (int l = 0; l < L; l++) {
                        m_Uin[l] = Uin[l][e + Offset, n];
                        m_Uout[l] = Uout[l][e + Offset, n];
                    }

                    double retIn;
                    double retOut;
                    InnerEdgeFlux(time, m_x, m_normal, m_Uin, m_Uout, e + jEdge, out retIn, out retOut);

                    Output_InCell[e + Offset, n] += retIn;
                    Output_OutCell[e + Offset, n] += retOut;
                }
            }
        }

        /// <summary>
        /// override this method to implement the Riemann flux at border edges
        /// </summary>
        protected abstract double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge);

        /// <summary>
        /// override this method to implement the Riemann flux a interior edges
        /// </summary>
        protected abstract void InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge, out double FluxIn, out double FluxOut);


        double[] m_Uin;

        double[] m_Uout;

        double[] m_normal;

        double[] m_x;

        double[] m_Flux;

        void AllocateHelperArrays(int D) {
            if (m_Uin == null) {
                int l = this.ArgumentOrdering.Count;
                if (this.ParameterOrdering != null)
                    l += this.ParameterOrdering.Count;

                m_Uin = new double[l];
                m_Uout = new double[l];

                m_x = new double[D];
                m_normal = new double[D];
                m_Flux = new double[D];
            }
        }


        #endregion

        #region IEquationComponent Members

        /// <summary>
        /// see <see cref="IEquationComponent.ArgumentOrdering"/>;
        /// </summary>
        abstract public IList<string> ArgumentOrdering { get; }

        /// <summary>
        /// see <see cref="IEquationComponent.ParameterOrdering"/>;
        /// </summary>
        virtual public IList<string> ParameterOrdering { get { return null; } }

        #endregion
    }
}
