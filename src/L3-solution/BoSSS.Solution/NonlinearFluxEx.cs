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

using BoSSS.Foundation;
using System.Collections.Generic;
using BoSSS.Platform;
using System.Collections;
using ilPSP;

namespace BoSSS.Solution.Utils {
    
    /// <summary>
    /// An abstract baseclass that helps with the implementation of 
    /// <see cref="INonlinearFluxEx"/>;
    /// Instead of implementing the (a bit complicated) vectorized functions in the interface,
    /// the user only has to override scalar versions;
    /// </summary>
    public abstract class NonlinearFluxEx : INonlinearFluxEx {

        /// <summary>
        /// not in use, returning null
        /// </summary>
        public virtual IList<string> ParameterOrdering { get { return null; } }

        /// <summary>
        /// override this method to implement the Riemann flux at border edges
        /// </summary>
        protected abstract double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, double[] UinMean, int jEdge);

        /// <summary>
        /// override this method to implement the Riemann flux a interior edges
        /// </summary>
        protected abstract double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] UinMean, double[] Uout, double[] UoutMean, int jEdge);

        /// <summary>
        /// override this method to implement the flux function.
        /// </summary>
        protected abstract void Flux(double time, double[] x, double[] U, double[] output, int jCell);

        double[] m_Uin;

        double[] m_Uout;

        double[] m_UinMean;

        double[] m_UoutMean;

        double[] m_normal;

        double[] m_x;

        double[] m_Flux;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="D">spatial dimension</param>
        void AllocateHelperArrays(int D) {
            if (m_Uin == null) {
                int l = this.ArgumentOrdering.Count;
                if (this.ParameterOrdering != null)
                    l += this.ParameterOrdering.Count;

                m_Uin = new double[l];
                m_Uout = new double[l];
                m_UinMean = new double[l];
                m_UoutMean = new double[l];

                m_x = new double[D];
                m_normal = new double[D];
                m_Flux = new double[D];
            }
        }

        #region INonlinearFlux Member

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// </summary>
        public void BorderEdgeFlux(double time, int jEdg,  MultidimensionalArray x,
                                   MultidimensionalArray normal, bool normalFlipped, 
                                   byte[] EdgeTags, int EdgeTagsOffset,
                                   MultidimensionalArray[] Uin, MultidimensionalArray[] UinMean,
                                   int IndexOffset, int Lenght, 
                                   MultidimensionalArray Output) {

            // preparation
            // -----------
            int D = normal.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = Uin[0].GetLength(1);
            int L = m_Uin.Length;
            double sign = normalFlipped ? -1.0 : 1.0;

            // loop
            // ----

            for (int e = 0; e < Lenght; e++) {
                for (int l = 0; l < L; l++) m_UinMean[l] = UinMean[l][e + IndexOffset];
                
                byte EdgeTag = EdgeTags[e + EdgeTagsOffset];

                for (int n = 0; n < NoOfNodes; n++) {
                    for (int d = 0; d < D; d++) m_x[d] = x[e + IndexOffset, n, d]; 
                    for (int l = 0; l < L; l++) m_Uin[l] = Uin[l][e + IndexOffset, n];
                    for (int d = 0; d < D; d++) m_normal[d] = normal[e + IndexOffset, n, d] * sign;

                    Output[e + IndexOffset, n] += BorderEdgeFlux(time, m_x, m_normal, EdgeTag, m_Uin, m_UinMean, jEdg + e);
                }
            }
        }

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// </summary>
        public void InnerEdgeFlux(double time, int jEdg,  MultidimensionalArray x,
                                  MultidimensionalArray normal,
                                  MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, 
                                  MultidimensionalArray[] UinMean, MultidimensionalArray[] UoutMean, 
                                  int Offset, int Lenght, 
                                  MultidimensionalArray Output) {

            // preparation
            // -----------
            int D = normal.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = Uin[0].GetLength(1);
            int L = m_Uin.Length;

            // loop
            // ----

            for (int e = 0; e < Lenght; e++) {
                for (int l = 0; l < L; l++) {
                    m_UinMean[l] = UinMean[l][e + Offset];
                    m_UoutMean[l] = UoutMean[l][e + Offset];
                }
                for (int n = 0; n < NoOfNodes; n++) {

                    for (int d = 0; d < D; d++) {
                        m_x[d] = x[e + Offset, n, d];
                        m_normal[d] = normal[e + Offset, n, d];
                    }

                    for (int l = 0; l < L; l++) {
                        m_Uin[l] = Uin[l][e + Offset, n];
                        m_Uout[l] = Uout[l][e + Offset, n];
                    }

                    Output[e + Offset, n] += InnerEdgeFlux(time, m_x, m_normal, m_Uin, m_UinMean, m_Uout, m_UoutMean, jEdg + e);
                }
            }
        }

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// <see cref="Flux(double,double[],double[],double[],int)"/>;
        /// </summary>
        public void Flux(double time, MultidimensionalArray x, MultidimensionalArray[] U, int Offset, int Length, MultidimensionalArray Output, int jCell) {
            // preparation
            // -----------
            int D;
            D = x.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = U[0].GetLength(1);
            int L = m_Uin.Length;

            // loop
            // ----

            for (int j = 0; j < Length; j++) {
                for (int n = 0; n < NoOfNodes; n++) {

                    for (int l = 0; l < L; l++) {
                        m_Uin[l] = U[l][j + Offset, n];
                    }
                    for (int d = 0; d < D; d++) {
                        m_x[d] = x[j + Offset, n, d];
                    }


                    Flux(time, m_x, m_Uin, m_Flux, j + jCell);

                    for (int d = 0; d < D; d++) Output[j + Offset, n, d] += m_Flux[d];
                }
            }

        }

        #endregion

        #region IEquationComponent Member

        /// <summary>
        /// <see cref="IEquationComponent.ArgumentOrdering"/>
        /// </summary>
        public abstract IList<string> ArgumentOrdering { get; }

        #endregion

    }
}
