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
    /// <see cref="INonlinearFlux"/>;
    /// Instead of implementing the (a bit complicated) vectorized functions in the interface,
    /// the user only has to override scalar versions;
    /// </summary>
    public abstract class NonlinearFlux : INonlinearFlux {

        /// <summary>
        /// not in use, returning null
        /// </summary>
        public virtual IList<string> ParameterOrdering { get { return null; } }

        /// <summary>
        /// override this method to implement the Riemann flux at border edges
        /// </summary>
        protected abstract double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge);

        /// <summary>
        /// override this method to implement the Riemann flux a interior edges
        /// </summary>
        protected abstract double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge);

        /// <summary>
        /// override this method to implement the flux function.
        /// </summary>
        protected abstract void Flux(double time, double[] x, double[] U, double[] output);

        double[] m_Uin;

        double[] m_Uout;

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
        public void BorderEdgeFlux(double time, int jEdge, MultidimensionalArray x, MultidimensionalArray normal, bool normalFlipped,  byte[] EdgeTags, int EdgeTagsOffset, MultidimensionalArray[] Uin, int Offset, int Lenght, MultidimensionalArray Output) {

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
                
                byte EdgeTag = EdgeTags[e + EdgeTagsOffset];

                for( int n = 0; n < NoOfNodes; n++) {

                    for (int d = 0; d < D; d++) {
                        m_x[d] = x[e + Offset, n, d];
                        m_normal[d] = normal[e + Offset, n, d] * sign;
                    }
                    for (int l = 0; l < L; l++) m_Uin[l] = Uin[l][e+Offset, n];

                    double fluxVal = BorderEdgeFlux(time, m_x, m_normal, EdgeTag, m_Uin, jEdge+e)*sign;
#if DEBUG
                    if(double.IsNaN(fluxVal) || double.IsInfinity(fluxVal))
                        throw new ArithmeticException("some flux is NAN or INF.");
#endif
                    Output[e + Offset, n] += fluxVal;
                }
            }
        }

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// </summary>
        public void InnerEdgeFlux(double time, int jEdge, MultidimensionalArray x, MultidimensionalArray normal, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, int Offset, int Lenght, MultidimensionalArray Output) {

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

                    double fluxVal = InnerEdgeFlux(time, m_x, m_normal, m_Uin, m_Uout, jEdge + e);
#if DEBUG
                    if(double.IsNaN(fluxVal) || double.IsInfinity(fluxVal))
                        throw new ArithmeticException("some flux is NAN or INF.");
#endif
                    Output[e + Offset, n] += fluxVal;
                }
            }
        }

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearFlux"/>;
        /// makes multiple calls to the serial version
        /// <see cref="Flux(double,double[],double[],double[])"/>;
        /// </summary>
        public void Flux(double time, MultidimensionalArray x, MultidimensionalArray[] U, int Offset, int Length, MultidimensionalArray Output) {
            
            // preparation
            // -----------
            int D;
            D = x.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = Output.GetLength(1);
            int L = m_Uin.Length;

            // loop
            // ----

            for (int e = 0; e < Length; e++) {
                for (int n = 0; n < NoOfNodes; n++) {

                    for (int l = 0; l < L; l++) {
                        m_Uin[l] = U[l][e+Offset, n];
                    }
                    for (int d = 0; d < D; d++) {
                        m_x[d] = x[e+Offset, n, d];
                    }
                    
                    Flux(time, m_x, m_Uin, m_Flux);

                    for(int d = 0; d < D; d++) {                       
#if DEBUG
                        if(double.IsNaN(m_Flux[d]) || double.IsInfinity(m_Flux[d]))
                            throw new ArithmeticException("some flux is NAN or INF.");
#endif
                        Output[e + Offset, n, d] += m_Flux[d];
                    }
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
