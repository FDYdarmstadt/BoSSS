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
using BoSSS.Foundation;
using ilPSP;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// An abstract base class that helps with the implementation of 
    /// <see cref="INonlinearSource"/>;
    /// Instead of implementing the (a bit complicated) vectorized functions in the interface,
    /// the user only has to override scalar versions;
    /// </summary>
    public abstract class NonlinearSource : INonlinearSource {

        /// <summary>
        /// not in use, returning null
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        /// <summary>
        /// override this method to implement the nonlinear source
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="j">cell index</param>
        /// <param name="U"></param>
        /// <returns></returns>
        abstract protected double Source(double time, int j, double[] x, double[] U);


        double[] m_U;

        double[] m_x;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="D">spatial dimension</param>
        void AllocateHelperArrays(int D) {
            if (m_U == null) {
                int l = this.ArgumentOrdering.Count;
                if (this.ParameterOrdering != null)
                    l += this.ParameterOrdering.Count;
                m_U = new double[l];
                m_x = new double[D];
            }
        }


        #region INonlinearSource Member

        /// <summary>
        /// vectorized version defined by <see cref="INonlinearSource"/>;
        /// makes multiple calls to the serial version
        /// <see cref="Source(double,int,double[],double[])"/>;
        /// </summary>
        public void Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int j0, int Lenght, MultidimensionalArray Output) {

            // preparation
            // -----------
            int D = x.GetLength(2);
            AllocateHelperArrays(D);
            int NoOfNodes = x.GetLength(1);
            int L = m_U.Length;

            // loop
            // ----

            for (int e = 0; e < Lenght; e++) {
                for (int n = 0; n < NoOfNodes; n++) {

                    for (int l = 0; l < L; l++) {
                        m_U[l] = U[l][e + IndexOffset, n];
                    }

                    for (int d = 0; d < D; d++) {
                        m_x[d] = x[e + IndexOffset, n, d];
                    }

                    Output[e + IndexOffset, n] += Source(time, j0 + e, m_x, m_U);
                }
            }

        }

        #endregion

        #region IEquationComponent Member

        /// <summary>
        /// <see cref="IEquationComponent.ArgumentOrdering"/>
        /// </summary>
        public abstract IList<string> ArgumentOrdering {
            get;
        }

        #endregion
    }
}


