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

using System.Linq;
using System.Collections.Generic;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Solution.CompressibleFlowCommon;

namespace CNS.Source {

    /// <summary>
    /// Gravity source term for the energy equation
    /// </summary>
    public class GravityEnergySource : INonlinearSource {

        /// <summary>
        /// Options
        /// </summary>
        private CNSControl control;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="control"></param>
        public GravityEnergySource(CNSControl control) {
            this.control = control;
        }

        #region INonlinearSource Members

        /// <summary>
        /// Computes the source term according to
        /// \f$ m_d / \text{Fr}^2\f$ 
        /// where d=1 in 2D and d=2 in 3D (cf.
        /// <see cref="ArgumentOrdering"/>).
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="U"></param>
        /// <param name="IndexOffset"></param>
        /// <param name="FirstCellInd"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        public void Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int FirstCellInd, int Lenght, MultidimensionalArray Output) {
            int noOfNodes = x.GetLength(1);
            double FroudeSquared = control.FroudeNumber * control.FroudeNumber;

            for (int i = 0; i < Lenght; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    Output[i + IndexOffset, j] += U[0][i + IndexOffset, j] / FroudeSquared;
                }
            }
        }

        #endregion

        #region IEquationComponent Members

        /// <summary>
        /// Always returns the last component of the momentum vector.
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return new string[] { CompressibleVariables.Momentum.Last() };
            }
        }

        /// <summary>
        /// None
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        #endregion
    }
}
