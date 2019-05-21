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
using BoSSS.Solution.CompressibleFlowCommon;
using ilPSP;

namespace CNS.Source {

    /// <summary>
    /// Source term of the momentum equation due to the volume force
    /// </summary>
    public class GravityMomentumSource : INonlinearSource {

        /// <summary>
        /// Options
        /// </summary>
        private CNSControl control;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="control"></param>
        public GravityMomentumSource(CNSControl control) {
            this.control = control;
        }

        #region INonlinearSource Members

        /// <summary>
        /// Evaluates the source term according to
        /// \f$ \rho / \text{Fr}^2\f$ 
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
        /// Only requires density
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                // Select density
                return new string[] { CompressibleVariables.Density };
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
