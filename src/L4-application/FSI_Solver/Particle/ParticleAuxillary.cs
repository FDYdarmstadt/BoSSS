/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    class ParticleAuxillary {

        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for onedimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        internal void SaveValueOfLastTimestep(List<double> variable) {
            variable.Insert(0, new double());
            variable[0] = 0;
            variable.RemoveAt(variable.Count - 1);
        }

        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for multidimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        internal void SaveMultidimValueOfLastTimestep(List<double[]> variable) {
            int Dim = variable[0].Length;
            variable.Insert(0, new double[Dim]);
            for (int d = 0; d < Dim; d++)
            {
                variable[0][d] = 0;
            }
            variable.RemoveAt(variable.Count - 1);
        }
    }
}
