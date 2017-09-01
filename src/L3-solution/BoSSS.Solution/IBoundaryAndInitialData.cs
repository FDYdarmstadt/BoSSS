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
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Control {

    /// <summary>
    /// Common interface for boundary or initial data values.
    /// </summary>
    public interface IBoundaryAndInitialData {

        /// <summary>
        /// Returns the value of a scalar, time-dependent field (e.g. pressure or some velocity component)
        /// at a specific position and time.
        /// </summary>
        /// <param name="X">global/physical coordinates.</param>
        /// <param name="t">Physical time.</param>
        /// <returns>Function value.</returns>
        double Evaluate(double[] X, double t);
    }

}
