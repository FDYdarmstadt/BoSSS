/* =======================================================================
Copyright 2025 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using System.Text;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// encapsulates the settings and computation of the additional Dong boundary energy term 
    /// </summary>
    public class DongBoundaryConditionTerm {

        double U0 = 1.0;

        double Delta = 1.0 / 20.0;

        /// <summary>
        /// Constructor for default settings
        /// </summary>
        public DongBoundaryConditionTerm() {
        }

        /// <summary>
        ///  custom constructor
        /// </summary>
        /// <param name="u0"></param>
        /// <param name="delta"></param>
        public DongBoundaryConditionTerm(double u0, double delta) {
            if (u0 < 0.0 || delta < 0.0) {
                throw new ArgumentException("Dong boundary condition term expecting pos values!");
            }
            U0 = u0;
            Delta = delta;
        }

        public double GetBoundaryTerm(double ndotu, double uAbs2) {
            double Sout = (U0 == 0.0 || Delta == 0.0) ? 0.0 : 0.5 * (1.0 - Math.Tanh(ndotu / (U0 * Delta)));
            return 0.5 * uAbs2 * Sout;
        }
    }
}
