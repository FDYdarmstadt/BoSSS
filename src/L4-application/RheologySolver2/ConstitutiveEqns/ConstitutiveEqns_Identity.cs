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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Volume integral of identity part of constitutive equations.
    /// </summary>
    public class ConstitutiveEqns_Identity : IVolumeForm, IEquationComponent {

        private int component; // equation index (0: xx, 1: xy, 2: yy)

        /// <summary>
        /// Initialize identity
        /// </summary>
        public ConstitutiveEqns_Identity(int component) {
            this.component = component;
        }

        /// <summary>
        /// Choosing the required terms (These Flags control, whether certain terms are evaluated during quadrature of the forms)
        /// </summary>
        public TermActivationFlags VolTerms
        {
            get { return TermActivationFlags.V | TermActivationFlags.UxV; }
        }

        /// <summary>
        /// Ordering of the dependencies as whole string
        /// </summary>
        static string[] allArg = new string[] { VariableNames.StressXX, VariableNames.StressXY, VariableNames.StressYY };

        /// <summary>
        /// Ordering of the dependencies
        /// </summary>
        public IList<string> ArgumentOrdering
        {
            get { return new string[] { allArg[component] }; }
        }

        /// <summary>
        /// Ordering of the parameters - null at identity part
        /// </summary>
        public IList<string> ParameterOrdering {get;}
 
        /// <summary>
        /// Calculating the integral of the volume part
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] T, double[,] Grad_T, double V, double[] GradV) {

            double res = 0.0;
            if (component == 1)
            {
                res = T[0];
            }
            else {
                res = T[0];
            }

            return res * V;
        }

    }
}
