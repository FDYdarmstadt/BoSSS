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


namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Implementation of the time derivative as a linearized source term in the low-Mach combustion solver.
    /// Based on the implicit Euler scheme.
    /// </summary>
    public class TimeDerivativeLinearSource : BoSSS.Solution.Utils.LinearSource {
        string[] m_ArgumentOrdering;
        //string[] m_ParameterOrdering;
        MaterialLaw EoS;
        bool m_conti;
        double rho;
        double dt;

        /// <summary>
        /// Ctor.
        /// </summary> 
        /// <param name="EoS">The material law</param>
        /// <param name="conti">Set conti: true for the continuity equation</param>
        /// <param name="ArgumentOrdering"></param>
        /// <param name="TimeStepSize"></param>
        public TimeDerivativeLinearSource(MaterialLaw EoS, double TimeStepSize, String[] ArgumentOrdering, bool conti = false) {
            m_ArgumentOrdering = ArgumentOrdering;
            //m_ParameterOrdering = ParameterOrdering;
            this.EoS = EoS;
            dt = TimeStepSize;
            m_conti = conti;
        }


        /// <summary>
        /// The argument
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }

        /// <summary>
        /// Paramaters used to compute the density
        /// </summary>
        public override IList<string> ParameterOrdering {
            get { return EoS.ParameterOrdering; }
        }



        protected override double Source(double[] x, double[] parameters, double[] U) {
            rho = EoS.GetDensity(parameters);
            if (m_conti)
                return rho / dt;
            else
                return rho * U[0] / dt;
        }
    }
}
