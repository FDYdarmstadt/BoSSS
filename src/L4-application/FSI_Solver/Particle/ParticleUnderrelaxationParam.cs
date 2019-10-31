﻿/* =======================================================================
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

namespace BoSSS.Application.FSI_Solver {
    public class ParticleUnderrelaxationParam {
        /// <summary>
        /// Collection of parameters for hydrodynamic underrelaxation.
        /// </summary>
        /// <param name="convergenceLimit">
        /// The convergence limit for the hydrodynamic forces.
        /// </param>
        /// <param name="relaxationFactor">
        /// The underrelaxation factor.
        /// </param>
        /// <param name="useAddaptiveUnderrelaxation">
        /// Set true if you want addaptive underrelaxation.
        /// </param>
        public ParticleUnderrelaxationParam(double convergenceLimit, UnderrelaxationMethod underrelaxationMethod = UnderrelaxationMethod.ProcentualRelaxation, double relaxationFactor = 0.1, bool useAddaptiveUnderrelaxation = true) {
            m_ConvergenceLimit = convergenceLimit;
            m_UnderrelaxationFactor = relaxationFactor;
            m_UseAdaptiveUnderrelaxation = useAddaptiveUnderrelaxation;
            m_Method = underrelaxationMethod;
        }

        public enum UnderrelaxationMethod {
            ProcentualRelaxation = 0,

            AitkenRelaxation = 1,

            //MinimalPolynomial = 2

        }

        public UnderrelaxationMethod m_Method = UnderrelaxationMethod.ProcentualRelaxation;

        /// <summary>
        /// The convergence limit for the hydrodynamic forces.
        /// </summary>
        public double m_ConvergenceLimit { get; }

        /// <summary>
        /// The underrelaxation factor.
        /// </summary>
        public double m_UnderrelaxationFactor { get; }

        /// <summary>
        /// True if the underrelaxation procedure is addaptive.
        /// </summary>
        public bool m_UseAdaptiveUnderrelaxation { get; }
    }
}