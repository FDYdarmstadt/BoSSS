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
        public ParticleUnderrelaxationParam(double convergenceLimit, double relaxationFactor, bool useAddaptiveUnderrelaxation) {
            ConvergenceLimit = convergenceLimit;
            UnderrelaxationFactor = relaxationFactor;
            UsaAddaptiveUnderrelaxation = useAddaptiveUnderrelaxation;
        }

        /// <summary>
        /// The convergence limit for the hydrodynamic forces.
        /// </summary>
        public double ConvergenceLimit { get; }

        /// <summary>
        /// The underrelaxation factor.
        /// </summary>
        public double UnderrelaxationFactor { get; }

        /// <summary>
        /// True if the underrelaxation procedure is addaptive.
        /// </summary>
        public bool UsaAddaptiveUnderrelaxation { get; }
    }
}
