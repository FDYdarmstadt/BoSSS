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

using FSI_Solver;
using ilPSP;
using System;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    [Serializable]
    public class InitializeMotion {

        /// <summary>
        /// The initialization of the particle motion. Depending on the given parameters the correct model is chosen.
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="particleDensity">
        /// The density of the particle.
        /// </param>
        /// <param name="isDry">
        /// Set true if no hydrodynamics should be included.
        /// </param>
        /// <param name="noRotation"></param>
        /// <param name="noTranslation"></param>
        /// <param name="underrelaxationParam">
        /// The underrelaxation parameters (convergence limit, prefactor and a bool whether to use addaptive underrelaxation) defined in <see cref="ParticleUnderrelaxationParam"/>.
        /// </param>
        /// <param name="addedDampingCoefficient">
        /// The added damping coefficient is a scaling factor for the model. If the value is smaller than zero no added damping is applied. Otherwise it should be between 0.5 and 1.5, for reference: Banks et.al. 2017.
        /// </param>
        public InitializeMotion(Vector gravity, double particleDensity = 0, bool isDry = false, bool noRotation = false, bool noTranslation = false, double addedDampingCoefficient = 0, bool calculateOnlyVelocity = false) {
            if (gravity.IsNullOrEmpty())
                gravity = new Vector(0, 0);
            this.gravity = new Vector(gravity);
            this.particleDensity = particleDensity == 0 ? 1 : particleDensity;
            this.isDry = isDry;
            this.noRotation = noRotation;
            this.noTranslation = noTranslation;
            this.addedDampingCoefficient = addedDampingCoefficient;
            this.calculateOnlyVelocity = calculateOnlyVelocity;
        }

        [DataMember]
        private readonly FSI_Auxillary Aux = new FSI_Auxillary();
        [DataMember]
        private readonly Vector gravity = new Vector(0,0);
        [DataMember]
        private readonly double particleDensity;
        [DataMember]
        private readonly bool isDry;
        [DataMember]
        private readonly bool noRotation;
        [DataMember]
        private readonly bool noTranslation;
        [DataMember]
        private readonly bool calculateOnlyVelocity;
        [DataMember]
        private readonly double addedDampingCoefficient;

        /// <summary>
        /// Initial check for the particle motion parameter.
        /// </summary>
        public void CheckInput() {
            if (addedDampingCoefficient != 0 && addedDampingCoefficient < 0.5 && addedDampingCoefficient > 1.5)
                throw new Exception("Error in control file: Added damping coefficient should be between 0.5 and 1.5! See for reference Banks et al.");
            if (addedDampingCoefficient != 0 && (noRotation || noTranslation))
                throw new Exception("Error in control file: The added damping model is designed to contain all possible motion types (translation and rotation).");

            Aux.TestArithmeticException(gravity, "gravity");
            Aux.TestArithmeticException(addedDampingCoefficient, "added damping coefficient");
        }

        /// <summary>
        /// Initialize the correct derived motion.cs.
        /// </summary>
        public Motion ParticleMotion {
            get {
                if (noRotation && noTranslation) {
                    if (calculateOnlyVelocity)
                        return new MotionFixedWithVelocity(gravity, particleDensity);
                    else
                        return new MotionFixed(gravity);
                }
                if (isDry) {
                    return noRotation ? new MotionDryNoRotation(gravity, particleDensity)
                        : noTranslation ? new MotionDryNoTranslation(gravity, particleDensity)
                        : new MotionDry(gravity, particleDensity);
                }
                if (addedDampingCoefficient != 0)
                    return new MotionAddedDamping(gravity, particleDensity, addedDampingCoefficient);
                else
                    return noRotation ? new MotionWetNoRotation(gravity, particleDensity)
                        : noTranslation ? new MotionWetNoTranslation(gravity, particleDensity)
                        : new Motion(gravity, particleDensity);
            }
        }
    }
}
