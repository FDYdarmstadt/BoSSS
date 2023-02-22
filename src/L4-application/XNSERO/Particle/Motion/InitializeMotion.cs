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
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    public class InitializeMotion {

        /// <summary>
        /// The initialization of the particle motion. Depending on the given parameters the correct model is chosen.
        /// </summary>
        /// <param name="particleDensity">
        /// The density of the particle.
        /// </param>
        /// <param name="isDry">
        /// Set true if no hydrodynamics should be included.
        /// </param>
        /// <param name="noRotation"></param>
        /// <param name="noTranslation"></param>
        public InitializeMotion(double particleDensity = 1, bool isDry = false, bool noRotation = false, bool noTranslation = false, bool calculateOnlyVelocity = false) {
            this.particleDensity = particleDensity;
            this.isDry = isDry;
            this.noRotation = noRotation;
            this.noTranslation = noTranslation;
            this.calculateOnlyVelocity = calculateOnlyVelocity;

            if (particleDensity == 0)
                throw new Exception("Zero particle density is not allowed.");
        }

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

        /// <summary>
        /// Initialize the correct derived motion.cs.
        /// </summary>
        public Motion ParticleMotion {
            get {
                if (noRotation && noTranslation) {
                    if (calculateOnlyVelocity)
                        return new MotionFixedWithVelocity(particleDensity);
                    else
                        return new MotionFixed();
                }
                else if (isDry) {
                    return noRotation ? new MotionDryNoRotation(particleDensity)
                        : noTranslation ? new MotionDryNoTranslation(particleDensity)
                        : new MotionDry(particleDensity);
                }
                else
                    return noRotation ? new MotionWetNoRotation(particleDensity)
                        : noTranslation ? new MotionWetNoTranslation(particleDensity)
                        : new Motion(particleDensity);
            }
        }
    }
}
