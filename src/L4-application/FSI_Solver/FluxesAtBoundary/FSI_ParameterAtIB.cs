/* =======================================================================
Copyright 2020 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using BoSSS.Application.FSI_Solver;
using ilPSP;

namespace FSI_Solver {
    public class FSI_ParameterAtIB {
        /// <summary>
        /// Coupling parameters at immersed boundary
        /// </summary>
        /// <param name="currentParticle">The particle which contains the current point</param>
        /// <param name="currentPoint">The current point</param>
        public FSI_ParameterAtIB(Particle currentParticle, Vector currentPoint) {
            m_CurrentParticle = currentParticle;
            m_CurrentParticle.CalculateRadialVector(currentPoint, out m_RadialVector, out m_radialLength);
        }

        private readonly Particle m_CurrentParticle;
        private Vector m_RadialVector;
        private readonly double m_radialLength;

        /// <summary>
        /// Veloctiy of the point X
        /// </summary>
        public Vector VelocityAtPointOnLevelSet() {
            return new Vector(m_CurrentParticle.Motion.GetTranslationalVelocity(0)[0] - m_radialLength * m_CurrentParticle.Motion.GetRotationalVelocity(0) * m_RadialVector[1],
                              m_CurrentParticle.Motion.GetTranslationalVelocity(0)[1] + m_radialLength * m_CurrentParticle.Motion.GetRotationalVelocity(0) * m_RadialVector[0]);
        }

        /// <summary>
        /// Current particle angle
        /// </summary>
        public double Angle() {
            return m_CurrentParticle.Motion.GetAngle(0);
        }

        /// <summary>
        /// Active tangential stress on the surface of the particle.
        /// </summary>
        public double ActiveStress() {
            return m_CurrentParticle.ActiveStress;
        }
    }
}
