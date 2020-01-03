using BoSSS.Application.FSI_Solver;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
