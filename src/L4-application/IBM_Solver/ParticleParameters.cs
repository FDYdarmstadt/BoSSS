using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IBM_Solver {
    public class ParticleParameters {
        public ParticleParameters(Vector pVelocity, Vector aVelocity) {
            m_pointVelocity = pVelocity;
            m_AngularVelocity = aVelocity;
        }

        private Vector m_pointVelocity;
        private Vector m_AngularVelocity;

        public Vector PointVelocity {
            get { return m_pointVelocity; }
        }
        public Vector AngularVelocity {
            get { return m_AngularVelocity; }
        }


    }
}
