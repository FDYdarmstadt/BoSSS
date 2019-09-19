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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {

    public class ParticleMotionInit {
        public ParticleMotionInit(double[] gravity = null,
            bool isDry = false,
            bool noRotation = false,
            bool noTranslation = false,
            ParticleUnderrelaxationParam underrelaxationParam = null,
            double addedDampingCoefficient = -1) {
            m_Gravity = gravity.IsNullOrEmpty() ? (new double[] { 0, 9.81 }) : gravity;
            m_IsDry = isDry;
            m_NoRotation = noRotation;
            m_NoTranslation = noTranslation;
            m_UnderrelaxationParam = underrelaxationParam;
            m_AddedDampingCoefficient = addedDampingCoefficient;
        }

        private readonly FSI_Auxillary Aux = new FSI_Auxillary();
        private readonly double[] m_Gravity;
        private readonly bool m_IsDry;
        private readonly bool m_NoRotation;
        private readonly bool m_NoTranslation;
        private readonly ParticleUnderrelaxationParam m_UnderrelaxationParam;
        private readonly double m_AddedDampingCoefficient;

        public void CheckInput() {
            if (m_IsDry && m_UnderrelaxationParam != null)
                throw new Exception("Error in control file: Cannot perform a dry simulation with full coupling between the particles and the (non-existing) fluid");

            if (m_AddedDampingCoefficient != -1 && m_AddedDampingCoefficient < 0.5 && m_AddedDampingCoefficient > 1.5)
                throw new Exception("Error in control file: Added damping coefficient should be between 0.5 and 1.5! See for reference Banks et al.");
            if (m_AddedDampingCoefficient != -1 && (m_NoRotation || m_NoTranslation))
                throw new Exception("Error in control file: The added damping model is designed to contain all possible motion types (translation and rotation).");

            Aux.TestArithmeticException(m_Gravity, "gravity");
            Aux.TestArithmeticException(m_AddedDampingCoefficient, "added damping coefficient");
        }

        public Motion_Wet GetParticleMotion() {
            if (m_NoRotation && m_NoTranslation)
                return new Motion_Fixed();
            if (m_IsDry) {
                return m_NoRotation ? new Motion_Dry_NoRotation(m_Gravity)
                    : m_NoTranslation ? new Motion_Dry_NoTranslation(m_Gravity)
                    : new Motion_Dry(m_Gravity);
            }
            if (m_AddedDampingCoefficient != -1)
                return new Motion_AddedDamping(m_Gravity, m_UnderrelaxationParam, m_AddedDampingCoefficient);
            else
                return m_NoRotation ? new Motion_Wet_NoRotation(m_Gravity)
                    : m_NoTranslation ? new Motion_Wet_NoTranslation(m_Gravity)
                    : new Motion_Wet(m_Gravity, m_UnderrelaxationParam);
        }
    }
}
