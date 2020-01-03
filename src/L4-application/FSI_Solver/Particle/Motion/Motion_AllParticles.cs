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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    internal class Motion_AllParticles {
        internal Motion_AllParticles(LevelSetTracker lsTrk) {
            m_LsTrk = lsTrk;
        }
        [DataMember]
        private static readonly int m_Dim = 2;
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorquePreviousIteration = new List<double[]>();
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorqueWithoutRelaxation = new List<double[]>();
        [DataMember]
        private readonly LevelSetTracker m_LsTrk;

        /// <summary>
        /// ...
        /// </summary>
        /// <param name="AllParticles"></param>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        /// <param name="underrelax"></param>
        internal void CalculateHydrodynamics(List<Particle> AllParticles, ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity, bool underrelax) {
            double[] hydrodynamics = new double[m_Dim * AllParticles.Count() + AllParticles.Count()];
            for (int p = 0; p < AllParticles.Count(); p++) {
                Particle currentParticle = AllParticles[p];
                CellMask cutCells = currentParticle.CutCells_P(m_LsTrk);
                int offset = p * (m_Dim + 1);
                double[] tempForces = currentParticle.Motion.CalculateHydrodynamicForces(hydrodynamicsIntegration, fluidDensity, cutCells);
                double tempTorque = currentParticle.Motion.CalculateHydrodynamicTorque(hydrodynamicsIntegration, cutCells);
                for (int d = 0; d < m_Dim; d++) {
                    hydrodynamics[offset + d] = tempForces[d];
                }
                hydrodynamics[offset + m_Dim] = tempTorque;
            }
            double[] relaxatedHydrodynamics = hydrodynamics.CloneAs();
            double omega = AllParticles[0].Motion.omega;
            if (underrelax)
                relaxatedHydrodynamics = HydrodynamicsPostprocessing(hydrodynamics, ref omega);
            AllParticles[0].Motion.omega = omega;
            for (int p = 0; p < AllParticles.Count(); p++) {
                Particle currentParticle = AllParticles[p];
                currentParticle.Motion.UpdateForcesAndTorque(p, relaxatedHydrodynamics);
            }
        }

        /// <summary>
        /// Post-processing of the hydrodynamics. If desired the underrelaxation is applied to the forces and torque.
        /// </summary>
        /// <param name="hydrodynamics"></param>
        private double[] HydrodynamicsPostprocessing(double[] hydrodynamics, ref double omega) {
            m_ForcesAndTorqueWithoutRelaxation.Insert(0, hydrodynamics.CloneAs());
            ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation(hydrodynamics, m_ForcesAndTorquePreviousIteration);
            double[] relaxatedHydrodynamics;
            if (m_ForcesAndTorquePreviousIteration.Count >= 3) {
                relaxatedHydrodynamics = Underrelaxation.AitkenUnderrelaxation(ref omega, m_ForcesAndTorqueWithoutRelaxation);
            }
            else {
                relaxatedHydrodynamics = Underrelaxation.StaticUnderrelaxation();
            }
            return relaxatedHydrodynamics;
        }

        /// <summary>
        /// ...
        /// </summary>
        /// <param name="AllParticles"></param>
        internal void SaveHydrodynamicOfPreviousIteration(List<Particle> AllParticles) {
            double[] hydrodynamics = new double[(m_Dim + 1) * AllParticles.Count()];
            for (int p = 0; p < AllParticles.Count(); p++) {
                Particle currentParticle = AllParticles[p];
                int offset = p * (m_Dim + 1);
                double[] tempForces = currentParticle.Motion.GetHydrodynamicForces(0);
                double tempTorque = currentParticle.Motion.GetHydrodynamicTorque(0);
                for (int d = 0; d < m_Dim; d++) {
                    hydrodynamics[offset + d] = tempForces[d];
                }
                hydrodynamics[offset + m_Dim] = tempTorque;
            }
            m_ForcesAndTorquePreviousIteration.Insert(0, hydrodynamics.CloneAs());
        }

        /// <summary>
        /// Residual for fully coupled system
        /// </summary>
        /// <param name="iterationCounter"></param>
        internal double CalculateParticleResidual(ref int iterationCounter) {
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            double residual = 0;
            double denom = 0;
            if (iterationCounter <= 2)
                residual = double.MaxValue;
            else {
                for (int i = 0; i < m_ForcesAndTorquePreviousIteration[0].Length; i++) {
                    residual += (m_ForcesAndTorquePreviousIteration[0][i] - m_ForcesAndTorquePreviousIteration[1][i]).Pow2();
                    denom += m_ForcesAndTorquePreviousIteration[0][i].Pow2();
                }
            }
            residual = Math.Sqrt(residual / denom);
            iterationCounter += 1;
            return residual;
        }
    }
}
