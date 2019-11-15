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
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    internal class Motion_AllParticles {
        internal Motion_AllParticles(LevelSetTracker lsTrk) {
            m_LsTrk = lsTrk;
        }
        [NonSerialized]
        internal readonly FSI_Auxillary Aux = new FSI_Auxillary();
        [DataMember]
        protected static int m_Dim = 2;
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorquePreviousIteration = new List<double[]>();
        [DataMember]
        private List<double[]> m_ForcesAndTorqueWithoutRelaxation = new List<double[]>();
        [DataMember]
        private double relaxationCoeff = 1;
        LevelSetTracker m_LsTrk;

        internal void CalculateHydrodynamics(List<Particle> AllParticles, ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity, bool underrelax) {
            double[] hydrodynamics = new double[m_Dim * AllParticles.Count() + AllParticles.Count()];
            for (int p = 0; p < AllParticles.Count(); p++) {
                Particle currentParticle = AllParticles[p];
                CellMask cutCells = currentParticle.CutCells_P(m_LsTrk);
                int offset = p * (m_Dim + 1);
                double[] tempForces = currentParticle.Motion.CalculateHydrodynamicForces(hydrodynamicsIntegration, fluidDensity, cutCells);
                double tempTorque = currentParticle.Motion.CalculateHydrodynamicTorque(hydrodynamicsIntegration, cutCells);
                StabilizeHydrodynamics(ref tempForces, ref tempTorque, currentParticle.Motion.MaxParticleLengthScale);
                for (int d = 0; d < m_Dim; d++) {
                    hydrodynamics[offset + d] = tempForces[d];
                }
                hydrodynamics[offset + m_Dim] = tempTorque;
            }
            if (underrelax)
                HydrodynamicsPostprocessing(hydrodynamics);
            for (int p = 0; p < AllParticles.Count(); p++) {
                Particle currentParticle = AllParticles[p];
                currentParticle.Motion.UpdateForcesAndTorque(p, hydrodynamics);
            }
        }

        private void StabilizeHydrodynamics(ref double[] tempForces, ref double tempTorque, double lengthScale) {
            double averageForcesAndTorque = Math.Abs(CalculateAverageForces(tempForces, tempTorque, lengthScale));
            for (int d = 0; d < m_Dim; d++) {
                if (Math.Abs(tempForces[d] * 1e6) < averageForcesAndTorque || Math.Abs(tempForces[d]) < 1e-10)
                    tempForces[d] = 0;
            }
            if (Math.Abs(tempTorque * 1e6) < averageForcesAndTorque || Math.Abs(tempTorque) < 1e-10)
                tempTorque = 0;
        }

        /// <summary>
        /// Does what it says.
        /// </summary>
        /// <param name="forces">
        /// The hydrodynamic forces.
        /// </param>
        /// <param name="torque">
        /// The hydrodynamic torque.
        /// </param>
        /// <param name="lengthScale">
        /// The average Lengthscale of the particle.
        /// </param>
        private double CalculateAverageForces(double[] forces, double torque, double lengthScale) {
            double averageForces = Math.Abs(torque) / lengthScale;
            for (int d = 0; d < forces.Length; d++) {
                averageForces += forces[d];
            }
            return averageForces /= 3;
        }

        /// <summary>
        /// Post-processing of the hydrodynamics. If desired the underrelaxation is applied to the forces and torque.
        /// </summary>
        /// <param name="tempForces"></param>
        /// <param name="tempTorque"></param>
        /// <param name="firstIteration"></param>
        private void HydrodynamicsPostprocessing(double[] hydrodynamics) {
            m_ForcesAndTorqueWithoutRelaxation.Insert(0, hydrodynamics.CloneAs());
            ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
            if (m_ForcesAndTorquePreviousIteration.Count >= 4) {
                Underrelaxation.ForcesAndTorque(ref hydrodynamics, m_ForcesAndTorquePreviousIteration, ref relaxationCoeff, m_ForcesAndTorqueWithoutRelaxation);
            }
            else {
                Underrelaxation.ForceAndTorque(ref hydrodynamics, m_ForcesAndTorquePreviousIteration[1]);
            }
        }

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
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="iterationCounter"></param>
        /// <param name="MaximumIterations"></param>
        /// <param name="Residual"></param>
        /// <param name="IterationCounter_Out"> </param>
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
