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
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    /// <summary>
    /// Contains methods to call the calculation of the hydrodynamic forces and torques. Also provides relaxations methods for hydrodynamics.
    /// </summary>
    internal class ParticleHydrodynamics {
        /// <summary>
        /// Enables calculation of the hydrodynamic forces and torques.
        /// </summary>
        /// <param name="SpatialDimension"></param>
        internal ParticleHydrodynamics(int SpatialDimension) {
            this.SpatialDimension = SpatialDimension;
            TorqueVectorDimension = SpatialDimension switch {
                2 => 1,
                3 => throw new NotImplementedException("XNSERO is currently only capable of 2D."),
                _ => throw new ArithmeticException("Error in spatial dimensions, only dim == 2 and dim == 3 are allowed!"),
            };
        }

        [DataMember]
        private readonly int SpatialDimension;
        [DataMember]
        private readonly int TorqueVectorDimension;
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorquePreviousIteration = new List<double[]>();
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorqueWithoutRelaxation = new List<double[]>();

        /// <summary>
        /// Calculation of hydrodynamic forces F and torque T for all particles. 
        /// The relaxation (Aitken-relaxation) is applied to the entire array of all F+T, leading to better convergence.
        /// </summary>
        /// <param name="Particles"></param>
        /// <param name="HydrodynamicIntegrator"></param>
        /// <param name="FluidSpecies"></param>
        /// <param name="Gravity"></param>
        /// <param name="dt"></param>
        internal void CalculateHydrodynamics(Particle[] Particles, CellMask AllCutCells, ParticleHydrodynamicsIntegration HydrodynamicIntegrator, string[] FluidSpecies, Vector Gravity, double dt) {
            using (new FuncTrace()) {
                double[] hydrodynamics = new double[(SpatialDimension + TorqueVectorDimension) * Particles.Length];
               
                // Calculate hydrodynamics.
                for (int p = 0; p < Particles.Length; p++) {
                    CellMask CutCells = Particles[p].ParticleCutCells(Particles[0].LsTrk, AllCutCells);
                    int offset = p * (SpatialDimension + TorqueVectorDimension);
                    double[] tempHydrodynamics = Particles[p].Motion.CalculateHydrodynamics(HydrodynamicIntegrator, CutCells, FluidSpecies, dt);
                    for (int d = 0; d < SpatialDimension + 1; d++) 
                        hydrodynamics[offset + d] = tempHydrodynamics[d];
                }

                hydrodynamics = hydrodynamics.MPISum();

                // Gravity, has to be added after MPISum()!
                for (int p = 0; p < Particles.Length; p++) {
                    Particle currentParticle = Particles[p];
                    int offset = p * (SpatialDimension + TorqueVectorDimension);
                    Vector gravity = currentParticle.Motion.GravityForce(Gravity);
                    for(int d = 0; d < SpatialDimension; d++) {
                        hydrodynamics[offset + d] += gravity[d];
                    }
                }

                // Relaxation
                if (hydrodynamics.Sum() > 1e-10){
                    double omega = Particles[0].Motion.RelaxationParameter;
                    hydrodynamics = HydrodynamicsRelaxation(hydrodynamics, ref omega);
                    Particles[0].Motion.RelaxationParameter = omega; 
                }

                // Write to single particles.
                for (int p = 0; p < Particles.Length; p++) {
                    Particle currentParticle = Particles[p];
                    currentParticle.Motion.UpdateForcesAndTorque(p, hydrodynamics);
                }
            }
        }

        /// <summary>
        /// Saves hydrodynamic forces and torques of the last time-step.
        /// </summary>
        /// <param name="AllParticles"></param>
        internal void SaveHydrodynamicOfPreviousTimestep(Particle[] AllParticles) {
            double[] hydrodynamics = new double[(SpatialDimension + 1) * AllParticles.Length];
            m_ForcesAndTorquePreviousIteration.Clear();
            m_ForcesAndTorqueWithoutRelaxation.Clear();
            for(int p = 0; p < AllParticles.Length; p++) {
                Particle currentParticle = AllParticles[p];
                int offset = p * (SpatialDimension + 1);
                double[] tempForces = currentParticle.Motion.GetHydrodynamicForces(0);
                double tempTorque = currentParticle.Motion.GetHydrodynamicTorque(0);
                for(int d = 0; d < SpatialDimension; d++) {
                    hydrodynamics[offset + d] = tempForces[d];
                }
                hydrodynamics[offset + SpatialDimension] = tempTorque;
            }
            m_ForcesAndTorquePreviousIteration.Insert(0, hydrodynamics.CloneAs());
        }

        /// <summary>
        /// Saves hydrodynamic forces and torques of the last iteration of the solver in the current time-step.
        /// </summary>
        /// <param name="AllParticles"></param>
        internal void SaveHydrodynamicOfPreviousIteration(Particle[] AllParticles) {
            double[] hydrodynamics = new double[(SpatialDimension + 1) * AllParticles.Length];
            if (m_ForcesAndTorquePreviousIteration.Count() > 4) {
                m_ForcesAndTorquePreviousIteration.RemoveAt(4);
            }
            if(m_ForcesAndTorqueWithoutRelaxation.Count() > 4) {
                m_ForcesAndTorqueWithoutRelaxation.RemoveAt(4);
            }
            for (int p = 0; p < AllParticles.Length; p++) {
                Particle currentParticle = AllParticles[p];
                int offset = p * (SpatialDimension + 1);
                double[] tempForces = currentParticle.Motion.GetHydrodynamicForces(0);
                double tempTorque = currentParticle.Motion.GetHydrodynamicTorque(0);
                for (int d = 0; d < SpatialDimension; d++) {
                    hydrodynamics[offset + d] = tempForces[d];
                }
                hydrodynamics[offset + SpatialDimension] = tempTorque;
            }
            m_ForcesAndTorquePreviousIteration.Insert(0, hydrodynamics.CloneAs());
        }

        /// <summary>
        /// Post-processing of the hydrodynamics. If desired the underrelaxation is applied to the forces and torque. Küttler+Wall 2008: "Fixed-point fluid–structure interaction solvers with dynamic relaxation"
        /// </summary>
        /// <param name="hydrodynamics"></param>
        private double[] HydrodynamicsRelaxation(double[] hydrodynamics, ref double omega) {
            m_ForcesAndTorqueWithoutRelaxation.Insert(0, hydrodynamics.CloneAs());
            return m_ForcesAndTorqueWithoutRelaxation.Count > 2 ? AitkenRelaxation(hydrodynamics, ref omega) : StaticUnderrelaxation(hydrodynamics);
        }

        private double[] StaticUnderrelaxation(double[] variable) {
            using (new FuncTrace()) {
                double[] returnVariable = variable.CloneAs();
                for (int d = 0; d < variable.Length; d++) {
                    if (variable[d] == 0)
                        continue;
                    returnVariable[d] = 0.5 * variable[d] + (1 - 0.5) * m_ForcesAndTorquePreviousIteration[1][d];
                }
                return returnVariable;
            }
        }

        private double[] AitkenRelaxation(double[] variable, ref double Omega) {
            using (new FuncTrace()) {
                double[][] residual = new double[variable.Length][];
                double[] residualDiff = new double[variable.Length];
                double residualScalar = 0;
                for (int i = 0; i < variable.Length; i++) {
                    residual[i] = new double[] { (variable[i] - m_ForcesAndTorquePreviousIteration[0][i]), (m_ForcesAndTorqueWithoutRelaxation[1][i] - m_ForcesAndTorquePreviousIteration[1][i]) };
                    residualDiff[i] = residual[i][0] - residual[i][1];
                    residualScalar += residual[i][1] * residualDiff[i];
                }
                Omega = -Omega * residualScalar / residualDiff.L2Norm().Pow2();
                double[] outVar = variable.CloneAs();
                for (int i = 0; i < variable.Length; i++) {
                    outVar[i] = Omega * (variable[i] - m_ForcesAndTorquePreviousIteration[0][i]) + m_ForcesAndTorquePreviousIteration[0][i];
                }
                return outVar;
            }
        }
    }
}
