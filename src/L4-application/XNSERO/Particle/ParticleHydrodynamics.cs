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
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    internal class ParticleHydrodynamics {
        internal ParticleHydrodynamics(int SpatialDimension) {
            this.SpatialDimension = SpatialDimension;
            if (SpatialDimension == 2)
                TorqueVectorDimension = 1;
            else if (SpatialDimension == 3)
                throw new NotImplementedException("FSI_Solver is currently only capable of 2D.");
            else
                throw new ArithmeticException("Error in spatial dimensions, only dim == 2 and dim == 3 are allowed!");
        }
        [DataMember]
        private readonly int SpatialDimension;
        [DataMember]
        private readonly int TorqueVectorDimension;
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorquePreviousIteration = new List<double[]>();
        [DataMember]
        private readonly List<double[]> m_ForcesAndTorqueWithoutRelaxation = new List<double[]>();

        public bool AlreadyUpdatedLevelSetParameter = false;
        public bool FromLevelSetParameterUpdate = false;

        /// <summary>
        /// Calculation of hydrodynamic forces F and torque T for all particles. 
        /// The relaxation (Aitken-relaxation) is applied to the entire array of all F+T, leading to better convergence.
        /// </summary>
        /// <param name="AllParticles"></param>
        /// <param name="HydrodynamicIntegrator"></param>
        /// <param name="FluidSpecies"></param>
        /// <param name="Gravity"></param>
        /// <param name="dt"></param>
        internal void CalculateHydrodynamics(Particle[] AllParticles, ParticleHydrodynamicsIntegration HydrodynamicIntegrator, string[] FluidSpecies, Vector Gravity, double dt) {
            using (new FuncTrace()) {
                double[] hydrodynamics = new double[(SpatialDimension + TorqueVectorDimension) * AllParticles.Count()];
               
                // Calculate hydrodynamics.
                CellMask allCutCells = AllParticles[0].LsTrk.Regions.GetCutCellMask();
                for (int p = 0; p < AllParticles.Count(); p++) {
                    Particle currentParticle = AllParticles[p];
                    CellMask CutCells = currentParticle.ParticleCutCells(AllParticles[0].LsTrk, allCutCells, 0);
                    int offset = p * (SpatialDimension + TorqueVectorDimension);
                    double[] tempHydrodynamics = currentParticle.Motion.CalculateHydrodynamics(HydrodynamicIntegrator, CutCells, FluidSpecies, dt);
                    for (int d = 0; d < SpatialDimension + 1; d++) 
                        hydrodynamics[offset + d] = tempHydrodynamics[d];
                }
                hydrodynamics = hydrodynamics.MPISum();

                // Gravity, has to be added after MPISum()!
                for (int p = 0; p < AllParticles.Count(); p++) {
                    Particle currentParticle = AllParticles[p];
                    int offset = p * (SpatialDimension + TorqueVectorDimension);
                    Vector gravity = currentParticle.Motion.GetGravityForces(Gravity);
                    for(int d = 0; d < SpatialDimension; d++) {
                        hydrodynamics[offset + d] += gravity[d];
                    }
                }

                // Relaxation
                double omega = AllParticles[0].Motion.omega;
                hydrodynamics = HydrodynamicsRelaxation(hydrodynamics, ref omega);
                AllParticles[0].Motion.omega = omega;

                // Write to single particles.
                for (int p = 0; p < AllParticles.Count(); p++) {
                    Particle currentParticle = AllParticles[p];
                    currentParticle.Motion.UpdateForcesAndTorque(p, hydrodynamics);
                }
            }
        }

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
        /// ...
        /// </summary>
        /// <param name="AllParticles"></param>
        internal void SaveHydrodynamicOfPreviousIteration(Particle[] AllParticles) {
            double[] hydrodynamics = new double[(SpatialDimension + 1) * AllParticles.Length];
            if (m_ForcesAndTorquePreviousIteration.Count() > 5)
                m_ForcesAndTorquePreviousIteration.RemoveAt(5);
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
        /// Post-processing of the hydrodynamics. If desired the underrelaxation is applied to the forces and torque.
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
                    if (variable[d] == 0)//ghost Particle
                        continue;
                    returnVariable[d] = 0.001 * variable[d] + (1 - 0.001) * m_ForcesAndTorquePreviousIteration[1][d];
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
                    if (variable[i] == 0) {// ghost particle
                        residualDiff[i] = 0;
                        continue;
                    }
                    residual[i] = new double[] { (variable[i] - m_ForcesAndTorquePreviousIteration[0][i]), (m_ForcesAndTorqueWithoutRelaxation[1][i] - m_ForcesAndTorquePreviousIteration[1][i]) };
                    residualDiff[i] = residual[i][0] - residual[i][1];
                    residualScalar += residual[i][1] * residualDiff[i];
                }
                Omega = -Omega * residualScalar / residualDiff.L2Norm().Pow2();
                double[] outVar = variable.CloneAs();
                for (int i = 0; i < variable.Length; i++) {
                    if (variable[i] == 0)// ghost particle
                        continue;
                    outVar[i] = Omega * (variable[i] - m_ForcesAndTorquePreviousIteration[0][i]) + m_ForcesAndTorquePreviousIteration[0][i];
                }
                return outVar;
            }
        }
    }
}
