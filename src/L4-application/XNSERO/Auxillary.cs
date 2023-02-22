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
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Application.XNSERO_Solver {
    /// <summary>
    /// Helper class for the FSISolver. Contains additional methods for testing.
    /// </summary>
    [Serializable]
    public class Auxillary {
        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for one dimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        internal void SaveValueOfLastTimestep(List<double> variable) {
            variable.Insert(0, new double());
            variable[0] = 0;
            variable.RemoveAt(variable.Count - 1);
        }

        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for multidimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        internal void SaveVectorOfLastTimestep(List<Vector> variable) {
            int dim = variable[0].Dim;
            variable.Insert(0, new Vector(dim));
            variable.RemoveAt(variable.Count - 1);
        }

        internal void TestArithmeticException(double[,] variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(double[] variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(Vector variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(double variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void ThrowIsNaNException(double[,] variable, string variableName) {
            for (int i = 0; i < variable.GetLength(0); i++) {
                for (int j = 0; j < variable.GetLength(1); j++) {
                    if (double.IsNaN(variable[i, j]))
                        throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
                }
            }
        }

        internal void ThrowIsNaNException(double[] variable, string variableName) {
            for (int i = 0; i < variable.Length; i++) {
                if (double.IsNaN(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
            }
        }

        internal void ThrowIsNaNException(Vector variable, string variableName) {
            for (int i = 0; i < variable.Dim; i++) {
                if (double.IsNaN(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
            }
        }

        internal void ThrowIsInfinityException(Vector variable, string variableName) {
            for (int i = 0; i < variable.Dim; i++) {
                if (double.IsInfinity(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
            }
        }

        internal void ThrowIsNaNException(double variable, string variableName) {
            if (double.IsNaN(variable))
                throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
        }

        internal void ThrowIsInfinityException(double[,] variable, string variableName) {
            for (int i = 0; i < variable.GetLength(0); i++) {
                for (int j = 0; j < variable.GetLength(1); j++) {
                    if (double.IsInfinity(variable[i, j]))
                        throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
                }
            }
        }

        internal void ThrowIsInfinityException(double[] variable, string variableName) {
            for (int i = 0; i < variable.Length; i++) {
                if (double.IsInfinity(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
            }
        }

        internal void ThrowIsInfinityException(double variable, string variableName) {
            if (double.IsInfinity(variable))
                throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
        }

        /// <summary>
        /// Check consistency of particle properties on all MPI-processes
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="GridData">
        /// IGridData
        /// </param>
        /// <param name="MPISize">
        /// No of processes
        /// </param>
        internal void ParticleStateMPICheck(Particle[] Particles, IGridData GridData, int MPISize, int TimeStepNo) {
            using (new FuncTrace()) {
                if ((TimeStepNo % 10) != 0)// Do this check every ten time-steps to save some time.
                    return;
                int D = GridData.SpatialDimension;
                int NoOfParticles = Particles.Length;
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                {
                    // verify that we have the same number of particles on each processor
                    int NoOfParticles_min = NoOfParticles.MPIMin();
                    int NoOfParticles_max = NoOfParticles.MPIMax();
                    if (NoOfParticles_min != NoOfParticles || NoOfParticles_max != NoOfParticles)
                        throw new ApplicationException("mismatch in number of MPI particles");

                    // now, compare those particles:
                    int matrixDim = Particles[0].Motion.AddedDampingTensor.GetLength(0);
                    int NoOfVars = (7 + D * 8 + matrixDim * matrixDim); // variables per particle; size can be increased if more values should be compared
                    double[] CheckSend = new double[NoOfParticles * NoOfVars];

                    for (int p = 0; p < NoOfParticles; p++) {
                        var P = Particles[p];

                        // scalar values
                        CheckSend[p * NoOfVars + 0] = P.Motion.GetAngle(0);
                        CheckSend[p * NoOfVars + 1] = P.Motion.GetAngle(1);
                        CheckSend[p * NoOfVars + 2] = P.Motion.GetRotationalVelocity(0);
                        CheckSend[p * NoOfVars + 3] = P.Motion.GetRotationalVelocity(1);
                        CheckSend[p * NoOfVars + 4] = P.Motion.GetRotationalAcceleration(0);
                        CheckSend[p * NoOfVars + 5] = P.Motion.GetRotationalAcceleration(1);
                        CheckSend[p * NoOfVars + 6] = P.Motion.Mass;

                        // vector values
                        int Offset = 7;
                        for (int d = 0; d < D; d++) {
                            CheckSend[p * NoOfVars + Offset + 0 * D + d] = P.Motion.GetPosition(0)[d];
                            CheckSend[p * NoOfVars + Offset + 1 * D + d] = P.Motion.GetPosition(1)[d];
                            CheckSend[p * NoOfVars + Offset + 2 * D + d] = P.Motion.GetTranslationalVelocity(0)[d];
                            CheckSend[p * NoOfVars + Offset + 3 * D + d] = P.Motion.GetTranslationalVelocity(1)[d];
                            CheckSend[p * NoOfVars + Offset + 4 * D + d] = P.Motion.GetTranslationalAcceleration(0)[d];
                            CheckSend[p * NoOfVars + Offset + 5 * D + d] = P.Motion.GetTranslationalAcceleration(1)[d];
                            CheckSend[p * NoOfVars + Offset + 6 * D + d] = P.Motion.GetHydrodynamicForces(0)[d];
                            CheckSend[p * NoOfVars + Offset + 7 * D + d] = P.Motion.GetHydrodynamicForces(1)[d];
                        }

                        // matrix values
                        Offset += D * 8;
                        for (int d1 = 0; d1 < matrixDim; d1++) {
                            for (int d2 = 0; d2 < matrixDim; d2++) {
                                CheckSend[p * NoOfVars + Offset + d1 + d1 * d2] = P.Motion.AddedDampingTensor[d1, d2];
                            }
                        }
                    }

                    double[] CheckReceive = new double[NoOfParticles * NoOfVars * MPISize];
                    unsafe {
                        fixed (double* pCheckSend = CheckSend, pCheckReceive = CheckReceive) {
                            csMPI.Raw.Allgather((IntPtr)pCheckSend, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                        }
                    }

                    for (int iP = 0; iP < NoOfParticles; iP++) {
                        for (int iVar = 0; iVar < NoOfVars; iVar++) {
                            // determine a tolerance...
                            int idx_l =
                                iP * NoOfVars // particle index offset
                               + iVar; // variable index offset
                            double VarTol = Math.Abs(CheckSend[idx_l]) * 1.0e-10;

                            // compare
                            for (int r = 0; r < MPISize; r++) {

                                int idx_g = CheckSend.Length * r // MPI index offset
                                    + idx_l;

                                if (Math.Abs(CheckReceive[idx_g] - CheckSend[idx_l]) > VarTol)
                                    throw new ApplicationException("Mismatch in particle state among MPI ranks. Index:  " + idx_l + " iP " + iP + " NoOfVars " + NoOfVars + " iVar " + iVar + " CheckReceive[idx_g] " + CheckReceive[idx_g] + " CheckSend[idx_l] " + CheckSend[idx_l] + " idx_g " + idx_g + "idx_l" + idx_l + " r " + r);
                            }
                        }
                    }
                }
            }
        }
    }
}
