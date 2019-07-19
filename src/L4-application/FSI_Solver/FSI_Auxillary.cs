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

using BoSSS.Application.FSI_Solver;
using BoSSS.Application.IBM_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FSI_Solver
{
    class FSI_Auxillary 
    {
        /// <summary>
        /// This method saves the list value at list position "0" to the next position.
        /// Use this method for onedimensional vars.
        /// </summary>
        /// <param name="variable">
        /// Name of the list.
        /// </param>
        internal void SaveValueOfLastTimestep(List<double> variable)
        {
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
        internal void SaveMultidimValueOfLastTimestep(List<double[]> variable)
        {
            int Dim = variable[0].Length;
            variable.Insert(0, new double[Dim]);
            for (int d = 0; d < Dim; d++)
            {
                variable[0][d] = 0;
            }
            variable.RemoveAt(variable.Count - 1);
        }

        internal void TestArithmeticException(double[,] variable, string variableName)
        {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(double[] variable, string variableName)
        {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(double variable, string variableName)
        {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void ThrowIsNaNException(double[,] variable, string variableName)
        {
            for (int i = 0; i < variable.GetLength(0); i++)
            {
                for (int j = 0; j < variable.GetLength(1); j++)
                {
                    if (double.IsNaN(variable[i, j]))
                        throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
                }
            }
        }

        internal void ThrowIsNaNException(double[] variable, string variableName)
        {
            for (int i = 0; i < variable.Length; i++)
            {
                if (double.IsNaN(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
            }
        }

        internal void ThrowIsNaNException(double variable, string variableName)
        {
            if (double.IsNaN(variable))
                throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
        }

        internal void ThrowIsInfinityException(double[,] variable, string variableName)
        {
            for (int i = 0; i < variable.GetLength(0); i++)
            {
                for (int j = 0; j < variable.GetLength(1); j++)
                {
                    if (double.IsInfinity(variable[i, j]))
                        throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
                }
            }
        }

        internal void ThrowIsInfinityException(double[] variable, string variableName)
        {
            for (int i = 0; i < variable.Length; i++)
            {
                if (double.IsInfinity(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
            }
        }

        internal void ThrowIsInfinityException(double variable, string variableName)
        {
            if (double.IsInfinity(variable))
                throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
        }

        /// <summary>
        /// Calculate componentwise sum of two vectors.
        /// </summary>
        /// <param name="Vector0">
        /// </param>
        /// <param name="Vector1"></param>
        internal double[] VectorSum(double[] Vector0, double[] Vector1)
        {
            int Dim = Vector0 != null ? Vector0.Length : Vector1.Length;
            if (Vector0 == null)
            {
                Vector0 = Vector1.CloneAs();
                for (int d = 0; d < Dim; d++)
                {
                    Vector0[d] = 0;
                }
            }
            if (Vector1 == null)
            {
                Vector1 = Vector0.CloneAs();
                for (int d = 0; d < Dim; d++)
                {
                    Vector1[d] = 0;
                }
            }
            double[] ResultVector = new double[Dim];
            if (Vector0.Length != Vector1.Length)
                throw new ArithmeticException("Mismatch in vector dimension");
            for (int d = 0; d < Dim; d++)
            {
                ResultVector[d] = Vector0[d] + Vector1[d];
            }
            TestArithmeticException(ResultVector, "result of vector summation");
            return ResultVector;
        }

        /// <summary>
        /// Calculate componentwise difference of two vectors.
        /// </summary>
        /// <param name="Vector0">
        /// </param>
        /// <param name="Vector1"></param>
        internal double[] VectorDiff(double[] Vector0, double[] Vector1)
        {
            int Dim = Vector0 != null ? Vector0.Length : Vector1.Length;
            if (Vector0 == null)
            {
                Vector0 = Vector1.CloneAs();
                for (int d = 0; d < Dim; d++)
                {
                    Vector0[d] = 0;
                }
            }
            if (Vector1 == null)
            {
                Vector1 = Vector0.CloneAs();
                for (int d = 0; d < Dim; d++)
                {
                    Vector1[d] = 0;
                }
            }
            double[] ResultVector = new double[Dim];
            if (Vector0.Length != Vector1.Length)
                throw new ArithmeticException("Mismatch in vector dimension");
            for (int d = 0; d < Dim; d++)
            {
                ResultVector[d] = Vector0[d] - Vector1[d];
            }
            TestArithmeticException(ResultVector, "result of vector difference");
            return ResultVector;
        }

        /// <summary>
        /// Calculate teh dot product of two vectors.
        /// </summary>
        /// <param name="Vector0">
        /// </param>
        /// <param name="Vector1"></param>
        internal double DotProduct(double[] Vector0, double[] Vector1)
        {
            int Dim = Vector0.Length;
            double DotProduct = new double();
            if (Vector0.Length != Vector1.Length)
                throw new ArithmeticException("Mismatch in vector dimension");
            for (int d = 0; d < Dim; d++)
            {
                DotProduct += Vector0[d] * Vector1[d];
            }
            TestArithmeticException(DotProduct, "dot product of two vectors");
            return DotProduct;
        }

        /// <summary>
        /// Quicksort algorithm
        /// </summary>
        /// <param name="Leftelement">
        /// The first element of the list to be sorted
        /// </param>
        /// <param name="RightElement">
        /// The last element of the list to be sorted
        /// </param>
        /// <param name="Data">
        /// The data to be sorted.
        /// </param>
        internal void Quicksort(int Leftelement, int RightElement, ref int[] Data)
        {
            if (Leftelement < RightElement)
            {
                int division = Quicksort_Divide(Leftelement, RightElement, ref Data);
                Quicksort(Leftelement, division - 1, ref Data);
                Quicksort(division + 1, RightElement, ref Data);
            }
        }

        /// <summary>
        /// Core routine of the quicksort algorithm, splits the list to be sorted
        /// </summary>
        /// <param name="Leftelement">
        /// The first element of the list to be sorted
        /// </param>
        /// <param name="RightElement">
        /// The last element of the list to be sorted
        /// </param>
        /// <param name="Data">
        /// The data to be sorted.
        /// </param>
        private int Quicksort_Divide(int Leftelement, int RightElement, ref int[] Data)
        {
            int i = Leftelement;
            int j = RightElement - 1;
            int pivot = Data[RightElement];

            do
            {
                while (Data[i] <= pivot && i < RightElement)
                    i += 1;

                while (Data[j] >= pivot && j > Leftelement)
                    j -= 1;

                if (i < j)
                {
                    int z = Data[i];
                    Data[i] = Data[j];
                    Data[j] = z;
                }

            } while (i < j);

            if (Data[i] > pivot)
            {
                int z = Data[i];
                Data[i] = Data[RightElement];
                Data[RightElement] = z;
            }
            return i; 
        }

        /// <summary>
        /// Binary search algorithm
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="Target">
        /// The element to be found
        /// </param>
        int BinarySearch(List<int> SortedList, int Target)
        {
            int L = 0;
            int R = SortedList.Count() - 1;
            int TargetIndex = 0;

            while (L <= R)
            {
                TargetIndex = (L + R) / 2;
                if (SortedList[TargetIndex] < Target)
                    L = TargetIndex + 1;
                else if (SortedList[TargetIndex] > Target)
                    R = TargetIndex - 1;
                else
                    break;
            }

            return TargetIndex;
        }

        /// <summary>
        /// Binary search algorithm, output is the target and its neighbours
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="Target">
        /// The element to be found
        /// </param>
        internal int[] BinarySearchWithNeighbors(List<int> SortedList, int Target)
        {
            int StartIndex = BinarySearch(SortedList, Target);
            if (StartIndex == 0)
                return LeftSearch(SortedList, StartIndex);
            else if (StartIndex == SortedList.Count() - 1)
                return RightSearch(SortedList, StartIndex);
            else
            {
                int[] Right = RightSearch(SortedList, StartIndex);
                int[] Left = LeftSearch(SortedList, StartIndex);
                return Left.Union(Right).ToArray();
            }

        }

        /// <summary>
        /// Search for left neighbours
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="StartIndex"></param>
        private int[] LeftSearch(List<int> SortedList, int StartIndex)
        {
            List<int> temp = new List<int> { StartIndex };
            while (temp.Last() + 1 < SortedList.Count() && SortedList[temp.Last() + 1] == SortedList[temp.Last()])
            {
                temp.Add(temp.Last() + 1);
            }
            return temp.ToArray();
        }

        /// <summary>
        /// Search for right neighbours
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="StartIndex"></param>
        private int[] RightSearch(List<int> SortedList, int StartIndex)
        {
            List<int> temp = new List<int> { StartIndex };
            while (temp.Last() - 1 > 0 && SortedList[temp.Last() - 1] == SortedList[temp.Last()])
            {
                temp.Add(temp.Last() - 1);
            }
            return temp.ToArray();
        }

        /// <summary>
        /// MPI exchange of the damping tensors
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        internal void ExchangeDampingTensors(List<Particle> Particles)
        {
            int NoOfParticles = Particles.Count;
            int NoOfVars = 3;
            double[] StateBuffer = new double[NoOfParticles * NoOfVars * NoOfParticles * NoOfVars];
            for (int p = 0; p < NoOfParticles; p++)
            {
                Particle P = Particles[p];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        StateBuffer[NoOfVars * NoOfVars * p + i + NoOfVars * j] = P.addedDampingTensor[i, j];
                    }
                }

            }
            double[] GlobalStateBuffer = StateBuffer.MPISum();
            for (int p = 0; p < NoOfParticles; p++)
            {
                var P = Particles[p];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        P.addedDampingTensor[i, j] = GlobalStateBuffer[NoOfVars * NoOfVars * p + i + NoOfVars * j];
                    }
                }
            }
        }

        /// <summary>
        /// Calls the calculation of the Acceleration and the velocity
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="dt">
        /// The time step
        /// </param>
        /// <param name="FullyCoupled">
        /// Do you use FSI_Fully_Coupled?
        /// </param>
        /// <param name="IterationCounter">
        /// No of iterations
        /// </param>
        /// <param name="IncludeHydrodynamics"></param>
        internal void CalculateParticleVelocity(List<Particle> Particles, double dt, bool FullyCoupled, int IterationCounter, int TimestepInt, bool IncludeHydrodynamics = true)
        {
            foreach (Particle p in Particles)
            {
                p.iteration_counter_P = IterationCounter;
                if (IterationCounter == 0 && FullyCoupled)
                {
                    Console.WriteLine("Predicting forces for the next timestep...");
                    if (p.UseAddedDamping)
                    {
                        p.UpdateDampingTensors();
                        //ExchangeDampingTensors(Particles);
                    }
                    p.PredictForceAndTorque(TimestepInt);
                }
                p.CalculateAcceleration(dt, FullyCoupled, IncludeHydrodynamics);
                p.UpdateParticleVelocity(dt);
            }
        }

        /// <summary>
        /// Calls the calculation of position
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="dt">
        /// The time step
        /// </param>
        internal void CalculateParticlePosition(List<Particle> Particles, double dt)
        {
            for (int p = 0; p < Particles.Count(); p++)
            {
                Particle CurrentParticle = Particles[p];
                CurrentParticle.CalculateParticlePosition(dt);
                CurrentParticle.CalculateParticleAngle(dt);
                CurrentParticle.CollisionTimestep = 0;
            }
        }

        /// <summary>
        /// Residual for fully coupled system
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="IterationCounter_In"></param>
        /// <param name="MaximumIterations"></param>
        /// <param name="Residual"></param>
        /// <param name="IterationCounter_Out"> </param>
        internal void CalculateParticleResidual(List<Particle> Particles, int IterationCounter_In, int MaximumIterations, out double Residual, out int IterationCounter_Out)
        {
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            if (IterationCounter_In == 0)
                Residual = 1e12;
            else
            {
                double[] ForcesNewSquared = new double[2];
                double TorqueNewSquared = new double();
                Residual = 0;
                foreach (Particle p in Particles)
                {
                    p.ForceTorqueResidual = Math.Sqrt((p.forcesPrevIteration[0] - p.hydrodynamicForces[0][0]).Pow2() + (p.forcesPrevIteration[1] - p.hydrodynamicForces[0][1]).Pow2() + (p.torquePrevIteration - p.hydrodynamicTorque[0]).Pow2());
                    for (int d = 0; d < 2; d++)
                        ForcesNewSquared[d] += p.hydrodynamicForces[0][d].Pow2();
                    TorqueNewSquared += p.hydrodynamicTorque[0].Pow2();
                    Residual += p.ForceTorqueResidual;
                }
            }
            Residual /= Particles.Count();
            if (IterationCounter_In != 0)
            {
                Console.WriteLine();
                Console.WriteLine("Forces and torque residual: " + Residual);
                Console.WriteLine("=======================================================");
                Console.WriteLine();
            }
            if (IterationCounter_In > MaximumIterations)
                throw new ApplicationException("No convergence in coupled iterative solver, number of iterations: " + IterationCounter_In);
            IterationCounter_Out = IterationCounter_In + 1;
        }

        /// <summary>
        /// Residual for fully coupled system
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="FluidViscosity"></param>
        /// <param name="phystime"></param>
        /// <param name="TimestepInt"></param>
        /// <param name="IterationCounter"> </param>
        /// /// <param name="Finalresult"></param>
        /// <param name="MPIangularVelocity"></param>
        /// <param name="Force"></param>
        internal void PrintResultToConsole(List<Particle> Particles, double phystime, int IterationCounter, out double MPIangularVelocity, out double[] Force)
        {
            if (Particles.Count() == 0)
            {
                MPIangularVelocity = 0;
                Force = new double[1];
                return;
            }
            double[] TranslationalMomentum = new double[2] { 0, 0 };
            double RotationalMomentum = 0;
            double[] totalKE = new double[3] { 0, 0, 0 };

            for (int p = 0; p < Particles.Count(); p++) 
            {
                Particle CurrentParticle = Particles[p];
                double[] SingleParticleMomentum = CurrentParticle.CalculateParticleMomentum();
                double[] SingleParticleKineticEnergy = CurrentParticle.CalculateParticleKineticEnergy();
                TranslationalMomentum[0] += SingleParticleMomentum[0];
                TranslationalMomentum[1] += SingleParticleMomentum[1];
                RotationalMomentum += SingleParticleMomentum[SingleParticleMomentum.Length - 1];
                totalKE[0] += SingleParticleKineticEnergy[0];
                totalKE[1] += SingleParticleKineticEnergy[1];
                totalKE[2] += SingleParticleKineticEnergy[SingleParticleMomentum.Length - 1];
            }

            Force = Particles[0].hydrodynamicForces[0];
            MPIangularVelocity = Particles[0].rotationalVelocity[0];

            StringBuilder OutputBuilder = new StringBuilder();

            OutputBuilder.AppendLine("Total kinetic energy in system:  " + (totalKE[0] + totalKE[1] + totalKE[2]));
            OutputBuilder.AppendLine("Total momentum in system:  " + Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
            OutputBuilder.AppendLine("Total kinetic energy in system:  " + (totalKE[0] + totalKE[1] + totalKE[2]));
            for (int p = 0; p < Particles.Count(); p++)
            {
                Particle CurrentParticle = Particles[p];
                // only print particles with some action
                if(CurrentParticle.IncludeTranslation || CurrentParticle.IncludeRotation)
                {
                    int PrintP = p + 1;
                    OutputBuilder.AppendLine("=======================================================");
                    OutputBuilder.AppendLine("Status report particle #" + PrintP + ", Time: " + phystime + ", Iteration #" + IterationCounter);
                    if (CurrentParticle.Collided)
                        OutputBuilder.AppendLine("The particle is collided");
                    OutputBuilder.AppendLine("-------------------------------------------------------");
                    OutputBuilder.AppendLine("Drag Force: " + CurrentParticle.hydrodynamicForces[0][0]);
                    OutputBuilder.AppendLine("Lift Force: " + CurrentParticle.hydrodynamicForces[0][1]);
                    OutputBuilder.AppendLine("Torqe: " + CurrentParticle.hydrodynamicTorque[0]);
                    OutputBuilder.AppendLine("Transl VelocityX: " + CurrentParticle.translationalVelocity[0][0]);
                    OutputBuilder.AppendLine("Transl VelocityY: " + CurrentParticle.translationalVelocity[0][1]);
                    OutputBuilder.AppendLine("Angular Velocity: " + CurrentParticle.rotationalVelocity[0]);
                }
            }
            Console.WriteLine(OutputBuilder.ToString());
        }

        /// <summary>
        /// Residual for fully coupled system
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="FluidViscosity"></param>
        /// <param name="phystime"></param>
        /// <param name="TimestepInt"></param>
        /// <param name="IterationCounter"> </param>
        /// /// <param name="Finalresult"></param>
        /// <param name="MPIangularVelocity"></param>
        /// <param name="Force"></param>
        internal void PrintResultToConsole(List<Particle> Particles, double FluidViscosity, double FluidDensity, double phystime, int TimestepInt, out double MPIangularVelocity, out double[] Force)
        {
            if (Particles.Count() == 0)
            {
                MPIangularVelocity = 0;
                Force = new double[1];
                return;
            }
            double[] TranslationalMomentum = new double[2] { 0, 0 };
            double RotationalMomentum = 0;
            double[] totalKE = new double[3] { 0, 0, 0 };
            double[] ParticleReynoldsNumber = new double[Particles.Count()];
            double[] ParticleStokesNumber = new double[Particles.Count()];

            for (int p = 0; p < Particles.Count(); p++)
            {
                Particle CurrentParticle = Particles[p];
                double[] SingleParticleMomentum = CurrentParticle.CalculateParticleMomentum();
                double[] SingleParticleKineticEnergy = CurrentParticle.CalculateParticleKineticEnergy();
                TranslationalMomentum[0] += SingleParticleMomentum[0];
                TranslationalMomentum[1] += SingleParticleMomentum[1];
                RotationalMomentum += SingleParticleMomentum[SingleParticleMomentum.Length - 1];
                totalKE[0] += SingleParticleKineticEnergy[0];
                totalKE[1] += SingleParticleKineticEnergy[1];
                totalKE[2] += SingleParticleKineticEnergy[SingleParticleMomentum.Length - 1];
                ParticleReynoldsNumber[Particles.IndexOf(CurrentParticle)] = CurrentParticle.ComputeParticleRe(FluidViscosity);
                ParticleStokesNumber[Particles.IndexOf(CurrentParticle)] = CurrentParticle.ComputeParticleSt(FluidViscosity, FluidDensity);
            }

            Force = Particles[0].hydrodynamicForces[0];
            MPIangularVelocity = Particles[0].rotationalVelocity[0];

            StringBuilder OutputBuilder = new StringBuilder();

            OutputBuilder.AppendLine("Total kinetic energy in system:  " + (totalKE[0] + totalKE[1] + totalKE[2]));
            OutputBuilder.AppendLine("Total momentum in system:  " + Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
            OutputBuilder.AppendLine("Total kinetic energy in system:  " + (totalKE[0] + totalKE[1] + totalKE[2]));
            for (int p = 0; p < Particles.Count(); p++)
            {
                Particle CurrentParticle = Particles[p];
                // only print particles with some action
                if (CurrentParticle.IncludeTranslation || CurrentParticle.IncludeRotation)
                {
                    int PrintP = p + 1;
                    OutputBuilder.AppendLine("=======================================================");
                    OutputBuilder.AppendLine("Final status report for timestep #" + TimestepInt + ", particle #" + PrintP + ", Time: " + phystime);
                    if (CurrentParticle.Collided)
                        OutputBuilder.AppendLine("The particle is collided");
                    OutputBuilder.AppendLine("-------------------------------------------------------");
                    OutputBuilder.AppendLine("Drag Force: " + CurrentParticle.hydrodynamicForces[0][0]);
                    OutputBuilder.AppendLine("Lift Force: " + CurrentParticle.hydrodynamicForces[0][1]);
                    OutputBuilder.AppendLine("Torqe: " + CurrentParticle.hydrodynamicTorque[0]);
                    OutputBuilder.AppendLine("Transl VelocityX: " + CurrentParticle.translationalVelocity[0][0]);
                    OutputBuilder.AppendLine("Transl VelocityY: " + CurrentParticle.translationalVelocity[0][1]);
                    OutputBuilder.AppendLine("Angular Velocity: " + CurrentParticle.rotationalVelocity[0]);
                    OutputBuilder.AppendLine("X-position: " + CurrentParticle.position[0][0]);
                    OutputBuilder.AppendLine("Y-position: " + CurrentParticle.position[0][1]);
                    OutputBuilder.AppendLine("Angle: " + CurrentParticle.angle[0]);
                    OutputBuilder.AppendLine();
                    OutputBuilder.AppendLine("Particle Reynolds number: " + ParticleReynoldsNumber[p]);
                    OutputBuilder.AppendLine("Particle Stokes number: " + ParticleStokesNumber[p]);
                    OutputBuilder.AppendLine();
                    OutputBuilder.AppendLine("=======================================================");
                    OutputBuilder.AppendLine();
                }
            }
            Console.WriteLine(OutputBuilder.ToString());
        }

        /// <summary>
        /// Does what it say
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        /// <param name="IterationCounter"></param>
        /// <param name="ForceTorqueConvergenceCriterion"></param>
        /// /// <param name="IsFullyCoupled"></param>
        internal void SaveOldParticleState(List<Particle> Particles, int IterationCounter, double ForceTorqueConvergenceCriterion, bool IsFullyCoupled)
        {
            foreach (Particle p in Particles)
            {
                p.iteration_counter_P = IterationCounter;
                // Save the old hydrondynamic forces, only necessary if no iteration is applied
                // ============================================================================
                if (IterationCounter == 0 && IsFullyCoupled == false)
                {
                    p.Aux.SaveMultidimValueOfLastTimestep(p.hydrodynamicForces);
                    p.Aux.SaveValueOfLastTimestep(p.hydrodynamicTorque);
                }
                // Save status for residual
                // ========================
                p.forceAndTorque_convergence = ForceTorqueConvergenceCriterion;
                p.forcesPrevIteration[0] = p.hydrodynamicForces[0][0];
                p.forcesPrevIteration[1] = p.hydrodynamicForces[0][1];
                p.torquePrevIteration = p.hydrodynamicTorque[0];
            }
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
        internal void ParticleState_MPICheck(List<Particle> Particles, IGridData GridData, int MPISize)
        {
            int D = GridData.SpatialDimension;
            int NoOfParticles = Particles.Count;
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            {
                // verify that we have the same number of particles on each processor
                int NoOfParticles_min = NoOfParticles.MPIMin();
                int NoOfParticles_max = NoOfParticles.MPIMax();
                if (NoOfParticles_min != NoOfParticles || NoOfParticles_max != NoOfParticles)
                    throw new ApplicationException("mismatch in number of MPI particles");

                // nor, compare those particles:
                int NoOfVars = (10 + D * 10); // variables per particle; size can be increased if more values should be compared
                double[] CheckSend = new double[NoOfParticles * NoOfVars];

                for (int p = 0; p < NoOfParticles; p++)
                {
                    var P = Particles[p];

                    // scalar values
                    CheckSend[p * NoOfVars + 0] = P.angle[0];
                    CheckSend[p * NoOfVars + 1] = P.angle[1];
                    CheckSend[p * NoOfVars + 2] = P.rotationalVelocity[0];
                    //CheckSend[iP*NoOfVars + 2] = P.Area_P;
                    CheckSend[p * NoOfVars + 3] = P.clearSmallValues ? 1.0 : 0.0;
                    CheckSend[p * NoOfVars + 4] = P.forceAndTorque_convergence;
                    CheckSend[p * NoOfVars + 5] = P.Mass_P;
                    CheckSend[p * NoOfVars + 6] = P.particleDensity;
                    CheckSend[p * NoOfVars + 7] = P.addedDampingTensor[0,0];
                    // todo: add more values here that might be relevant for the particle state;

                    // vector values
                    for (int d = 0; d < D; d++)
                    {
                        int Offset = 10;
                        CheckSend[p * NoOfVars + Offset + 0 * D + d] = P.position[0][d];
                        CheckSend[p * NoOfVars + Offset + 1 * D + d] = P.position[1][d];
                        CheckSend[p * NoOfVars + Offset + 2 * D + d] = P.translationalAcceleration[0][d];
                        CheckSend[p * NoOfVars + Offset + 3 * D + d] = P.translationalAcceleration[1][d];
                        CheckSend[p * NoOfVars + Offset + 4 * D + d] = P.translationalVelocity[0][d];
                        CheckSend[p * NoOfVars + Offset + 5 * D + d] = P.translationalVelocity[1][d];
                        CheckSend[p * NoOfVars + Offset + 6 * D + d] = P.hydrodynamicForces[0][d];
                        CheckSend[p * NoOfVars + Offset + 7 * D + d] = P.hydrodynamicForces[1][d];
                        // todo: add more vector values here that might be relevant for the particle state;
                    }
                }

                double[] CheckReceive = new double[NoOfParticles * NoOfVars * MPISize];
                unsafe
                {
                    fixed (double* pCheckSend = CheckSend, pCheckReceive = CheckReceive)
                    {
                        csMPI.Raw.Allgather((IntPtr)pCheckSend, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                    }
                }

                for (int iP = 0; iP < NoOfParticles; iP++)
                {
                    for (int iVar = 0; iVar < NoOfVars; iVar++)
                    {
                        // determine a tolerance...
                        int idx_l =
                            iP * NoOfVars // particle index offset
                           + iVar; // variable index offset
                        double VarTol = Math.Abs(CheckSend[idx_l]) * 1.0e-10;

                        // compare
                        for (int r = 0; r < MPISize; r++)
                        {

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
