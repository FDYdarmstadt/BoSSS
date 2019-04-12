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
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FSI_Solver
{
    class FSI_Auxillary 
    {
        internal void Quicksort(int Leftelement, int RightElement, ref int[] Data)
        {
            if (Leftelement < RightElement)
            {
                int division = Divide(Leftelement, RightElement, ref Data);
                Quicksort(Leftelement, division - 1, ref Data);
                Quicksort(division + 1, RightElement, ref Data);
            }
        }

        private int Divide(int Leftelement, int RightElement, ref int[] Data)
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

        int BinarySearch(List<int> SortedList, int Target)
        {
            int L = 0;
            int R = SortedList.Count();
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

        int[] LeftSearch(List<int> SortedList, int StartIndex)
        {
            List<int> temp = new List<int> { StartIndex };
            while (temp.Last() + 1 < SortedList.Count() && SortedList[temp.Last() + 1] == SortedList[temp.Last()])
            {
                temp.Add(temp.Last() + 1);
            }
            return temp.ToArray();
        }

        int[] RightSearch(List<int> SortedList, int StartIndex)
        {
            List<int> temp = new List<int> { StartIndex };
            while (temp.Last() - 1 > 0 && SortedList[temp.Last() - 1] == SortedList[temp.Last()])
            {
                temp.Add(temp.Last() - 1);
            }
            return temp.ToArray();
        }

        internal void ExchangeDampingTensors(List<Particle> Particles)
        {
            // Sum forces and moments over all MPI processors
            // ==============================================
            {
                // step 1: collect all variables that we need to sum up
                int NoOfParticles = Particles.Count;
                int NoOfVars = 3; //only for 2D at the moment
                double[] StateBuffer = new double[NoOfParticles * NoOfVars * NoOfParticles * NoOfVars];
                for (int p = 0; p < NoOfParticles; p++)
                {
                    Particle P = Particles[p];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            StateBuffer[NoOfVars * NoOfVars * p + i + NoOfVars * j] = P.AddedDampingTensor[i, j];
                        }
                    }

                }
                // step 2: sum over MPI processors
                // note: we want to sum all variables by a single MPI call, which is way more efficient
                // B. Deußen: a single call of MPISum() would only consider the first entry of StateBuffer, thus I implemented the loop over all entries
                //double[,] GlobalStateBuffer = new double[NoOfParticles * NoOfVars, NoOfParticles * NoOfVars];
                //for (int i = 0; i < NoOfParticles * NoOfVars; i++)
                //{
                //    for (int j = 0; j < NoOfParticles * NoOfVars; j++)
                //    {
                //        GlobalStateBuffer[i, j] = StateBuffer[i, j].MPISum();
                //    }

                //}
                double[] GlobalStateBuffer = StateBuffer.MPISum();
                // step 3: write sum variables back 
                for (int p = 0; p < NoOfParticles; p++)
                {
                    var P = Particles[p];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            P.AddedDampingTensor[i, j] = GlobalStateBuffer[NoOfVars * NoOfVars * p + i + NoOfVars * j];
                        }
                    }
                }
            }
        }

        internal void UpdateParticleAccelerationAndDamping(List<Particle> Particles, int IterationCounter, double dt, bool LieSplittingFullyCoupled)
        {
            for (int p = 0; p < Particles.Count(); p++)
            {
                if (IterationCounter == 0 && LieSplittingFullyCoupled)
                {
                    if (Particles[p].neglectAddedDamping == false)
                    {
                        Particles[p].UpdateDampingTensors();
                        //ExchangeDampingTensors(Particles);
                    }
                    Particles[p].PredictAcceleration();
                }
                else
                {
                    Particles[p].CalculateAcceleration(dt);
                }
            }
        }

        internal void CalculateParticleResidual(List<Particle> Particles, double[] ForcesOldSquared, double TorqueOldSquared, int IterationCounter, int MaximumIterations, out double Residual, out int _IterationCounter)
        {
            if (IterationCounter == 0)
                Residual = 1e12;
            else
            {
                double[] ForcesNewSquared = new double[2];
                double TorqueNewSquared = new double();
                foreach (Particle p in Particles)
                {
                    for (int d = 0; d < 2; d++)
                        ForcesNewSquared[d] += p.HydrodynamicForces[0][d].Pow2();
                    TorqueNewSquared += p.HydrodynamicTorque[0].Pow2();
                }
                Residual = Math.Sqrt((Math.Sqrt(ForcesNewSquared[0]) - Math.Sqrt(ForcesOldSquared[0])).Pow2() + (Math.Sqrt(ForcesNewSquared[1]) - Math.Sqrt(ForcesOldSquared[1])).Pow2() + (Math.Sqrt(TorqueNewSquared) - Math.Sqrt(TorqueOldSquared)).Pow2());
            }
            int PrintIteration = IterationCounter + 1;
            Console.WriteLine("Fully coupled system, number of iterations:  " + PrintIteration);
            Console.WriteLine("Forces and torque residual: " + Residual);
            Console.WriteLine();
            if (IterationCounter > MaximumIterations)
                throw new ApplicationException("no convergence in coupled iterative solver, number of iterations: " + IterationCounter);
            _IterationCounter = IterationCounter + 1;
        }

        void PrintResultToConsole(List<Particle> Particles, double phystime, double dt, out double MPIangularVelocity, out double[] force)
        {
            double[] TranslationalMomentum = new double[2] { 0, 0 };
            double RotationalMomentum = 0;
            double[] totalKE = new double[3] { 0, 0, 0 };

            foreach (Particle p in Particles)
            {
                double[] SingleParticleMomentum = p.CalculateParticleMomentum(dt);
                double[] SingleParticleKineticEnergy = p.CalculateParticleKineticEnergy(dt);
                TranslationalMomentum[0] += SingleParticleMomentum[0];
                TranslationalMomentum[1] += SingleParticleMomentum[1];
                RotationalMomentum += SingleParticleMomentum[SingleParticleMomentum.Length - 1];
                totalKE[0] += SingleParticleKineticEnergy[0];
                totalKE[1] += SingleParticleKineticEnergy[1];
                totalKE[2] += SingleParticleKineticEnergy[SingleParticleMomentum.Length - 1];
            }

            Console.WriteLine("Total momentum in system:  " + Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
            Console.WriteLine("Total kinetic energy in system:  " + (totalKE[0] + totalKE[1] + totalKE[2]));

            force = Particles[0].HydrodynamicForces[0];
            MPIangularVelocity = Particles[0].RotationalVelocity[0];

            /*
            if ((base.MPIRank == 0) && (Log_DragAndLift != null)) {
                double drag = force[0];
                double lift = force[1];
                //string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", TimestepNo, phystime, m_Particles[0].positionAtIteration[0][0], m_Particles[0].positionAtIteration[0][1], m_Particles[0].particleAnglePerIteration[0], m_Particles[0].transVelocityAtIteration[0][0], m_Particles[0].transVelocityAtIteration[0][1], 0.0, (totalKE[0] + totalKE[1] + totalKE[2]), Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", phystime, m_Particles[0].Position[0][0], m_Particles[0].Position[0][1], m_Particles[0].Angle[0], m_Particles[0].TranslationalVelocity[0][0], m_Particles[0].TranslationalVelocity[0][1], 0.0, (totalKE[0] + totalKE[1] + totalKE[2]), Math.Sqrt(TranslationalMomentum[0].Pow2() + TranslationalMomentum[1].Pow2()));
                Log_DragAndLift.WriteLine(line);
                Log_DragAndLift.Flush();
            }
            */
            for (int p = 0; p < Particles.Count(); p++)
            {
                Particle CurrentP = Particles[p];
                int PrintP = p + 1;
                Console.WriteLine("=======================================================");
                Console.WriteLine("Status report particle #" + PrintP + ",Time: " + phystime);
                Console.WriteLine("-------------------------------------------------------");
                Console.WriteLine("Drag Force:   {0}", CurrentP.HydrodynamicForces[0][0]);
                Console.WriteLine("Lift Force:   {0}", CurrentP.HydrodynamicForces[0][1]);
                Console.WriteLine("Torqe:   {0}", CurrentP.HydrodynamicTorque[0]);
                Console.WriteLine("Transl VelocityX:   {0}", CurrentP.TranslationalVelocity[0][0]);
                Console.WriteLine("Transl VelocityY:   {0}", CurrentP.TranslationalVelocity[0][1]);
                Console.WriteLine("Angular Velocity:   {0}", CurrentP.RotationalVelocity[0]);
                Console.WriteLine("X-position:   {0}", CurrentP.Position[0][0]);
                Console.WriteLine("Y-position:   {0}", CurrentP.Position[0][1]);
                Console.WriteLine("Angle:   {0}", CurrentP.Angle[0]);
                Console.WriteLine();
                Console.WriteLine("=======================================================");
                Console.WriteLine();
            }
        }
    }
}
