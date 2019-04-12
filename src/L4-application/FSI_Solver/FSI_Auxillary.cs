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
    }
}
