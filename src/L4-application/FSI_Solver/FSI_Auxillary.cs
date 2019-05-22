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
using BoSSS.Solution.XdgTimestepping;
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

        internal void UpdateParticleAccelerationAndDamping(Particle Particle, int IterationCounter, double dt, bool LieSplittingFullyCoupled) {

            if (IterationCounter == 0 && LieSplittingFullyCoupled) {
                if (Particle.neglectAddedDamping == false) {
                    Particle.UpdateDampingTensors();
                    //ExchangeDampingTensors(Particles);
                }
                Particle.PredictAcceleration();
            } else {
                Particle.CalculateAcceleration(dt, LieSplittingFullyCoupled);
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

        internal void PrintResultToConsole(List<Particle> Particles, double phystime, double dt, out double MPIangularVelocity, out double[] force)
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

        internal void SaveOldParticleState(List<Particle> Particles, int IterationCounter, int SpatialDim, double ForceTorqueConvergenceCriterion, bool IsFullyCoupled, out double[] ForcesOldSquared, out double TorqueOldSquared)
        {
            ForcesOldSquared = new double[SpatialDim];
            TorqueOldSquared = 0;
            foreach (Particle p in Particles)
            {
                p.iteration_counter_P = IterationCounter;
                // Save the old hydrondynamic forces, only necessary if no iteration is applied
                // ============================================================================
                if (IterationCounter == 0 && IsFullyCoupled == false)
                {
                    p.Aux.SaveMultidimValueOfLastTimestep(p.HydrodynamicForces);
                    p.Aux.SaveValueOfLastTimestep(p.HydrodynamicTorque);
                }
                // Save status for residual
                // ========================
                p.ForceAndTorque_convergence = ForceTorqueConvergenceCriterion;
                ForcesOldSquared[0] += p.HydrodynamicForces[0][0].Pow2();
                ForcesOldSquared[1] += p.HydrodynamicForces[0][1].Pow2();
                TorqueOldSquared += p.HydrodynamicTorque[0].Pow2();
                p.ForcesPrevIteration[0] = p.HydrodynamicForces[0][0];
                p.ForcesPrevIteration[1] = p.HydrodynamicForces[0][1];
                p.TorquePrevIteration = p.HydrodynamicTorque[0];
            }
        }

        // Initial check: is the motion state of the particles equal on all MPI processors?
        // ================================================================================
        internal void ParticleState_MPICheck(List<Particle> Particles, IGridData GridData, int MPISize)
        {
            int D = GridData.SpatialDimension;
            int NoOfParticles = Particles.Count;

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
                    CheckSend[p * NoOfVars + 0] = P.Angle[0];
                    CheckSend[p * NoOfVars + 1] = P.Angle[1];
                    //CheckSend[iP*NoOfVars + 2] = P.Area_P;
                    CheckSend[p * NoOfVars + 3] = P.ClearSmallValues ? 1.0 : 0.0;
                    CheckSend[p * NoOfVars + 4] = P.ForceAndTorque_convergence;
                    CheckSend[p * NoOfVars + 5] = P.Mass_P;
                    CheckSend[p * NoOfVars + 6] = P.particleDensity;
                    // todo: add more values here that might be relevant for the particle state;

                    // vector values
                    for (int d = 0; d < D; d++)
                    {
                        int Offset = 10;
                        CheckSend[p * NoOfVars + Offset + 0 * D + d] = P.Position[0][d];
                        CheckSend[p * NoOfVars + Offset + 1 * D + d] = P.Position[1][d];
                        CheckSend[p * NoOfVars + Offset + 2 * D + d] = P.TranslationalAcceleration[0][d];
                        CheckSend[p * NoOfVars + Offset + 3 * D + d] = P.TranslationalAcceleration[1][d];
                        CheckSend[p * NoOfVars + Offset + 4 * D + d] = P.TranslationalVelocity[0][d];
                        CheckSend[p * NoOfVars + Offset + 5 * D + d] = P.TranslationalVelocity[1][d];
                        CheckSend[p * NoOfVars + Offset + 6 * D + d] = P.HydrodynamicForces[0][d];
                        CheckSend[p * NoOfVars + Offset + 7 * D + d] = P.HydrodynamicForces[1][d];
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
                                throw new ApplicationException("Mismatch in particle state among MPI ranks. Index:  " + idx_l);
                        }
                    }
                }
            }
        }

        internal void Collision_MPICommunication(List<Particle> Particles, Particle CurrentParticle, int MPISize, bool WallCollision = false)
        {
            int NoOfVars = 3;
            double[] BoolSend = new double[1];
            bool NoCurrentCollision = true;
            if (CurrentParticle.m_collidedWithWall[0] && WallCollision)
            {
                BoolSend[0] = -1;
            }
            else
            {
                for (int p = 0; p < Particles.Count(); p++)
                {
                    if (CurrentParticle.m_collidedWithParticle[p])
                    {
                        BoolSend[0] = p + 1;
                    }

                }
            }

            double[] BoolReceive = new double[MPISize];
            unsafe
            {
                fixed (double* pCheckSend = BoolSend, pCheckReceive = BoolReceive)
                {
                    csMPI.Raw.Allgather((IntPtr)pCheckSend, BoolSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, BoolSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                }
            }
            for (int i = 0; i < BoolReceive.Length; i++)
            {
                if (BoolReceive[i] != 0)
                {
                    double[] CheckSend = new double[NoOfVars];
                    CheckSend[0] = CurrentParticle.RotationalVelocity[0];
                    CheckSend[1] = CurrentParticle.TranslationalVelocity[0][0];
                    CheckSend[2] = CurrentParticle.TranslationalVelocity[0][1];

                    double[] CheckReceive = new double[NoOfVars * MPISize];
                    unsafe
                    {
                        fixed (double* pCheckSend = CheckSend, pCheckReceive = CheckReceive)
                        {
                            csMPI.Raw.Allgather((IntPtr)pCheckSend, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                        }
                    }
                    CurrentParticle.RotationalVelocity[0] = CheckReceive[0 + i * 3];
                    CurrentParticle.TranslationalVelocity[0][0] = CheckReceive[1 + i * 3];
                    CurrentParticle.TranslationalVelocity[0][1] = CheckReceive[2 + i * 3];
                    if (!WallCollision)
                    {
                        int p = Convert.ToInt32(BoolReceive[i]);
                        CurrentParticle.m_collidedWithParticle[p - 1] = true;
                        CurrentParticle.skipForceIntegration = true;
                    }
                    NoCurrentCollision = false;
                }
            }
            if (NoCurrentCollision)
            {
                //CurrentParticle.skipForceIntegration = false;
                CurrentParticle.m_collidedWithWall[0] = false;
                for (int p = 0; p < CurrentParticle.m_collidedWithParticle.Length; p++)
                {
                    CurrentParticle.m_collidedWithParticle[p] = false;
                }
            }
        }

        internal void GJK_DistanceAlgorithm(Particle p0, Particle p1, double[] Point0_old, double[] Point1_old, int SpatialDim, out double Min_Distance, out double[] ClosestPoint0, out double[] ClosestPoint1)
        {
            ClosestPoint0 = new double[SpatialDim];
            ClosestPoint1 = new double[SpatialDim];
            Initialize_GJK(SpatialDim, Point0_old, Point1_old, out double[] v0, out List<double[]> Simplex);
            double[] v = v0.CloneAs();
            double[] SupportPoint = new double[SpatialDim];
            
            for (int i = 0; i < 1000; i++)
            {
                double[] vt = v.CloneAs();
                for (int d = 0; d < SpatialDim; d++)
                {
                    vt[d] = -v[d];
                }
                CalculateSupportPoint(p0, SpatialDim, vt, out ClosestPoint0);
                CalculateSupportPoint(p1, SpatialDim, v, out ClosestPoint1);
                for (int d = 0; d < SpatialDim; d++)
                {
                    SupportPoint[d] = ClosestPoint0[d] - ClosestPoint1[d];
                }
                double test = Math.Sqrt((v[0] - SupportPoint[0]).Pow2() + (v[1] - SupportPoint[1]).Pow2());
                if (Math.Sqrt((v[0] - SupportPoint[0]).Pow2() + (v[1] - SupportPoint[1]).Pow2()) <= 1e-12)
                    break;
                Simplex.Add(SupportPoint.CloneAs());
                DistanceAlgorithm(Simplex, out v);
            }
            Min_Distance = Math.Sqrt(v[0].Pow2() + v[1].Pow2());
        }
        private void Initialize_GJK(int SpatialDim, double[] Point0_old, double[] Point1_old, out double[] v0, out List<double[]> Simplex)
        {
            Simplex = new List<double[]>();
            v0 = new double[SpatialDim];
            for (int d = 0; d< SpatialDim; d++)
            {
                v0[d] = Point0_old[d] - Point1_old[d];
            }
            Simplex.Add(v0.CloneAs());
        }
        private void CalculateSupportPoint(Particle _Particle, int SpatialDim, double[] Vector, out double[] SupportPoint)
        {
            double VectorLength = Math.Sqrt(Vector[0].Pow2() + Vector[1].Pow2());
            double CosT = Vector[0] / VectorLength;
            double SinT = Vector[1] / VectorLength;
            _Particle.GetSupportPoint(SpatialDim, CosT, SinT, out SupportPoint);
        }

        private void DistanceAlgorithm(List<double[]> Simplex, out double[] v)
        {
            v = new double[2];
            
            List<double[]> DotProd_Simplex = new List<double[]>();
            for (int s1 = 0; s1 < Simplex.Count(); s1++)
            {
                DotProd_Simplex.Add(new double[Simplex.Count()]);
                for (int s2 = s1; s2 < Simplex.Count(); s2++)
                {
                    DotProd_Simplex[s1][s2] = Simplex[s1][0] * Simplex[s2][0] + Simplex[s1][1] * Simplex[s2][1];
                }
            }
            if(Simplex.Count() == 1)
            {
                v = Simplex[0];
            }
            else if (Simplex.Count() == 2)
            {
                if (DotProd_Simplex[0][0] - DotProd_Simplex[0][1] <= 0)
                {
                    v = Simplex[0].CloneAs();
                    Simplex.RemoveAt(1);
                }

                else if (DotProd_Simplex[1][1] - DotProd_Simplex[0][1] <= 0)
                {
                    v = Simplex[1].CloneAs();
                    Simplex.RemoveAt(0);
                }
                else
                {
                    double a0 = (Simplex[1][1] - Simplex[0][1]) / (Simplex[1][0] - Simplex[0][0]);
                    double b = Simplex[1][1] - Simplex[1][0] * a0;
                    double a1 = -(a0 + 1 / a0);
                    v[0] = b / a1;
                    v[1] = -v[0] / a0;
                }
            }
            else if (Simplex.Count() == 3)
            {
                bool Return = false;
                for (int s1 = 0; s1 < Simplex.Count(); s1++)
                {
                    int s2 = s1 == 2 ? 2 : 1;
                    int s3 = s1 == 0 ? 0 : 1;
                    double test1 = DotProd_Simplex[s1][s1] - DotProd_Simplex[0][s2];
                    double test2 = DotProd_Simplex[s1][s1] - DotProd_Simplex[s3][2];
                    if (DotProd_Simplex[s1][s1] - DotProd_Simplex[0][s2] <= 0 && DotProd_Simplex[s1][s1] - DotProd_Simplex[s3][2] <= 0)
                    {
                        v = Simplex[s1].CloneAs();
                        Simplex.Clear();
                        Simplex.Add(v.CloneAs());
                        Return = true;
                    }
                }
                if (!Return)
                {
                    for (int s1 = Simplex.Count() - 1; s1 >= 0; s1--)
                    {
                        int s2 = s1 == 2 ? 1 : 2;
                        int s3 = s1 == 0 ? 2 : 0;
                        int s4 = s1 == 0 ? 1 : 0;
                        double temp1 = (Simplex[s2][0] - Simplex[s1][0]) * (Simplex[s2][1] - Simplex[s4][1]);
                        double temp2 = (Simplex[s2][1] - Simplex[s1][1]) * (Simplex[s2][0] - Simplex[s4][0]);
                        double CrossProd = Simplex[s2][0] * (-temp1 + temp2) * (Simplex[s2][1] - Simplex[s4][1]) + Simplex[s2][1] * (temp1 - temp2) * (Simplex[s2][0] - Simplex[s4][0]);
                        CrossProd *= 1;
                        double test1 = DotProd_Simplex[s4][s4] - DotProd_Simplex[s4][s2];
                        double test2 = DotProd_Simplex[s2][s2] - DotProd_Simplex[s4][s2];
                        if (DotProd_Simplex[s4][s4] - DotProd_Simplex[s4][s2] >= 0 && DotProd_Simplex[s2][s2] - DotProd_Simplex[s4][s2] >= 0 && CrossProd >= 0)
                        {
                            double a0 = (Simplex[s2][1] - Simplex[s4][1]) / (Simplex[s2][0] - Simplex[s4][0]);
                            double b = Simplex[s2][1] - Simplex[s2][0] * a0;
                            double a1 = -(a0 + 1 / a0);
                            v[0] = b / a1;
                            v[1] = -v[0] / a0;

                            double[] test = new double[2];
                            test[1] = Simplex[s2][1] - Simplex[s4][1];
                            test[0] = (Simplex[s2][0] - Simplex[s4][0]);
                            double[] test12 = new double[2];
                            test12[0] = -test[1];
                            test12[1] = test[0];
                            double[] tempSimplex1 = Simplex[s2].CloneAs();
                            double[] tempSimplex2 = Simplex[s4].CloneAs();
                            Simplex.Clear();
                            Simplex.Add(tempSimplex1.CloneAs());
                            Simplex.Add(tempSimplex2.CloneAs());
                            break;
                        }
                    }
                }
            }   
        }
    }
}
