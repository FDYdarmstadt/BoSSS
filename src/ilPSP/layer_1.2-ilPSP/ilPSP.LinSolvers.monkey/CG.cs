/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using System;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;
using MPI.Wrappers;
using System.Diagnostics;

namespace ilPSP.LinSolvers.monkey {

    /*
    /// <summary>
    /// only for debug purposes - compares the result of two <see cref="Device"/>'s
    /// </summary>
    public class DebugSolver : Solver {

        public DebugSolver(Device dev, Device __ComparisonDev)
            : base(dev) {
                ComparisonDev = __ComparisonDev;
        }

        Device ComparisonDev;

        MatrixBase comp_Matrix;

        public override void DefineMatrix(MsrMatrix M) {
            comp_Matrix = ComparisonDev.CreateMatrix(M, M.RowPartiton);
            base.DefineMatrix(M);
        }

        protected override void CallSolver(VectorBase x, VectorBase rhs, ref SolverResult stats) {

            VectorBase comp_x = ComparisonDev.CreateVector(x.Part);
            VectorBase comp_rhs = ComparisonDev.CreateVector(rhs.Part);

            Random rnd = new Random(0);
            for (int i = 0; i < comp_x.Part.LocalLength; i++) {
                double rndval = rnd.NextDouble();
                comp_x[i] = rndval;
                x[i] = rndval;
            }

            for (int i = 0; i < comp_rhs.Part.LocalLength; i++) {
                double rndval = rnd.NextDouble();
                //comp_rhs[i] = rndval;
                //rhs[i] = rndval;
            }


            VectorBase.CommVector comp_x_comm = comp_x.CreateCommVector(comp_Matrix);
            VectorBase.CommVector x_comm = x.CreateCommVector(m_Matrix);

            x.Lock();
            rhs.Lock();
            m_Matrix.Lock();
            comp_x.Lock();
            comp_rhs.Lock();
            comp_Matrix.Lock();


            m_Matrix.SpMV_Expert(1.0, x_comm, 0.0, rhs);
            comp_Matrix.SpMV_Expert(1.0, comp_x_comm, 0.0, comp_rhs);

            x.Unlock();
            rhs.Unlock();
            m_Matrix.Unlock();
            comp_x.Unlock();
            comp_rhs.Unlock();
            comp_Matrix.Unlock();

            //int rank = Enviroment.MPIEnv.MPI_Rank;
            //if (rank == 0)
            //    Debugger.Break();


            // compare receive buffers
            // =======================
            {
                SortedList<int, double[]> recv_x = x_comm.RecvBuffers;
                SortedList<int, double[]> recv_x_comp = x_comm.RecvBuffers;

                // test wheter the processor lists are equal
                if (recv_x.Keys.Count != recv_x_comp.Count)
                    throw new ApplicationException();
                foreach (int p in recv_x.Keys)
                    if (!recv_x_comp.Keys.Contains(p))
                        throw new ApplicationException();
                foreach (int p in recv_x_comp.Keys)
                    if (!recv_x.Keys.Contains(p))
                        throw new ApplicationException();
                
                // compare the lists
                double L1diff = 0;
                foreach (int proc in recv_x.Keys) {
                    double[] rcv_a = recv_x[proc];
                    double[] rcv_b = recv_x_comp[proc];

                    if (rcv_a.Length != rcv_b.Length)
                        throw new ApplicationException();

                    for (int i = 0; i < rcv_a.Length; i++)
                        L1diff += Math.Abs(rcv_a[i] - rcv_b[i]);
                }
                if (L1diff != 0)
                    Console.WriteLine("L1 err of rcv buffer: " + L1diff);
            }

            // compare results
            // ===============
            {
                double L1errX = 0;
                double L1err = 0;
                for (int i = 0; i < comp_x.Part.LocalLength; i++) {
                    L1errX += Math.Abs(comp_x[i] - x[i]);

                    double cval = rhs[i];
                    //double cval = ((CUDA.CudaMatrix)m_Matrix).ext_result[i];
                    L1err += Math.Abs(comp_rhs[i] - cval);
                }
                double L1errglob = 0;
                unsafe {
                    csMPI.Raw.Allreduce(&L1err, &L1errglob, 1, MPI_Datatype.DOUBLE, MPI_Op.MPI_SUM, csMPI.Raw.MPI_COMM_WORLD);
                }
                //Console.WriteLine("local difference in L1 - norm: " + L1err + "\t MPI pro: " + Enviroment.MPIEnv.MPI_Rank);
                if (Enviroment.MPIEnv.MPI_Rank == 0) {
                    Console.WriteLine("Total difference in L1 - norm: " + L1errglob);
                }

                if (L1errX != 0.0)
                    throw new ApplicationException();
            }
        }
    }*/

    /*
    /// <summary>
    /// only for debug purposes - compares the result of two <see cref="Device"/>'s
    /// </summary>
    public class MThrSolver : Solver {
        [DllImport("Kernel32.dll")]
        static extern bool QueryPerformanceCounter(out long lpPerformanceCount);

        public MThrSolver(Device dev)
            : base(dev) {
        }

        protected override void CallSolver(VectorBase x, VectorBase rhs, ref SolverResult stats) {


            Random rnd = new Random(0);
            for (int i = 0; i < x.Part.LocalLength; i++) {
                double rndval = rnd.NextDouble();
                x[i] = rndval;
            }

            for (int i = 0; i < rhs.Part.LocalLength; i++) {
                double rndval = rnd.NextDouble();
                rhs[i] = rndval;
            }


            VectorBase.CommVector x_comm = x.CreateCommVector(m_Matrix);

            x.Lock();
            rhs.Lock();
            m_Matrix.Lock();

            long[] parRunTime = new long[4];

            for (int i = 0; i < parRunTime.Length; i++) {

                int numThr = i % parRunTime.Length;

                ilPSP.Threading.Paralleism.NumThreads = numThr + 1;
                //DateTime st = DateTime.Now;
                long st;
                QueryPerformanceCounter(out st );
                m_Matrix.SpMV_Expert(1.0, x_comm, 0.0, rhs);

                long en;
                QueryPerformanceCounter(out en);


                parRunTime[numThr] = en - st;
                Console.Write((numThr + 1));
                Console.Write(" threads: ");
                Console.Write(parRunTime[numThr]);
                Console.Write(" msec; ");

                if (numThr > 0) {
                    double speedup = 1.0 / ((double)parRunTime[numThr] / (double)parRunTime[0]);
                    Console.Write("speedup = ");
                    Console.Write(speedup);
                }

                Console.WriteLine();
                

            }

            

            x.Unlock();
            rhs.Unlock();
            m_Matrix.Unlock();


        }
    }
     * */
    
    /// <summary>
    /// a very simple conjugated gradient solver
    /// </summary>
    public class CG : Solver {

        //[DllImport("Kernel32.dll")]
        //static extern bool QueryPerformanceCounter(out long lpPerformanceCount);
        
        
        /// <summary>
        /// implementation of the CG algorithm
        /// </summary>
        /// <param name="x"></param>
        /// <param name="rhs"></param>
        /// <param name="stats"></param>
        protected override void CallSolver(VectorBase x, VectorBase rhs, ref SolverResult stats) {

            VectorBase P = Device.CreateVector(x.Part);
            VectorBase.CommVector commP = P.CreateCommVector(m_Matrix);

            VectorBase R = rhs; // rhs is only needed once, so we can use it to store residuals
            VectorBase V = Device.CreateVector(x.Part);

            // lock objects
            // ============
            x.Lock();
            P.Lock();
            R.Lock();
            V.Lock();

            m_Matrix.Lock();


            // compute P0, R0
            // ==============

            // we only need to multiply x once by the Matrix, so we don't want to create
            // a seperate VectorBase.CommVector - object for x;
            // Instead, we're temporatily exchangeing the roles of x and P;
            P.Swap(x); 
            x.CopyFrom(rhs);                               // x = rhs
            m_Matrix.SpMV_Expert(-1.0, commP, 1.0, x);     // x = rhs - M*x
            P.Swap(x);                                     
            R.CopyFrom(P);

            double alpha = R.TwoNormSquare();
            double alpha_0 = alpha;
            double ResNorm;

            if (m_ConvergenceType == ConvergenceTypes.Absolute)
                ResNorm = Math.Sqrt(alpha);
            else if (m_ConvergenceType == ConvergenceTypes.Relative)
                ResNorm = 1.0;
            else {
                throw new NotImplementedException("unknown convergence type: " + m_ConvergenceType.ToString());
            }

            //long total = 0;
            //long gemv = 0;
            //long rest = 0;
            //long st, en;

            // iterate
            // =======
            stats.Converged = false;
            stats.NoOfIterations = 1; // one iteration has already been performed (P0, R0)
            for (int n = m_MaxIterations - 2; n >= 0; n--) {

                if (ResNorm <= m_Tolerance && stats.NoOfIterations >= base.m_MinIterations) {
                    stats.Converged = true;
                    break;
                }

                if (Math.Abs(alpha) <= double.Epsilon)
                    // numerical breakdown
                    break;

                m_Matrix.SpMV_Expert(1.0, commP, 0, V);
                double lambda = alpha / V.InnerProd(P);


                x.Acc(lambda, P);
                
                R.Acc(-lambda, V);
                
                double alpha_neu = R.TwoNormSquare();

                // compute residual norm
                if (m_ConvergenceType == ConvergenceTypes.Absolute)
                    ResNorm = Math.Sqrt(alpha);
                else
                    ResNorm = Math.Sqrt(alpha/alpha_0);

                
                P.Scale(alpha_neu / alpha);
                P.Acc(1.0, R);
                                
                alpha = alpha_neu;
                stats.NoOfIterations++;
                //QueryPerformanceCounter(out st);
                //rest += (st - en);           
            }

            //Console.WriteLine("CG: R" + stats.NoOfIterations + " = " + ResNorm);
            
            // unlock objects
            // ==============
            x.Unlock();
            P.Unlock();
            R.Unlock();
            V.Unlock();

            m_Matrix.Unlock();
            
            commP.Dispose();
            P.Dispose();
            V.Dispose();
        }
    }
}
