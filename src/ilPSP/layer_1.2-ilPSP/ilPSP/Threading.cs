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
using System.Threading;
using System.Runtime.InteropServices;

namespace ilPSP.Threading {


    /// <summary>
    /// a general parallel section, 
    /// used by <see cref="Paralleism.Run"/>;
    /// </summary>
    /// <param name="ThreadRank">current thread rank</param>
    /// <param name="NoOfThreads">number of running threads</param>
    public delegate void RunParallel(int ThreadRank, int NoOfThreads);
    
    /// <summary>
    /// a parallel FOR-loop, 
    /// used by <see cref="Paralleism.For"/>;
    /// </summary>
    /// <param name="i0">start index, inclusive</param>
    /// <param name="iE">end index, exclusive</param>
    public delegate void ForParallel(int i0, int iE);

    /// <summary>
    /// a parallel FOR-loop, 
    /// used by <see cref="Paralleism.ReduceFor"/>;
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="i0">start index, inclusive</param>
    /// <param name="iE">end index, exclusive</param>
    /// <returns>
    /// the result of the threac
    /// </returns>
    public delegate T ReduceForParallel<T>(int i0, int iE) where T : struct;

    /// <summary>
    /// implements the reduction of a per-thread-result into a global result
    /// (e.g. summing up over all threads), 
    /// used by <see cref="Paralleism.ReduceFor"/>;
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="TotalRes"></param>
    /// <param name="ThreadRes"></param>
    public delegate void ReduceOp<T>(ref T TotalRes, ref T ThreadRes) where T : struct;
        
    /// <summary>
    /// Helpers to realize quasi-OpenMP - parallelism;<br/>
    /// This class should substitute the .NET 4 class 'System.Threading.Parallel' - class which currently (20sep10) isnt available in Mono.
    /// </summary>
    public static class Paralleism {
        
        static int m_NumThreads = -1;

        /// <summary>
        /// the currently predefined number of threads;
        /// The default value is 'Number-of-Procesors' over 'Number-Of-MPI-Processes on current machine';
        /// </summary>
        static public int NumThreads {
            get {
                if (m_NumThreads < 0) {
                    m_NumThreads = System.Environment.ProcessorCount / ilPSP.Environment.MPIEnv.ProcessesOnMySMP;
                }
                return m_NumThreads;
            }
            set {
                if (value < 1) {
                    throw new ArgumentException("number of threads must be greater or equal to 1.");
                }
                m_NumThreads = value;
            }
        }


        static RunParallel m_CurrentWorker = null;

        static ManualResetEvent[] m_sync;

        static Exception[] m_exc;

        /// <summary>
        /// A parallel FOR-loop;
        /// runs <paramref name="operation"/> in parallel, with <see cref="NumThreads"/> number of threads;
        /// </summary>
        /// <param name="operation"></param>
        /// <param name="i0">
        /// start index
        /// </param>
        /// <param name="Len">
        /// loop lenght
        /// </param>
        /// <remarks>
        /// if at least one thread throws an exception, an exception it hrown.
        /// </remarks>
        public static void For(int i0, int Len, ForParallel operation) {
            Run(delegate(int i, int N) {
                int thr_i0 = (Len * i) / N + i0;
                int thr_iE = (Len * (i + 1)) / N + i0;
                operation(thr_i0, thr_iE);
            });
        }

        /// <summary>
        /// A parallel FOR-loop;
        /// runs <paramref name="operation"/> in parallel, with <see cref="NumThreads"/> number of threads;
        /// </summary>
        /// <param name="operation"></param>
        /// <param name="ReduceOp"></param>
        /// <param name="i0">
        /// start index
        /// </param>
        /// <param name="Len">
        /// loop lenght
        /// </param>
        /// <remarks>
        /// if at least one thread throws an exception, an exception is thrown.
        /// </remarks>
        public static T ReduceFor<T>(int i0, int Len, ReduceForParallel<T> operation, ReduceOp<T> ReduceOp) 
            where T : struct {
            if (m_CurrentWorker != null) {
                throw new ApplicationException("allready in parallel section");
            }

            int N = NumThreads;
            T[] thr_res = new T[N];

            Run(delegate(int rnk, int Sz) {
                int thr_i0 = (Len * rnk) / Sz + i0;
                int thr_iE = (Len * (rnk + 1)) / Sz + i0;
                thr_res[rnk] = operation(thr_i0, thr_iE);
            });


            T retval = thr_res[0];// default(T);
            for (int i = 1; i < N; i++) {
                ReduceOp(ref retval, ref thr_res[i]);
            }
            return retval;
        }

        //[DllImport("Kernel32.dll")]
        //static extern bool QueryPerformanceCounter(out long lpPerformanceCount);


        /// <summary>
        /// A general parallel section;
        /// runs <paramref name="operation"/> in parallel, with <see cref="NumThreads"/> number of threads;
        /// </summary>
        /// <param name="operation">operation to perform in parallel</param>
        /// <remarks>
        /// if at least one thread throws an exception, an exception it hrown.
        /// </remarks>
        public static void Run(RunParallel operation) {
            //long start;
            //QueryPerformanceCounter(out start);


            if (m_CurrentWorker != null) {
                throw new ApplicationException("already in parallel section");
            }

            m_CurrentWorker = operation;
            int N = NumThreads;
            if (m_sync == null || m_sync.Length != (N - 1)) {
                m_sync = new ManualResetEvent[N - 1];
                m_exc = new Exception[N];
                //ThreadRuntime = new long[N];
            }
            for (int i = 0; i < N - 1; i++) {
                m_sync[i] = new ManualResetEvent(false);
                m_exc[i] = null;
            }

            //System.Threading.ThreadPool.SetMinThreads(N-1, 0);
            //System.Threading.ThreadPool.SetMaxThreads(N-1, 0);
            for (int i = 1; i < N; i++) {
                System.Threading.ThreadPool.QueueUserWorkItem(Worker, new int[] { i, N });
            }

            try {
                Worker(new int[] { 0, N }); // run the 0-th pice of work in the current thread !
            } catch (Exception _e) {
                m_exc[0] = _e;
            }
            //long sync_st;
            //QueryPerformanceCounter(out sync_st);
            if (N > 1)
                WaitHandle.WaitAll(m_sync);
            //long sync_end;
            //QueryPerformanceCounter(out sync_end);

            m_CurrentWorker = null;
            for (int i = 0; i < N; i++) {
                if (m_exc[i] != null)
                    throw new ApplicationException("at least one thread (" + i + " of " + N + " throwed an exception;", m_exc[i]);
            }

            //long end;
            //QueryPerformanceCounter(out end);
            //{
            //    long runtime = end - start;
            //    long synctime = sync_end - sync_st;
            //    double perc = (double)synctime / (double)runtime;
            //    Console.WriteLine(N + " thr, total runtime: " + runtime + ", synctime: " + synctime + " = " + perc);

            //    /*for (int i = 0; i < N; i++) {
            //        Console.Write("thread #");
            //        Console.Write(i);
            //        Console.Write(": ");
            //        Console.Write(ThreadRuntime[i]);
            //        double perc = (double)ThreadRuntime[i] / (double)runtime;
            //        Console.Write(", = ");
            //        Console.Write(perc);
            //        Console.Write(" of total.");
            //        Console.WriteLine();
            //    }*/
            //}
        }

        //static long[] ThreadRuntime;

        static void Worker(object p) {
            
            //long st;
            //QueryPerformanceCounter(out st);

            int[] nums = (int[])p;
            int rnk = nums[0];
            int sz = nums[1];
            
            try {
                m_CurrentWorker(rnk, sz);
            } catch (Exception e) {
                m_exc[rnk] = e;
            }
            if( rnk > 0)
                m_sync[rnk-1].Set();

            //long en;
            //QueryPerformanceCounter(out en);

            //ThreadRuntime[rnk] = en - st;

        }
    }
}
