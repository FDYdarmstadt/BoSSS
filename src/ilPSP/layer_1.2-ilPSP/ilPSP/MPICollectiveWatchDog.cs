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
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.Threading;
using MPI.Wrappers;
using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;

namespace ilPSP {

    /// <summary>
    /// A watchdog that should ensure proper use of MPI-collective methods.
    /// </summary>
    public static class MPICollectiveWatchDog {

        /// <summary>
        /// Use this for debugging only! This will impact the performance.
        /// For large parallel runs and if you are only interested to watch for specific methods.
        /// Will be fired at DEBUG and RELEASE alike.
        /// </summary>
        public static void WatchAtRelease(MPI_Comm comm, double WaitTimeSeconds = 10, bool quitAppIfExpired = false, int token = 0) {
            WatchInternal(comm, WaitTimeSeconds, quitAppIfExpired, token);
        }


        /// <summary>
        /// This is a debug/testing method: it
        /// terminates the application if the method is not called by all MPI processes in 
        /// the MPI world communicator.
        /// </summary>
        [Conditional("DEBUG")]
        public static void Watch(int token = 0) {
            Watch(csMPI.Raw._COMM.WORLD, token: token);
        }

        static volatile int Fire = 100;
        static volatile int Fired = 100;

        /// <summary>
        /// This is a debug/testing method: it
        /// terminates the application if the method is not called by all MPI processes in  
        /// the MPI communicator <paramref name="comm"/>
        /// </summary>
        [Conditional("DEBUG")]
        [MethodImpl(MethodImplOptions.NoInlining)]
        public static void Watch(MPI_Comm comm, double WaitTimeSeconds = 10, bool quitAppIfExpired = false, int token = 0) {
            WatchInternal(comm, WaitTimeSeconds, quitAppIfExpired, token);
        }

        private static void WatchInternal(MPI_Comm comm, double WaitTimeSeconds, bool quitAppIfExpired, int token) {
            if (!Tracing.Tracer.InstrumentationSwitch)
                // if tracing is off for performance reasons,
                // also this method should be turned of.
                return;

            // init
            // ====
            string m_functionName;
            DateTime entryTime = DateTime.Now;
            {
                // ==================================================================================================
                // find name of calling function: i.e., the first method on stack which is NOT a member of this class
                // ==================================================================================================
                m_functionName = "undetermined";
                for (int i = 0; i < 5000; i++) {
                    StackFrame fr = new StackFrame(i, true);
                    if (fr == null) {
                        m_functionName = $"undetermined (StackIndex = {i})";
                        break;
                    }
                    System.Reflection.MethodBase m = fr.GetMethod();
                    if (m == null) {
                        m_functionName = $"undetermined (StackIndex = {i})";
                        break;
                    }
                    string nn = m.DeclaringType.FullName + "." + m.Name;
                    if(m.DeclaringType != typeof(MPICollectiveWatchDog)) {
                        m_functionName = nn;
                        break;
                    }
                }
            }

            var st = new StackTrace();
            int rank;
            csMPI.Raw.Comm_Rank(comm, out rank);

            // start watchdog ...
            // ==================
            Fire = 100;
            Fired = 200;
            Thread wDog = (new Thread(delegate() {
                Stopwatch stw = new Stopwatch();
                stw.Start();


                while (true) {
                    if (Fire == 0) {
                        Fired = 0;
                        if (stw.Elapsed.TotalSeconds > WaitTimeSeconds)
                            Console.WriteLine("Returning from out-of-sync after " + stw.Elapsed.TotalSeconds + " seconds.");
                        return;
                    } else if (stw.Elapsed.TotalSeconds > WaitTimeSeconds) {
                        Console.Error.WriteLine("WARNING: MPI out of sync on rank " + rank + " for more than " + WaitTimeSeconds + " seconds.");
                        Console.Error.WriteLine("Method '" + m_functionName + "' was called at " + entryTime +  " on process " + rank + " but not on all other processes.");
                        Console.Error.WriteLine("Call stack:");
                        Console.Error.WriteLine(st.ToString());

                        if(quitAppIfExpired)
                            System.Environment.Exit(-666);
                        // dbg_launch();
                        WaitTimeSeconds += 10.0;
                    }
                }
            }));
            wDog.Start();


            // wait for all processors to arrive here 
            // ======================================
            int tokenMin = token.MPIMin();
            int tokenMax = token.MPIMax();
            if (tokenMin != tokenMax)
                throw new ApplicationException($"synchronization error: token minimum {tokenMin}, maximum {tokenMax}");
            csMPI.Raw.Barrier(comm); // must be in the 'main' thread
                                     // MPI is not necessarily thread-safe!

            bool b = ilPSP.Environment.StdoutOnlyOnRank0;
            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            //if (token != 0)
            //    Console.WriteLine($"rank {rank}: reached {token} in {m_functionName}");
            //ilPSP.Environment.StdoutOnlyOnRank0 = b;

            // ok
            Fire = 0; // quit watchdog
            while (Fired != 0) ;
            
        }
    }
}

