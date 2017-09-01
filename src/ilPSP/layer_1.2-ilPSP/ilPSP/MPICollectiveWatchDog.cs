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
        /// This is a debug/testing method: it
        /// terminates the application if the method is not called by all MPI processes in 
        /// the MPI world communicator.
        /// </summary>
        [Conditional("DEBUG")]
        public static void Watch() {
            Watch(csMPI.Raw._COMM.WORLD);
        }

        static volatile int Fire = 100;
        static volatile int Fired = 100;

        /// <summary>
        /// This is a debug/testing method: it
        /// terminates the application if the method is not called by all MPI processes in  
        /// the MPI communicator <paramref name="comm"/>
        /// </summary>
        /// <param name="comm"></param>
        [Conditional("DEBUG")]
        [MethodImpl(MethodImplOptions.NoInlining)]
        public static void Watch(MPI_Comm comm) {
            if(!Tracing.Tracer.InstrumentationSwitch)
                // if tracing is off for performance reasons,
                // also this method should be turned of.
                return;

            // init
            // ====
            string m_functionName;
            DateTime entryTime = DateTime.Now;
            {
                StackFrame fr = new StackFrame(1, true);
                _MethodBase m = fr.GetMethod();
                m_functionName = m.DeclaringType.FullName + "." + m.Name;
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
                double WaitTime = 10.0;

                while (true) {
                    if (Fire == 0) {
                        Fired = 0;
                        if (stw.Elapsed.TotalSeconds > WaitTime)
                            Console.WriteLine("Returning from out-of-sync after " + stw.Elapsed.TotalSeconds + " seconds.");
                        return;
                    } else if (stw.Elapsed.TotalSeconds > WaitTime) {
                        Console.Error.WriteLine("WARNING: MPI out of sync on rank " + rank + " for more than " + WaitTime + " seconds.");
                        Console.Error.WriteLine("Method '" + m_functionName + "' was called at " + entryTime +  " on process " + rank + " but not on all other processes.");
                        Console.Error.WriteLine("Call stack:");
                        Console.Error.WriteLine(st.ToString());
                        //System.Environment.Exit(-666);
                        //Debugger.Launch();
                        WaitTime += 10.0;
                    }
                }
            }));
            wDog.Start();


            // wait for all processors to arrive here 
            // ======================================
            csMPI.Raw.Barrier(comm); // must be in the 'main' thread
                                 // MPI is not necessarily thread-safe!

            // ok
            Fire = 0; // quit watchdog
            while (Fired != 0) ;
            
        }
    }
}

