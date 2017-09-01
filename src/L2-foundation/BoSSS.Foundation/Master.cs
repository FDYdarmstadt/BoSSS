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
using System.IO;
using System.Net;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP.Tracing;

namespace BoSSS.Foundation.Comm {
    
    /// <summary>
    /// Communication "DatabaseDriver"
    /// </summary>
    public class Master {

        /// <summary>
        /// the MPI communicator on which the BoSSS context, and the grid, ... 'live' on. Currently,
        /// this is always <see cref="IMPI_CommConstants.WORLD"/>.
        /// </summary>
        MPI_Comm Comm {
            get { return csMPI.Raw._COMM.WORLD; }
        }

        /// <summary>
        /// 
        /// </summary>
        internal Master() {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out m_Size);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out m_MyRank);

            m_hostname = Dns.GetHostName();
            IPHostEntry e = Dns.GetHostEntry(m_hostname);
            if ( e != null) {
                string hostname = e.HostName;
                if (hostname != null)
                    m_hostname = hostname;
            }  
            m_hostname = m_hostname.ToLowerInvariant();
                        
            //Console.WriteLine("BoSSS: Hello from " + (m_MyRank + 1) + " of " + m_Size + " on host '" + m_hostname + "'");

            SMPEvaluation();
        }

        /// <summary>
        /// detects how the MPI nodes are distributed over compute nodes (SMP nodes)
        /// </summary>
        private void SMPEvaluation() {
            //int ht = m_Context.IOMaster.tracer.EnterFunction("BoSSS.Foundation.Comm.DatabaseDriver.SMPEvaluation");
            using (new FuncTrace()) {
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                // define SMP rank;
                // for each MPI process, the SMP node index
                // index: MPI rank; content: SMP rank;
                int[] SMPRank = null;
                int NoOfSMPs = -1;
                {
                    // we are using the computer name to determine
                    // which MPI processes run on the same physical machine


                    // send host name to proc 0.
                    SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
                    if (MyRank > 0)
                        sms.SetCommPath(0);
                    sms.CommitCommPaths();

                    if (MyRank > 0)
                        sms.Transmitt(0, m_hostname);

                    int recvRnk;
                    string nmn;
                    sms.GetNext(out recvRnk, out nmn);
                    if (MyRank == 0) {
                        // receiving names form all processors

                        List<string> hosts_unique = new List<string>();
                        hosts_unique.Add(m_hostname);

                        string[] hosts = new string[Size];
                        SMPRank = new int[Size];
                        ArrayTools.SetAll(SMPRank, int.MinValue);
                        hosts[0] = m_hostname;
                        SMPRank[0] = hosts_unique.IndexOf(m_hostname);

                        while (nmn != null) {

                            if (hosts[recvRnk] != null)
                                throw new ApplicationException("should not happen.");
                            hosts[recvRnk] = nmn;

                            int smpRnk = hosts_unique.IndexOf(nmn);
                            if (smpRnk < 0) {
                                hosts_unique.Add(nmn);
                                smpRnk = hosts_unique.Count - 1;
                            }

                            SMPRank[recvRnk] = smpRnk;

                            sms.GetNext(out recvRnk, out nmn);
                        }
                        NoOfSMPs = hosts_unique.Count;

                        for (int i = 0; i < Size; i++) {
                            if (hosts[i] == null || SMPRank[i] < 0)
                                throw new ApplicationException("fatal error in algorithm.");
                        }

                    } else {
                        // don't receive anything

                        if (nmn != null)
                            // fatal error in algorithm
                            throw new ApplicationException("ha?");
                    }

                    sms.Dispose();
                }

                m_SMPSize = NoOfSMPs.MPIBroadcast(0, csMPI.Raw._COMM.WORLD);
                m_SMPRanks = SMPRank.MPIBroadcast(0, csMPI.Raw._COMM.WORLD);


                {
                    // number of MPI processes per SMP rank; index: SMP rank
                    m_MPIProcessesPerSMP = new int[m_SMPSize];
                    int[] _MPIProcessesPerSMP = new int[m_SMPSize];
                    _MPIProcessesPerSMP[m_SMPRanks[m_MyRank]]++;

                    unsafe {
                        fixed (int* pSnd = &_MPIProcessesPerSMP[0], pRcv = &m_MPIProcessesPerSMP[0]) {
                            csMPI.Raw.Allreduce((IntPtr)pSnd, (IntPtr)pRcv, m_SMPSize, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                        }
                    }
                }

                //m_Context.IOMaster.tracer.LeaveFunction(ht);
            }
        }

        string m_hostname;

        /// <summary>
        /// name of the compute node of the actual MPI process;
        /// </summary>
        public string Hostname { get { return m_hostname; } }





        /// <summary>
        /// rank of current process
        /// </summary>
        int m_MyRank;

        /// <summary>
        /// number of processors
        /// </summary>
        int m_Size;

        /// <summary>
        /// rank of current MPI process
        /// </summary>
        public int MyRank { get { return m_MyRank; } }
       
        /// <summary>
        /// number of MPI processes
        /// </summary>
        public int Size { get { return m_Size; } }

        /// <summary>
        /// <see cref="SMPSize"/>
        /// </summary>
        int m_SMPSize;

        /// <summary>
        /// <see cref="SMPRank"/>
        /// </summary>
        int[] m_SMPRanks;
        
        /// <summary>
        /// number of compute nodes (Symmetric MultiProcessing - nodes) over which 
        /// the MPI processes are distributed
        /// </summary>
        public int SMPSize { get { return m_SMPSize; }}

        /// <summary>
        /// for each MPI rank, a number that identifies the compute node (Symmetric MultiProcessing - node)
        /// on which the MPI process is allocated
        /// </summary>
        /// <param name="rank">
        /// an MPI rank, within the current communicator
        /// </param>
        public int SMPRank(int rank) {
            return m_SMPRanks[rank];
        }

        /// <summary>
        /// see <see cref="MPIProcessesPerSMP"/>;
        /// </summary>
        int[] m_MPIProcessesPerSMP;

        /// <summary>
        /// returns the number of MPI processes which are assigned to the compute node (Symmetric MultiProcessing - node)
        /// <paramref name="_SMPrank"/>.
        /// </summary>
        /// <param name="_SMPrank"></param>
        /// <returns></returns>
        public int MPIProcessesPerSMP(int _SMPrank) {
            return m_MPIProcessesPerSMP[_SMPrank];
        }
    }
}
