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
using System.Net;
using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.Runtime.Serialization;
using MPI.Wrappers;
using ilPSP.Utils;
using System.Diagnostics;

namespace ilPSP {

    /// <summary>
    /// the issue of this class is to answer 
    /// </summary>
    public class MPIEnvironment {

        MPI_Comm m_mpi_comm;

        /// <summary>
        /// MPI communicator handle for this environment
        /// </summary>
        public MPI_Comm Mpi_comm {
            get { return m_mpi_comm; }
        }


        /// <summary>
        /// ctor
        /// </summary>
        public MPIEnvironment() {
            m_mpi_comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out m_MPISize);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out m_MPI_Rank);
            
            
            m_hostname = Dns.GetHostName();
            IPHostEntry e = Dns.GetHostEntry(m_hostname);
            if (e != null) {
                string hostname = e.HostName;
                if (hostname != null)
                    m_hostname = hostname;
            }
            m_hostname = m_hostname.ToLowerInvariant();

            m_HostnameForRank = m_hostname.MPIAllGatherO(csMPI.Raw._COMM.WORLD);

            {
                

                m_AllHostNames = new Dictionary<string, int[]>();
                for(int iRnk = 0; iRnk < MPI_Size; iRnk++) {
                    var hn = m_HostnameForRank[iRnk];
                    
                    int[] Ranks;
                    if(!m_AllHostNames.TryGetValue(hn, out Ranks)) {
                        Ranks = new int[0];
                        m_AllHostNames.Add(hn, Ranks);
                    }

                    iRnk.AddToArray(ref Ranks);
                    m_AllHostNames[hn] = Ranks;
                }
            }

            SMPEvaluation();
        }

        /// <summary>
        /// rank of current process
        /// </summary>
        int m_MPI_Rank;

        /// <summary>
        /// number of processors
        /// </summary>
        int m_MPISize;

        /// <summary>
        /// detects how the MPI nodes are distributed over compute nodes (SMP nodes)
        /// </summary>
        private void SMPEvaluation() {


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
                if (MPI_Rank > 0)
                    sms.SetCommPath(0);
                sms.CommitCommPaths();

                if (MPI_Rank > 0)
                    sms.Transmit(0, m_hostname);

                int recvRnk;
                string nmn;
                sms.GetNext(out recvRnk, out nmn);
                if (MPI_Rank == 0) {
                    // receiving names form all processors

                    List<string> hosts_unique = new List<string>();
                    hosts_unique.Add(m_hostname);

                    string[] hosts = new string[MPI_Size];
                    SMPRank = new int[MPI_Size];
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

                    for (int i = 0; i < MPI_Size; i++) {
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

            m_NoOfSMPs = MPIExtensions.MPIBroadcast(NoOfSMPs, 0, csMPI.Raw._COMM.WORLD);
            m_RankOfSMPs = MPIExtensions.MPIBroadcast(SMPRank, 0, csMPI.Raw._COMM.WORLD);


            {
                // number of MPI processes per SMP rank; index: SMP rank
                m_MPIProcessesPerSMP = new int[m_NoOfSMPs];
                int[] _MPIProcessesPerSMP = new int[m_NoOfSMPs];
                _MPIProcessesPerSMP[m_RankOfSMPs[m_MPI_Rank]]++;

                unsafe {
                    fixed (int* pSnd = &_MPIProcessesPerSMP[0], pRcv = &m_MPIProcessesPerSMP[0]) {
                        csMPI.Raw.Allreduce((IntPtr)pSnd, (IntPtr)pRcv, m_NoOfSMPs, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    }
                }
            }

            // define SMP-local MPI ranks
            // ==========================
            {
                m_ProcessRankOnSMP = 0;
                int MySMP = MPIrank_2_SMPrank(m_MPI_Rank);
                for (int r = 0; r < m_MPI_Rank; r++) {
                    if (m_RankOfSMPs[r] == MySMP)
                        m_ProcessRankOnSMP++;
                }
            }

            //m_Context.IOMaster.tracer.LeaveFunction(ht);
        }

        string m_hostname;

        /// <summary>
        /// name of the compute node of the actual MPI process;
        /// </summary>
        public string Hostname { get { return m_hostname; } }


        string[] m_HostnameForRank;

        /// <summary>
        /// Hostname for each MPI rank
        /// - index: MPI rank within the world communicator
        /// </summary>
        public IReadOnlyList<string> HostnameForRank {
            get {
                return m_HostnameForRank.CloneAs();
            }
        }

        Dictionary<string, int[]> m_AllHostNames;

        /// <summary>
        /// Hostnames and respective MPI ranks
        /// </summary>
        public IReadOnlyDictionary<string,int[]> AllHostNames {
            get {
                

                return m_AllHostNames;
            }
        }

        /// <summary>
        /// rank of current MPI process, within the world communicator
        /// </summary>
        public int MPI_Rank { get { return m_MPI_Rank; } }

        /// <summary>
        /// number of MPI processes, within the world communicator
        /// </summary>
        public int MPI_Size { get { return m_MPISize; } }

        /// <summary>
        /// <see cref="NoOfSMPs"/>
        /// </summary>
        int m_NoOfSMPs;

        /// <summary>
        /// <see cref="MPIrank_2_SMPrank"/>
        /// </summary>
        int[] m_RankOfSMPs;

        /// <summary>
        /// number of compute nodes (Symmetric MultiProcessing - nodes) over which 
        /// the MPI processes are distributed
        /// </summary>
        public int NoOfSMPs { get { return m_NoOfSMPs; } }

        /// <summary>
        /// the SMP rank of the current MPI process
        /// ("on which SMP am I running?");
        /// </summary>
        public int SMPrank { get { return MPIrank_2_SMPrank(MPI_Rank); } }

        /// <summary>
        /// for each MPI rank (within the world communicator), a number that identifies the compute node (Symmetric MultiProcessing - node)
        /// on which the MPI process is allocated
        /// </summary>
        /// <param name="MPI_rank">
        /// an MPI rank, within the current communicator
        /// </param>
        public int MPIrank_2_SMPrank(int MPI_rank) {
            return m_RankOfSMPs[MPI_rank];
        }

        /// <summary>
        /// see <see cref="MPIProcessesPerSMP"/>;<br/>
        /// index: SMP-node rank; <br/>
        /// content: number of MPI processes on corresponding SMP node;
        /// </summary>
        int[] m_MPIProcessesPerSMP;

        /// <summary>
        /// returns the number of MPI processes (in the world communicator) which are asigned to one compute node (Symmetric MultiProcessing - node)
        /// <paramref name="_RankOfSMP"/>.
        /// </summary>
        /// <param name="_RankOfSMP">
        /// rank of SMP-node (equal for all MPI processes on one compute node).
        /// </param>
        /// <returns></returns>
        public int MPIProcessesPerSMP(int _RankOfSMP) {
            return m_MPIProcessesPerSMP[_RankOfSMP];
        }

        /// <summary>
        /// total number of MPI processes on the actual Compute Node (Symmetric Multi Processing - node);
        /// </summary>
        public int ProcessesOnMySMP {
            get { return MPIProcessesPerSMP(MPIrank_2_SMPrank(this.MPI_Rank)); }
        }
        
        private int m_ProcessRankOnSMP;

        /// <summary>
        /// "SMP-local process-rank": a SMP-local rank, i.e. a number between 0 (including) and <see cref="ProcessesOnMySMP"/> (excluding); <br/>
        /// </summary>
        public int ProcessRankOnSMP {
            get { return m_ProcessRankOnSMP; }
        }

        

    }
}
