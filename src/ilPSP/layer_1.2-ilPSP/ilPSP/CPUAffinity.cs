using ilPSP.Tracing;
using MPI.Wrappers;
using NUnit.Framework.Constraints;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;
using System.Text;

namespace ilPSP.Utils {
    /// <summary>
    /// Returns the CPU affinity, 
    /// supporting systems with more than 64 processors (unlike <see cref="Process.ProcessorAffinity"/>).
    /// </summary>
    public class CPUAffinity {


        /// <summary>
        /// Returns the list of CPU's to which the current process is assigned to;
        /// Driver which calls either the respective Linux or Windows API functions.
        /// </summary>
        public static IEnumerable<int> GetAffinity() {
            IEnumerable<int> ret;
            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                ret = CPUAffinityWindows.GetAffinity();
            } else if (System.Environment.OSVersion.Platform == PlatformID.Unix) {
                ret = CPUAffinityLinux.GetAffinity();
            } else {
                throw new NotSupportedException("Not implemented for system: " + System.Environment.OSVersion.Platform);
            }

            int NumProcs = TotalNumberOfCPUs;
            foreach (int iCPU in ret) {
                if (iCPU >= NumProcs) {
                    Console.WriteLine($"Got affinity to CPU #{iCPU}, but system reports at maximum {NumProcs} CPUs.");
                }
            }
            return ret;
        }

        static public int TotalNumberOfCPUs {
            get {
                if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                    return CPUAffinityWindows.TotalNumberOfCPUs;
                } else if (System.Environment.OSVersion.Platform == PlatformID.Unix) {
                    return CPUAffinityLinux.TotalNumberOfCPUs;
                } else {
                    return System.Environment.ProcessorCount;
                }
            }
        }


        /// <summary>
        /// Deals with a special kink of handing CPU indices to Intel OpenMP on Windows:
        /// - On Windows, E.g. a system with 96 CPUs is organized into two groups of 48 CPUs.
        /// - Normally, the CPUs in the first group are 0 to 47, second group is 48 to 95, etc.
        /// - OpenMP, however, starts each group at a multiple of 64: fist group is 0 to 47, second is 64 to 111, ...
        /// On other systems, there is nothing to do.
        /// </summary>
        public static IEnumerable<int> ToOpenMpCPUindices(IEnumerable<int> cpuList) {
            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                int NumberOfCPUsPerGroup = CPUAffinityWindows.NumberOfCPUsPerGroup;

                var ret = new List<int>();
                foreach(int iCPU in cpuList) {
                    int iGroup = iCPU / NumberOfCPUsPerGroup;
                    int iWithinGroup = iCPU % NumberOfCPUsPerGroup;

                    ret.Add(iGroup*64 + iWithinGroup); // 
                }

                return ret.ToArray();

            } else {
                return cpuList;
            }

        }







        /// <summary>
        /// Configuration of OpenMP Environment variable `OMP_PLACES` to a given CPU affinity.
        /// </summary>
        /// <param name="CPUlist">e.g., return value from <see cref="GetAffinity"/></param>
        /// <returns></returns>
        public static int SetOMP_PLACESFromCPUList(IEnumerable<int> CPUlist) {
            using (var tr = new FuncTrace()) {
                /*
                var GlobalCPUlist = CpuListOnSMP(CPUlist);

                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIrank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPIsize);

                int SMPsize = ilPSP.Environment.MPIEnv.ProcessesOnMySMP; // number of MPI ranks on compute node
                int SMPrank = ilPSP.Environment.MPIEnv.ProcessRankOnSMP;
                string OMP_PLACES;
                int MaxNumOMPThreads;
                if (CPUlist.SetEquals(GlobalCPUlist)) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // the same CPU list for all ranks on the SMP node
                    // give to each process its
                    // dedicated portion of CPUs
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    //
                    // (It seems, OpenMP is not smart enough to negotiate thread ownership on one SMP node;
                    // some locking was observed, i.e., it seems that multiple ranks on one SMP node try to grab the same CPU.)
                    //


                    var subGroup = CPUlist.ToArray().GetSubVector(SMPrank*iThreads, iThreads);
                    var sanSubGroup = SanitzeGroup(subGroup);


                    OMP_PLACES = sanSubGroup.ToConcatString("{", ",", "}");
                    MaxNumOMPThreads= sanSubGroup.Count();

                    /*

                    if (CPUlist.Count < iThreads*SMPsize) {
                        throw new NotSupportedException($"Less CPU's reserved by MS HPC ({CPUlist.Count}) than required; number of threads: ({SMPsize}*{iThreads} = {SMPsize*iThreads})");
                    }
                    if (allInOneGroup) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // all CPU's in one CPU group: all processors may use all CPU's
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                        OMP_PLACES = CPUlist.ToConcatString("{", ",", "}");
                        MaxNumOMPThreads = CPUlist.Count();

                    } else {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // CPU's span over different groups: give to each process its
                        // dedicated portion of CPUs
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        var subGroup = CPUlist.GetSubVector(SMPrank*iThreads, iThreads);
                        var sanSubGroup = SanitzeGroup(subGroup);


                        OMP_PLACES = sanSubGroup.ToConcatString("{", ",", "}");
                        MaxNumOMPThreads= sanSubGroup.Count();
                    }
                     * /
                } else {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // MS HPC gave us different groups for each process
                    // use the entire group for this process and hope that Windows is smart enough
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    if (CPUlist.Count() < iThreads) {
                        throw new NotSupportedException($"Less CPU's reserved by MS HPC ({CPUlist.Count()}) than required; number of threads: ({iThreads})");
                    }


                    var sanSubGroup = SanitzeGroup(CPUlist);
                    OMP_PLACES = sanSubGroup.ToConcatString("{", ",", "}");
                    MaxNumOMPThreads= sanSubGroup.Count();
                }
                */

                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIrank);
                int SMPrank = ilPSP.Environment.MPIEnv.ProcessRankOnSMP;

                var sanSubGroup = SanitzeGroup(ToOpenMpCPUindices(CPUlist));
                int MaxNumOMPThreads = sanSubGroup.Count();

                var OMP_PLACES = sanSubGroup.ToConcatString("{", ",", "}");
                tr.Info($"R{MPIrank}, SMP rank {SMPrank}: setting OMP_PLACES = {OMP_PLACES}, MaxNumOMPThreads = {MaxNumOMPThreads}");

                System.Environment.SetEnvironmentVariable("OMP_PLACES", OMP_PLACES);
                System.Environment.SetEnvironmentVariable("OMP_PROC_BIND", "spread");
                return MaxNumOMPThreads;
            }
        }

        /// <summary>
        /// Configuration of OpenMP Environment variable `OMP_PLACES` to a given CPU affinity.
        /// </summary>
        /// <param name="iThreads">number of threads on MPI rank</param>
        /// <param name="CPUlist">e.g., return value from <see cref="GetAffinity"/></param>
        /// <returns></returns>
        public static int SetKMP_AFFINITYFromCPUList(int iThreads, IEnumerable<int> CPUlist) {
            using (var tr = new FuncTrace()) {
                
                
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIrank);
                int SMPrank = ilPSP.Environment.MPIEnv.ProcessRankOnSMP;

                var sanSubGroup = ToOpenMpCPUindices(CPUlist);
                int MaxNumOMPThreads = sanSubGroup.Count();

                string KMP_AFFINITY_proclist = sanSubGroup.ToConcatString("[", ",", "]");
                    


                string KMP_AFFINITY = $"verbose,proclist={KMP_AFFINITY_proclist},explicit";

                tr.Info($"R{MPIrank}, SMP rank {SMPrank}: setting KMP_AFFINITY = {KMP_AFFINITY}, MaxNumOMPThreads = {MaxNumOMPThreads}");

                System.Environment.SetEnvironmentVariable("KMP_AFFINITY", KMP_AFFINITY);
                return MaxNumOMPThreads;
            }
        }

        /// <summary>
        /// Collects the CPU indices of all processes running on the same node.
        /// </summary>
        /// <param name="__CPUlist">
        /// list of reserved CPU's indices for the current MPI process
        /// </param>
        /// <param name="allequal">
        /// if true, all MPI ranks on the current Compute Node are assigned to the same CPUs.
        /// </param>
        /// <param name="allequal">
        /// if true, all MPI ranks on the current Compute Node have exclusive sets of CPUs.
        /// </param>        /// <returns>
        /// list of CPU indices for all MPI processes running on the current Symmetric-Multi-Threading (SMT) Node (aka. the current Computer)
        /// </returns>
        public static int[] CpuListOnSMP(IEnumerable<int> __CPUlist, out bool disjoint, out bool allequal) {
            int mySMPRank = ilPSP.Environment.MPIEnv.SMPrank;
            int SizeOnSMP = ilPSP.Environment.MPIEnv.ProcessesOnMySMP;

            int[] CPUList = __CPUlist.ToSet().ToArray();
            int[] SMPRank = new int[CPUList.Length]; 
            SMPRank.SetAll(mySMPRank);

            var GlobalCPUList = CPUList.MPIAllGather();
            var GlobalSMPRank = SMPRank.MPIAllGather();

            var r = new HashSet<int>();
            disjoint = true;
            for(int i = 0; i < GlobalCPUList.Length; i++) {
                if (GlobalSMPRank[i] == mySMPRank) { // only add if on this SMP (compute Node)


                    bool added = r.Add(GlobalCPUList[i]);

                    if(added == false)
                        // CPU was already added from some other rank;
                        disjoint = false;
                }
            }

            allequal = GlobalCPUList.SetEquals(CPUList);

            if(SizeOnSMP == 1) {
                Debug.Assert(allequal == true);
                Debug.Assert(disjoint== true);
                disjoint = false;
            }


            var rr = r.ToArray();
            Array.Sort(rr);
            return rr;
        }


        /// <summary>
        /// make sure the selected cores for OMP_PLACES are in one CPU group; 
        /// it does not seem to work to specify OMP_PLACES across different CPU groups
        /// </summary>
        static IEnumerable<int> SanitzeGroup(IEnumerable<int> CPUlist) {
            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                var CPUsPerGroup = new Dictionary<int, List<int>>();
                foreach (int iCPU in CPUlist) {
                    int iGroup = iCPU / 64;

                    if (!CPUsPerGroup.TryGetValue(iGroup, out var CPUpg)) {
                        CPUpg = new List<int>();
                        CPUsPerGroup.Add(iGroup, CPUpg);
                    }

                    CPUpg.Add(iCPU);
                }

                // select the largest possible group for OpenMP
                int[] ret = new int[0];
                foreach (var kv in CPUsPerGroup) {
                    if (kv.Value.Count > ret.Length) {
                        kv.Value.Sort();
                        ret = kv.Value.ToArray();
                    }
                }

                return ret;
            } else {
                return CPUlist;
            }
        }


    }
}
