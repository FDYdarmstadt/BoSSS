﻿using log4net.Core;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;

namespace ilPSP.Utils {


    /// <summary>
    /// Returns the CPU affinity, 
    /// supporting systems with more than 64 processors (unlike <see cref="Process.ProcessorAffinity"/>).
    /// Therefore, windows processor groups have to be considered:
    ///
    /// - In Windows, CPUs (aka. processor cores, processors) are organized into Processor Groups.
    /// - each processor group contains at most 64 cores; thus, on any ordinary desktop, there is only one CPU group.
    ///   - according to Microsoft documentation, the number of CPU groups is as low as possible; 
    ///   - it also seems that the groups are symmetric, i.e., all groups contain the same amount of processors.
    ///     E.g., a system with 96 logical cores is organized into 2 groups with 48 cores for each group.
    ///   - it also seems that a singe process is always assigned to one CPU group.
    ///     At least, on Windows Server 2019, setting the affinity via the task manager seemingly only allows CPUs from one group.
    ///     (If one adds CPUs in e.g. group 1, all CPUs from the original group 0 get de-selected for the respective process.)
    /// - the C#-way of examining the processor affinity (<see cref="Process.ProcessorAffinity"/>) only returns the affinity within the current group;
    ///   there is no information whether we are dealing with group 0 or 1 or ...; thus, the CPU indices are only relative within the group
    /// - we must set `OMP_PLACES` environment variable to assign OpenMP-threads (PARDISO, BLAS, ...) to correct core,
    ///   but we need global CPU indices, i.e. 0-64 for group 0, 64-127 for group 1, etc.; 
    ///   therefore, we have to invoke the Win32-API to get the right processor group, to compute the CPU index offset.
    /// </summary>
    public static class CPUAffinityWindows {

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool GetProcessGroupAffinity(IntPtr hProcess, out ushort GroupCount, [Out] ushort[] GroupArray);


        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool GetThreadGroupAffinity(IntPtr hThread, out GROUP_AFFINITY lpGroupAffinity);


        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr GetCurrentThread();


        [StructLayout(LayoutKind.Sequential)]
        struct GROUP_AFFINITY {
            public UIntPtr Mask;
            public ushort Group;
            public ushort Reserved1;
            public ushort Reserved2;
            public ushort Reserved3;
        }

        public static void HelloGroup() {
            GROUP_AFFINITY _groupAffinity;
            if (GetThreadGroupAffinity(GetCurrentThread(), out _groupAffinity)) {
                Console.Error.WriteLine($"Group aff is {_groupAffinity.Group}, mask = {_groupAffinity.Mask:X}");
            } else {
                Console.Error.WriteLine($"Group aff is err");
            }
        }

        /// <summary>
        /// (Windows version) Returns the list of CPU's to which the current process is assigned to.
        /// </summary>
        public static IEnumerable<int> GetAffinity() {
            Process currentProcess = Process.GetCurrentProcess();
            IntPtr processHandle = currentProcess.Handle;
            
            ushort groupCount = 0;
            GetProcessGroupAffinity(processHandle, out groupCount, null);
            if (groupCount != 1) {
                //Console.WriteLine($"Process associated to more than one processor group ({groupCount}) -- i don't know what to do about it (tell Florian)!");
                //throw new NotSupportedException("Process associated to more than one processor group -- i don't know what to do about it (tell Florian)!");
            }
            
            
            ushort[] groups = new ushort[groupCount];
            if (!GetProcessGroupAffinity(processHandle, out groupCount, groups)) {
                Console.Error.WriteLine("Failed to get processor group affinity.");
                int errorCode = Marshal.GetLastWin32Error();
                throw new Win32Exception(errorCode);
            }

                //Console.WriteLine("Groups are " + groups.ToConcatString("", ", ", ""));

            ushort group = groups[0];

            /*
            ilPSP.Environment.ParallelFor(0, 1024, delegate (int i) {

                GROUP_AFFINITY _groupAffinity;
                if (GetThreadGroupAffinity(GetCurrentThread(), out _groupAffinity)) {
                Console.WriteLine($"Group aff in {i} = {_groupAffinity.Group}, mask = {_groupAffinity.Mask:X}");
                } else {
                    Console.WriteLine($"Group aff in {i} = err");
                }
            });
            */

            GROUP_AFFINITY groupAffinity;
            if (GetThreadGroupAffinity(GetCurrentThread(), out groupAffinity)) {
                return CheckCpuAffinity(groupAffinity.Mask, group, System.Environment.ProcessorCount / groupCount);
            } else {
                int errorCode = Marshal.GetLastWin32Error();
                throw new Win32Exception(errorCode);
            }
        }


        /// <summary>
        /// When MS HPC server is used, it seems to be necessary to define the `OMP_PLACES` environment variable.
        /// Otherwise, it seems that many OpenMP-threads (PARDISO, BLAS) seem to concentrate on the same cores.
        /// While C#-threads seem to respect the core affinity set by the HPC server, the OpenMP-threads don't care.
        /// 
        /// Luckily, MS HPC typically defines the Environment variable `CCP_AFFINITY`;
        /// We are going to use this 
        /// </summary>
        public static int SetOMP_PLACESfromCCPVar(int NumThreads) {

            // Check if CCP_AFFINITY is defined and defined on all ranks
            // =========================================================

            string CCP_AFFINITY = System.Environment.GetEnvironmentVariable("CCP_AFFINITY");
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIrank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPIsize);

            bool CCP_AFFINITY_DEFINED = (CCP_AFFINITY != null);
            bool glCCP_AFFINITY_DEFINED = CCP_AFFINITY_DEFINED.MPIOr();

            if (glCCP_AFFINITY_DEFINED != CCP_AFFINITY_DEFINED) {
                string errMsg = $"`CCP_AFFINITY` defined on some ranks, but not on all; defined on {MPIrank}? {CCP_AFFINITY_DEFINED}, globally? {glCCP_AFFINITY_DEFINED}";
                Console.Error.WriteLine(errMsg);
                throw new ApplicationException(errMsg);
            }

            if (glCCP_AFFINITY_DEFINED == false)
                // make all processors on system available for OpenMP
                return NumThreads;

            Console.WriteLine($"R{MPIrank}, rank on node {ilPSP.Environment.MPIEnv.ProcessRankOnSMP}: CCP_AFFINITY = {CCP_AFFINITY}");


            // decode the variable
            // ===================


            var affGroup = CCP_AFFINITY.Split(new string[] { "," }, StringSplitOptions.RemoveEmptyEntries);

            var CPUlist = new List<int>();
            int iGroup = 0;
            var groupOccupied = new List<bool>();
            foreach (string aff in affGroup) {
                //
                // note: at least in our MKL version, it seems that the indices for OMP_PLACES always start at 0 for group 0 and 64 for group 1; Even if the system has e.g. 48 processors per group.
                //

                var groupCPUs = CheckCpuAffinity(new UIntPtr(Convert.ToUInt64(aff, 16)), iGroup, 64); // always 64 procs per group!!!
                groupOccupied.Add(groupCPUs.Count() > 0);
                CPUlist.AddRange(groupCPUs);
                iGroup++;

            }
            CPUlist.Sort();

            //bool allInOneGroup = groupOccupied.Where(bg => bg).Count() == 1;

            
            return CPUAffinity.SetOMP_PLACESFromCPUList(NumThreads, CPUlist.ToArray());
        }


      

        static IEnumerable<int> CheckCpuAffinity(UIntPtr mask, int iProcessorGroup, int procsPerGroup) {
            var res = new List<int>();
            ulong bitmask = (ulong)mask;
            for (int cpu = 0; cpu < 64; cpu++)  // Assuming a maximum of 64 CPUs per group
            {
                ulong cpuBit = 1UL << cpu;
                if ((bitmask & cpuBit) != 0) {
                    //Console.WriteLine($"  CPU {cpu + iProcessorGroup*64} is available in this group.");
                    res.Add(cpu + iProcessorGroup*procsPerGroup);
                }
            }
            return res.ToArray();
        }
    }
}

