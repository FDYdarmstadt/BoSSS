using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace ilPSP.Utils {


    /// <summary>
    /// Returns the CPU affinity, 
    /// supporting systems with more than 64 processors (unlike <see cref="Process.ProcessorAffinity"/>).
    /// </summary>
    public class CPUAffinityLinux {
        [DllImport("libc.so.6", SetLastError = true)]
        unsafe private static extern int sched_getaffinity(int pid, IntPtr cpusetsize, void* mask);



        /// <summary>
        /// (Linux version) Returns the list of CPU's to which the current process is assigned to.
        /// </summary>
        public static IEnumerable<int> GetAffinity() {

            const int CPU_SETSIZE = 1024;
            const int ULONG_SIZE = 8; // Size of ulong in bytes
            int arraySize = CPU_SETSIZE / (8 * ULONG_SIZE); // Calculate array size for 1024 bits

            //cpu_set_t mask = new cpu_set_t(arraySize);
            var bits = new ulong[arraySize];
            unsafe {
                fixed (ulong* pbits = bits) {
                    int result = sched_getaffinity(0, new IntPtr(arraySize * ULONG_SIZE), pbits);

                    if (result == -1) {
                        throw new ApplicationException("Error retrieving CPU affinity.");
                    }
                }
            }

            var ret = new List<int>();
            for (int i = 0; i < CPU_SETSIZE; i++) {
                int idx = i / (8 * ULONG_SIZE);
                int bit = i % (8 * ULONG_SIZE);

                if ((bits[idx] & (1UL << bit)) != 0) {
                    //Console.WriteLine($"CPU {i} is in the affinity set.");
                    ret.Add(i);
                }
            }
            return ret.ToArray();
        }


        // Define the sysconf call
        [DllImport("libc.so.6", SetLastError = true)]
        private static extern long sysconf(int name);

        // Constant for _SC_NPROCESSORS_ONLN from sysconf.h
        private const int _SC_NPROCESSORS_ONLN = 84;

        /// <summary>
        /// The total number of CPUs in a system; 
        /// This might be larger than the number reported from <see cref="System.Environment.ProcessorCount"/>,
        /// since SLURM seems to mess with this value
        /// </summary>
        public static int TotalNumberOfCPUs {
            get {
                try {
                    long processorCount = sysconf(_SC_NPROCESSORS_ONLN);
                    if (processorCount == -1) {
                        throw new InvalidOperationException("Failed to get processor count");
                    }

                    Console.WriteLine("Linux No Of CPUS: " + processorCount);
                    return (int)processorCount;
                } catch (Exception ex) {
                    Console.WriteLine("Error: " + ex.Message);
                    return System.Environment.ProcessorCount;
                }
            }
        }
    }
}
    
