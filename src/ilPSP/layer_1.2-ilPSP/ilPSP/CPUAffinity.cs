using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace ilPSP.Utils {
    /// <summary>
    /// Returns the CPU affinity, 
    /// supporting systems with more than 64 processors (unlike <see cref="Process.ProcessorAffinity"/>).
    /// </summary>
    public class CPUAffinity {


        /// <summary>
        /// (Linux version) Returns the list of CPU's to which the current process is assigned to.
        /// </summary>
        public static IEnumerable<int> GetAffinity() {

            //Console.WriteLine("bootstrapping necessary.");
            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                return CPUAffinityWindows.GetAffinity();

            } else if (System.Environment.OSVersion.Platform == PlatformID.Unix) {
                return CPUAffinityLinux.GetAffinity();
            } else {
                throw new NotSupportedException("Not implemented for system: " + System.Environment.OSVersion.Platform);
            }

        }
    }
}
