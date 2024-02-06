using System;
using System.Collections.Generic;
using System.Text;

namespace ilPSP {
    
    /// <summary>
    /// This attribute tells the test runner with how many threads (<see cref="ilPSP.Environment.InitThreading"/>) a test should be executed,
    /// **if submitted via the batch processor** 
    /// </summary>
    public class NUnitNumThreads : Attribute {

        /// <summary>
        /// Wäh
        /// </summary>
        public NUnitNumThreads(int __numthreads) {
            this.NumThreads = __numthreads;
        }

        /// <summary>
        /// number of threads for some test
        /// </summary>
        readonly public int NumThreads;
    }
}
