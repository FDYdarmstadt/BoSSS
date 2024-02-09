using System;
using System.Collections.Generic;
using System.Text;

namespace ilPSP {
    public enum Parallelism {
        /// <summary>
        /// This wrapper will use the sequential version of the third party library
        /// </summary>
        SEQ,
        /// <summary>
        /// This wrapper will use the MPI parallel version of the third party library
        /// </summary>
        MPI,
        /// <summary>
        /// This wrapper will use the OMP parallel version of the third party library
        /// </summary>
        OMP
    }
}
