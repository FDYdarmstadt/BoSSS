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
using System.Runtime.InteropServices;
using MPI.Wrappers;

namespace ilPSP.LinSolvers.HYPRE.Wrappers
{
    /// <summary>
    /// Create a solver object
    /// </summary>
    class ParCSRGMRES 
    {
        
        /// <summary>
        /// Create a solver object
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRGMRESCreate")]
        public static extern int CreateGMRES(MPI_Comm MPI_Comm, out T_Solver solver);

        /// <summary>
        /// Destroy a solver object
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRGMRESDestroy")]
        public static extern int Destroy(T_Solver solver);
    }

}
