/* =======================================================================
Copyright 2026 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Selects the parallel runtime used by the dedicated PETSc native library.
    /// </summary>
    public enum PETScNativeVariant {

        /// <summary>
        /// Distributed-memory OpenMPI build. This is the default and also supports one rank.
        /// </summary>
        Mpi,

        /// <summary>
        /// Shared-memory OpenMP build using PETSc MPIUNI. BoSSS must run with one MPI rank.
        /// </summary>
        OpenMp,

        /// <summary>
        /// Hybrid OpenMPI and OpenMP build.
        /// </summary>
        MpiOpenMp
    }
}

