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
        /// Selects the MPI-enabled build for an MPI-initialized BoSSS process.
        /// </summary>
        Auto = 0,

        // Seq = 1, // we don't have this :)

        /// <summary>
        /// Shared-memory OpenMP build using PETSc MPIUNI. This is only valid
        /// when no real MPI implementation has been initialized.
        /// </summary>
        OpenMp = 2,

        /// <summary>
        /// Distributed-memory OpenMPI build. This also supports one rank.
        /// </summary>
        Mpi = 3,

        /// <summary>
        /// Hybrid OpenMPI and OpenMP build.
        /// </summary>
        MpiOpenMp = 4        
    }

    /// <summary>Resolves the PETSc shared-library variant for a BoSSS execution.</summary>
    public static class PETScNativeVariantResolver {

        /// <summary>
        /// Returns the explicitly requested variant or resolves <see cref="PETScNativeVariant.Auto"/> for a process
        /// without initialized MPI. Use the MPI-aware overload for BoSSS jobs.
        /// </summary>
        public static PETScNativeVariant Resolve(PETScNativeVariant requestedVariant, int mpiSize, int numberOfThreads) {
            return Resolve(requestedVariant, mpiSize, numberOfThreads, mpiIsInitialized: false);
        }

        /// <summary>
        /// Returns the explicitly requested variant or resolves <see cref="PETScNativeVariant.Auto"/>.
        /// An initialized real-MPI process must never load PETSc MPIUNI, because its MPI symbols and communicator
        /// constants conflict with the process MPI implementation even for one rank.
        /// </summary>
        public static PETScNativeVariant Resolve(PETScNativeVariant requestedVariant, int mpiSize, int numberOfThreads, bool mpiIsInitialized) {
            if(mpiSize < 1)
                throw new System.ArgumentOutOfRangeException(nameof(mpiSize), mpiSize, "At least one MPI rank is required.");
            if(numberOfThreads < 1)
                throw new System.ArgumentOutOfRangeException(nameof(numberOfThreads), numberOfThreads, "At least one thread is required.");

            if(requestedVariant == PETScNativeVariant.OpenMp && mpiIsInitialized)
                throw new System.NotSupportedException(
                    "The PETSc OpenMP-only runtime uses MPIUNI and cannot be loaded after real MPI is initialized. "
                    + "Select the MPI or MPI+OpenMP PETSc runtime.");
            if(requestedVariant != PETScNativeVariant.Auto)
                return requestedVariant;

            if(!mpiIsInitialized && mpiSize == 1)
                return PETScNativeVariant.OpenMp;
            return numberOfThreads == 1 ? PETScNativeVariant.Mpi : PETScNativeVariant.MpiOpenMp;
        }
    }
}
