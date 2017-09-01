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
using System.Runtime.InteropServices;
using MPI.Wrappers;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {

    /// <summary>
    /// bindings to the HYYPRE_IJVector - Interface 
    /// </summary>
    class IJVector {

        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorCreate")]
        static extern int Create4(uint MPI_Comm, int jlower, int jupper, out T_IJVector vector);

        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorCreate")]
        static extern int Create8(ulong MPI_Comm, int jlower, int jupper, out T_IJVector vector);

        /// <summary>
        /// Create a vector object
        /// </summary>
        public static int Create(MPI_Comm MPI_Comm, int jlower, int jupper, out T_IJVector vector) {
            ulong _com8;
            uint _com4;
            // we need to convert the MPI comm in ilPSP (which is a FORTRAN MPI comm)
            // to a C-MPI comm: can be either 4 or 8 bytes!
            int sz = csMPI.Raw.MPI_Comm_f2c(MPI_Comm, out _com4, out _com8);
            switch (sz) {
                case 4: return Create4(_com4, jlower, jupper, out vector);
                case 8: return Create8(_com8, jlower, jupper, out vector);
                default: throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Destroy a vector object 
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorDestroy")]
        public static extern int Destroy(T_IJVector vector);

        /// <summary>
        /// Prepare a vector object for setting coefficient values
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorInitialize")]
        public static extern int Initialize(T_IJVector vector);

        /// <summary>
        /// (Optional) Sets the maximum number of elements that are expected to be set
        /// (or added) on other processors from this processor This routine can significantly
        /// improve the efficiency of matrix construction, and should always be
        /// utilized if possible 
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorSetMaxOffProcElmts")]
        public static extern int SetMaxOffProcElmts(T_IJVector vector, int max_off_proc_elmts);

        /// <summary>
        /// Sets values in vector
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorSetValues")]
        public static extern int SetValues(T_IJVector vector, int nvalues, int[] indices, double[] values);

        /// <summary>
        /// Adds to values in vector
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorAddToValues")]
        public static extern int AddToValues(T_IJVector vector, int nvalues, int[] indices, double[] values);

        /// <summary>
        /// Finalize the construction of the vector before using
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorAssemble")]
        public static extern int Assemble(T_IJVector vector);

        /// <summary>
        /// Gets values in vector
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorGetValues")]
        public static extern int GetValues(T_IJVector vector, int nvalues, int[] indices, double[] values);

        /// <summary>
        /// Set the storage type of the vector object to be constructed
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorSetObjectType")]
        public static extern int SetObjectType(T_IJVector vector, int type);

        /// <summary>
        /// Get the storage type of the constructed vector object
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorGetObjectType")]
        public static extern int GetObjectType(T_IJVector vector, out int type);

        /// <summary>
        /// Returns range of the part of the vector owned by this processor
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorGetLocalRange")]
        public static extern int GetLocalRange(T_IJVector vector, out int jlower, out int jupper);

        /// <summary>
        /// Get a reference to the constructed vector object 
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJVectorGetObject")]
        public static extern int GetObject(T_IJVector vector, out T_ParCRS_vector mtx_object);
    }
}
