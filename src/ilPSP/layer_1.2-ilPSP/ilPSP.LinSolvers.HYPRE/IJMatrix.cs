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
    /// bindings to the HYPRE_IJMatrix - Interface 
    /// </summary>
    internal class IJMatrix {


        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixCreate")]
        extern static int Create4(uint MPI_Comm, int ilower, int iupper, int jlower, int jupper, out T_IJMatrix matrix);

        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixCreate")]
        extern static int Create8(ulong MPI_Comm, int ilower, int iupper, int jlower, int jupper, out T_IJMatrix matrix);

        /// <summary>
        /// creates a new IJ matrix object
        /// </summary>
        public static int Create(MPI_Comm MPI_Comm, int ilower, int iupper, int jlower, int jupper, out T_IJMatrix matrix) {
            ulong _com8;
            uint _com4;
            // we need to convert the MPI comm in ilPSP (which is a FORTRAN MPI comm)
            // to a C-MPI comm: can be either 4 or 8 bytes!
            int sz = csMPI.Raw.MPI_Comm_f2c(MPI_Comm, out _com4, out _com8);
            switch (sz) {
                case 4: return Create4(_com4, ilower, iupper, jlower, jupper, out matrix);
                case 8: return Create8(_com8, ilower, iupper, jlower, jupper, out matrix);
                default: throw new NotImplementedException();
            }
        }
        /// <summary>
        /// frees a IJ matrix object
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixDestroy")]
        extern public static int Destroy(T_IJMatrix matrix);

        /// <summary>
        /// Prepare a matrix object for setting coefficient values
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixInitialize")]
        extern public static int Initialize(T_IJMatrix matrix);


        /// <summary>
        /// Sets values for nrows rows or partial rows of the matrix
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixSetValues")]
        extern public static int SetValues(T_IJMatrix matrix, int nrows, int[] ncols, int[] rows, int[] cols, double[] values);

        /// <summary>
        /// Adds to values for nrows rows or partial rows of the matrix
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixAddToValues")]
        extern public static int AddToValues(T_IJMatrix matrix, int nrows, int[] ncols, int[] rows, int[] cols, double[] values);

        /// <summary>
        /// Finalize the construction of the matrix before using
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixAssemble")]
        extern public static int Assemble(T_IJMatrix matrix);

        /// <summary>
        /// Gets number of nonzeros elements for nrows rows specified in rows and
        /// returns them in ncols, which needs to be allocated by the user
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixGetRowCounts")]
        extern public static int GetRowCounts(T_IJMatrix matrix, int nrows, int[] rows, int[] ncols);

        /// <summary>
        /// Gets values for nrows rows or partial rows of the matrix
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixGetValues")]
        extern public static int GetValues(T_IJMatrix matrix, int nrows, int[] ncols, int[] rows, int[] cols, double[] values);
 
        /// <summary>
        /// Set the storage type of the matrix object to be constructed
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixSetObjectType")]
        extern public static int SetObjectType(T_IJMatrix matrix, int type);

        /// <summary>
        /// Get the storage type of the constructed matrix object;
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixGetObjectType")]
        extern public static int GetObjectType(T_IJMatrix matrix, out int type);

        /// <summary>
        /// Gets range of rows owned by this processor and range of column partitioning for this processor
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixGetLocalRange")]
        extern public static int GetLocalRange(T_IJMatrix matrix, ref int ilower, ref int iupper, ref int jlower, ref int jupper);

        /// <summary>
        /// Get a reference to the constructed matrix object;
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixGetObject")]
        extern public static int GetObject(T_IJMatrix matrix, out T_ParCSR_matrix mtx_object);

        /// <summary>
        /// (Optional) Set the max number of nonzeros to expect in each row
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixSetRowSizes")]
        extern public static int SetRowSizes(T_IJMatrix matrix, int[] sizes);
 
        /// <summary>
        /// (Optional) Set the max number of nonzeros to expect in each row of the
        /// diagonal and off-diagonal blocks 
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixSetDiagOffdSizes")]
        extern public static int SetDiagOffdSizes(T_IJMatrix matrix, int[] diag_sizes, int[] offdiag_sizes);

        /// <summary>
        /// (Optional) Sets the maximum number of elements that are expected to be set
        /// (or added) on other processors from this processor This routine can significantly
        /// improve the efficiency of matrix construction, and should always be
        /// utilized if possible;
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_IJMatrixSetMaxOffProcElmts")]
        extern public static int SetMaxOffProcElmts(T_IJMatrix matrix, int max_off_proc_elmts);
    }
}
