using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
using ilPSP;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using static ilPSP.LinSolvers.ILU.Wrapper_MKL;

namespace ilPSP.LinSolvers.ILU {
    public class ILUSolver : ISparseSolver, IDisposable {

        public void DefineMatrix(IMutableMatrixEx inputMatrix) {
            m_Matrix= inputMatrix;
            wrapper = new Wrapper_MKL();
            out_Matrix = new MsrMatrix(inputMatrix.RowPartitioning);

            if (tmpMatrix == null)
                tmpMatrix = new Matrix(m_Matrix, this.UseDoublePrecision);
            Debug.Assert(tmpMatrix.Symmetric == false);

            SetParameters();

            if (m_Matrix.RowPartitioning.MpiSize == 1)
                LocalFactorization();
            else
                throw new NotImplementedException("ILU is only implemented for local use only! Note: There is no parallel version available in the Intel MKL right now.");
            tmpMatrix.Dispose();
            tmpMatrix = null;
        }

        IMutableMatrixEx m_Matrix;
        MsrMatrix out_Matrix;
        Matrix tmpMatrix;
        Wrapper_MKL wrapper;
        bool UseDoublePrecision=true;
        ILUParameters m_parms= new ILUParameters();
        Matrix ILUfactorization;
        IntPtr mklHandleForCSRMatrix;

        public SolverResult Solve<Tunknowns, TRHS>(Tunknowns X, TRHS B)
            where Tunknowns : IList<double>
            where TRHS : IList<double>
            {
            using (new FuncTrace()) {
                double[] _x;
                double[] _rhs;
                if (X.GetType() == typeof(double[])) {
                    _x = X as double[];
                } else {
                    _x = new double[X.Count];
                }
                if (B.GetType() == typeof(double[]))
                    _rhs = (double[])((ICloneable)B).Clone();
                else {
                    int L = m_Matrix.RowPartitioning.LocalLength;
                    _rhs = new double[L];
                    for (int i = 0; i < L; i++) _rhs[i] = B[i];
                }

                var solverrun = new Stopwatch();
                solverrun.Start();
                if (m_Matrix.RowPartitioning.MpiSize == 1) {
                    double[] _y = new double[_x.Length];
                    LocalSubstitution(_y, _rhs, sparse_fill_mode_t.SPARSE_FILL_MODE_LOWER); // Ly=b
                    LocalSubstitution(_x, _y, sparse_fill_mode_t.SPARSE_FILL_MODE_UPPER); // Ux=y
                } else
                    throw new NotImplementedException("ILU is only implemented for local use only! Note: There is no parallel version available in the Intel MKL right now.");
                solverrun.Stop();

                var sresult = new SolverResult() {
                    RunTime = solverrun.Elapsed,
                    Converged = true,
                    NoOfIterations = 1,
                };

                return sresult;
            }
        }

        private void SetParameters() {
            // Note: Fortran library so first index 1 in arrays!
            m_parms.iparm[30] = 0; //if 0, then computation will be stopped if zero on diagonal occur, otherwise a value defined in double parameters will be inserted
            //m_parms.dparm[30] = //defines a threshold for insertion, default is 1.0e-16, if value < dpar[31], then value = dpar[32]
            m_parms.dparm[31] = 1.0; //replace zero diagonal with this value, default is 1.0e-10
        }


        private unsafe void LocalSubstitution(double[] sol, double[] rhs, sparse_fill_mode_t Mode) {
            fixed (double* p_sol = sol, p_rhs = rhs) {
                var spM = CreateInternalMatrixFormat();

                //double* outA;
                //sparse_index_base_t outindexing;
                //int outNrows, outNcols;
                //int* outp_ia, outp_ia2, outp_ja;
                //wrapper.Export(spM, out outindexing, out outNrows, out outNcols, out outp_ia, out outp_ia2, out outp_ja, out outA);


                sparse_operation_t operation;
                double alpha = 1.0;
                sparse_diag_type_t Diag;
                switch (Mode) {
                    case sparse_fill_mode_t.SPARSE_FILL_MODE_LOWER:
                        Diag = sparse_diag_type_t.SPARSE_DIAG_UNIT;
                    operation = sparse_operation_t.SPARSE_OPERATION_NON_TRANSPOSE;
                        break;
                    case sparse_fill_mode_t.SPARSE_FILL_MODE_UPPER:
                        Diag = sparse_diag_type_t.SPARSE_DIAG_NON_UNIT;
                        operation = sparse_operation_t.SPARSE_OPERATION_TRANSPOSE;
                        break;
                    default:
                        throw new NotSupportedException();
                }
                var Mdescr = new MATRIX_DESCR() {
                    type = sparse_matrix_type_t.SPARSE_MATRIX_TYPE_TRIANGULAR,
                    mode = Mode,
                    diag = Diag
                };
                //int* operation, double* alpha, SPARSE_MATRIX_T A, MATRIX_DESCR descr, double* x, double* y
                int err = wrapper.Substitute(operation, alpha, spM, Mdescr, p_rhs, p_sol);
                if (err != 0)
                    throw new Exception("PadaBuuum! "+ MKLstatus2string(err));
            }
        }

        unsafe private IntPtr CreateInternalMatrixFormat() {
            int nonZ = Convert.ToInt32(m_Matrix.GetTotalNoOfNonZeros());
            int Nrows = m_Matrix.RowPartitioning.LocalLength;
            int Ncols = m_Matrix.ColPartition.LocalLength;
            Debug.Assert(Nrows == Ncols);
            IntPtr spM = new IntPtr();
            fixed (int* p_ia = ILUfactorization.ia, p_ja = ILUfactorization.ja) {
                double* A = (double*)ILUfactorization.aPtr;
                var indexing = sparse_index_base_t.SPARSE_INDEX_BASE_ONE;
                //Wrapper_MKL.SPARSE_MATRIX_T spM;
                //GCHandle pin_ia = GCHandle.Alloc(ILUfactorization.ia, GCHandleType.Pinned);
                //GCHandle pin_ja = GCHandle.Alloc(ILUfactorization.ja, GCHandleType.Pinned);
                //GCHandle pin_A = GCHandle.Alloc(*A, GCHandleType.Pinned);
                
                //var spM = Marshal.AllocHGlobal(nonZ*10 * sizeof(double));
                //_d_crcsr(SPARSE_MATRIX_T A, int* indexing, int* rows, int* cols, int* rows_start, int* rows_end, int* col_indx, double* values)
                int err = wrapper.CreateCSRMatrix(&spM, indexing, Nrows, Ncols, p_ia, p_ia + 1, p_ja, A);
                if (err != 0)
                    throw new Exception("Buuum! " + MKLstatus2string(err));
                Console.WriteLine("status of csr creation: "+MKLstatus2string(err));
            }
            //IntPtr* source, sparse_index_base_t indexing, int rows, int cols, int* rows_start, int* rows_end, int* col_indx, double* values
            //double* outA;
            //sparse_index_base_t outindexing;
            //int outNrows, outNcols;
            //int* outp_ia, outp_ia2, outp_ja;
            //wrapper.Export(spM, out outindexing, out outNrows, out outNcols, out outp_ia, out outp_ia2, out outp_ja, out outA);

            //int stat1 = wrapper.Optimize(spM);
            //Console.WriteLine("1:" + MKLstatus2string(stat1));
            return spM;
        }

        private unsafe void LocalFactorization() {
            int nonZ = Convert.ToInt32(m_Matrix.GetTotalNoOfNonZeros());
            var bilu = new double[nonZ];
            var ia = tmpMatrix.ia;
            var ja = tmpMatrix.ja;

            fixed (double* dparm = m_parms.dparm, p_bilu = bilu) {
                fixed (int* p_ia = ia, p_ja = ja, iparm = m_parms.iparm) {
                    int N = Convert.ToInt32(m_Matrix.NoOfRows);
                    double* a = (double*)(tmpMatrix.aPtr);
                    //double* bilu = (double*)Marshal.AllocHGlobal(nonZ * sizeof(double));
                    int error = 0;
                    //int* n, double* a, int* ia, int* ja, double* bilu0, int* ipar, double* dpar, int* ierr);
                    wrapper.ILU0(&N, a, p_ia, p_ja, p_bilu, iparm, dparm, &error);
                    Console.WriteLine(ILUerror2string(error));

                    ILUfactorization = new Matrix((IntPtr)p_bilu, ia, ja, N);
                    //TranslateMatrixBack(bilu,ia,ja,nonZ);
                    //Marshal.FreeHGlobal((IntPtr)bilu);
                }
            }
            Debug.Assert(ILUfactorization.aPtr != null);
        }

        public void Dispose() {
            m_Matrix = null;
            tmpMatrix.Dispose();
            tmpMatrix = null;
            ILUfactorization.Dispose();
            ILUfactorization=null;
        }

        private unsafe void TranslateMatrixBack(double* M, int* rowoffsetinarray, int* colidx, int NoNZ) {
            out_Matrix.Clear();
            int r = -1;
            for(int i =0;i< NoNZ; i++) {
                int c = colidx[i]-1;
                if (i >= rowoffsetinarray[r+1]-1)
                    r++;
                out_Matrix[r, c] = M[i];
            }
        }

        //public MsrMatrix ILUfactorization {
        //    get { return out_Matrix;  }
        //}
    }



    class ILUParameters {
        public int[] iparm = new int[128];
        public double[] dparm = new double[128];
    }

}
