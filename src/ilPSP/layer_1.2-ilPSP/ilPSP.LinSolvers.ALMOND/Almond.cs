using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MPI.Wrappers;
using System.Diagnostics;
using ilPSP;
using ilPSP.Utils;
using ilPSP.LinSolvers;

namespace ilPSP.LinSolvers.ALMOND {

    public enum KrylovType {
        PCG = 0, 
        GMRES = 1
    }
    
    
    public class Almond : ilPSP.LinSolvers.ISparseSolver {

        public Almond() {
            int size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            if(size > 1)
                throw new NotSupportedException("not supported for MPI parallelism.");
        }


        static bool KernelsInitialized = false;

        IntPtr MatrixHandle = IntPtr.Zero;
        IntPtr SolverHandle = IntPtr.Zero;

        MsrMatrix Matrix;

        public void DefineMatrix(IMutableMatrixEx M) {
            if(M.RowPartitioning.TotalLength != M.ColPartition.TotalLength)
                throw new ArgumentException("Matrix must be symmetric.");
            this.Matrix = M.ToMsrMatrix();
        }


        void Init() {
            if(SolverHandle != IntPtr.Zero)
                return;

            int ierr;
            if(KernelsInitialized == false) {
                int occaPlatId = 0, occaDevId = 0;
                Wrappers.ALMOND_Init("OpenMP", ref occaPlatId, ref occaDevId, out ierr);
                Wrappers.ALMOND_buildAlmondKernels_d("C:\\Users\\florian\\Documents\\ALMOND", out ierr);
                CheckErr(ierr);
                KernelsInitialized = true;
            }

            var M = this.Matrix;

            int nnz = M.GetTotalNoOfNonZeros();
            int[] ia = new int[nnz];
            int[] ja = new int[nnz];
            double[] aa = new double[nnz];
            int c = 0;
            int i0 = M.RowPartitioning.i0, L = M.RowPartitioning.LocalLength;
            int LR;
            int[] col = null;
            double[] val = null;
            for (int i = 0; i < L; i++) {
                int iRow = i0 + i;
                LR = M.GetRow(iRow, ref col, ref val);
                for(int j = 0; j < LR; j++) {
                    ia[c] = iRow;
                    ja[c] = col[j];
                    aa[c] = val[j];

                    c++;
                }
            }
            Debug.Assert(c == nnz);

            int N = M.RowPartitioning.TotalLength;
            Wrappers.ALMOND_DefineMatrix_d(ref N, ref N, ref nnz, ia, ja, aa, out MatrixHandle, out ierr);
            CheckErr(ierr);
            Wrappers.ALMOND_CreateSolver_d(ref MatrixHandle, out SolverHandle, out ierr);
            CheckErr(ierr);
            
        }

        private static void CheckErr(int ierr) {
            if(ierr != 0)
                throw new ApplicationException("ALMOND error: " + ierr);
        }

        /// <summary>
        /// implementation of <see cref="MaxIterations"/>
        /// </summary>
        protected int m_MaxIterations = 10000;

        /// <summary>
        /// upper limit for the number of iterations of the sparse solver;
        /// </summary>
        public int MaxIterations {
            get {
                return m_MaxIterations;
            }
            set {
                if(value < 0)
                    throw new ArgumentOutOfRangeException();
                m_MaxIterations = value;
            }
        }

        /// <summary>
        /// see <see cref="Tolerance"/>;
        /// </summary>
        protected double m_Tolerance = 1.0e-9;

        /// <summary>
        /// Residual threshold for the iterative solver termination
        /// </summary>
        public double Tolerance {
            get {
                return m_Tolerance;
            }
            set {
                if(value < 0.0)
                    throw new ArgumentOutOfRangeException();
                m_Tolerance = value;
            }
        }

        KrylovType m_KrylovType = KrylovType.PCG;

        public KrylovType KrylovType {
            get {
                return m_KrylovType;
            }
            set {
                m_KrylovType = value;
            }
        }



        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns _x, Trhs _rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> 
        {
            Init();

            if(_x.Count != this.Matrix.RowPartitioning.LocalLength)
                throw new ArgumentException("Mismatch in length of X vector.");
            if(_rhs.Count != this.Matrix.RowPartitioning.LocalLength)
                throw new ArgumentException("Mismatch in lenght of RHS vector.");
            
            double[] x = _x as double[];
            if(x == null)
                x = _x.ToArray();
            double[] rhs = _rhs as double[];
            if(x == null)
                rhs = _rhs.ToArray();

            SolverResult R = new SolverResult();
            Stopwatch stw = new Stopwatch();
            stw.Start();

            int ktype = (int) this.m_KrylovType;
            int maxI = this.m_MaxIterations;
            int ierr = 0;
            int Nrows = this.Matrix.RowPartitioning.LocalLength;
            int NoIter;
            Wrappers.ALMOND_UseSolver_d(ref Nrows, ref ktype, ref maxI, ref this.m_Tolerance, ref this.SolverHandle, x, rhs, out NoIter, out ierr);
            CheckErr(ierr);

            stw.Stop();
            R.RunTime = stw.Elapsed;
            R.Converged = NoIter >= 0;
            R.NoOfIterations = Math.Abs(NoIter);

            if(!object.ReferenceEquals(x, _x)) {
                _x.SetV(x);
            }
            if(!object.ReferenceEquals(rhs, _rhs)) {
                _rhs.SetV(rhs);
            }
            return R;
        }

        public void Dispose() {
            if(SolverHandle != IntPtr.Zero) {
                int ierr;
                Wrappers.ALMOND_DeleteSolver_d(ref SolverHandle, out ierr);
                SolverHandle = IntPtr.Zero;
            }
            if(MatrixHandle != IntPtr.Zero) {
                int ierr;
                Wrappers.ALMOND_DeleteMatrix_d(ref MatrixHandle, out ierr);
                MatrixHandle = IntPtr.Zero;
            }
        }
    }
}
