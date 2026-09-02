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
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

namespace ilPSP.LinSolvers.UMFPACK {

    /// <summary>
    /// Sparse direct solver using UMFPACK.
    /// </summary>
    public class UMFPACKSolver : IExternalDirectSolver, IDisposable {

        const int UMFPACK_At = 1;
        const int UMFPACK_WARNING_singular_matrix = 1;

        /// <summary>
        /// ctor
        /// </summary>
        public UMFPACKSolver() {
        }

        Parallelism m_Parallelism = Parallelism.OMP;

        /// <summary>
        /// Level of parallelism used to select the native UMFPACK library.
        /// </summary>
        public Parallelism Parallelism {
            get {
                return m_Parallelism;
            }
            set {
                if(m_UnsafeUMFPACK != null)
                    throw new NotSupportedException("Cannot be changed after init.");
                m_Parallelism = value;
            }
        }

        /// <summary>
        /// Communicator on which the solver is defined.
        /// </summary>
        public MPI_Comm MpiComm {
            get;
            private set;
        }

        UnsafeUMFPACK wrapper {
            get {
                if(m_UnsafeUMFPACK == null)
                    m_UnsafeUMFPACK = new UnsafeUMFPACK(this.Parallelism);
                return m_UnsafeUMFPACK;
            }
        }

        UnsafeUMFPACK m_UnsafeUMFPACK;
        IMutableMatrixEx m_OrgMatrix;
        Matrix m_UMFPACKMatrix;
        int m_MpiRank = -1;
        int m_OrgMatrixLength = -1;

        unsafe void* m_Symbolic = null;
        unsafe void* m_Numeric = null;

        bool m_CacheFactorization = true;

        /// <summary>
        /// Cache the UMFPACK symbolic and numeric factorization across repeated solves with the same matrix.
        /// </summary>
        public bool CacheFactorization {
            get {
                return m_CacheFactorization;
            }
            set {
                if(!value && m_CacheFactorization)
                    UMFPACKDispose();
                m_CacheFactorization = value;
            }
        }

        int OrgLength {
            get {
                if(m_OrgMatrix != null)
                    m_OrgMatrixLength = m_OrgMatrix.RowPartitioning.LocalLength;
                if(m_OrgMatrixLength == -1 && m_OrgMatrix == null)
                    throw new Exception("Input Matrix deleted!");
                return m_OrgMatrixLength;
            }
        }

        /// <summary>
        /// <see cref="ISparseSolver.DefineMatrix"/>
        /// </summary>
        public void DefineMatrix(IMutableMatrixEx M) {
            if(m_OrgMatrix != null)
                throw new ApplicationException("matrix is already defined. 'DefineMatrix'-method can be invoked only once in the lifetime of this object.");
            if(M.NoOfCols != M.NoOfRows)
                throw new ArgumentException("Expecting quadratic matrix.");
            if(M.NoOfCols <= 0 || M.NoOfRows <= 0)
                throw new ArgumentException("Matrix must have a non-zero size.");

            m_MpiRank = M.RowPartitioning.MpiRank;
            m_OrgMatrix = M;
            MpiComm = M.MPI_Comm;
            MPICollectiveWatchDog.Watch(this.MpiComm);
        }

        /// <summary>
        /// <see cref="ISparseSolver.Solve{Tunknowns, Trhs}"/>
        /// </summary>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            return Solve<double[], Tunknowns, Trhs>(1.0, null, x, rhs);
        }

        /// <summary>
        /// <see cref="ISparseSolverExt.Solve{Tdiag, Tunknowns, Trhs}"/>
        /// </summary>
        public SolverResult Solve<Tdiag, Tunknowns, Trhs>(double Scale, Tdiag d, Tunknowns x, Trhs rhs)
            where Tdiag : IList<double>
            where Tunknowns : IList<double>
            where Trhs : IList<double> {

            using(var tr = new FuncTrace()) {
                CheckNestedParallism();

                int LOrg = OrgLength;

                if(x.Count != LOrg)
                    throw new ArgumentException("length of x must be equal to matrix size.");
                if(rhs.Count != LOrg)
                    throw new ArgumentException("length of rhs must be equal to matrix size.");
                if(Math.Abs(Scale) <= double.Epsilon)
                    throw new ArgumentException("scale to small.");

                double[] _x;
                double[] _rhs;
                if(x.GetType() == typeof(double[])) {
                    _x = x as double[];
                } else {
                    _x = new double[x.Count];
                }
                if(rhs.GetType() == typeof(double[]))
                    _rhs = (double[])((ICloneable)rhs).Clone();
                else {
                    _rhs = new double[LOrg];
                    for(int i = 0; i < LOrg; i++)
                        _rhs[i] = rhs[i];
                }

                IMutableMatrixEx Mtx;
                bool throwAwayMtx;
                if(d != null || Scale != 1.0) {
                    MsrMatrix _Mtx = new MsrMatrix(m_OrgMatrix);
                    Mtx = _Mtx;
                    if(Scale != 1.0)
                        _Mtx.Scale(Scale);
                    if(d != null) {
                        int dLen = d.Count, L = _Mtx.RowPartitioning.LocalLength, i0 = (int)(_Mtx.RowPartitioning.i0);
                        if((L % dLen) != 0)
                            throw new ApplicationException("wrong length of 'd'.");

                        for(int l = 0; l < L; l++)
                            _Mtx[i0 + l, i0 + l] += d[l % dLen];
                    }
                    throwAwayMtx = true;
                } else {
                    Mtx = m_OrgMatrix;
                    throwAwayMtx = false;
                }

                if(throwAwayMtx)
                    ClearUMFPACKMatrix();

                SolverResult r = new SolverResult();
                using(new BlockTrace("CALL SOLVER", tr)) {
                    r.Converged = true;
                    r.NoOfIterations = 1;
                    Stopwatch st = new Stopwatch();
                    st.Start();

                    double[] gath_x, gath_b;
                    GatherOnProc0(_x, _rhs, out gath_x, out gath_b);
                    Debug.Assert(m_MpiRank == 0 || gath_x == null);
                    Debug.Assert(m_MpiRank == 0 || gath_b == null);

                    r.Converged = UMFPACKInitAndSolve(Mtx, gath_x, gath_b);

                    if(!m_CacheFactorization || throwAwayMtx)
                        UMFPACKDispose();
                    if(throwAwayMtx)
                        ClearUMFPACKMatrix();

                    csMPI.Raw.Barrier(MpiComm);
                    unsafe {
                        bool snd = r.Converged;
                        csMPI.Raw.Bcast((IntPtr)(&snd), 1, csMPI.Raw._DATATYPE.BYTE, 0, MpiComm);
                        r.Converged = snd;
                    }

                    ScatterFromProc0(gath_x, _x);
                    csMPI.Raw.Barrier(MpiComm);

                    if(x.GetType() != typeof(double[])) {
                        for(int i = 0; i < LOrg; i++)
                            x[i] = _x[i];
                    }

                    st.Stop();
                    r.RunTime = st.Elapsed;
                }

                return r;
            }
        }

        void CheckNestedParallism() {
            if(ilPSP.Environment.InParallelSection && Parallelism == Parallelism.OMP)
                throw new ApplicationException("Nested parallelism detected: UMFPACK with OpenMP parallelism is called from within a parallel region. This is not supported.");
        }

        bool UMFPACKInitAndSolve(IMutableMatrixEx Mtx, double[] _x, double[] _b) {
            using(var tr = new FuncTrace()) {
                if(m_UMFPACKMatrix == null)
                    m_UMFPACKMatrix = new Matrix(Mtx);

                int rank;
                csMPI.Raw.Comm_Rank(MpiComm, out rank);

                if(rank == 0) {
                    unsafe {
                        fixed(double* px = _x, pb = _b) {
                            fixed(long* Ap = m_UMFPACKMatrix.Ap, Ai = m_UMFPACKMatrix.Ai) {
                                double* Ax = (double*)m_UMFPACKMatrix.AxPtr;

                                if(m_Symbolic == null) {
                                    using(new BlockTrace("UMFPACK_symbolic", tr)) {
                                        void* symbolic = null;
                                        int status = wrapper.Symbolic(m_UMFPACKMatrix.n, m_UMFPACKMatrix.n, Ap, Ai, Ax, &symbolic, null, null);
                                        if(status < 0 || status == UMFPACK_WARNING_singular_matrix)
                                            return HandleStatus(status);
                                        m_Symbolic = symbolic;
                                    }
                                }

                                if(m_Numeric == null) {
                                    using(new BlockTrace("UMFPACK_numeric", tr)) {
                                        void* numeric = null;
                                        int status = wrapper.Numeric(Ap, Ai, Ax, m_Symbolic, &numeric, null, null);
                                        if(status < 0 || status == UMFPACK_WARNING_singular_matrix)
                                            return HandleStatus(status);
                                        m_Numeric = numeric;
                                    }
                                }

                                using(new BlockTrace("UMFPACK_solve", tr)) {
                                    // The ilPSP matrix is assembled as CSR(A). UMFPACK expects CSC, and the same
                                    // arrays are CSC(A^T). Solving A^T x = b for the interpreted matrix therefore
                                    // solves the original system A x = b.
                                    int status = wrapper.Solve(UMFPACK_At, Ap, Ai, Ax, px, pb, m_Numeric, null, null);
                                    if(status < 0 || status == UMFPACK_WARNING_singular_matrix)
                                        return HandleStatus(status);
                                }
                            }
                        }
                    }
                    SendErrorToOtherRanks(0);
                } else {
                    int[] error = new int[] { -1 };
                    unsafe {
                        fixed(int* Buffer = error) {
                            MPI_Status status;
                            csMPI.Raw.Recv((IntPtr)Buffer, 1, csMPI.Raw._DATATYPE.INT, 0, 112, MpiComm, out status);
                        }
                    }
                    if(error[0] != 0)
                        wrapper.UMFPACKerror2string(error[0]);
                }

                return true;
            }
        }

        bool HandleStatus(int status) {
            UMFPACKDispose();
            SendErrorToOtherRanks(status);
            Console.Error.WriteLine("UMFPACK ERROR: " + wrapper.UMFPACKerror2string(status));
            return false;
        }

        void SendErrorToOtherRanks(int error) {
            csMPI.Raw.Comm_Size(MpiComm, out int size);
            for(int iProc = 1; iProc < size; iProc++) {
                int[] buffer = new int[] { error };
                unsafe {
                    fixed(int* Buffer = buffer) {
                        csMPI.Raw.Send((IntPtr)Buffer, 1, csMPI.Raw._DATATYPE.INT, iProc, 112, MpiComm);
                    }
                }
            }
        }

        void ScatterFromProc0(double[] P0__x_IN, double[] Px_x_OUT) {
            int size, rank;
            csMPI.Raw.Comm_Size(this.MpiComm, out size);
            csMPI.Raw.Comm_Rank(this.MpiComm, out rank);
            Debug.Assert(size == m_OrgMatrix.RowPartitioning.MpiSize);
            Debug.Assert(size == m_OrgMatrix.ColPartitioning.MpiSize);
            Debug.Assert(rank == m_OrgMatrix.RowPartitioning.MpiRank);
            Debug.Assert(rank == m_OrgMatrix.ColPartitioning.MpiRank);

            if(rank == 0)
                Debug.Assert(P0__x_IN.Length == m_OrgMatrix.RowPartitioning.TotalLength);
            else
                Debug.Assert(P0__x_IN == null);
            Debug.Assert(Px_x_OUT.Length == m_OrgMatrix.RowPartitioning.LocalLength);

            if(size > 1) {
                int[] SendCounts = new int[size];
                int[] Displ = new int[size];
                for(int r = 0; r < size; r++) {
                    SendCounts[r] = m_OrgMatrix.RowPartitioning.GetLocalLength(r);
                    Displ[r] = checked((int)(m_OrgMatrix.RowPartitioning.GetI0Offest(r)));
                }

                unsafe {
                    fixed(int* pSendCounts = SendCounts, pDispl = Displ) {
                        fixed(double* pSend = P0__x_IN, pRecv = Px_x_OUT) {
                            csMPI.Raw.Scatterv((IntPtr)pSend, (IntPtr)pSendCounts, (IntPtr)pDispl,
                                csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pRecv, SendCounts[rank],
                                csMPI.Raw._DATATYPE.DOUBLE, 0, this.MpiComm);
                        }
                    }
                }
            } else {
                if(!object.ReferenceEquals(P0__x_IN, Px_x_OUT)) {
                    int L = P0__x_IN.Length;
                    if(P0__x_IN.Length != Px_x_OUT.Length)
                        throw new ApplicationException("internal error.");
                    Array.Copy(P0__x_IN, Px_x_OUT, L);
                }
            }
        }

        void GatherOnProc0(double[] xIN, double[] bIN, out double[] gath_x_OUT, out double[] gath_b_OUT) {
            int size, rank;
            csMPI.Raw.Comm_Size(this.MpiComm, out size);
            csMPI.Raw.Comm_Rank(this.MpiComm, out rank);
            Debug.Assert(size == m_OrgMatrix.RowPartitioning.MpiSize);
            Debug.Assert(size == m_OrgMatrix.ColPartitioning.MpiSize);
            Debug.Assert(rank == m_OrgMatrix.RowPartitioning.MpiRank);
            Debug.Assert(rank == m_OrgMatrix.ColPartitioning.MpiRank);

            if(size > 1) {
                int[] RcvCounts = new int[size];
                for(int r = 0; r < size; r++)
                    RcvCounts[r] = m_OrgMatrix.RowPartitioning.GetLocalLength(r);

                gath_x_OUT = MPIExtensions.MPIGatherv(xIN, RcvCounts, 0, this.MpiComm);
                gath_b_OUT = MPIExtensions.MPIGatherv(bIN, RcvCounts, 0, this.MpiComm);

                Debug.Assert((rank == 0) != (gath_x_OUT == null));
                Debug.Assert((rank == 0) || (gath_b_OUT == null));
            } else {
                gath_x_OUT = xIN;
                gath_b_OUT = bIN;
            }
        }

        /// <summary>
        /// Returns the internally stored matrix.
        /// </summary>
        public ISparseMatrix GetMatrix() {
            return m_OrgMatrix;
        }

        /// <summary>
        /// Frees UMFPACK and matrix storage.
        /// </summary>
        public void Dispose() {
            UMFPACKDispose();
            ClearUMFPACKMatrix();
        }

        /// <summary>
        /// calls <see cref="Dispose"/>
        /// </summary>
        ~UMFPACKSolver() {
            Dispose();
        }

        void UMFPACKDispose() {
            unsafe {
                if(m_Numeric != null) {
                    void* numeric = m_Numeric;
                    wrapper.FreeNumeric(&numeric);
                    m_Numeric = null;
                }
                if(m_Symbolic != null) {
                    void* symbolic = m_Symbolic;
                    wrapper.FreeSymbolic(&symbolic);
                    m_Symbolic = null;
                }
            }
        }

        void ClearUMFPACKMatrix() {
            if(m_UMFPACKMatrix != null) {
                m_UMFPACKMatrix.Dispose();
                m_UMFPACKMatrix = null;
            }
        }
    }
}
