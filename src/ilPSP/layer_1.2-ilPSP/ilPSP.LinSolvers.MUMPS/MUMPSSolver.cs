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

using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.Diagnostics;
using ilPSP.Utils;

namespace ilPSP.LinSolvers.MUMPS {
    public class MUMPSSolver : ISparseSolver, IDisposable {

        MUMPS_csharp.DMUMPS_STRUC_CS mumps_par = new MUMPS_csharp.DMUMPS_STRUC_CS();

        Matrix m_MumpsMatrix;

        IMutableMatrixEx m_OrgMatrix;

        bool verbose;

        MPI_Comm m_MPI_Comm;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="verbose"></param>
        /// <param name="MPI"></param>
        public MUMPSSolver(bool verbose = false, bool MPI = true) {
            this.verbose = verbose;
            if (MPI == true) {
                this.m_MPI_Comm = csMPI.Raw._COMM.WORLD;
            } else {
                this.m_MPI_Comm = csMPI.Raw._COMM.SELF;
            }
            
        }

        public void DefineMatrix(IMutableMatrixEx M) {
            m_MumpsMatrix = new Matrix(M);
            m_OrgMatrix = M;

        }

        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            return Solve<double[], Tunknowns, Trhs>(1.0, null, x, rhs);
        }

        public SolverResult Solve<Tdiag, Tunknowns, Trhs>(double Scale, Tdiag d, Tunknowns x, Trhs rhs)
            where Tdiag : System.Collections.Generic.IList<double>
            where Tunknowns : System.Collections.Generic.IList<double>
            where Trhs : System.Collections.Generic.IList<double> {

            using (var tr = new FuncTrace()) {

                // check input arguments
                // =====================

                if (x.Count != m_OrgMatrix.RowPartitioning.LocalLength)
                    throw new ArgumentException("length of x must be equal to matrix size.");
                if (rhs.Count != m_OrgMatrix.RowPartitioning.LocalLength)
                    throw new ArgumentException("length of rhs must be equal to matrix size.");
                if (Math.Abs(Scale) <= double.Epsilon)
                    throw new ArgumentException("scale to small.");

                // prepare vectors
                // ===============
                double[] _x;
                double[] _rhs;
                if (x.GetType() == typeof(double[])) {
                    _x = x as double[];
                } else {
                    _x = new double[x.Count];
                }
                if (rhs.GetType() == typeof(double[]))
                    _rhs = (double[])((ICloneable)rhs).Clone();
                else {
                    int L = m_MumpsMatrix.RowPart.LocalLength;
                    _rhs = new double[L];
                    for (int i = 0; i < L; i++) _rhs[i] = rhs[i];
                }

                // define matrix
                // =============

                IMutableMatrixEx Mtx;
                //bool throwAwayMtx;
                if (d != null || Scale != 1.0) {
                    // not very efficient, but _I_dont_care_ !

                    MsrMatrix _Mtx = new MsrMatrix(m_OrgMatrix);
                    Mtx = _Mtx;
                    _Mtx.Scale(Scale);

                    int dLen = d.Count, L = _Mtx.RowPartitioning.LocalLength, i0 = (int)(_Mtx.RowPartitioning.i0);
                    if ((L % dLen) != 0)
                        throw new ApplicationException("wrong length of 'd'.");

                    for (int l = 0; l < L; l++) {
                        _Mtx[i0 + l, i0 + l] += d[l % dLen];
                    }

                    //throwAwayMtx = true;

                } else {
                    //throwAwayMtx = false;
                    Mtx = m_OrgMatrix;
                }

                // call solver
                // ===========

                SolverResult r = new SolverResult();
                using (new BlockTrace("CALL SOLVER", tr)) {
                    r.Converged = true;
                    r.NoOfIterations = 1;
                    Stopwatch st = new Stopwatch();
                    st.Reset();
                    st.Start();

                    double[] gath_x, gath_b;
                    GatherOnProc0(_x, _rhs, out gath_x, out gath_b);

                    r.NoOfIterations = 1;
                    int rank; int size;
                    csMPI.Raw.Comm_Rank(this.m_MPI_Comm, out rank);
                    csMPI.Raw.Comm_Size(this.m_MPI_Comm, out size);

                    r.Converged = MUMPSInitAndSolve(m_OrgMatrix, gath_b);

                    csMPI.Raw.Barrier(this.m_MPI_Comm);
                    unsafe
                    {
                        bool snd = r.Converged;
                        csMPI.Raw.Bcast((IntPtr)(&snd), 1, csMPI.Raw._DATATYPE.BYTE, 0, this.m_MPI_Comm);
                        r.Converged = snd;
                    }

                    ScatterFromProc0(_x, gath_b);

                    if (x.GetType() == typeof(double[])) {
                        // do nothing
                    }
                    else {
                        x.SetV(_x);
                    }
                    

                    csMPI.Raw.Barrier(this.m_MPI_Comm);

                    st.Stop();

                    r.RunTime = st.Elapsed;
                }
                return r;

            }
        }

        bool MUMPSInitAndSolve(IMutableMatrixEx Mtx, double[] _b) {
            using (new FuncTrace()) {
                if (m_MumpsMatrix == null)
                    m_MumpsMatrix = new Matrix(Mtx);


                // Definition of Struc variables
                mumps_par.icntl = new int[40];
                mumps_par.keep = new int[500];
                mumps_par.cntl = new double[15];
                mumps_par.dkeep = new double[130];
                mumps_par.keep8 = new long[150];
                mumps_par.info = new int[40];
                mumps_par.infog = new int[40];
                mumps_par.rinfo = new double[40]; 
                mumps_par.rinfog = new double[40];

                int rank, size;
                csMPI.Raw.Comm_Rank(this.m_MPI_Comm, out rank);
                csMPI.Raw.Comm_Size(this.m_MPI_Comm, out size);



                //Initialize MUMPS
                mumps_par.par = 1; mumps_par.job = -1; mumps_par.comm_fortran = this.m_MPI_Comm.m1; mumps_par.sym = m_MumpsMatrix.Symmetric;

                MUMPS_csharp.mumps_cs(ref mumps_par);

                if (rank == 0) {
                    mumps_par.n = m_MumpsMatrix.n;
                    mumps_par.rhs = _b;
                }

                // Matrix Input Format
                // 0: assembled, 1: elemental
                mumps_par.icntl[4] = 0;

                // Defines the strategy for the distribution of the input matrx (only for assembled matrix)
                mumps_par.icntl[17] = 3;

                // if solution vector should be distributed set 1, otherwise set to 0 for rhs
                mumps_par.icntl[20] = 0;
                
                mumps_par.nz_loc = m_MumpsMatrix.nz_loc;
                mumps_par.irn_loc = m_MumpsMatrix.irn_loc;
                mumps_par.jcn_loc = m_MumpsMatrix.jrn_loc;


                if (!verbose) {
                    // NOTE: Comments are FORTRAN-style indexing, copied from the reference manual
                    // thus, our indices are always FORTRANINDEX - 1
                    //
                    // ICNTL(2) is the output stream for diagnostic printing, statistics, and warning messages.
                    mumps_par.icntl[1] = -1;
                    // ICNTL(3) is the output stream for global information, collected on the host.
                    mumps_par.icntl[2] = -1;
                    // ICNTL(4) is the level of printing for error, warning, and diagnostic messages
                    mumps_par.icntl[3] = 0;
                }
                
                // Call MUMPS
                mumps_par.job = 1;
                MUMPS_csharp.mumps_cs(ref mumps_par);

                mumps_par.a_loc = m_MumpsMatrix.a_loc;

                // corresponds to the maximum size of the working memory in MegaBytes that MUMPS can allocate per working processor
                if (rank == 0) {
                   // mumps_par.icntl[22] = 10000;
                    //Console.WriteLine("Estimated memory needed for factorization over all processes: " + mumps_par.infog[16]);
                }

                // corresponds to the percentage increase in the estimated working space: Default value: 20 (which corresponds to a 20 % increase).
                mumps_par.icntl[13] = 200;

                mumps_par.job = 2;
                MUMPS_csharp.mumps_cs(ref mumps_par);

                if (mumps_par.icntl[20] == 1) {
                    mumps_par.isol_loc = new int[mumps_par.info[22]];
                    mumps_par.lsol_loc = mumps_par.info[22];
                    mumps_par.sol_loc = new double[mumps_par.info[22]];
                }


                switch (mumps_par.info[0]) {
                    // no error, no warning
                    //=====================
                    case 0: break; 

                    // throw memory errors
                    //====================
                    case -8:
                    case -9:
                    case -11:
                    case -12:
                    case -14:
                    case -15:
                        throw new ApplicationException("A MUMPS memory error occured on proc with rank: " + rank +
                            ". Error Code:  " + mumps_par.info[0] +
                            "  (For further information see MUMPS handbook or contact your local MUMPS support)");

                    // throw singular Matrix
                    //======================
                    case -10:
                        throw new ApplicationException("MUMPS encountered a numerically singular Matrix");

                    // throw all other errors and warnings
                    // Just in case, those should be written to console anyway, since ICNTL(1) is set to default output
                    //=================================================================================================
                    default:
                        if (mumps_par.info[0] < 0) {
                            string ErrorString = String.Format("MUMPS: An Error occurred Info(1) ={0}, Info(2) = {1}", mumps_par.info[0], mumps_par.info[1]);
                            throw new ApplicationException(ErrorString);
                        }
                        else { // i.e. mumps_par.info[0] > 0 which are warnings
                            string WarningString = String.Format("MUMPS: An Error occurred Info(1) ={0}, Info(2) = {1}", mumps_par.info[0], mumps_par.info[1]);
                            Console.WriteLine(WarningString);
                            break;
                        }
                }
                               
                mumps_par.job = 3;
                MUMPS_csharp.mumps_cs(ref mumps_par);
                //if (rank != 0)
                //mumps_par.irn = new int[] { };
                //mumps_par.irn.SaveToTextFile("MUMPS_irn.txt");
                mumps_par.job = -2;
                MUMPS_csharp.mumps_cs(ref mumps_par);
                // }       

                return true;
            }
        }

        /// <summary>
        /// Scatters solution vector from process 0 to other MPI processors.
        /// </summary>
        /// <param name="__x">
        /// input; long vector on proc 0
        /// </param>
        /// <param name="_x">
        /// output;
        /// </param>
        private void ScatterFromProc0(double[] __x, double[] _x) {
            int size, rank;
            csMPI.Raw.Comm_Size(this.m_MPI_Comm, out size);
            csMPI.Raw.Comm_Rank(this.m_MPI_Comm, out rank);

            if (size > 1) {
                // distribute solution to other processors
                // +++++++++++++++++++++++++++++++++++++++

                if (rank == 0) {
                    Array.Copy(_x, 0, __x, 0, m_OrgMatrix.RowPartitioning.LocalLength);

                    unsafe
                    {
                        fixed (double* px = &_x[0])
                        {
                            for (int targ_rank = 1; targ_rank < size; targ_rank++) {
                                csMPI.Raw.Send((IntPtr)(px + m_OrgMatrix.RowPartitioning.GetI0Offest(targ_rank)), m_OrgMatrix.RowPartitioning.GetLocalLength(targ_rank), csMPI.Raw._DATATYPE.DOUBLE, targ_rank, 4444 + targ_rank, this.m_MPI_Comm);
                            }
                        }
                    }

                } else {
                    unsafe
                    {
                        fixed (double* px = &__x[0])
                        {
                            MPI_Status _st;
                            csMPI.Raw.Recv((IntPtr)px, m_OrgMatrix.RowPartitioning.LocalLength, csMPI.Raw._DATATYPE.DOUBLE, 0, 4444 + rank, this.m_MPI_Comm, out _st);
                        }
                    }
                }
            } else {
                if (!object.ReferenceEquals(__x, _x)) {
                    int L = __x.Length;
                    if (__x.Length != _x.Length)
                        throw new ApplicationException("internal error.");
                    Array.Copy(_x, __x, L);
                }

            }
        }

        /// <summary>
        /// Gathers RHS and solution vector on process 0 for more than one MPI process.
        /// </summary>
        /// <param name="__b">local part of rhs</param>
        /// <param name="__x">local part of solution/initial guess</param>
        /// <param name="_x">gathered rhs</param>
        /// <param name="_b">gathered solution vectors</param>
        private void GatherOnProc0(double[] __x, double[] __b, out double[] _x, out double[] _b) {
            int size, rank;
            csMPI.Raw.Comm_Size(this.m_MPI_Comm, out size);
            csMPI.Raw.Comm_Rank(this.m_MPI_Comm, out rank);

            if (size > 1) {
                // gather rhs on processor 0
                // +++++++++++++++++++++++++

                if (rank == 0) {
                    _x = new double[m_OrgMatrix.RowPartitioning.TotalLength];
                    _b = new double[m_OrgMatrix.RowPartitioning.TotalLength];

                    Array.Copy(__b, 0, _b, 0, m_OrgMatrix.RowPartitioning.LocalLength);

                    unsafe
                    {
                        fixed (double* pb = &_b[0])
                        {
                            for (int rcv_rank = 1; rcv_rank < size; rcv_rank++) {
                                MPI_Status stat;
                                csMPI.Raw.Recv((IntPtr)(pb + m_OrgMatrix.RowPartitioning.GetI0Offest(rcv_rank)), m_OrgMatrix.RowPartitioning.GetLocalLength(rcv_rank), csMPI.Raw._DATATYPE.DOUBLE, rcv_rank, 342346 + rcv_rank, this.m_MPI_Comm, out stat);
                            }
                        }
                    }

                } else {
                    // send my part to P0
                    unsafe
                    {
                        fixed (double* pb = &__b[0])
                        {
                            csMPI.Raw.Send((IntPtr)pb, m_OrgMatrix.RowPartitioning.LocalLength, csMPI.Raw._DATATYPE.DOUBLE, 0, 342346 + rank, this.m_MPI_Comm);
                        }
                    }

                    _x = null;
                    _b = null;
                }
            } else {
                _x = __x;
                _b = __b;
            }

        }

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing) {
            if (!disposedValue) {
                if (disposing) {
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }



        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~MUMPSSolver() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose() {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }
}
