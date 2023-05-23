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
using System.Collections.Generic;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Tracing;
using System.IO;
using ilPSP.Utils;

namespace ilPSP.LinSolvers.PARDISO {

    
    /// <summary>
    /// a direct solver using the PARDISO library from the Intel MKL
    /// </summary>
    /// <remarks>
    /// **IMPORTANT: Licensing issues:**
    /// PARDISO does not ship with free license, neither source nor 
    /// binaries compiled from it can be shipped with this software;
    /// 
    /// PARDISO is ether distributed with the INTEL MKL library, or it may be downloaded
    /// from http://www.pardiso-project.org/;
    /// </remarks>
    public class PARDISOSolver : ISparseSolverExt, IDisposable {

        /// <summary>
        /// ctor
        /// </summary>
        public PARDISOSolver() {
            
        }

        Parallelism m_Parallelism = Parallelism.OMP;

        /// <summary>
        /// Level of parallelism, which should be used for this solver instance 
        /// </summary>
        public Parallelism Parallelism {
            get {
                return m_Parallelism;

            }
            set {
                if(m_wrapper != null) {
                    throw new NotSupportedException("Cannot be changed after init.");
                }
                m_Parallelism = value;
            }
        }


        /// <summary>
        /// total time, over entire application instance lifetime, spent in 
        /// symbolic factorization (PARDISO phase 11)
        /// </summary>
        public static Stopwatch Phase_11 = new Stopwatch();
        
        
        /// <summary>
        /// total time, over entire application instance lifetime, spent in 
        /// numerical factorization (PARDISO phase 22)
        /// </summary>
        public static Stopwatch Phase_22 = new Stopwatch();
        
        /// <summary>
        /// total time, over entire application instance lifetime, spent in 
        /// back substitution and iterative refinement (PARDISO phase 33)
        /// </summary>
        public static Stopwatch Phase_33 = new Stopwatch();


        /// <summary>
        /// interface to native PARDISO library
        /// </summary>
        MetaWrapper wrapper {
            get {
                System.Environment.SetEnvironmentVariable("PARDISOLICMESSAGE", "1");
                if (m_wrapper == null) {
                    if ((this.Version == PARDISO.Version.v5) && this.LicenseCode != null && this.LicenseCode.Length > 0) {
                        string appDir = AppDomain.CurrentDomain.BaseDirectory;

                        StreamWriter licfile = new StreamWriter(Path.Combine(appDir, "pardiso.lic"));
                        licfile.WriteLine(this.LicenseCode);
                        licfile.Close();

                    }
                    if (this.Version == PARDISO.Version.v5)
                        System.Environment.SetEnvironmentVariable("PARDISOLICMESSAGE", "1");
                    
                    m_wrapper = new MetaWrapper(this.Version, this.Parallelism);
                }
                return m_wrapper;
            }
        }

        MetaWrapper m_wrapper = null;

        /// <summary>
        /// calls <see cref="Dispose"/>
        /// </summary>
        ~PARDISOSolver() {
            Dispose();
        }

        /// <summary>
        /// original matrix, as defined by <see cref="DefineMatrix"/>;
        /// must be stored in order to support the <see cref="ISparseSolverExt"/>-interface.
        /// </summary>
        /// <remarks>
        /// It is obvious that storing the matrix twice (i.e. <see cref="m_PardisoMatrix"/> and <see cref="m_OrgMatrix"/>) is not very 
        /// efficient in terms of memory use.
        /// However: if the matrix is quite large, than the size of the matrix negligible in comparison to
        /// the size of the factorization, so this does not play a role;
        /// on the other hand, is the matrix is small, the computer should have enough memory to handle both. 
        /// </remarks>
        IMutableMatrixEx m_OrgMatrix;

        /// <summary>
        /// PARDISO matrix
        /// </summary>
        Matrix m_PardisoMatrix;


        /// <summary>
        /// Communicator on which the solver is defined.
        /// </summary>
        public MPI_Comm MpiComm {
            get;
            private set;
        }

       
        
        /// <summary>
        /// converts <paramref name="M"/> into suitable structures for PARDISO;
        /// </summary>
        /// <param name="M"></param>
        public void DefineMatrix(IMutableMatrixEx M) {
            if (m_OrgMatrix != null)
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

        Version m_Version = Version.MKL;

        string m_LicenseCode = "";
        int m_OrgMatrixLength = -1;

        private int OrgLength {
            get {
                if (m_OrgMatrix != null) {
                    m_OrgMatrixLength = m_OrgMatrix.RowPartitioning.LocalLength;
                }
                if (m_OrgMatrixLength == -1 && m_OrgMatrix == null)
                    throw new Exception("Input Matrix deleted!");
                return m_OrgMatrixLength;
            }
        }

        /// <summary>
        /// PARDISO license code; if this is set, a file called 'pardiso.lic' will be 
        /// created in the application directory, containing this code.
        /// </summary>
        public string LicenseCode {
            get {
                return m_LicenseCode;
            }
            set {
                m_LicenseCode = value;
            }
        }

        /// <summary>
        /// Switch Pardiso version.
        /// </summary>
        public Version Version {
            get {
                return m_Version;
            }
            set {
                if(m_wrapper != null) {
                    throw new NotSupportedException("Cannot be changed after init.");
                }

                m_Version = value;
            }
        }

        /// <summary>
        /// not appropriate for direct solver;
        /// setting has no effect; returns always 0.0;
        /// </summary>
        public double Tolerance {
            get {
                return 0.0;
            }
            set { }
        }
        
        /// <summary>
        /// not appropriate for direct solver;
        /// setting has no effect; returns always 1;
        /// </summary>
        public int MaxIterations {
            get {
                return 1;
            }
            set { }
        }

        /// <summary>
        /// not appropriate for direct solver;
        /// setting has no effect;
        /// returns always <see cref="ConvergenceTypes.Other"/>
        /// </summary>
        public ConvergenceTypes ConvergenceType {
            get { return ConvergenceTypes.Other; }
            set { }
        }

        bool m_UseDoublePrecision = true;

        /// <summary>
        /// Use double (8 byte) or single (4 byte) fp numbers.
        /// </summary>
        public bool UseDoublePrecision {
            get {
                return m_UseDoublePrecision;
            }
            set {
                //if (!value)
                //    throw new ApplicationException("blocked for user protection.");
                m_UseDoublePrecision = value;
            }
        }


        /// <summary>
        /// solves the equation system 
        /// M*<paramref name="x"/>=<paramref name="rhs"/>,
        /// where M denotes the 
        /// that was provided by the constructor of this solver;
        /// </summary>
        /// <typeparam name="Tunknowns"></typeparam>
        /// <typeparam name="Trhs"></typeparam>
        /// <param name="x">vector for the unknowns</param>
        /// <param name="rhs">right hand side</param>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            return Solve<double[], Tunknowns, Trhs>(1.0, null, x, rhs);
        }

        /// <summary>
        /// solves the equation system 
        /// (diag(<paramref name="d"/>) + <paramref name="Scale"/>*M)*<paramref name="x"/>=<paramref name="rhs"/>,
        /// where diag(<paramref name="d"/>) denotes a diagonal matrix with the diagonal vector <paramref name="d"/>
        /// and M denotes the
        /// that was provided to the constructor of this solver;
        /// </summary>
        /// <typeparam name="Tunknowns"></typeparam>
        /// <typeparam name="Trhs"></typeparam>
        /// <typeparam name="Tdiag"></typeparam>
        /// <param name="Scale">
        /// scaling of the matrix (does not apply to <paramref name="d"/>);
        /// must be unequal to 0;
        /// this scaling does not alter the matrix-it's the same after the return of this method.
        /// </param>
        /// <param name="d">
        /// diagonal vector; can be null;
        /// length must be equal to the Nupdate-length of the matrix mapping; 
        /// may be altered during execution;
        /// </param>
        /// <param name="x">on exit, the approximate solution to the equation system;
        /// </param>
        /// <param name="rhs">
        /// right hand side of the equation system; may be overwritten;
        /// </param>
        public SolverResult Solve<Tdiag, Tunknowns, Trhs>(double Scale, Tdiag d, Tunknowns x, Trhs rhs)
            where Tdiag : System.Collections.Generic.IList<double>
            where Tunknowns : System.Collections.Generic.IList<double>
            where Trhs : System.Collections.Generic.IList<double> // 
        {

            //
            using (var tr = new FuncTrace()) {

                

                // check input arguments
                // =====================

                int LOrg = OrgLength;

                if (x.Count != LOrg)
                    throw new ArgumentException("length of x must be equal to matrix size.");
                if (rhs.Count != LOrg)
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
                    int L = LOrg;
                    _rhs = new double[L];
                    for (int i = 0; i < L; i++) _rhs[i] = rhs[i];
                }

                
                // define matrix
                // =============

                IMutableMatrixEx Mtx;
                bool throwAwayMtx;
                if(d != null || Scale != 1.0) {
                    MsrMatrix _Mtx = new MsrMatrix(m_OrgMatrix);
                    Mtx = _Mtx;
                    if(Scale != 1.0) {
                        // not very efficient, but _I_dont_care_ !
                        _Mtx.Scale(Scale);
                    }
                    if(d != null) {
                        int dLen = d.Count, L = _Mtx.RowPartitioning.LocalLength, i0 = (int)(_Mtx.RowPartitioning.i0);
                        if((L % dLen) != 0)
                            throw new ApplicationException("wrong length of 'd'.");

                        for(int l = 0; l < L; l++) {
                            _Mtx[i0 + l, i0 + l] += d[l % dLen];
                        }
                    }
                    throwAwayMtx = true;
                } else {
                    throwAwayMtx = false;
                    Mtx = m_OrgMatrix;
                }

                // call solver
                // ===========
                
                SolverResult r = new SolverResult();
                using(new BlockTrace("CALL SOLVER", tr)) {
                    r.Converged = true;
                    r.NoOfIterations = 1;
                    Stopwatch st = new Stopwatch();
                    st.Reset();
                    st.Start();

                    double[] gath_x, gath_b;
                    GatherOnProc0(_x, _rhs, out gath_x, out gath_b);
                    Debug.Assert(m_MpiRank == 0 || gath_x == null);
                    Debug.Assert(m_MpiRank == 0 || gath_b == null);

                    r.NoOfIterations = 1;
                    r.Converged = PARDISOInitAndSolve(Mtx, gath_x, gath_b);

                    if(!m_CacheFactorization || throwAwayMtx) {
                        PARDISODispose();
                    }

                    csMPI.Raw.Barrier(MpiComm);
                    unsafe {
                        bool snd = r.Converged;
                        csMPI.Raw.Bcast((IntPtr)(&snd), 1, csMPI.Raw._DATATYPE.BYTE, 0, MpiComm);
                        r.Converged = snd;
                    }

                    ScatterFromProc0(gath_x, _x);

                    csMPI.Raw.Barrier(MpiComm);

                    // write back / return
                    if(x.GetType() != typeof(double[])) {
                        int NN = LOrg;
                        for(int i = 0; i < NN; i++)
                            x[i] = _x[i];
                    }

                    // measure time
                    st.Stop();
                    r.RunTime = st.Elapsed;
                }
                return r;
            }
        }


        class PARDISOinternals {
            /// <summary>
            /// PARDISO auxiliary variable
            /// </summary>
            public int[] m_parm = new int[164];

            /// <summary>
            /// PARDISO auxiliary variable
            /// </summary>
            public int[] m_pt = new int[128];


            public double[] m_dparam = new double[100]; // eigentlich werden nur 64 gebraucht, zur sicherheit...


            public bool m_PardisoInitialized = false;

            //public int mtype = 11;        /* Real unsymmetric matrix */

            public int maxfct = 1;
            
            public int mnum = 1;
            
            public int msglvl = 0;
        }

        /// <summary>
        /// turn printing of solver statistics to stdout on/off;
        /// </summary>
        public bool PrintStats {
            get {
                return (m_PardInt.msglvl == 1);
            }
            set {
                if (value)
                    m_PardInt.msglvl = 1;
                else
                    m_PardInt.msglvl = 0;
            }
        }

        /// <summary>
        /// Activates pivoting strategy for highly indefinite symmetric matrices.
        /// Applies 1x1 and 2x2 Bunch-Kaufman pivoting during factorization.
        /// </summary>
        public bool SymmIndefPivot {
            get {
                return (m_PardInt.m_parm[20] == 3);
            }
            set {
                if (value) {
                    m_PardInt.m_parm[20] = 3; // pivoting for symm indefenite matrices}
                    //m_PardInt.m_parm[9] = 8;
                    m_PardInt.m_parm[10] = 1;
                    m_PardInt.m_parm[12] = 1;
                }
            }
        }

        /// <summary>
        /// From PARDISO manual: Actual matrix for the solution phase. 
        /// With this scalar you can define the matrix that you would like to factorize. 
        /// The value must be: 1 &lt; <see cref="mnum"/> &lt; <see cref="maxfct"/>
        /// in most applications this value is equal to 1. 
        /// </summary>
        public int mnum {
            get { return m_PardInt.mnum; }
            set {
                if (m_PardInt.m_PardisoInitialized && value != m_PardInt.mnum)
                    PARDISODispose();
                m_PardInt.mnum = value;
            }
        }

        /// <summary>
        /// From PARDISO manual: Actual matrix for the solution phase. 
        /// With this scalar you can define the matrix that you would like to factorize. 
        /// The value must be: 1 &lt; <see cref="mnum"/> &lt; <see cref="maxfct"/>
        /// in most applications this value is equal to 1. 
        /// </summary>
        public int maxfct {
            get { return m_PardInt.maxfct; }
            set {
                if (m_PardInt.m_PardisoInitialized && value != m_PardInt.maxfct)
                    PARDISODispose();
                m_PardInt.maxfct = value;
            }
        }
                
        PARDISOinternals m_PardInt = new PARDISOinternals();


        private int GetMType() {
            if (m_PardisoMatrix.Symmetric)
                return -2; // real indefinite symmetric
            else
                return 11; // real unsymmetric
        }


        /// <summary>
        /// Scatters solution vector from process 0 to other MPI processors.
        /// </summary>
        /// <param name="P0__x_IN">
        /// input; long vector on proc 0
        /// </param>
        /// <param name="Px_x_OUT">
        /// output;
        /// </param>
        private void ScatterFromProc0(double[] P0__x_IN, double[] Px_x_OUT) {
            int size, rank;
            csMPI.Raw.Comm_Size(this.MpiComm, out size);
            csMPI.Raw.Comm_Rank(this.MpiComm, out rank);
            Debug.Assert(size == m_OrgMatrix.RowPartitioning.MpiSize);
            Debug.Assert(size == m_OrgMatrix.ColPartition.MpiSize);
            Debug.Assert(rank == m_OrgMatrix.RowPartitioning.MpiRank);
            Debug.Assert(rank == m_OrgMatrix.ColPartition.MpiRank);

            if(rank == 0)
                Debug.Assert(P0__x_IN.Length == m_OrgMatrix.RowPartitioning.TotalLength);
            else
                Debug.Assert(P0__x_IN == null);
            Debug.Assert(Px_x_OUT.Length == m_OrgMatrix.RowPartitioning.LocalLength);
            
            if (size > 1) {
                // distribute solution to other processors
                // +++++++++++++++++++++++++++++++++++++++

                int[] SendCounts = new int[size];
                int[] Displ = new int[size];
                for(int r = 0; r < size; r++) {
                    SendCounts[r] = m_OrgMatrix.RowPartitioning.GetLocalLength(r);
                    Displ[r] = checked((int)(m_OrgMatrix.RowPartitioning.GetI0Offest(r)));
                }


                unsafe
                {
                    fixed(int* pSendCounts = SendCounts, pDispl = Displ) {
                        fixed (double* pSend = P0__x_IN, pRecv = Px_x_OUT) {
                            csMPI.Raw.Scatterv(
                            (IntPtr)pSend,
                            (IntPtr)pSendCounts,
                            (IntPtr)pDispl,
                            csMPI.Raw._DATATYPE.DOUBLE,
                            (IntPtr)pRecv,
                            SendCounts[rank],
                            csMPI.Raw._DATATYPE.DOUBLE,
                            0,
                            this.MpiComm);
                        }
                    }

                }


            } else {
                if (!object.ReferenceEquals(P0__x_IN, Px_x_OUT)) {
                    int L = P0__x_IN.Length;
                    if (P0__x_IN.Length != Px_x_OUT.Length)
                        throw new ApplicationException("internal error.");
                    Array.Copy(P0__x_IN, Px_x_OUT, L);
                }

            }
        }

        /// <summary>
        /// Gathers RHS and solution vector on process 0 for more than one MPI process.
        /// </summary>
        /// <param name="bIN">local part of rhs</param>
        /// <param name="xIN">local part of solution/initial guess</param>
        /// <param name="gath_x_OUT">gathered rhs</param>
        /// <param name="gath_b_OUT">gathered solution vectors</param>
        private void GatherOnProc0( double[] xIN, double[] bIN, out double[] gath_x_OUT, out double[] gath_b_OUT) {
            int size, rank;
            csMPI.Raw.Comm_Size(this.MpiComm, out size);
            csMPI.Raw.Comm_Rank(this.MpiComm, out rank);
            Debug.Assert(size == m_OrgMatrix.RowPartitioning.MpiSize);
            Debug.Assert(size == m_OrgMatrix.ColPartition.MpiSize);
            Debug.Assert(rank == m_OrgMatrix.RowPartitioning.MpiRank);
            Debug.Assert(rank == m_OrgMatrix.ColPartition.MpiRank);


            // This is bad style: 
            // Driver for the driver for the driver...

            if (size > 1) {
                // gather rhs on processor 0
                // +++++++++++++++++++++++++

                int[] RcvCounts = new int[size];
                for(int r = 0; r < size; r++) {
                    RcvCounts[r] = m_OrgMatrix.RowPartitioning.GetLocalLength(r);
                }

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
        /// only executed on proc 0
        /// </summary>
        bool PARDISOInitAndSolve(IMutableMatrixEx Mtx, double[] _x, double[] _b) {
            using (var tr = new FuncTrace()) {
                //InitAndSolve.Start();
                              
                if (m_PardisoMatrix == null)
                    m_PardisoMatrix = new Matrix(Mtx, this.UseDoublePrecision);


                int rank;
                csMPI.Raw.Comm_Rank(MpiComm, out rank);
                if (rank == 0) {
#if DEBUG
                    long aSize = m_PardisoMatrix.ja.LongLength * (UseDoublePrecision ? sizeof(double) : sizeof(float));
                    //double[] a_DClone = null;
                    //float[] a_SClone = null;
                    //if(UseDoublePrecision)
                    //    a_DClone = m_PardisoMatrix.a_D.CloneAs();
                    //else
                    //    a_SClone = m_PardisoMatrix.a_S.CloneAs();
                    IntPtr aClone = Marshal.AllocHGlobal((IntPtr)aSize);
                    unsafe
                    {
                        int* psrc = (int*)m_PardisoMatrix.aPtr;
                        int* pdst = (int*)aClone;
                        for(long l = 0; l < aSize; l+=sizeof(int))
                        {
                            *pdst = *psrc;
                            pdst++;
                            psrc++;
                        }
                    }

                    int[] iaClone = m_PardisoMatrix.ia.CloneAs();
                    int[] jaClone = m_PardisoMatrix.ja.CloneAs();
#endif
                   
                    unsafe {
                        fixed (double* px = _x, pb = _b, dparam = m_PardInt.m_dparam) {
                            fixed (int* ia = m_PardisoMatrix.ia, ja = m_PardisoMatrix.ja, iparm = m_PardInt.m_parm, __pt = m_PardInt.m_pt) {
                                int n = m_PardisoMatrix.n;
                                //int nnz = ia[n];
                                int mtype = GetMType();

                                int nrhs = 1;                       /* Number of right hand sides. */


                                /* Internal solver memory pointer pt,                  */
                                /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
                                /* or void *pt[64] should be OK on both architectures  */
                                void* pt = (void*)__pt;


                                /* Pardiso control parameters. */
                                int maxfct = m_PardInt.maxfct, mnum = m_PardInt.mnum, msglvl = m_PardInt.msglvl;
                                int phase, error;

                                /* Auxiliary variables. */
                                //int      i;

                                double ddum;              /* Double dummy */
                                int idum;              /* Integer dummy. */

                                void* a = (void*)(m_PardisoMatrix.aPtr);
                                double* b, x;
                                if (UseDoublePrecision) {
                                    b = pb;
                                    x = px;
                                } else {

                                    b = (double*)Marshal.AllocHGlobal(n * sizeof(float));
                                    x = (double*)Marshal.AllocHGlobal(n * sizeof(float));

                                    float* bS = (float*)b;
                                    float* xS = (float*)x;

                                    for (int i = 0; i < n; i++) {
                                        bS[i] = (float)(pb[i]);
                                        xS[i] = (float)(px[i]);
                                    }
                                    //SingleCalls.Start();

                                }
                                //Inner.Start();


                                if (!m_PardInt.m_PardisoInitialized) {



                                    /* -------------------------------------------------------------------- */
                                    /* ..  Setup Pardiso control parameters und initialize the solvers      */
                                    /*     internal adress pointers. This is only necessary for the FIRST   */
                                    /*     call of the PARDISO solver.                                      */
                                    /* ---------------------------------------------------------------------*/

                                    if (this.Version == PARDISO.Version.MKL) {
                                        // Do nothing NUM_THREADS is set by environment variables in MKLPardiso
                                        // iparm[2] stays set to zero
                                        // If you want to control this manually, do something like this:
                                        //System.Environment.SetEnvironmentVariable("OMP_NUM_THREADS", "4");
                                        //System.Environment.SetEnvironmentVariable("MKL_NUM_THREADS", "4");
                                        //If this is true: MKL determines the number of OMP threads automatically, based on proc information
                                        //System.Environment.SetEnvironmentVariable("MKL_DYNAMIC", "false");

                                    } else if (this.Version == PARDISO.Version.v5) {
                                        // Get Value for IPARM(3) from the Environment Variable
                                        int NumOfProcs = Convert.ToInt32(System.Environment.GetEnvironmentVariable("OMP_NUM_THREADS"));
                                        iparm[2] = NumOfProcs;
#if DEBUG
                                            Console.WriteLine("IPARM(3) - NumberOfProcessors set to {0}", iparm[2]);
#endif

                                    }

                                    using (new BlockTrace("PARDISOINIT", tr)) {
                                        wrapper.PARDISOINIT(pt, &mtype, iparm, dparam);
                                    }

                                    //Console.WriteLine("init: IPARAM(22) = {0}, IPARAM(23) = {1}", iparm[21], iparm[22]);

                                    iparm[27] = this.UseDoublePrecision ? 0 : 1; // set single or double precision



                                    maxfct = 1;         /* Maximum number of numerical factorizations.  */
                                    mnum = 1;         /* Which factorization to use. */

                                    msglvl = 0;         /* Print statistical information  */
                                    error = 0;         /* Initialize error flag */


                                    /* -------------------------------------------------------------------- */
                                    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
                                    /*     all memory that is necessary for the factorization.              */
                                    /* -------------------------------------------------------------------- */
                                    using (new BlockTrace("PARDISO_phase11", tr)) {
                                        phase = 11;
                                        iparm[59] = 0; // in-core (1 == out-of-core)                            

                                        //Console.Write("calling pardiso, phase 11... ");
                                        Phase_11.Start();
                                        wrapper.PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                                                          &n, a, ia, ja, &idum, &nrhs,
                                                          iparm, &msglvl, &ddum, &ddum, &error, dparam);
                                        Phase_11.Stop();
                                        SetMemoryConsumption();
                                        //Console.WriteLine("11: IPARAM(22) = {0}, IPARAM(23) = {1}", iparm[21], iparm[22]);
                                    }
                                    if (error != 0) {
                                        PARDISODispose();
                                        SendErrorToOtherRanks(error);
                                        Console.Error.WriteLine("PARDISO ERROR: " + wrapper.PARDISOerror2string(error));
                                        return false;
                                    }
                                    //Console.Write("\nReordering completed ... ");
                                    //Console.Write("\nNumber of nonzeros in factors  = %d", iparm[17]);
                                    //Console.Write("\nNumber of factorization MFLOPS = %d", iparm[18]);

                                    /* -------------------------------------------------------------------- */
                                    /* ..  Numerical factorization.                                         */
                                    /* -------------------------------------------------------------------- */
                                    using (new BlockTrace("PARDISO_phase22", tr)) {
                                        phase = 22;

                                        Phase_22.Start();
                                        wrapper.PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                                                          &n, a, ia, ja, &idum, &nrhs,
                                                          iparm, &msglvl, &ddum, &ddum, &error, dparam);
                                        Phase_22.Stop();
                                    }

                                    if (error != 0) {
                                        // some error occured: release mem, dispose objects...
                                        PARDISODispose();
                                        SendErrorToOtherRanks(error);
                                        Console.Error.WriteLine("PARDISO ERROR: " + wrapper.PARDISOerror2string(error));
                                        //InitAndSolve.Stop();
                                        return false;
                                    }
                                }

                                /* -------------------------------------------------------------------- */
                                /* ..  Back substitution and iterative refinement.                      */
                                /* -------------------------------------------------------------------- */
                                phase = 33;

                                iparm[7] = 0;       /* Max numbers of iterative refinement steps, 0 == 2 steps */

                                //m_foo.mkl_serv_mkl_set_num_threads(num_procs);
                                using (new BlockTrace("PARDISO_phase33", tr)) {
                                    Phase_33.Start();
                                    wrapper.PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                                                      &n, a, ia, ja, &idum, &nrhs,
                                                      iparm, &msglvl, b, x, &error, dparam);
                                    Phase_33.Stop();
                                }
                                if (error != 0) {
                                    // some error occurred: release mem, dispose objects...
                                    PARDISODispose();
                                    SendErrorToOtherRanks(error);
                                    Console.Error.WriteLine("PARDISO ERROR: " + wrapper.PARDISOerror2string(error));
                                    //InitAndSolve.Stop();
                                    return false;
                                }

                                //Inner.Stop();
                                if (UseDoublePrecision) {

                                } else {
                                    //SingleCalls.Stop();
                                    float* bS = (float*)b;
                                    float* xS = (float*)x;
                                    for (int i = 0; i < n; i++) {
                                        pb[i] = (bS[i]);
                                        px[i] = (xS[i]);
                                    }

                                    Marshal.FreeHGlobal((IntPtr)b);
                                    Marshal.FreeHGlobal((IntPtr)x);
                                }
                            }
                        }
                    }



#if DEBUG
                    //if(UseDoublePrecision)
                    //    Debug.Assert(ArrayTools.ListEquals<double>(a_DClone, m_PardisoMatrix.a_D), "PARDISO changed the matrix.");
                    //else
                    //    Debug.Assert(ArrayTools.ListEquals<float>(a_SClone, m_PardisoMatrix.a_S), "PARDISO changed the matrix.");

                    unsafe
                    {
                        int* psrc = (int*)m_PardisoMatrix.aPtr;
                        int* pdst = (int*)aClone;
                        for (long l = 0; l < aSize; l += sizeof(int))
                        {
                            Debug.Assert( *pdst == *psrc, "PARDISO changed the matrix.");
                            pdst++;
                            psrc++;
                        }
                    }


                    Debug.Assert(ArrayTools.ListEquals<int>(iaClone, m_PardisoMatrix.ia), "PARDISO changed the matrix.");
                    Debug.Assert(ArrayTools.ListEquals<int>(jaClone, m_PardisoMatrix.ja), "PARDISO changed the matrix.");
#endif
                    SendErrorToOtherRanks(0); // Everything is fine
                } else {
                    /*
                     procs != 0 are waiting for the error code
                     if error != 0 an error is thrown on all ranks of MPI_Comm (error code is distributed by rank0)
                     see SendErrorToOtherRanks
                     This makes pardiso exceptions catchable
                    */
                    int[] error = new int[] { -1 };
                    unsafe {
                        fixed (int* Buffer = error) {
                            MPI_Status status;
                            MPI.Wrappers.csMPI.Raw.Recv((IntPtr)Buffer, 1, csMPI.Raw._DATATYPE.INT, 0, 112, MpiComm, out status);
                        }
                    }
                    if (error[0] != 0)
                        wrapper.PARDISOerror2string(error[0]);
                }


                m_PardInt.m_PardisoInitialized = true;
                //InitAndSolve.Stop();
                return true;
            }
        }

        // Broadcasts the error code over MPI_Comm.
        // Receive is at the end of PardisoInitAndSolver
        private void SendErrorToOtherRanks(int error) {
            csMPI.Raw.Comm_Size(MpiComm, out int size);
            for (int iProc = 1; iProc < size; iProc++) {
                int[] buffer = new int[] { error };
                unsafe {
                    fixed (int* Buffer = buffer) {
                        csMPI.Raw.Send((IntPtr)Buffer, 1, csMPI.Raw._DATATYPE.INT, iProc, 112, MpiComm);
                    }
                }
            }
        }

        int m_MpiRank = -1;

        int m_PeakmemoryFac = -1;
        int m_PeakmemorySol = -1;
        int m_PermanentMemoryFac = -1;

        void SetMemoryConsumption() {
            m_PeakmemoryFac = Math.Max(m_PardInt.m_parm[14], m_PeakmemoryFac);
            m_PermanentMemoryFac = Math.Max(m_PardInt.m_parm[15], m_PermanentMemoryFac);
            m_PeakmemorySol = Math.Max(m_PardInt.m_parm[16], m_PeakmemorySol);
        }

        /// <summary>
        /// Used memory for factorization in bytes
        /// </summary>
        public long UsedMemory() {
            return m_PermanentMemoryFac*1024;
        }

        /// <summary>
        /// Returns the Peak Memory in MB
        /// </summary>
        /// <returns></returns>
        public int PeakMemory() {
            return Convert.ToInt32(Math.Max(m_PeakmemoryFac, m_PeakmemorySol + m_PermanentMemoryFac)/1024);
        }

        /// <summary>
        /// disposal of things cached by Pardiso
        /// </summary>
        private void PARDISODispose() {
           //if (!m_PardInt.m_PardisoInitialized)
           //    throw new ApplicationException("Pardiso is not initialized");

            if (m_MpiRank == 0) {
                if (m_PardisoMatrix != null) {
                    unsafe {
                        fixed (double* dparam = m_PardInt.m_dparam) {
                            fixed (int* ia = m_PardisoMatrix.ia, ja = m_PardisoMatrix.ja, iparm = m_PardInt.m_parm, __pt = m_PardInt.m_pt) {


                                /* Internal solver memory pointer pt,                  */
                                /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
                                /* or void *pt[64] should be OK on both architectures  */
                                void* pt = (void*)__pt;
                                int n = m_PardisoMatrix.n;
                                int maxfct = m_PardInt.maxfct, mnum = m_PardInt.mnum, msglvl = m_PardInt.msglvl;
                                int mtype = GetMType();
                                int phase, error;

                                double ddum; // dummy
                                int idum;              /* Integer dummy. */

                                int nrhs = 1;          /* Number of right hand sides. */


                                /* -------------------------------------------------------------------- */
                                /* ..  Termination and release of memory.                               */
                                /* -------------------------------------------------------------------- */
                                phase = -1;                 /* Release internal memory. */

                                //if (m_wrapper.PARDISO == null)
                                //    m_wrapper = new Wrapper_MKL(); // maybe the library has already been unloaded by the GC

                                wrapper.PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                                                  &n, &ddum, ia, ja, &idum, &nrhs,
                                                  iparm, &msglvl, &ddum, &ddum, &error, dparam);
                            }
                        }
                    }

                }

                Array.Clear(m_PardInt.m_parm, 0, m_PardInt.m_parm.Length);
                Array.Clear(m_PardInt.m_pt, 0, m_PardInt.m_pt.Length);
            }
            m_PardisoMatrix.Dispose();
            m_PardisoMatrix = null;
            m_PardInt.m_PardisoInitialized = false;
        }
        
        /// <summary>
        /// <see cref="CacheFactorization"/>
        /// </summary>
        bool m_CacheFactorization = true;

        /// <summary>
        /// If the solver should be used more than once (with the same matrix),
        /// setting this variable to true stores the matrix factorization;
        /// This occupies a large bit of memory but is strongly recommended,
        /// it may speed up subsequent calls to <see cref="Solve{Tunknowns, Trhs}(Tunknowns, Trhs)"/>
        /// by a factor of 10 or higher!
        /// <br/>
        /// default is true;
        /// </summary>
        public bool CacheFactorization {
            get { return m_CacheFactorization; }
            set {
                if (value == false && m_PardInt.m_PardisoInitialized)
                    PARDISODispose();
                m_CacheFactorization = value;
            }
        }

        #region IDisposable Members

        /// <summary>
        /// releases unmanaged resources (internal memory of PARDISO)
        /// </summary>
        public void Dispose() {
            if (m_PardInt.m_PardisoInitialized)
                PARDISODispose();
            m_OrgMatrix = null;
        }

        #endregion

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public ISparseMatrix GetMatrix() {
            return m_OrgMatrix;
        }
    }
}
