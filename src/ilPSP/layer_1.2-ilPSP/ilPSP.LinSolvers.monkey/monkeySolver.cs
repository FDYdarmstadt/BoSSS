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
using MPI.Wrappers;
using ilPSP.Tracing;
using log4net;
using System.Diagnostics;

namespace ilPSP.LinSolvers.monkey {

    
    /// <summary>
    /// Base class for all monkey solvers;
    /// </summary>
    abstract public class Solver : ISparseSolverExt {

        static ILog Logger = LogManager.GetLogger(typeof(Solver));
        
        /// <summary>
        /// constructor
        /// </summary>
        protected Solver() {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
        }

        DeviceType m_DevType = DeviceType.Auto;

        /// <summary>
        /// encodes the driver (subclass of <see cref="Device"/>) that should be used in this
        /// solver; default is <see cref="DeviceType.Auto"/>
        /// </summary>
        public DeviceType DevType {
            get { return m_DevType; }
            set {
                if (m_Device != null)
                    throw new ApplicationException("unable to change the Device after first use of solver (Device already created).");
                m_DevType = value;
            }
        }

        MatrixType m_MatrixType = MatrixType.Auto;

        /// <summary>
        /// matrix storage format, ignored if running on CPU, but effective when running on CUDA
        /// </summary>
        public MatrixType MatrixType {
            get { return m_MatrixType; }
            set {
                if (m_Matrix != null)
                    throw new ApplicationException("unable to change the MatrixType after matrix is defined (DefineMatrix(...) allready called).");
                m_MatrixType = value;
            }
        }
        
        #region ISparseSolver Members

        /// <summary>
        /// cache for different device types (we want e.g. all CUDA -- solvers to use the same device)
        /// </summary>
        /// <remarks>
        /// this is static to ensure that only one devise per application will be created; Therefore, only one CUDA context will be created.
        /// </remarks>
        static SortedDictionary<DeviceType, Device> m_DeviceS = new SortedDictionary<DeviceType, Device>();

        static Device GetOrCereateDevice(DeviceType DevType) {
            Device Dev = null;
            m_DeviceS.TryGetValue(DevType, out Dev);
            if (Dev == null) {
                switch(DevType) {
                    case DeviceType.Cuda: Dev = new CUDA.CudaDevice(new CUDA.CudaEnviroment(Environment.MPIEnv)); break;
                    case DeviceType.OpenCL: Dev = new CL.clDevice(new CL.clEnvironment(Environment.MPIEnv)); break;
                    case DeviceType.CPU: Dev = new CPU.ReferenceDevice(); break;
                    case DeviceType.MultiThreadCPU: Dev = new mtCPU.MtDevice(); break;
                    case DeviceType.Auto: {
                        // try cuda at frist:
                        try {
                            Dev = GetOrCereateDevice(DeviceType.Cuda);
                        } catch (Exception) {
                            Dev = null;
                        }
                        if (Dev != null) break;

                        // try OpenCL next:
                        //try {
                        //    Dev = GetOrCereateDevice(DeviceType.OpenCL);
                        //} catch (Exception) {
                        //    Dev = null;
                        //}
                        //if (Dev != null) break;

                        // fall back to CPU:
                        Dev = GetOrCereateDevice(DeviceType.CPU);

                        break;
                        }
                    default:
                        throw new NotImplementedException("monkey device type: " + DevType.ToString() + " missing in factory.");
                }
                m_DeviceS.Add(DevType, Dev);
            }
            return Dev;
        }


        /// <summary>
        /// see <see cref="Device"/>;
        /// </summary>
        Device m_Device;

        /// <summary>
        /// the factory which will create vector/matrix objects that either live on main CPU or on GPU
        /// </summary>
        protected Device Device {
            get {
                if (m_Device == null) {
                    m_Device = GetOrCereateDevice(this.DevType);
                }
                return m_Device;
            }
        }

        /// <summary>
        /// common interface (for all main memory/accelerator devices) to the solver matrix.
        /// </summary>
        protected MatrixBase m_Matrix;

        /// <summary>
        /// sets the matrix of the solver
        /// </summary>
        virtual public void DefineMatrix(IMutableMatrixEx _M) {
            using (new FuncTrace()) {
                MsrMatrix M = _M as MsrMatrix;
                if (M == null)
                    // provisorium, geht sicher besser
                    M = _M.ToMsrMatrix();


                if (M.RowPartitioning.TotalLength != M.NoOfCols)
                    throw new ArgumentException("Matrix has to be quadratic.", "_M");
                m_Matrix = Device.CreateMatrix(M, m_MatrixType);
            }
        }

        /// <summary>
        /// see <see cref="ISparseSolverExt.GetMatrix"/>
        /// </summary>
        /// <returns></returns>
        public ISparseMatrix GetMatrix() {
            return m_Matrix;
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
                m_Tolerance = value;
            }
        }

        /// <summary>
        /// implementation of <see cref="ConvergenceType"/>;
        /// </summary>
        protected ConvergenceTypes m_ConvergenceType = ConvergenceTypes.Relative;

        /// <summary>
        /// different types of residual computation, 
        /// see <see cref="ConvergenceTypes"/>;
        /// </summary>
        public ConvergenceTypes ConvergenceType {
            get {
                return m_ConvergenceType;
            }
            set {
                m_ConvergenceType = value;
            }
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
                m_MaxIterations = value;
            }
        }

        /// <summary>
        /// implementation of <see cref="MaxIterations"/>
        /// </summary>
        protected int m_MinIterations = 0;

        /// <summary>
        /// minimum number of iterations that should be performed; default is 0;<br/>
        /// The only case when this is usefull is maybe when a solver is used as preconditioner,
        /// where a fixed number of iterations should be performed, i.e. <see cref="MinIterations"/>== <see cref="MaxIterations"/>.
        /// </summary>
        public int MinIterations {
            get {
                return m_MinIterations;
            }
            set {
                m_MinIterations = value;
            }
        }

        /// <summary>
        /// internal implementation of the solver
        /// </summary>
        /// <param name="x">on exit, (hopefully) the solution to the equation</param>
        /// <param name="rhs">right-hand-side</param>
        /// <param name="stats">
        /// an implementor should at least set <see cref="SolverResult.Converged"/>
        /// and <see cref="SolverResult.NoOfIterations"/>
        /// </param>
        abstract protected void CallSolver(VectorBase x, VectorBase rhs, ref SolverResult stats);
        
        /// <summary>
        /// see <see cref="ISparseSolver.Solve{Tunknowns, Trhs}(Tunknowns,Trhs)"/>;
        /// </summary>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            return Solve<double[], Tunknowns, Trhs>(1.0, null, x, rhs);
        }

        /// <summary>
        /// see <see cref="ISparseSolverExt.Solve{Tdiag, Tunknowns, Trhs}(double,Tdiag,Tunknowns,Trhs)"/>;
        /// </summary>
        public SolverResult Solve<Tdiag, Tunknowns, Trhs>(double Scale, Tdiag d, Tunknowns x, Trhs rhs)
            where Tdiag : IList<double>
            where Tunknowns : IList<double>
            where Trhs : IList<double> {

            using (var tr = new ilPSP.Tracing.FuncTrace()) {

                SolverResult res = new SolverResult();

                Stopwatch st = new Stopwatch();
                st.Reset();
                st.Start();

                // modify diagonal
                // ===============

                // truly, we're not solving (diag(d) + Scale*M)*x = rhs,
                // but ((1.0/Scale)*diag(d) + M) = (1.0/Scale)*rhs

                double ooScale = 1.0 / Scale;
                int N = int.MinValue;
                if (d != null) {

                    int i0 = (int)m_Matrix.RowPartitioning.i0;

                    N = m_Matrix.RowPartitioning.LocalLength;
                    int Nd = d.Count;
                    if (d.Count > N || N % Nd != 0) {
                        throw new ArgumentException("length must be equal to or a factor of the number of rows stored on this processor", "d");
                    }

                    int ix = 0;
                    for (int i = 0; i < N; i++) {
                        double vadd = d[ix];
                        ix++;
                        if (ix >= Nd) ix = 0;

                        if (vadd != 0.0) {
                            int iglob = i + i0;
                            double v = m_Matrix.GetDiagonalElement(iglob);
                            v += ooScale * vadd;
                            m_Matrix.SetDiagonalElement(iglob, v);
                        }
                    }
                }

                // pass values to monkey
                // =====================
                bool shallow, dummy2;
                VectorBase X = Device.CreateVector(m_Matrix.ColPartition, x, out shallow);
                VectorBase Rhs = Device.CreateVector(m_Matrix.RowPartitioning, rhs, out dummy2);

                // scale rhs
                // =========

                if (ooScale != 1.0) {
                    Rhs.Lock();
                    Rhs.Scale(ooScale);
                    Rhs.Unlock();
                }

                // call Solver
                // ===========

                CallSolver(X, Rhs, ref res);

                if(res.Converged != true)
                    Logger.Warn("Solver did NOT CONVERGE: " + res.ToString());


                // return
                // ======

                if (d != null) {
                    int ix = 0;
                    int Nd = d.Count;

                    int i0 = (int)m_Matrix.RowPartitioning.i0;

                    for (int i = 0; i < N; i++) {
                        double vadd = d[ix];
                        ix++;
                        if (ix >= Nd) ix = 0;

                        if (vadd != 0.0) {
                            int iglob = i + i0;
                            double v = m_Matrix.GetDiagonalElement(iglob);
                            v -= ooScale * vadd;
                            m_Matrix.SetDiagonalElement(iglob, v);
                        }
                    }
                }

                if (!shallow)
                    X.GetValues(x, 0, 0, m_Matrix.ColPartition.LocalLength);
                X.Dispose();
                Rhs.Dispose();

                st.Stop();
                res.RunTime = st.Elapsed;
                return res;

            }
        }

        #endregion

        #region IDisposable Members

        /// <summary>
        /// %
        /// </summary>
        virtual public void Dispose() {
            m_Matrix.Dispose();
        }


        #endregion
    }
}
