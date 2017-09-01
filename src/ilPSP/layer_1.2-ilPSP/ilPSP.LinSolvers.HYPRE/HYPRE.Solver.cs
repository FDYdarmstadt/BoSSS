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
using System.IO;
using MPI.Wrappers;
using System.Diagnostics;
using ilPSP.Tracing;

namespace ilPSP.LinSolvers.HYPRE {

    
    /// <summary>
    /// baseclass for all HYPRE solvers
    /// </summary>
    public abstract class Solver : ISparseSolver, ISparseSolverExt, IDisposable {

        /// <summary>
        /// If the solver should be used as a preconditioner (e.g. like <see cref="Euclid"/> or <see cref="BoomerAMG"/>)
        /// this member must be initialized as the function address
        /// of the HYPRE-solver function
        /// </summary>
        internal IntPtr m_NativeSolverFuncPtr;

        /// <summary>
        /// If the solver should be used as a preconditioner (e.g. like <see cref="Euclid"/> or <see cref="BoomerAMG"/>)
        /// this member must be initialized as the function address
        /// of the HYPRE-solver function
        /// </summary>
        internal IntPtr m_NativeSetupFuncPtr;
        
        /// <summary>
        /// creates a HYPRE solver
        /// </summary>
        public Solver() {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            CreateSolver();
        }

        /// <summary>
        /// calls <see cref="Dispose"/>;
        /// </summary>
        ~Solver() {
            Dispose();
        }

        /// <summary>
        /// pointer/handle to the solver object
        /// </summary>
        internal Wrappers.T_Solver m_Solver;


        /// <summary>
        /// allocates and initializes matrix, and allocates
        /// memory for unknowns and right hand side.
        /// </summary>
        /// <param name="M"></param>
        public void DefineMatrix(IMutableMatrixEx M) {
            if (m_Matrix != null)
                throw new ApplicationException("matrix is allready defined. 'DefineMatrix'-method can be invoked only once in the lifetime of this object.");

            m_Matrix = new IJMatrix(M);
            //m_Unknowns = new IJVector(M.RowPartiton);
            //m_Rhs = new IJVector(M.RowPartiton);
            //if (m_MatrixToFile != null)
            //    DumpMatrix(M);
        }

        /// <summary>
        /// <see cref="ISparseSolverExt.GetMatrix"/>
        /// </summary>
        /// <returns></returns>
        public ISparseMatrix GetMatrix() {
            return m_Matrix;
        }

        /// <summary>
        /// override to implement the solver init;
        /// called by Constructor
        /// </summary>
        protected abstract void CreateSolver();


        /*
        public string m_MatrixToFile;


        void DumpMatrix(MsrMatrix M) {
            StreamWriter fs = new StreamWriter(m_MatrixToFile);

            fs.Write(M.RowMap.GlobalCount);
            fs.Write(" ");
            fs.Write(M.ColMap.GlobalCount);
            fs.WriteLine();

            fs.Write(M.RowMap.i0Offset);
            fs.Write(" ");
            fs.Write(M.RowMap.LocalLength);
            fs.Write(" ");
            fs.Write(M.ColMap.i0Offset);
            fs.Write(" ");
            fs.Write(M.ColMap.LocalLength);
            fs.Write(" ");
            fs.WriteLine();

            fs.WriteLine(M.GetMaxNoOfOffDiagonalNonZerosPerRow() + 1);


            for (int i = 0; i < M.RowMap.LocalLength; i++) {
                long iRow = M.RowMap.Local2GlobalIndex(i);
                fs.Write(iRow);
                fs.Write(" ");

                MsrMatrix.MSREntry[] row = M.GetRow(i);
                fs.Write(row.Length);
                fs.Write(" ");
                for (int j = 0; j < row.Length; j++) {
                    long jCol = M.ColMap.Local2GlobalIndex(row[j].ColIndex);
                    fs.Write(jCol);
                    fs.Write(" ");
                    fs.Write(row[j].val.ToString(System.Globalization.NumberFormatInfo.InvariantInfo));
                    fs.Write(" ");
                }
                fs.WriteLine();




            }




            fs.Flush();
            fs.Close();

        }


        */



        /// <summary>
        /// the 
        /// </summary>
        internal  IJMatrix m_Matrix;

        ///// <summary>
        ///// memory for unknowns/solution
        ///// </summary>
        //public IJVector m_Unknowns;

        ///// <summary>
        ///// Right-hand-side of the equation system.
        ///// </summary>
        //internal  IJVector m_Rhs;
        
        ///// <summary>
        ///// get/set relative tolerance (convergence criterion)
        ///// </summary>
        //abstract public double Tolerance { get; set; }

        ///// <summary>
        ///// get/set maximum number of solver iterations
        ///// </summary>
        //abstract public int MaxIterations { get; set; }
        
        /// <summary>
        /// solves the equation system 
        /// M*<paramref name="x"/>=<paramref name="rhs"/>,
        /// where M denotes the matrix of the equation system;
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
        /// where diag(<paramref name="d"/>) denotes a diagonal matrix with the diagonal vector <paramref name="d"/>;
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
            where Tdiag : IList<double>
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
                using (new FuncTrace()) {

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


                    // scale rhs
                    // =========

                    if (ooScale != 1.0) {
                        for (int i = 0; i < N; i++)
                            rhs[i] *= ooScale;
                    }


                    // pass values to HYPRE
                    // ====================
                    IJVector Unknowns = new IJVector(m_Matrix.ColPartition);
                    IJVector Rhs = new IJVector(m_Matrix.RowPartitioning);

                    Unknowns.SetValues<Tunknowns>(x);
                    Rhs.SetValues<Trhs>(rhs);

                    CallSolver(out res.NoOfIterations, out res.Converged, Unknowns, Rhs);


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



                    Unknowns.GetValues<Tunknowns>(x);

                    Unknowns.Dispose();
                    Rhs.Dispose();

                    st.Stop();
                    res.RunTime = st.Elapsed;
                    return res;
                }
        }
        
        /// <summary>
        /// calls the corresponding HYPRE solver function
        /// </summary>
        /// <param name="Converged">true if converged</param>
        /// <param name="NoOfIter">no of iterations done by solver</param>
        /// <param name="Unknowns">on input, the starting value for the iterative process;
        /// on output, an approximpate solution to the linear problem.</param>
        /// <param name="Rhs">
        /// on input, the right-hand-side of the problem.
        /// </param>
        protected abstract void CallSolver(out int NoOfIter, out bool Converged, IJVector Unknowns, IJVector Rhs);
        
        /// <summary>
        /// frees internal memory (Matrix, RHS and solution vector);
        /// </summary>
        public virtual void Dispose() {
            if (m_Matrix != null) {
                m_Matrix.Dispose();
                //m_Rhs.Dispose();
                //m_Unknowns.Dispose();

                m_Matrix = null;
                //m_Rhs = null;
                //m_Unknowns = null;
            }
        }
        
        /// <summary>
        /// <see cref="ConvergenceType"/>
        /// </summary>
        protected ConvergenceTypes _ConvType;
        
        /// <summary>
        /// Choose between absolute and relative residual convergence criterion;
        /// Remark: <see cref="ilPSP.LinSolvers.ConvergenceTypes.Other"/> is not
        /// accepted by set and results in an exception.
        /// </summary>
        public ConvergenceTypes ConvergenceType {
            get {
                return _ConvType;
            }
            set {
                if (value == ConvergenceTypes.Other)
                    throw new ArgumentException("ConvergenceTypes.Other is not accepted as an input");
                _ConvType = value;
            }
        }


    }
}
