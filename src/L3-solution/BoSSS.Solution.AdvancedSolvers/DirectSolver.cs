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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.LinSolvers.MUMPS;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Platform;
using System.IO;
using ilPSP.Tracing;
using ilPSP.Connectors.Matlab;
using System.Diagnostics;
using ilPSP.LinSolvers.monkey;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Sparse direct solver. 
    /// This class is a wrapper around either 
    /// - PARDISO (<see cref="PARDISOSolver"/>) or 
    /// - MUMPS (<see cref="MUMPSSolver"/>) or
    /// - LAPACK.
    /// </summary>
    public class DirectSolver : ISubsystemSolver, ISolverWithCallback {

        /// <summary>
        /// 
        /// </summary>
        [Serializable]
        public class Config : ISolverFactory {

            /// <summary>
            /// Switch between PARDISO and MUMPS.
            /// </summary>
            public _whichSolver WhichSolver = _whichSolver.PARDISO;

            /// <summary>
            /// 
            /// </summary>
            public string Name => "Sparse direct solver " + WhichSolver.ToString();

            /// <summary>
            /// 
            /// </summary>
            public string Shortname => WhichSolver.ToString();

            /// <summary>
            /// 
            /// </summary>
            public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new DirectSolver();
                instance.m_config = this;
                instance.Init(level);
                return instance;
            }

            /// <summary>
            /// 
            /// </summary>
            public bool Equals(ISolverFactory other) {
                return EqualsImpl(other);
            }

            /// <summary>
            /// 
            /// </summary>
            public override bool Equals(object obj) {
                return EqualsImpl(obj);
            }


            public override int GetHashCode() {
                return (int)(this.WhichSolver);
            }

            private bool EqualsImpl(object o) {
                var other = o as Config;
                if (other == null)
                    return false;

                return (this.WhichSolver == other.WhichSolver);
            }


            /// <summary>
            /// If set to true, the solution returned by the direct solver is tested by computing the residual norm.
            /// Currently, the default is true, since the direct solvers seem unreliable.
            /// </summary>
            public bool TestSolution {
                get;
                set;
            } = true;

            /// <summary>
            /// - true: use double precision floats internally (default)
            /// - false: use single precision; might be useful for solution guesses in preconditioners; only supported for <see cref="WhichSolver"/>==<see cref="_whichSolver.PARDISO"/>.
            /// </summary>
            public bool UseDoublePrecision {
                get;
                set;
            } = true;
        }


        Config m_config = new Config();

        /// <summary>
        /// Solver configuration
        /// </summary>
        public Config config {
            get {
                return m_config;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public enum _whichSolver {

            /// <summary>
            /// Using <see cref="PARDISOSolver"/>.
            /// </summary>
            PARDISO,

            /// <summary>
            /// Using <see cref="MUMPSSolver"/>.
            /// </summary>
            MUMPS,

            /// <summary>
            /// Conversion to dense matrix, solution 
            /// via LU-decomposition from LAPACK, see also <see cref="IMatrixExtensions.Solve{T, W}(T, W)"/>.
            /// Only suitable for small systems (less than 10000 DOF).
            /// </summary>
            Lapack,

            /// <summary>
            /// MATLAB 'backslash' solver, see <see cref="ilPSP.Connectors.Matlab.Extensions.SolveMATLAB{T1, T2}(IMutableMatrixEx, T1, T2, string)"/> 
            /// </summary>
            Matlab,

        }

        void InitImpl(IOperatorMappingPair op) {
            using(var tr = new FuncTrace()) {
                if(object.ReferenceEquals(op, MatrixNMapping))
                    return; // already initialized
                else
                    this.Dispose();

                var Mtx = op.OperatorMatrix;
                var MgMap = op.DgMapping;
                MatrixNMapping = op;

                if(!Mtx.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if(!Mtx.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");

                Mtx.CheckForNanOrInfM(typeof(DirectSolver) + ", matrix definition: ");

                if(m_Solver != null) {
                    m_Solver.Dispose();
                    m_Solver = null;
                }

                m_Mtx = Mtx;
            }
        }
    

        IOperatorMappingPair MatrixNMapping;

        class MatlabSolverWrapper : ISparseSolver {

            MsrMatrix Mtx;

            public void DefineMatrix(IMutableMatrixEx M) {
                Mtx = M.ToMsrMatrix();
            }

            public void Dispose() {
                Mtx = null;
            }

            public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
                where Tunknowns : IList<double>
                where Trhs : IList<double> {
                var StartTime = DateTime.Now;

                Mtx.SolveMATLAB(x, rhs);

                return new SolverResult() {
                    Converged = true,
                    NoOfIterations = 1,
                    RunTime = DateTime.Now - StartTime
                };
            }
        }

        class DenseSolverWrapper : ISparseSolver {

            MultidimensionalArray FullMatrix;

            public void DefineMatrix(IMutableMatrixEx M) {
                FullMatrix = M.ToFullMatrixOnProc0();
            }

            public void Dispose() {
                FullMatrix = null;
            }

            public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
                where Tunknowns : IList<double>
                where Trhs : IList<double> {
                var StartTime = DateTime.Now;

                double[] Int_x = x as double[];
                bool writeBack = false;
                if(Int_x == null) {
                    Int_x = new double[x.Count];
                    writeBack = true;
                }

                double[] Int_rhs = rhs as double[];
                if(Int_rhs == null)
                    Int_rhs = rhs.ToArray();

                FullMatrix.Solve(Int_x, Int_rhs);

                if(writeBack)
                    x.SetV(Int_x);

                return new SolverResult() {
                    Converged = true,
                    NoOfIterations = 1,
                    RunTime = DateTime.Now - StartTime
                };
            }
        }



        ISparseSolver GetSolver(IMutableMatrixEx Mtx) {
            ISparseSolver solver;

            bool RunSerial = Mtx.MPI_Comm == csMPI.Raw._COMM.SELF;

            if(config.WhichSolver != _whichSolver.PARDISO) {
                if (!config.UseDoublePrecision)
                    throw new NotSupportedException($"{config.WhichSolver} is not supported in single precision");
            }


            switch (config.WhichSolver) {
                case _whichSolver.PARDISO:
                bool CachingOn = false;
                if (ActivateCaching != null) {
                    CachingOn = ActivateCaching.Invoke(m_ThisLevelIterations, (MatrixNMapping as MultigridOperator)?.LevelIndex ?? -1);
                }

                    
                solver = new PARDISOSolver() {
                    CacheFactorization = CachingOn,
                    UseDoublePrecision = config.UseDoublePrecision,
                    Parallelism = RunSerial ? Parallelism.SEQ : Parallelism.OMP
                };
                break;

                case _whichSolver.MUMPS:
                solver = new MUMPSSolver() {
                    Parallelism = RunSerial ? Parallelism.SEQ : Parallelism.MPI
                };
                break;

                case _whichSolver.Matlab:
                solver = new MatlabSolverWrapper();
                break;

                case _whichSolver.Lapack:
                solver = new DenseSolverWrapper();
                break;

                default:
                throw new NotImplementedException();

            }

            solver.DefineMatrix(Mtx);

            return solver;
        }

        BlockMsrMatrix m_Mtx;
        int IterCnt = 1;


        long m_UsedMemoryInLastCall = 0; 

        /// <summary>
        /// %
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using(var tr = new FuncTrace()) {
                B.CheckForNanOrInfV(true, true, true, typeof(DirectSolver).Name + ", RHS on entry: ");
                
                double[] Residual = this.config.TestSolution ? B.ToArray() : null;

                string SolverName = "NotSet";
                

                {
                    if(m_Solver == null)
                        m_Solver = GetSolver(m_Mtx);
                    SolverName = m_Solver.GetType().FullName;
                    //Console.Write("Direct solver run {0}, using {1} ... ", IterCnt, solver.GetType().Name);
                    IterCnt++;
                    m_Solver.Solve(X, B);
                    this.Converged = true; 
                    m_ThisLevelIterations++;
                    //Console.WriteLine("done.");

                    if(m_Solver is PARDISOSolver pslv) {
                        m_UsedMemoryInLastCall = pslv.UsedMemory();
                    }
                }

                X.CheckForNanOrInfV(true, true, true, typeof(DirectSolver).Name + ", solution after solver call: ");


                if(Residual != null) {
                    double MatrixInfNorm = m_Mtx.InfNorm();

                    double RhsNorm_loc = Residual.L2NormPow2();

                    m_Mtx.SpMV(-1.0, X, 1.0, Residual);

                    double[] normsGlobPow2 = (new[] { RhsNorm_loc, Residual.L2NormPow2(), B.L2NormPow2(), X.L2NormPow2() }).MPISum(m_Mtx.MPI_Comm);
                    double RhsNorm = normsGlobPow2[0].Sqrt();
                    double ResidualNorm = normsGlobPow2[1].Sqrt();
                    double SolutionNorm = normsGlobPow2[2].Sqrt();
                    double Denom = Math.Max(MatrixInfNorm, Math.Max(RhsNorm, Math.Max(SolutionNorm, Math.Sqrt(BLAS.MachineEps))));
                    double RelResidualNorm = ResidualNorm / Denom;

                    //Console.WriteLine("done: Abs.: {0}, Rel.: {1}", ResidualNorm, RelResidualNorm);

                    if(RelResidualNorm > 1.0e-10) {

                        //Console.WriteLine("High residual from direct solver: abs {0}, rel {1}", ResidualNorm , ResidualNorm / SolutionNorm);
#if TEST
                        m_Mtx.SaveToTextFileSparse("Mtx.txt");
                        X.SaveToTextFile("X.txt");
                        B.SaveToTextFile("B.txt");
#endif

                        string ErrMsg;
                        using(var stw = new StringWriter()) {
                            stw.WriteLine("High residual from direct solver (using {0}).", SolverName);
                            stw.WriteLine("    L2 Norm of RHS:         " + RhsNorm);
                            stw.WriteLine("    L2 Norm of Solution:    " + SolutionNorm);
                            stw.WriteLine("    L2 Norm of Residual:    " + ResidualNorm);
                            stw.WriteLine("    Relative Residual norm: " + RelResidualNorm);
                            stw.WriteLine("    Matrix Inf norm:        " + MatrixInfNorm);
#if TEST
                            stw.WriteLine("Dumping text versions of Matrix, Solution and RHS.");
#endif
                            ErrMsg = stw.ToString();
                        }
                        Console.Error.WriteLine(ErrMsg);
                        this.Converged = false;
                    }
                }

                if(this.IterationCallback != null) {
                    double[] _xl = X.ToArray();
                    double[] _bl = B.ToArray();
                    m_Mtx.SpMV(-1.0, _xl, 1.0, _bl);
                    this.IterationCallback(1, _xl, _bl, this.MatrixNMapping as MultigridOperator);
                }
            }
        }

        ISparseSolver m_Solver = null;

        int m_ThisLevelIterations;

        

        /// <summary>
        /// Instruction for delayed caching of the factorization of block solver.
        /// Useful if memory peaks in linear solver tend to burst the memory.
        /// - 1st int: number of iterations
        /// - 2nd int: multigrid level
        /// </summary>
        public Func<int, int, bool> ActivateCaching {
            private get;
            set;
        }

        

        public int IterationsInNested {
            get { return 0; }
        }

        public int ThisLevelIterations {
            get { return m_ThisLevelIterations; }
        }

        public bool Converged {
            get;
            private set;
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        public void ResetStat() {
            m_ThisLevelIterations = 0;
        }

        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        /// <summary>
        /// Release internal memory
        /// </summary>
        public void Dispose() {
            if(m_Solver != null)
                m_Solver.Dispose();
            m_Solver = null;
            this.m_Mtx = null;
        }

        public long UsedMemory() {
            return m_UsedMemoryInLastCall;
        }

        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        public void Init(MultigridOperator op) {
            InitImpl(op);
        }
    }



}
