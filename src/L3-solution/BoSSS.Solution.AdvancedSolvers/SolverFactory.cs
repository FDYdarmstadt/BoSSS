using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.Control;
using BoSSS.Solution.Queries;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Instantiation of nonlinear solvers
    /// </summary>
    public class SolverFactory {

        /// <summary>
        /// This <see cref="SolverFactory"/> enables creation of nonlinear solver objects. The configuration <see cref="ISolverFactory"/> and <see cref="NonLinearSolverConfig"/> can be set in Control-Files (defined in <see cref="AppControl"/>).
        /// </summary>
        public SolverFactory(NonLinearSolverConfig nc, ISolverFactory lc) {
            m_lc = lc;
            m_nc = nc;
        }

        /// <summary>
        /// A <see cref="SolverFactory"/> which supports queries. Provides additional solver specific output: e.g. Used Multigridlevels, etc.
        /// </summary>
        public SolverFactory(NonLinearSolverConfig nc, ISolverFactory lc, QueryHandler qh) {
            m_lc = lc;
            m_nc = nc;
            m_qh = qh;
        }

        /// <summary>
        /// Enables additional Output via queries
        /// </summary>
        private QueryHandler m_qh;


        ///// <summary>
        ///// Is get and set by <see cref="Selfmade_nonlinsolver"/> and used by <see cref="GenerateNonLinear"/>.
        ///// </summary>
        //private NonlinearSolver m_nonlinsolver;

        /// <summary>
        /// Internal linear solver configuration. Shall always be != null. And will not be edited within this method!
        /// </summary>
        private ISolverFactory m_lc;

        /// <summary>
        /// Internal nonlinear solver configuration. Shall always be != null. And will not be edited within this method!
        /// </summary>
        private NonLinearSolverConfig m_nc;

        /// <summary>
        /// Internal nonlinear solver configuration. Shall always be != null. And will not be edited within this method!
        /// </summary>
        public NonLinearSolverConfig Config {
            get {
                return m_nc;
            }
        }


        private MultigridOperator.ChangeOfBasisConfig[][] m_MGchangeofBasis;
        
        private IEnumerable<AggregationGridBasis[]> m_MGBasis;

        /// <summary>
        /// This will return which are configured according to <see cref="Config"/>
        /// </summary>
        public void GenerateNonLin(out NonlinearSolver nonlinSolver, OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig) {

            m_MGchangeofBasis = ts_MultigridOperatorConfig;
            m_MGBasis = ts_MultigridBasis;

            nonlinSolver = null;
            nonlinSolver = GenerateNonLin_body(ts_AssembleMatrixCallback, ts_MultigridBasis, ts_MultigridOperatorConfig); 

            Debug.Assert(nonlinSolver != null);
            return;
        }

        /// <summary>
        /// This one is the method-body of <see cref="GenerateNonLin"/> and shall not be called from the outside. The parameters are mainly handed over to the NonLinearSolver object, which lives in <see cref="AdvancedSolvers.NonlinearSolver"/>.
        /// </summary>
        private NonlinearSolver GenerateNonLin_body(OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig) {

            var lc = m_lc;
            var nc = m_nc;

            NonlinearSolver nonlinSolver;

            switch(nc.SolverCode) {
                case NonLinearSolverCode.Picard:

                nonlinSolver = new FixpointIterator(
                    ts_AssembleMatrixCallback,
                    ts_MultigridBasis,
                    MultigridOperatorConfig) {
                    MaxIter = nc.MaxSolverIterations,
                    MinIter = nc.MinSolverIterations,
                    m_LinearSolver = lc,
                    ConvCrit = nc.ConvergenceCriterion,
                    UnderRelax = nc.UnderRelax,
                };
                break;

                //Besides NonLinearSolverConfig Newton needs also LinearSolverConfig
                //Newton uses MUMPS as linearsolver by default
                case NonLinearSolverCode.Newton: {
                    nonlinSolver = new Newton(
                    ts_AssembleMatrixCallback,
                    ts_MultigridBasis,
                    MultigridOperatorConfig) {
                        //maxKrylovDim = lc.MaxKrylovDim,
                        MaxIter = nc.MaxSolverIterations,
                        MinIter = nc.MinSolverIterations,
                        Globalization = nc.Globalization,
                        ApproxJac = Newton.ApproxInvJacobianOptions.ExternalSolver,
                        PrecondConfig = lc,
                        ConvCrit = nc.ConvergenceCriterion,
                        constant_newton_it = nc.constantNewtonIterations,
                        HomotopyStepLongFail = nc.HomotopyStepLongFail,
                    };
                }
                break;
                                
                default:
                throw new NotImplementedException();
            }


            SetNonLinItCallback(nonlinSolver);

            //Console.WriteLine("nonlinear solver code: {0}", nc.SolverCode.ToString());

            return nonlinSolver;
        }

        private void SetNonLinItCallback(NonlinearSolver nonlinsolver) {

            int _caseselect = 0;
            string _type = "NLinSolver";
            DefaultItCallback = GenerateDefaultCallback<NonlinearSolver>(_type, nonlinsolver, _caseselect);
            if (m_nc.verbose) {
                nonlinsolver.IterationCallback += DefaultItCallback;
            }
            //nonlinsolver.IterationCallback += CustomizedCallback;
        }

        private int[] m_Iterations;
        private void FirstLineinCallBack() {
            if (m_Iterations == null) {
                string FirstLine = "NLinS-, Precond-, LinS-, total-Iterations : Type, SolverName, InfResi, Multigridlevel";
                m_Iterations = new int[4];
                Console.WriteLine(FirstLine);
            }
        }

        private Action<int, double[], double[], MultigridOperator> GenerateDefaultCallback<T>(string type, T solverwithcallback, int caseselect) {

            string name = String.Join(".", solverwithcallback.ToString().Split('.'), 3, 1);

            return delegate (int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
                FirstLineinCallBack();
                double res = currentRes.L2NormPow2().MPISum().Sqrt();
                m_Iterations[caseselect] = iterIndex;
                m_Iterations[3]++;
                string Its = "";
                Array.ForEach<int>(m_Iterations, i => Its += i.ToString() + ",");
                Console.WriteLine("{0} : {1}, {2}, {3}, {4}", Its, type, name, res, Mgop.LevelIndex);
            };
        }

        //private Action<int, double[], double[], MultigridOperator> CustomizedCallback;

        private Action<int, double[], double[], MultigridOperator> DefaultItCallback;

    }
}
