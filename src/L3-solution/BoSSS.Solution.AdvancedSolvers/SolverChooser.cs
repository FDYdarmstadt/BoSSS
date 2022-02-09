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
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.LinSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP;
using ilPSP.Utils;
using System.IO;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;
using BoSSS.Foundation.IO;
using static BoSSS.Solution.AdvancedSolvers.MultigridOperator;
using BoSSS.Solution.Queries;

namespace BoSSS.Solution {

    public class SolverFactory {

        /// <summary>
        /// This <see cref="SolverFactory"/> enables creation of linear and nonlinear solver objects. The configuration <see cref="LinearSolverConfig"/> and <see cref="NonLinearSolverConfig"/> can be set in Control-Files (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        public SolverFactory(NonLinearSolverConfig nc, LinearSolverConfig lc) {
            m_lc = lc;
            m_nc = nc;
        }

        /// <summary>
        /// A <see cref="SolverFactory"/> which supports queries. Provides additional solver specific output: e.g. Used Multigridlevels, etc.
        /// </summary>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        /// <param name="qh"></param>
        public SolverFactory(NonLinearSolverConfig nc, LinearSolverConfig lc, QueryHandler qh) {
            m_lc = lc;
            m_nc = nc;
            m_qh = qh;
        }

        #region private get/set
        private Action<int, double[], double[], MultigridOperator> CustomizedCallback;
        private Action<int, double[], double[], MultigridOperator> DefaultItCallback;
        private int[] m_Iterations;


        private MultigridOperator.ChangeOfBasisConfig[][] m_MGchangeofBasis;
        private IEnumerable<AggregationGridBasis[]> m_MGBasis;

        /// <summary>
        /// Is get and set by <see cref="Selfmade_linsolver"/> and used by <see cref="GenerateLinear"/>.
        /// </summary>
        private ISolverSmootherTemplate m_linsolver;

        /// <summary>
        /// Is get and set by <see cref="Selfmade_nonlinsolver"/> and used by <see cref="GenerateNonLinear"/>.
        /// </summary>
        private NonlinearSolver m_nonlinsolver;

        /// <summary>
        /// Is get and set by <see cref="Selfmade_precond"/> and used by <see cref="GenerateLinear"/>.
        /// </summary>
        private ISolverSmootherTemplate m_precond;

        /// <summary>
        /// Internal linear solver configuration. Shall always be != null. And will not be edited within this method!
        /// </summary>
        private LinearSolverConfig m_lc;

        /// <summary>
        /// Internal nonlinear solver configuration. Shall always be != null. And will not be edited within this method!
        /// </summary>
        private NonLinearSolverConfig m_nc;

        /// <summary>
        /// Enables additional Output via queries
        /// </summary>
        private QueryHandler m_qh;

        // Related to MGCallback
        private int m_MG_Counter = 0;
        private double[] ProlRes = new double[10];
        private double[] RestRes = new double[10];
        private double m_ResOfPreviousSolver = 0;
        private int[] m_LocalDOF = null;

        private int m_MaxMGLevel = -1;

        /// <summary>
        /// Maximal used Mg Level. Note: Is available during Solver run but variable in Solver Chooser.
        /// Therefore, only valid after lin. Solver generation terminates
        /// </summary>
        private int MaxUsedMGLevel {
            get { return m_MaxMGLevel; }
            set { m_MaxMGLevel = Math.Max(value, m_MaxMGLevel); }
        }

        private bool SetQuery(string desciption, double value, bool setonce) {
            if (m_qh != null) {
                m_qh.ValueQuery(desciption, value, !setonce);
            }
            return m_qh != null;
        }

        private int[] _LocalDOF {
            get {
                if (m_LocalDOF == null)
                    m_LocalDOF = GetLocalDOF();
                return m_LocalDOF;
            }
        }

        private int _NoOfBlocks {
            get {
                return (int)Math.Max(1, Math.Round(_LocalDOF[0] / (double)m_lc.TargetBlockSize));
            }
        }

        #endregion

        #region public get/set
        public int[] GetIterationcounter {
            get {
                return m_Iterations;
            }
        }

        /// <summary>
        /// For developers, who want full control over solvers: In <see cref="selfmade_linsolver"/> you can insert your own config of linear solver.
        /// Clear() will clear the selfmade stuff and enables solver creation from linear and nonlinear config again.
        /// </summary>
        public ISolverSmootherTemplate Selfmade_linsolver {
            set {
                m_linsolver = value;
            }
            get {
                return m_linsolver;
            }
        }

        /// <summary>
        /// For developers, who want full control over solvers: In <see cref="selfmade_nonlinsolver"/> you can insert your own config of nonlinear solver,
        /// which will overwrite the output of <see cref="GenerateNonLinear"/> with the overgiven solver.
        /// Clear() will clear the selfmade stuff and enables solver creation from linear and nonlinear config again.
        /// Note: The overgiven solver has to be completely defined (precond!=null and linsolve!=null) 
        /// </summary>
        public NonlinearSolver Selfmade_nonlinsolver {
            set {
                m_nonlinsolver = value;
            }
            get {
                return m_nonlinsolver;
            }
        }

        /// <summary>
        /// For developers, who want full control over solvers: In <see cref="Selfmade_precond"/> you can insert your own preconditioner.
        /// Clear() will clear the selfmade stuff and enables solver creation from linear and nonlinear config again.
        /// </summary>
        public ISolverSmootherTemplate Selfmade_precond {
            set {
                m_precond = value;
            }
            get {
                return m_precond;
            }
        }

        /// <summary>
        /// Get the active linear config of the solver chooser
        /// </summary>
        public LinearSolverConfig GetLinearConfig {
            get { return m_lc; }
        }

        /// <summary>
        /// Get the active non linear solver config of the solver chooser
        /// </summary>
        public NonLinearSolverConfig GetNonLinearConfig {
            get { return m_nc; }
        }

        /// <summary>
        /// 
        /// </summary>
        public bool IsNewtonGmresType {
            get {
                bool LinIsGmres = false;
                switch (m_lc.SolverCode) {
                    case LinearSolverCode.exp_gmres_levelpmg:
                    case LinearSolverCode.exp_gmres_schwarz_pmg:
                    case LinearSolverCode.exp_softgmres:
                    case LinearSolverCode.exp_gmres_AS:
                    case LinearSolverCode.exp_gmres_AS_MG:
                    LinIsGmres = true;
                    break;
                    default:
                    break;
                }

                return m_nc.SolverCode == NonLinearSolverCode.Newton && LinIsGmres;
            }
        }
        #endregion

        #region Generate Lin/NonLin Solver - public access
        /// <summary>
        /// This will return <code>linear</code> and <code>nonlinear</code> solver objects, which are configured according to <see cref="LinearSolverConfig"/> and <see cref="NonLinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>).
        /// </summary>
        public void GenerateNonLin(out NonlinearSolver nonlinSolver, out ISolverSmootherTemplate linsolver, OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig, AggregationGridData[] ts_MGS) {

            if (m_nonlinsolver != null)
                m_nc.SolverCode = NonLinearSolverCode.selfmade;
            if (m_linsolver != null)
                m_lc.SolverCode = LinearSolverCode.selfmade;

            m_MGchangeofBasis = ts_MultigridOperatorConfig;
            m_MGBasis = ts_MultigridBasis;

            linsolver = null;
            nonlinSolver = null;
            ISolverSmootherTemplate precondonly = null;
            ISolverSmootherTemplate[] pretmp;

            linsolver = GenerateLinear_body(ts_MultigridBasis, ts_MultigridOperatorConfig, out pretmp);
            if (pretmp != null) {
                precondonly = pretmp[0];
            }
            Debug.Assert(linsolver != null);

            nonlinSolver = GenerateNonLin_body(ts_AssembleMatrixCallback, ts_MultigridBasis, ts_MultigridOperatorConfig, linsolver, precondonly); // preconditioner may be used as standalone, then linsolver pointer is redirected

                Debug.Assert(nonlinSolver != null);
            return;
        }

        /// <summary>
        /// This will return a linear solver object, which is configured according to <see cref="LinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>). 
        /// </summary>
        public void GenerateLinear(out ISolverSmootherTemplate linSolve, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig) {
            if (m_linsolver != null)
                m_lc.SolverCode = LinearSolverCode.selfmade;

            m_MGchangeofBasis = ts_MultigridOperatorConfig;
            m_MGBasis = ts_MultigridBasis;

            ISolverSmootherTemplate[] preconds;
            linSolve = GenerateLinear_body(ts_MultigridBasis, ts_MultigridOperatorConfig, out preconds);

            Debug.Assert(linSolve != null);
        }

        /// <summary>
        /// This will return a linear solver object, which is configured according to <see cref="LinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="linSolve">the linear solver object</param>
        /// <param name="ts_MultigridBasis"></param>
        /// <param name="ts_MultigridOperatorConfig"></param>
        /// <param name="IterationCallbacks">insert callback specified by user here</param>
        public void GenerateLinear(out ISolverSmootherTemplate linSolve, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig, List<Action<int, double[], double[], MultigridOperator>> IterationCallbacks) {
            IterationCallbacks.ForEach(ICB => this.CustomizedCallback += ICB);
            GenerateLinear(out linSolve, ts_MultigridBasis, ts_MultigridOperatorConfig);
        }
        #endregion

        #region Generate Lin/NonLin Solver - private method bodies
        /// <summary>
        /// This one is the method-body of <see cref="GenerateNonLin"/> and shall not be called from the outside. The parameters are mainly handed over to the NonLinearSolver object, which lives in <see cref="AdvancedSolvers.NonlinearSolver"/>.
        /// </summary>
        private NonlinearSolver GenerateNonLin_body(OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig, ISolverSmootherTemplate linsolver, ISolverSmootherTemplate precondonly) {

            var lc = m_lc;
            var nc = m_nc;

            NonlinearSolver nonlinSolver;

            switch (nc.SolverCode) {
                case NonLinearSolverCode.Picard:

                    nonlinSolver = new FixpointIterator(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                        MaxIter = nc.MaxSolverIterations,
                        MinIter = nc.MinSolverIterations,
                        m_LinearSolver = linsolver,
                        ConvCrit = nc.ConvergenceCriterion,
                        UnderRelax = nc.UnderRelax,
                    };
                    break;

                //Besides NonLinearSolverConfig Newton needs also LinearSolverConfig
                //Newton uses MUMPS as linearsolver by default
                case NonLinearSolverCode.Newton:

                    if (IsNewtonGmresType) {
                        nonlinSolver = new Newton(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                            //maxKrylovDim = lc.MaxKrylovDim,
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            ApproxJac = Newton.ApproxInvJacobianOptions.MatrixFreeGMRES,
                            Precond = precondonly,
                            ConvCrit = nc.ConvergenceCriterion,
                            printLambda = false,
                            Globalization = nc.Globalization,
                        };
                        linsolver = precondonly; // put out the solver, which is actually used!
                    } else {
                        nonlinSolver = new Newton(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                            //maxKrylovDim = lc.MaxKrylovDim,
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            printLambda = nc.verbose,
                            Globalization = nc.Globalization,
                            ApproxJac = Newton.ApproxInvJacobianOptions.ExternalSolver,
                            Precond = linsolver,
                            ConvCrit = nc.ConvergenceCriterion,
                            constant_newton_it = nc.constantNewtonIterations,
                            HomotopyStepLongFail = nc.HomotopyStepLongFail,
                        };
                    }
                    break;


                case NonLinearSolverCode.NLSolverSequence:

                    var myFixPoint = new FixpointIterator(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                        MaxIter = nc.MinSolverIterations,
                        MinIter = nc.MinSolverIterations,
                        m_LinearSolver = linsolver,
                        ConvCrit = nc.ConvergenceCriterion,
                        UnderRelax = nc.UnderRelax,
                    };
                    SetNonLinItCallback(myFixPoint);
                    // this is normal Newton
                    var myNewton = new Newton(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                        //maxKrylovDim = lc.MaxKrylovDim,
                        MaxIter = nc.MaxSolverIterations,
                        MinIter = nc.MinSolverIterations,
                        Globalization = nc.Globalization,
                        ApproxJac = Newton.ApproxInvJacobianOptions.ExternalSolver,
                        Precond = linsolver,
                        ConvCrit = nc.ConvergenceCriterion,
                        constant_newton_it = nc.constantNewtonIterations
                    };
                    SetNonLinItCallback(myNewton);
                    nonlinSolver = new NLSolverSequence(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                        m_NLSequence = new NonlinearSolver[] { myFixPoint, myNewton }
                    };
                    break;

                case NonLinearSolverCode.selfmade:
                    Console.WriteLine("INFO: Selfmade Nonlinear Solver used!");
                    nonlinSolver = m_nonlinsolver;
                    break;
                default:
                    throw new NotImplementedException();
            }

            if(nc.SolverCode != NonLinearSolverCode.selfmade)
                Check_NonLinearSolver(nonlinSolver);

            SetNonLinItCallback(nonlinSolver);

            //Console.WriteLine("nonlinear solver code: {0}", nc.SolverCode.ToString());

            return nonlinSolver;
        }



        private ISolverSmootherTemplate[] GeneratePrecond(IEnumerable<AggregationGridBasis[]> MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig) {

            var lc = m_lc;
            var precond = new ISolverSmootherTemplate[1];

            // some solver require special parameters, like ...
            
            //int SpaceDim = MultigridSequence[0].SpatialDimension;
            int MaxMGDepth = MultigridBasis.Count();

            switch (lc.SolverCode) {

                //no preconditioner used ...
                case LinearSolverCode.classic_mumps:
                case LinearSolverCode.classic_pardiso:
                case LinearSolverCode.classic_cg:
                case LinearSolverCode.exp_softgmres:
                case LinearSolverCode.exp_Kcycle_schwarz:
                case LinearSolverCode.exp_Kcycle_ILU:
                case LinearSolverCode.exp_decomposedMG_OrthoScheme:
                case LinearSolverCode.exp_Kcycle_schwarz_4Rheology:
                case LinearSolverCode.exp_AS:
                case LinearSolverCode.exp_AS_MG:
                case LinearSolverCode.exp_pTG:
                case LinearSolverCode.automatic:
                precond[0] = null;
                break;
                case LinearSolverCode.exp_gmres_ILU:
                precond[0] = new HypreILU();
                break;
                case LinearSolverCode.exp_gmres_AS_MG:
                var dirSolver = new DirectSolver() {
                    WhichSolver = DirectSolver._whichSolver.PARDISO,
                };
                var Smoother = new Schwarz() {
                    m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                        NoOfPartsOnCurrentProcess = _NoOfBlocks,
                    },
                    Overlap = 1,
                    EnableOverlapScaling = false,
                    
                };
                precond[0] = ClassicMGwithSmoother(3, dirSolver, Smoother);
                (precond[0] as ClassicMultigrid).TerminationCriterion = (iIter, r0, ri) => iIter <= 10;
                break;
                case LinearSolverCode.exp_gmres_MG_gmres:
                    var dirSolver2 = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                    };
                    var addSchwarz = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = 4,
                        },
                        Overlap = 1,
                        EnableOverlapScaling = false,

                    };
                    var smoother2 = new SoftGMRES() {
                        MaxKrylovDim = 1,
                        Precond = addSchwarz
                    };
                    precond[0] = ClassicMGwithSmoother(1, dirSolver2, smoother2);
                    (precond[0] as ClassicMultigrid).TerminationCriterion = (int iter, double r0, double r) => iter < 1;
                    break;
                /*
                case LinearSolverCode.exp_gmres_localPrec:
                    precond[0] = new LocalizedOperatorPrec() {
                        m_dt = lc.exp_localPrec_Min_dt,
                        m_muA = lc.exp_localPrec_muA,
                    };
                    break;
                */
                case LinearSolverCode.exp_gmres_AS:

                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = _NoOfBlocks,
                        },
                        Overlap = 1,
                        EnableOverlapScaling = true,
                        //TwoGrid = new Schwarz.LowOrderCoarseSolver() { isMultiplicative = false }
                    };
                    break;

                case LinearSolverCode.exp_softpcg_schwarz:

                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = _NoOfBlocks,
                        },
                        Overlap = 1,
                        EnableOverlapScaling = true,
                        //TwoGrid = new Schwarz.LowOrderCoarseSolver() { isMultiplicative = false }
                    };
                break;

                case LinearSolverCode.exp_softpcg_schwarz_directcoarse:

                    precond[0] = new Schwarz() {
                        CoarseSolver = new Schwarz.ClassicMG(),
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = _NoOfBlocks
                        },
                        Overlap = 1,
                    };
                    break;

                case LinearSolverCode.exp_gmres_levelpmg:
                    precond[0] = new LevelPmg() { UseHiOrderSmoothing = true, OrderOfCoarseSystem = m_lc.pMaxOfCoarseSolver, FullSolveOfCutcells = true};
                    SetQuery("XdgCellsToLowBlock", ((LevelPmg)precond[0]).FullSolveOfCutcells ? 1 : 0, true);
                    break;

                case LinearSolverCode.exp_gmres_schwarz_pmg:
                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = _NoOfBlocks
                        },
                        Overlap = 1,
                        EnableOverlapScaling = false,
                        UsePMGinBlocks = true,

                    };
                    break;

                case LinearSolverCode.exp_softpcg_jacobi_mg:

                    ISolverSmootherTemplate[] _prechain = new ISolverSmootherTemplate[] {
                        new SoftPCG() {
                             NoOfIterations = 5,
                        }
                    };

                    ISolverSmootherTemplate[] _postchain = new ISolverSmootherTemplate[]{
                        new SoftPCG() {
                             NoOfIterations = 5,
                        }
                    };

                    ISolverSmootherTemplate toppre = new BlockJacobi() {
                        NoOfIterations = 3,
                        omega = 0.5
                    };

                    ISolverSmootherTemplate toppst = new BlockJacobi() {
                        NoOfIterations = 3,
                        omega = 0.5
                    };

                    precond[0] = My_MG_Precond(MaxMGDepth, _prechain, _postchain, toppre, toppst);
                    break;

                case LinearSolverCode.exp_OrthoS_pMG:
                    precond = new ISolverSmootherTemplate[]{
                        //new LevelPmg() {
                        //    UseHiOrderSmoothing =true,
                        //    CoarseLowOrder=1
                        //},

                        new Schwarz() {
                            FixedNoOfIterations = 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //NoOfPartsPerProcess = LocalNoOfSchwarzBlocks
                                NoOfPartsOnCurrentProcess = 4
                            },
                            Overlap = 1, // overlap seems to help; more overlap seems to help more
                            EnableOverlapScaling = true,
                            UsePMGinBlocks = false,
                            CoarseSolveOfCutcells = true,
                            CoarseLowOrder = m_lc.pMaxOfCoarseSolver
                        }
                    };
                    break;

                case LinearSolverCode.selfmade:
                    Console.WriteLine("INFO: Selfmade Preconditioner is used!");
                    precond[0] = m_precond;
                    break;
                case LinearSolverCode.exp_another_Kcycle:
                    precond[0] = null;
                    //precond[0]= expKcycleSchwarz(MaxMGDepth, LocalDOF, X => m_lc.TargetBlockSize);
                    break; 

                default:
                    throw new NotImplementedException($"Preconditioner for solver code {lc.SolverCode} not available.");
            }


            foreach (var pre in precond) {
                if (pre == null)
                    continue;
                if(lc.verbose)
                    GenerateInfoMessageAtSetup(pre);
                SetLinItCallback(pre, true);
            }

            return precond;
        }


        /// <summary>
        /// This one is the method-body of <see cref="GenerateLinear(out ISolverSmootherTemplate, IEnumerable{AggregationGridBasis[]}, ChangeOfBasisConfig[][], List{Action{int, double[], double[], MultigridOperator}})"/> and shall not be called from the outside. 
        /// Some Solver acquire additional information, thus the timestepper is passed as well.
        /// </summary>
        private ISolverSmootherTemplate GenerateLinear_body(IEnumerable<AggregationGridBasis[]> MultigridBasis, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig, out ISolverSmootherTemplate[] precond) {

            var lc = m_lc;
            //Console.WriteLine("linear solver code : {0}", lc.SolverCode.ToString());
            ISolverSmootherTemplate templinearSolve = null;
            precond = GeneratePrecond(MultigridBasis, MultigridOperatorConfig);

            // some solver require special parameters, like ...
            int NoCellsLoc = MultigridBasis.First()[0].AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int SpaceDim = MultigridBasis.First()[0].AggGrid.ParentGrid.SpatialDimension;
            int MaxMGDepth = MultigridBasis.Count();

            if (MaxMGDepth < 1)
                throw new ArgumentException("ERROR: At least one multigrid levels is required.");

            switch (lc.SolverCode) {
                case LinearSolverCode.automatic:
                if(m_nc != null) {
                    templinearSolve = AutomaticSolver( SpaceDim, NoCellsLoc,
                        precond);
                }
                break;

                case LinearSolverCode.classic_mumps:
                templinearSolve = new DirectSolver() {
                    WhichSolver = DirectSolver._whichSolver.MUMPS,
                    SolverVersion = Parallelism.MPI,
                };
                break;

                case LinearSolverCode.classic_pardiso:
                templinearSolve = new DirectSolver() {
                    WhichSolver = DirectSolver._whichSolver.PARDISO,
                    SolverVersion = Parallelism.OMP,
                };
                break;

                case LinearSolverCode.exp_AS:

                templinearSolve = new Schwarz() {
                    FixedNoOfIterations = m_lc.MaxSolverIterations,
                    m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                        NoOfPartsOnCurrentProcess = _NoOfBlocks,
                    },
                    Overlap = 1, // overlap seems to help; more overlap seems to help more
                    EnableOverlapScaling = true,
                    //UsePMGinBlocks = false,
                    //CoarseSolveOfCutcells = true,
                    //CoarseLowOrder = m_lc.pMaxOfCoarseSolver
                };
                break;

                case LinearSolverCode.exp_AS_MG:

                if(lc.NoOfMultigridLevels < 2)
                    throw new ApplicationException("At least 2 Multigridlevels are required");

                var smoother = new Schwarz() {
                    FixedNoOfIterations = m_lc.MaxSolverIterations,
                    m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                        NoOfPartsOnCurrentProcess = _NoOfBlocks,
                    },
                    Overlap = 1, // overlap seems to help; more overlap seems to help more
                    EnableOverlapScaling = true,
                };
                var coarsesolver = new DirectSolver() {
                    WhichSolver = DirectSolver._whichSolver.PARDISO,
                    SolverVersion = Parallelism.SEQ,
                };
                templinearSolve = ClassicMGwithSmoother(MaxMGDepth-1, coarsesolver, smoother);

                break;

                case LinearSolverCode.classic_cg:
                templinearSolve = new MonkeySolver() {
                    WhichSolver = MonkeySolver._whichSolver.CG,
                    LinConfig = lc
                };
                break;

                case LinearSolverCode.exp_softpcg_schwarz:
                case LinearSolverCode.exp_softpcg_schwarz_directcoarse:
                case LinearSolverCode.exp_softpcg_jacobi_mg:

                templinearSolve = new SoftPCG() {
                    //m_MaxIterations = lc.MaxSolverIterations,
                    //m_Tolerance = lc.ConvergenceCriterion,
                    Precond = precond[0]
                };
                break;

                case LinearSolverCode.exp_gmres_levelpmg:
                case LinearSolverCode.exp_gmres_schwarz_pmg:
                case LinearSolverCode.exp_softgmres:
                case LinearSolverCode.exp_gmres_AS:
                case LinearSolverCode.exp_gmres_AS_MG:
                
                templinearSolve = new SoftGMRES() {
                    //m_Tolerance = lc.ConvergenceCriterion,
                    //m_MaxIterations = lc.MaxSolverIterations,
                    MaxKrylovDim = lc.MaxKrylovDim,
                    Precond = precond[0]
                };
                break;
                case LinearSolverCode.exp_gmres_MG_gmres:
                    templinearSolve = new FlexGMRES() {
                        //m_Tolerance = lc.ConvergenceCriterion,
                        //m_MaxIterations = lc.MaxSolverIterations,
                        MaxKrylovDim = lc.MaxKrylovDim,
                        PrecondS = precond
                    };
                    break;
                case LinearSolverCode.exp_gmres_ILU:
                templinearSolve = new SoftGMRES() {
                    //m_Tolerance = lc.ConvergenceCriterion,
                    //m_MaxIterations = lc.MaxSolverIterations,
                    MaxKrylovDim = lc.MaxKrylovDim,
                    Precond = precond[0]
                };
                break;
                case LinearSolverCode.exp_Kcycle_schwarz: {
                    Func<int, int> SblkSizeFunc = delegate (int iLevel) { return m_lc.TargetBlockSize; };
                    templinearSolve = KcycleMultiSchwarz(MaxMGDepth, SblkSizeFunc);
                    break;
                }
                case LinearSolverCode.exp_Kcycle_ILU: {
                    Func<int, int> SblkSizeFunc = delegate (int iLevel) { return 10000; };
                    templinearSolve = KcycleMultiILU(MaxMGDepth, SblkSizeFunc);
                    break;
                }

                case LinearSolverCode.exp_Kcycle_schwarz_4Rheology:
                templinearSolve = KcycleMultiSchwarz_4Rheology();
                break;

                case LinearSolverCode.exp_decomposedMG_OrthoScheme:

                ISolverSmootherTemplate[] subsmoother = new ISolverSmootherTemplate[]{
                        new SoftPCG() {
                             NoOfIterations=5
                             //m_MaxIterations = 5,
                             //m_MinIterations = 5,
                        }
                    };

                ISolverSmootherTemplate[] topsmoother = new ISolverSmootherTemplate[]{
                        new BlockJacobi()
                        {
                            NoOfIterations = 5,
                            omega = 0.5
                        }
                    };

                templinearSolve = MakeOrthoNormMGDecomp(MaxMGDepth, new SolverSquence() { SolverChain = subsmoother }, new SolverSquence() { SolverChain = topsmoother });
                break;

                case LinearSolverCode.exp_OrthoS_pMG:

                templinearSolve = new OrthonormalizationScheme() {
                    PrecondS = precond,
                    MaxKrylovDim = lc.MaxKrylovDim,
                    MaxIter = lc.MaxSolverIterations,
                    Tolerance = lc.ConvergenceCriterion,
                    Restarted = false
                };
                break;
                case LinearSolverCode.exp_another_Kcycle:
                    templinearSolve = new FlexGMRES() {
                        PrecondS = new ISolverSmootherTemplate[] { expKcycleSchwarz(MaxMGDepth, X => m_lc.TargetBlockSize) },
                        MaxKrylovDim = 50,
                        TerminationCriterion = (int iter, double r0, double r) => iter <= 10,
                    };
                    break;
                case LinearSolverCode.exp_pTG:
                    templinearSolve = new LevelPmg() {
                        OrderOfCoarseSystem = m_lc.pMaxOfCoarseSolver,
                        UseHiOrderSmoothing = true,
                    };
                    break;
                case LinearSolverCode.selfmade:
                Console.WriteLine("INFO: Selfmade LinearSolver is used!");
                templinearSolve = m_linsolver;
                break;

                default:
                throw new NotImplementedException($"Linear solver for code {lc.SolverCode} not available.");
            }
            Debug.Assert(templinearSolve != null);
            SetLinItCallback(templinearSolve,false);

            // set Queries here which apply to all solvers:
            SetQuery("UsedMGDepth", MaxUsedMGLevel + 1, true);

            if (lc.verbose)
                Console.WriteLine("linear solver : {0}", templinearSolve.ToString().Split('.').Last());

            if(templinearSolve is IProgrammableTermination pt) {
                 // Sets default termination criterion if there is none yet
                bool LinearConvergence(int iter, double r0_l2, double r_l2) {
                    if (iter <= lc.MinSolverIterations)
                        return true;

                    if (iter > lc.MaxSolverIterations)
                        return false;

                    if (r_l2 < lc.ConvergenceCriterion)
                        return false; // terminate

                    return true; // keep running
                };
                pt.TerminationCriterion = LinearConvergence;
            }

            Check_linsolver(templinearSolve);
            return templinearSolve;
        }
        #endregion

        #region auxiliary methods
        private void GenerateInfoMessageAtSetup(ISolverSmootherTemplate solver) {
            string solvername = solver.ToString().Split('.').Last();
            Console.WriteLine("preconditioner : {0}", solvername);
            List<string> solverinfotxt = new List<string>();
            if (solver.GetType() == typeof(LevelPmg)) {
                solverinfotxt.Add("entries upto p=" + ((LevelPmg)solver).OrderOfCoarseSystem + " are assigned to low order blocks");
            }
            if (solver.GetType() == typeof(Schwarz)) {
                if (((Schwarz)solver).UsePMGinBlocks) {
                    solverinfotxt.Add("pmg used in Schwarz Blocks");
                    solverinfotxt.Add("entries upto p=" + ((Schwarz)solver).CoarseLowOrder + " are assigned to low order blocks");
                }
                solverinfotxt.Add("Additive Schwarz w. direct coarse, No of blocks: " + _NoOfBlocks);
            }
            foreach (string line in solverinfotxt)
                Console.WriteLine(solvername + " INFO:\t" + line);
        }

        private int getMaxDG (int iLevel, int iVar) {
            //This workaround takes into account, the wierd structure of the ChangeOfBasis-thing
            //prohibits out of bounds exception

            var MGChangeOfBasis = m_MGchangeofBasis;

            int tLevel = iLevel < MGChangeOfBasis.Length ? iLevel : MGChangeOfBasis.Length - 1;
            int tVar = iVar < MGChangeOfBasis[tLevel].Length ? iVar : MGChangeOfBasis[tLevel].Length - 1;
            return MGChangeOfBasis[tLevel][tVar].DegreeS[0];
        }

        /// <summary>
        /// Returns DOF per multi grid level on this process.
        /// If you are interested in the DOF of a subsystem, up to a polynomial degree of <paramref name="pOfLowOrderSystem"/>. Set a value for <paramref name="pOfLowOrderSystem"/>.
        /// </summary>
        /// <param name="pOfLowOrderSystem"></param>
        /// <returns></returns>
        private int[] GetLocalDOF(int pOfLowOrderSystem = -1) {
            var MGChangeOfBasis = m_MGchangeofBasis;
            var MultigridBasis = m_MGBasis;
            var MGBasis = MultigridBasis.ToArray();
            int MGDepth = MGBasis.Length;
            int[] LocalDOF = new int[MGDepth];
            int D = ((AggregationGridBasis)MultigridBasis.First()[0]).AggGrid.ParentGrid.SpatialDimension;

            bool IsSinglePhase = MGBasis[0][0] is AggregationGridBasis || ((XdgAggregationBasis)MGBasis[0][0]).UsedSpecies.Length == 1;
            bool IsXdg = MGBasis[0][0] is XdgAggregationBasis;

            if (IsSinglePhase) {
                LocalDOF = GetLocalDOF_DG(MultigridBasis, MGChangeOfBasis, pOfLowOrderSystem);
            } else if(IsXdg) {
                // Don't you dare to comment this out!!!
                LocalDOF = GetLocalDOF_XDG(MultigridBasis, pOfLowOrderSystem);
            } else {
                throw new NotSupportedException("the type of "+ MGBasis[0][0].GetType()+" is not supported");
            }
            return LocalDOF;
        }

        /// <summary>
        /// Gets the DOF of system with multiple species.
        /// If you are interested in the DOF of a subsystem, up to a polynomial degree of <paramref name="pOfLowOrderSystem"/>. Set a value for <paramref name="pOfLowOrderSystem"/>.
        /// </summary>
        /// <param name="MGBasisEnum"></param>
        /// <param name="pOfLowOrderSystem"></param>
        /// <returns></returns>
        private int[] GetLocalDOF_XDG(IEnumerable<AggregationGridBasis[]> MGBasisEnum, int pOfLowOrderSystem) {
            var MGBasis = MGBasisEnum.ToArray();
            int MGDepth = MGBasis.Length;
            var AggBasis = (XdgAggregationBasis[][])MGBasis;
            int D = AggBasis[0][0].AggGrid.SpatialDimension;
            int[] LocalDOF = new int[MGDepth];

            int NoOfFakeCells = 0;//number of single species cells that would have same localDOF

            for (int iLevel = 0; iLevel < MGDepth; iLevel++) {
                for (int iCell = 0; iCell < AggBasis[iLevel][0].AggCellsSpecies.Length; iCell++) {
                    for (int iVar = 0; iVar < AggBasis[iLevel].Length; iVar++) {
                        NoOfFakeCells += AggBasis[iLevel][iVar].AggCellsSpecies[iCell].Length;
                    }
                }
            }

            // This loop can be merged with GetLocalDOF_DG in further refactoring
            for (int iLevel = 0; iLevel < MGDepth; iLevel++) {
                for (int iVar = 0; iVar < AggBasis[iLevel].Length; iVar++) {
                    int p = -1;
                    if (pOfLowOrderSystem <= -1) {
                        p = getMaxDG(iLevel, iVar);
                    } else {
                        p = iVar == D ? pOfLowOrderSystem - 1 : pOfLowOrderSystem;
                    }
                    LocalDOF[iLevel] += DOFofDegree(p, D)* NoOfFakeCells;
                }
            }


            return LocalDOF;
        }



        /// <summary>
        /// extra DOF in cells with multiple species are not considered by this method.
        /// One can expect only reliable results for single phase problems.
        /// If you are only interested in the DOF of a subsystem, up to a polynomial degree of <paramref name="pOfLowOrderSystem"/>. Set a value for <paramref name="pOfLowOrderSystem"/>.
        /// </summary>
        /// <param name="MultigridBasis"></param>
        /// <param name="MGChangeOfBasis"></param>
        /// <param name="pOfLowOrderSystem">maximal order of the low order system</param>
        /// <returns></returns>
        private int[] GetLocalDOF_DG(IEnumerable<AggregationGridBasis[]> MultigridBasis, ChangeOfBasisConfig[][] MGChangeOfBasis, int pOfLowOrderSystem) {
            int NoOfLevels = MultigridBasis.Count();
            int[] DOFperCell = new int[NoOfLevels];
            int[] LocalDOF = new int[NoOfLevels];
            var MGBasisAtLevel = MultigridBasis.ToArray();
            int[] NoOFCellsAtLEvel = MGBasisAtLevel.Length.ForLoop(b=> MGBasisAtLevel[b].First().AggGrid.iLogicalCells.NoOfLocalUpdatedCells);
            int counter = 0;
            int D = ((AggregationGridBasis)MultigridBasis.First()[0]).AggGrid.ParentGrid.SpatialDimension;

            for (int iLevel = 0; iLevel < DOFperCell.Length; iLevel++) {
                counter = iLevel;
                if (iLevel >= MGChangeOfBasis.Length)
                    counter = MGChangeOfBasis.Length - 1;
                foreach (var cob in MGChangeOfBasis[counter]) {
                    for (int iVar = 0; iVar < cob.VarIndex.Length; iVar++) {
                        int p = -1;
                        if (pOfLowOrderSystem <= -1)
                            p = cob.DegreeS[iVar];
                        else
                            p = iVar == D ? pOfLowOrderSystem - 1 : pOfLowOrderSystem;
                        DOFperCell[iLevel]+=DOFofDegree(p,D);
                    }
                }
                LocalDOF[iLevel] = NoOFCellsAtLEvel[iLevel] * DOFperCell[iLevel];
            }
            return LocalDOF;
        }

        private static int DOFofDegree(int p, int D) {
            int toacc = 0;
            switch (D) {
                case 1:
                    toacc = p + 1 + p + 1;
                    break;
                case 2:
                    toacc = (p * p + 3 * p + 2) / 2;
                    break;
                case 3:
                    toacc = (p * p * p + 6 * p * p + 11 * p + 6) / 6;
                    break;
                default:
                    throw new Exception("wtf?Spacialdim=1,2,3 expected");
            }
            return toacc;
        }

        /// <summary>
        /// Determine max MG level from target blocksize
        /// </summary>
        /// <param name="DirectKickIn"></param>
        /// <param name="MSLength"></param>
        /// <param name="LocalDOF"></param>
        /// <returns></returns>
        private int GetMaxMGLevel(int DirectKickIn, int MSLength, int[] LocalDOF) {
            int maxlevel = -1;
            for (int iLevel = 0; iLevel < MSLength; iLevel++) {
                int SysSize = LocalDOF[iLevel].MPISum();
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == MSLength - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    maxlevel = iLevel;
                }
            }
            return maxlevel;
        }


        private int MPIsize {
            get {
                int MPIsize;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MPIsize);
                return MPIsize;
            }
        }

        private int MPIrank {
            get {
                int MPIrank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);
                return MPIrank;
            }
        }

        #endregion

        #region callback methods


        private void SetLinItCallback(ISolverSmootherTemplate solverwithoutcallback, bool IsLinPrecond) {

            ISolverWithCallback _solverwithcallback = solverwithoutcallback as ISolverWithCallback;
            if(_solverwithcallback == null)
                return;

            string _type = null;
            int _caseselect = -1;
            if (IsLinPrecond) {
                _type = "Precond";
                _caseselect = 1;
            } else {
                _type = "LinSolver";
                _caseselect = 2;
            }

            DefaultItCallback = GenerateDefaultCallback<ISolverWithCallback>(_type, _solverwithcallback, _caseselect);
            if (m_lc.verbose) {
                _solverwithcallback.IterationCallback += DefaultItCallback;
            }
            _solverwithcallback.IterationCallback += CustomizedCallback;
        }

        private void SetNonLinItCallback(NonlinearSolver nonlinsolver) {

            int _caseselect = 0;
            string _type = "NLinSolver";
            DefaultItCallback = GenerateDefaultCallback<NonlinearSolver>(_type, nonlinsolver, _caseselect);
            if (m_nc.verbose) {
                nonlinsolver.IterationCallback += DefaultItCallback;
            }
            nonlinsolver.IterationCallback += CustomizedCallback;
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

        private void FirstLineinCallBack() {
            if (m_Iterations == null) {
                string FirstLine = "NLinS-, Precond-, LinS-, total-Iterations : Type, SolverName, InfResi, Multigridlevel";
                m_Iterations = new int[4];
                Console.WriteLine(FirstLine);
            }
        }

        private void MultigridCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            int currentMGLevel = Mgop.LevelIndex;

            if (m_MG_Counter - currentMGLevel == +1) {
                double residualNormAf = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm after coarse grid correction
                ProlRes[currentMGLevel] = residualNormAf;
                Console.WriteLine("after Prolongation {0}<-{1}: {2}", currentMGLevel, currentMGLevel + 1, ProlRes[currentMGLevel] - ProlRes[currentMGLevel + 1]);
            }

            if (m_MG_Counter - currentMGLevel == 0 && currentMGLevel == MaxUsedMGLevel) {
                double residualNormAf = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm after coarse grid correction
                ProlRes[currentMGLevel] = residualNormAf;
                Console.WriteLine("after Prolongation {0}<-{1}: {2}", currentMGLevel, currentMGLevel + 1, ProlRes[currentMGLevel] - RestRes[currentMGLevel]);
            }

            if (m_MG_Counter - currentMGLevel == 0 && currentMGLevel == 0) {
                double residualNormB4 = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm before coarse grid correction
                RestRes[currentMGLevel] = residualNormB4;
                Console.WriteLine("before Restriction {0}->{1}: {2}", currentMGLevel, currentMGLevel + 1, RestRes[currentMGLevel] - ProlRes[currentMGLevel]);
            }

            if (m_MG_Counter - currentMGLevel == -1) {
                double residualNormB4 = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm before coarse grid correction
                RestRes[currentMGLevel] = residualNormB4;
                Console.WriteLine("before Restriction {0}->{1}: {2}", currentMGLevel, currentMGLevel + 1, RestRes[currentMGLevel] - RestRes[currentMGLevel - 1]);
            }
            m_MG_Counter = currentMGLevel;
            MaxUsedMGLevel = currentMGLevel;
        }
        #endregion

        #region special configs
        /// <summary>
        /// Obsolete. Has to be revisited ...
        /// 
        /// </summary>
        private ISolverSmootherTemplate AutomaticSolver( int Dim, int NoCellsLoc, ISolverSmootherTemplate[] PreCond) {

            NonLinearSolverConfig nc = null;
            if (m_nc != null) {
                nc = m_nc;
            } else {
                throw new ArgumentException("No NonlinearSolver specified");
            }

            var D = Dim;

            //int pV = Control.FieldOptions["VelocityX"].Degree;
            int pV = _LocalDOF[0];
            int pP = pV - 1; //Control.FieldOptions["Pressure"].Degree;


            // Detecting variables for solver determination 

            var cellsLoc = NoCellsLoc;
            var cellsGlo = NoCellsLoc.MPISum();

            ISolverSmootherTemplate tempsolve = null;

                switch (D) {
                    case 1:
                        tempsolve = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS
                        };
                        break;
                    case 2:
                        tempsolve = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS
                        };
                        break;
                    case 3:
                        var dofsPerCell3D = _LocalDOF[0] / NoCellsLoc;
                        var dofsLoc = dofsPerCell3D * cellsLoc;
                        var dofsGlo = dofsPerCell3D * cellsGlo;

                        if (dofsGlo > 10000) {

                            if (m_lc.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigrid levels are required");

                            tempsolve = new SoftGMRES() {
                                MaxKrylovDim = m_lc.MaxKrylovDim,
                                //m_Tolerance = lc.ConvergenceCriterion,
                                TerminationCriterion = ((iter, r0_l2, r_l2) => r_l2 > m_lc.ConvergenceCriterion && iter < m_lc.MaxSolverIterations),
                                Precond = new Schwarz() {
                                    m_BlockingStrategy = new Schwarz.SimpleBlocking() {
                                        NoOfPartsPerProcess = (int)Math.Ceiling(dofsLoc / 6500.0),
                                    },
                                    Overlap = 1,
                                },
                            };
                        } else {
                            tempsolve = new DirectSolver() {
                                WhichSolver = DirectSolver._whichSolver.MUMPS,
                            };
                        }
                        break;

                    default:
                        throw new NotImplementedException(String.Format("WTF?! You are sure, that your problem has {0} dimensions", D));

            }

            return tempsolve;
            //Timestepper.

            // Wenn Gesamtproblem in 2D < 100000 DoFs -> Direct Solver
            // Wenn Gesamtproblem in 3D < 10000 DoFs -> Direct Solver 

            // Block Solve 3D ca. 6000 DoFs per Process -> Adjust Blocks per Process
            // Coarse Solve ca. 5000 bis 10000 DoFs. -> Adjust Multigrid Levels

        }


        /// <summary>
        /// experimental. Is connected to Decomposed MG OrthoScheme. Can be deleted if not used anymore ...
        /// </summary>
        private ISolverSmootherTemplate ClassicMGwithSmoother(int MGlevels, ISolverSmootherTemplate coarseSolver, ISolverSmootherTemplate smoother=null)
        {
            ISolverSmootherTemplate solver;
            ISolverSmootherTemplate thislevelsmoother = null;
            if (smoother != null) {
                thislevelsmoother = smoother.CloneAs();
            }
            if (MGlevels > 0)
            {
                solver = new ClassicMultigrid() {                  
                    CoarserLevelSolver = ClassicMGwithSmoother(--MGlevels, coarseSolver, smoother),
                    PreSmoother = thislevelsmoother,
                    PostSmoother = thislevelsmoother,
                };
            }
            else
            {
                solver = coarseSolver;
            }
            return solver;
        }

        /// <summary>
        /// experimental. Is connected to Decomposed MG OrthoScheme. Can be deleted if not used anymore ...
        /// </summary>
        /// <param name="MGlevels"></param>
        /// <param name="lc"></param>
        /// <param name="coarseSolver"></param>
        /// <returns></returns>
        private ISolverSmootherTemplate BareMGSquence(int MGlevels, ISolverSmootherTemplate coarseSolver, ISolverSmootherTemplate smoother, ISolverSmootherTemplate topsmoother)
        {
            ISolverSmootherTemplate solver;
            if (MGlevels > 0)
            {
                solver = new ClassicMultigrid() {
                    CoarserLevelSolver = ClassicMGwithSmoother(MGlevels - 1, coarseSolver, smoother),
                    PreSmoother= topsmoother.CloneAs(),
                    PostSmoother= topsmoother.CloneAs()
                };
            }
            else
            {
                solver = coarseSolver;
            }
            return solver;
        }
            
        private int MGOcc = 0;
        private Stopwatch MGTimer = new Stopwatch();
        private void MultigridAnalysis(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop)
        {
            if (Mgop.LevelIndex == 0 && Mgop.GridData.MpiRank == 0)
            {
                string file = "MG-Analysis.txt";
                if (MGOcc == 0) { MultigridAnalysisHeader(file); MGTimer.Start(); }
                if (iterIndex == 0) { MGOcc++; MGTimer.Restart(); }
                StreamWriter sw = new StreamWriter(file, true);
                // IMPLEMENT A WAY TO KNOW IF THIS IS A SEPERATE MG OCCURRENCE
                var time = MGTimer.Elapsed;
                sw.Write(MGOcc);
                sw.Write(",");
                sw.Write(iterIndex);
                sw.Write(",");
                sw.Write(time);
                sw.Write(",");
                sw.Write(currentRes.L2Norm());
                sw.WriteLine();
                sw.Flush();
                sw.Close();
            }
        }

        private void MultigridAnalysisHeader(string file)
        {
            StreamWriter sw = new StreamWriter(file, false);
            sw.Write("MG-Occurrence");
            sw.Write(",");
            sw.Write("MG-Iteration");
            sw.Write(",");
            sw.Write("Time");
            sw.Write(",");
            sw.Write("Residual");
            sw.WriteLine();
            sw.Flush();
            sw.Close();
        }

        private ISolverSmootherTemplate My_MG_Precond(int MGLength, ISolverSmootherTemplate[] prechain, ISolverSmootherTemplate[] postchain, ISolverSmootherTemplate toplevelpre, ISolverSmootherTemplate toplevelpst)
        {

            int MGDepth = Math.Min(MGLength, m_lc.NoOfMultigridLevels);

            bool isLinPrecond = true;
            

            int DirectKickIn = m_lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            ISolverSmootherTemplate[] MultigridChain = new ISolverSmootherTemplate[MGDepth];
            for (int iLevel = 0; iLevel < MGDepth; iLevel++) {
                MaxUsedMGLevel = iLevel;
                int SysSize = _LocalDOF[iLevel].MPISum();
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == MGDepth - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    MultigridChain[iLevel] = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.MUMPS,
                        TestSolution = false
                    };
                } else {

                    ISolverSmootherTemplate[] newpostchain = new ISolverSmootherTemplate[postchain.Length];
                    ISolverSmootherTemplate[] newprechain = new ISolverSmootherTemplate[prechain.Length];

                    ClassicMultigrid MgLevel = new ClassicMultigrid() {
                        TerminationCriterion = ((iter, r0_l2, r_l2) => iter <= 1) // termination controlled by top level PCG
                    };

                    ((ISolverWithCallback)MgLevel).IterationCallback += MultigridCallback;


                    MultigridChain[iLevel] = MgLevel;

                    ISolverSmootherTemplate pre, pst;
                    if (iLevel > 0) {
                        

                        for(int i=0;i< prechain.Length;i++) {
                            newprechain[i] = prechain[i].CloneAs();
                            SetLinItCallback(newprechain[i], isLinPrecond);
                        }

                        for (int i = 0; i <postchain.Length; i++) {
                            newpostchain[i] = postchain[i].CloneAs();
                            SetLinItCallback(newpostchain[i], isLinPrecond);
                        }

                        pre = new SolverSquence() { SolverChain = newprechain };
                        pst = new SolverSquence() { SolverChain = newpostchain };

                    } else {

                        pre = toplevelpre;
                        pst = toplevelpst;

                        SetLinItCallback(toplevelpre,  isLinPrecond);
                        SetLinItCallback(toplevelpst,  isLinPrecond);
                    }

                    MgLevel.PreSmoother = pre;
                    MgLevel.PostSmoother = pst;
                }

                if (iLevel > 0) {
                    ((ClassicMultigrid)(MultigridChain[iLevel - 1])).CoarserLevelSolver = MultigridChain[iLevel];
                }

                if (useDirect) {
                    Console.WriteLine("MG: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }
            } 

            return MultigridChain[0];
        }

        private ISolverSmootherTemplate SpecialMultilevelSchwarz(LinearSolverConfig _lc, int[] _LocalDOF, int MSLength, MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig) {


            // my tests show that the ideal block size may be around 10'000
            int DirectKickIn = _lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            //MultigridOperator Current = op;
            ISolverSmootherTemplate[] MultigridChain = new ISolverSmootherTemplate[MSLength];
            for (int iLevel = 0; iLevel < MSLength; iLevel++) {
                MaxUsedMGLevel = iLevel;
                int SysSize = _LocalDOF[iLevel].MPISum();
                //int SysSize = Current.Mapping.TotalLength;
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == MSLength - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    MultigridChain[iLevel] = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.MUMPS,
                        TestSolution = false
                    };
                } else {

                    ClassicMultigrid MgLevel = new ClassicMultigrid() {
                        TerminationCriterion = ((iter, r0_l2, r_l2) => iter <= 1) // termination controlled by top level PCG
                    };

                    MultigridChain[iLevel] = MgLevel;

                    ISolverSmootherTemplate pre, pst;
                    if (iLevel > 0) {

                        Schwarz swz1 = new Schwarz() {
                            FixedNoOfIterations = 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsOnCurrentProcess = NoOfBlocks
                            },
                            Overlap = 0 // overlap does **NOT** seem to help
                        };
                        
                        SoftPCG pcg1 = new SoftPCG() {
                            TerminationCriterion = ((iter, r0_l2, r_l2) => iter <= 5)
                        };

                        SoftPCG pcg2 = new SoftPCG() {
                            TerminationCriterion = ((iter, r0_l2, r_l2) => iter <= 5)
                        };

                        SetLinItCallback(swz1, IsLinPrecond: true);
                        SetLinItCallback(pcg1, IsLinPrecond: true);
                        SetLinItCallback(pcg2, IsLinPrecond: true);

                        var preChain = new ISolverSmootherTemplate[] { swz1, pcg1 };
                        var pstChain = new ISolverSmootherTemplate[] { swz1, pcg2 };

                        pre = new SolverSquence() { SolverChain = preChain };
                        pst = new SolverSquence() { SolverChain = pstChain };
                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++
                        // top level - use only iterative (non-direct) solvers
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++

                        pre = new BlockJacobi() {
                            NoOfIterations = 3,
                            omega = 0.5
                        };

                        pst = new BlockJacobi() {
                            NoOfIterations = 3,
                            omega = 0.5
                        };

                        SetLinItCallback(pre, IsLinPrecond: true);
                        SetLinItCallback(pst, IsLinPrecond: true);
                    }

                    MgLevel.PreSmoother = pre;
                    MgLevel.PostSmoother = pst;
                }

                if (iLevel > 0) {
                    ((ClassicMultigrid)(MultigridChain[iLevel - 1])).CoarserLevelSolver = MultigridChain[iLevel];
                }

                if (useDirect) {
                    Console.WriteLine("MG: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }
                //Current = Current.CoarserLevel;
            } // end of level loop

            return MultigridChain[0];
        }


        private ISolverSmootherTemplate MakeOrthoNormMGDecomp( int MSLength, ISolverSmootherTemplate subsmootherchain, ISolverSmootherTemplate toplevelsmootherchain)
        {

            int MGDepth = Math.Min(MSLength, m_lc.NoOfMultigridLevels);

            List<ISolverSmootherTemplate> MG_list = new List<ISolverSmootherTemplate>();

            int DirectKickIn = m_lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            MaxUsedMGLevel = GetMaxMGLevel(DirectKickIn, MGDepth, _LocalDOF);
            SetLinItCallback(subsmootherchain, true);

            foreach(var solversmoother in ((SolverSquence)subsmootherchain).SolverChain)
                SetLinItCallback(solversmoother, true);

            foreach (var toplevelsmoother in ((SolverSquence)toplevelsmootherchain).SolverChain)
                SetLinItCallback(toplevelsmoother, true);

            for (
                
                int iDepth = MaxUsedMGLevel; iDepth >= 0; iDepth--)
            {
                ISolverSmootherTemplate solvertoinject;

                if (iDepth == 0) {
                    solvertoinject = toplevelsmootherchain.CloneAs();
                }
                else if (iDepth == MaxUsedMGLevel) {
                    solvertoinject = new DirectSolver()
                    {
                        WhichSolver = DirectSolver._whichSolver.MUMPS,
                        TestSolution = false
                    };
                    SetLinItCallback(solvertoinject, true);
                } else {

                    solvertoinject = subsmootherchain.CloneAs();
                }
               
                ISolverSmootherTemplate MG = BareMGSquence(iDepth, solvertoinject, subsmootherchain, toplevelsmootherchain);
                MG_list.Add(MG);
            }

            ISolverSmootherTemplate orthosolve = new OrthonormalizationScheme()
            {
                PrecondS = MG_list.ToArray(),
                MaxKrylovDim = m_lc.MaxKrylovDim,
                MaxIter = 30,
                Tolerance = m_lc.ConvergenceCriterion
            };
            SetLinItCallback(orthosolve, true);

            return orthosolve;
        }



        /// <summary>
        /// I want to ride my k-cycle ...
        /// </summary>
        ISolverSmootherTemplate KcycleMultiSchwarz(int MaxMGDepth, Func<int,int> SchwarzblockSize) {

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();
            int maxDG = getMaxDG(0, 0);
            bool UsePmg = false; //enables larger Schwarz Blocks

            var LocalDOF4directSolver = _LocalDOF;
            // if we use lvlpmg in Sblocks, we can have less and larger blocks ...
            if (m_lc.pMaxOfCoarseSolver < maxDG && UsePmg) {
                LocalDOF4directSolver = GetLocalDOF(m_lc.pMaxOfCoarseSolver);
            }

            int DirectKickIn = m_lc.TargetBlockSize; // 10'000 DOF seemed to be optimal at lowest lvl
            int LocSysSizeZeroLvl = _LocalDOF[0];

            //Precompute Total Systemsize
            var TotalSizeOnLevel = new int[MaxMGDepth];
            for (int iLevel = 0; iLevel < MaxMGDepth; iLevel++) {
                TotalSizeOnLevel[iLevel] = _LocalDOF[iLevel].MPISum();
            }

            for (int iLevel = 0; iLevel < MaxMGDepth; iLevel++) {
                MaxUsedMGLevel = iLevel;
                var SysSize = TotalSizeOnLevel[iLevel];
                var PrevSize = TotalSizeOnLevel[Math.Max(iLevel - 1, 0)];
                double SizeFraction = (double)LocalDOF4directSolver[iLevel] / (double)SchwarzblockSize(iLevel);
                Console.WriteLine("DOF on L{0}: {1}",iLevel,SysSize);
                if (SizeFraction < 1 && iLevel == 0) {
                    Console.WriteLine($"WARNING: local system size ({LocalDOF4directSolver[iLevel]}) < Schwarz-Block size ({SchwarzblockSize(iLevel)});");
                    Console.WriteLine($"resetting local number of Schwarz-Blocks to 1.");
                }
                int LocalNoOfSchwarzBlocks = Math.Max(1, (int)Math.Floor(SizeFraction));
                int TotalNoOfSchwarzBlocks = LocalNoOfSchwarzBlocks.MPISum();
                SetQuery("GlobalSblocks at Lvl" + iLevel, TotalNoOfSchwarzBlocks, true);
                SetQuery("SblockSize at Lvl"+ iLevel, SchwarzblockSize(iLevel), true);

                bool useDirect = false;
                // It has to be ensured, that directKickin takes place on all ranks at same level
                // therefore only global defined criterion have to be used here !!!
                useDirect |= (SysSize < DirectKickIn);
                //useDirect |= (double)PrevSize / (double)SysSize < 1.5 && SysSize < 50000; // degenerated MG-Agglomeration, because too few candidates
                useDirect |= iLevel == m_lc.NoOfMultigridLevels - 1;
                useDirect |= TotalNoOfSchwarzBlocks < MPIsize;
                useDirect |= iLevel == MaxMGDepth - 1; // otherwise error, due to missing coarse solver
                useDirect = useDirect.MPIOr();

                if (useDirect)
                    Console.WriteLine("KcycleMultiSchwarz: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("KcycleMultiSchwarz: lv {0}, no of blocks {1} : ", iLevel, TotalNoOfSchwarzBlocks);

                Func<int, int, int, bool> SmootherCaching = delegate (int Iter, int MgLevel, int iBlock) {
                    //return Iter >= ((MaxMGLevel - MgLevel) * 3) * (1);
                    return true;
                };

                Func<int, int, bool> CoarseCaching = delegate (int Iter, int MgLevel) {
                    //return Iter >= MaxMGDepth+1;
                    return true;
                };

                ISolverSmootherTemplate levelSolver;

                if (useDirect) {
                    levelSolver = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                        SolverVersion = Parallelism.OMP,
                        TestSolution = false,
                        ActivateCaching = CoarseCaching,
                    };

                    //levelSolver = new Schwarz() {
                    //    FixedNoOfIterations = 1,
                    //    CoarseSolver = null,
                    //    m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                    //        NoOfPartsOnCurrentProcess = LocalNoOfSchwarzBlocks,
                    //    },
                    //    ActivateCachingOfBlockMatrix = delayedCaching,
                    //    Overlap = 1, // overlap seems to help; more overlap seems to help more
                    //    EnableOverlapScaling = true,
                    //};

                    //levelSolver = new DirectSolver() {
                    //    WhichSolver = DirectSolver._whichSolver.MUMPS,
                    //    SolverVersion = Parallelism.SEQ,
                    //    TestSolution = false
                    //};
                } else {


                    var smoother1 = new Schwarz() {
                        FixedNoOfIterations = 1,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = LocalNoOfSchwarzBlocks,
                        },
                        ActivateCachingOfBlockMatrix = SmootherCaching,
                        Overlap = 1, // overlap seems to help; more overlap seems to help more
                        EnableOverlapScaling = true,
                    };

                    levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother1,
                        m_omega = 1,
                    };

                    if (iLevel > 0) {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    } else {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= m_lc.MaxSolverIterations && r > r0 * m_lc.ConvergenceCriterion + m_lc.ConvergenceCriterion;
                    }

                    /*
                    ((OrthonormalizationMultigrid)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " + iter + " : " + renorm);
                        };
                    */

                    // Extended Multigrid Analysis
                    //((OrthonormalizationMultigrid)levelSolver).IterationCallback += MultigridAnalysis;                    
                }
                SolverChain.Add(levelSolver);

                if (iLevel > 0) {

                    ((OrthonormalizationMultigrid)(SolverChain[iLevel - 1])).CoarserLevelSolver = levelSolver;

                }

                if (useDirect) {
                    Console.WriteLine("INFO: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }
            }

            return SolverChain[0];
        }

        /// <summary>
        /// I want to ride my k-cycle ...
        /// </summary>
        ISolverSmootherTemplate KcycleMultiILU(int MaxMGDepth, Func<int,int> SchwarzblockSize) {

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();
            int maxDG = getMaxDG(0, 0);
            bool UsePmg = false; //enables larger Schwarz Blocks

            var LocalDOF4directSolver = _LocalDOF;
            // if we use lvlpmg in Sblocks, we can have less and larger blocks ...
            if (m_lc.pMaxOfCoarseSolver < maxDG && UsePmg) {
                LocalDOF4directSolver = GetLocalDOF(m_lc.pMaxOfCoarseSolver);
            }

            int DirectKickIn = m_lc.TargetBlockSize; // 10'000 DOF seemed to be optimal at lowest lvl
            int LocSysSizeZeroLvl = _LocalDOF[0];
            

            for (int iLevel = 0; iLevel < MaxMGDepth; iLevel++) {
                MaxUsedMGLevel = iLevel;
                double SizeFraction = (double)LocalDOF4directSolver[iLevel] / (double)SchwarzblockSize(iLevel);
                int SysSize = _LocalDOF[iLevel].MPISum();
                Console.WriteLine("DOF on L{0}: {1}",iLevel,SysSize);
                if (SizeFraction < 1 && iLevel == 0) {
                    Console.WriteLine($"WARNING: local system size ({LocalDOF4directSolver[iLevel]}) < Schwarz-Block size ({SchwarzblockSize(iLevel)});");
                    Console.WriteLine($"resetting local number of Schwarz-Blocks to 1.");
                }

                bool useDirect = false;
                // It has to be ensured, that directKickin takes place on all ranks at same level
                // therefore only global criterion have to be used here !!!
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == m_lc.NoOfMultigridLevels - 1;
                useDirect = useDirect.MPIOr();

                if(useDirect)
                    Console.WriteLine("KcycleMultiILU: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("KcycleMultiILU: lv {0}, ", iLevel);

                ISolverSmootherTemplate levelSolver;
                if (useDirect) {
                    levelSolver = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                        TestSolution = false
                    };

                } else {

                    //var smoother1 = new DirectSolver() {
                    //    TestSolution = true
                    //};

                    var smoother1 = new HypreILU() {
                        LocalPreconditioning = true
                    };


                    levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother1,
                        m_omega = 1,
                    };

                    if (iLevel > 0) {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    } else {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= m_lc.MaxSolverIterations && r>r0*m_lc.ConvergenceCriterion;
                    }

                    /*
                    ((OrthonormalizationMultigrid)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " + iter + " : " + renorm);
                        };
                    */

                    // Extended Multigrid Analysis
                    //((OrthonormalizationMultigrid)levelSolver).IterationCallback += MultigridAnalysis;                    

                }
                SolverChain.Add(levelSolver);

                if (iLevel > 0) {

                    ((OrthonormalizationMultigrid)(SolverChain[iLevel - 1])).CoarserLevelSolver = levelSolver;

                }

                if (useDirect) {
                    Console.WriteLine("INFO: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }
            }

            return SolverChain[0];
        }


        ISolverSmootherTemplate expKcycleSchwarz(int MaxMGDepth, Func<int, int> SchwarzblockSize) {

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();

            int DirectKickIn = m_lc.TargetBlockSize; // 10'000 DOF seemed to be optimal at lowest lvl
            int LocSysSizeZeroLvl = _LocalDOF[0];
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize ({0}) < smallest blocksize of MG ({1}).", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }
            if (MaxMGDepth < 1)
                throw new ArgumentException("ERROR: At least two multigrid levels are required.");

            for (int iLevel = 0; iLevel < MaxMGDepth; iLevel++) {
                MaxUsedMGLevel = iLevel;
                double SizeFraction = (double)_LocalDOF[iLevel] / (double)SchwarzblockSize(iLevel);
                int SysSize = _LocalDOF[iLevel].MPISum();
                if (SizeFraction < 1 && iLevel == 0)
                    Console.WriteLine("WARNING: Schwarzblock size ({0}) exceeds local system size ({1}); \n resetting local number of Schwarzblocks to 1.", SchwarzblockSize, _LocalDOF[iLevel]);
                int LocalNoOfSchwarzBlocks = Math.Max(1, (int)Math.Ceiling(SizeFraction));
                int TotalNoOfSchwarzBlocks = LocalNoOfSchwarzBlocks.MPISum();
                SetQuery("GlobalSblocks at Lvl" + iLevel, TotalNoOfSchwarzBlocks, true);
                SetQuery("SblockSize at Lvl" + iLevel, SchwarzblockSize(iLevel), true);

                bool useDirect = false;
                // It has to be ensured, that directKickin takes place on all ranks at same level
                // therefore only global criterion have to be used here !!!
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == m_lc.NoOfMultigridLevels - 1;
                useDirect |= TotalNoOfSchwarzBlocks < MPIsize;

                if (useDirect && iLevel == 0)
                    Console.WriteLine("WARNING: You are using the direct solver. Recommendations: \n\tRaise the Number of Multigridlevels\n\tLower the target blocksize");

                if (useDirect)
                    Console.WriteLine("KcycleMultiSchwarz: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("KcycleMultiSchwarz: lv {0}, no of blocks {1} : ", iLevel, TotalNoOfSchwarzBlocks);

                ISolverSmootherTemplate levelSolver;

                if (useDirect) {
                    levelSolver = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                        TestSolution = false
                    };
                    //levelSolver = new SoftGMRES() {
                    //    MaxKrylovDim = 100,
                    //    TerminationCriterion = (int iter, double r0, double r) => iter <= 100,
                    //};

                } else {


                    var solve1 = new Schwarz() {
                        FixedNoOfIterations = 1,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = LocalNoOfSchwarzBlocks
                            //NoOfPartsOnCurrentProcess = 4
                        },
                        Overlap = 1, // overlap seems to help; more overlap seems to help more
                        EnableOverlapScaling = true,
                    };

                    //var solve1 = new HypreILU();

                    //var solve1 = new BlockJacobi() { omega = 0.5 };

                    //var smoother1 = new SolverSquence() { SolverChain = new ISolverSmootherTemplate[] { solve1, solve2 } };
                    //var smoother2 = new SolverSquence() { SolverChain = new ISolverSmootherTemplate[] { solve1, solve2 } };

                    //if (iLevel == 0) SetQuery("KcycleSchwarz:XdgCellsToLowBlock", ((Schwarz)smoother1).CoarseSolveOfCutcells ? 1 : 0, true);
                    //if (iLevel == 0) SetQuery("KcycleSchwarz:OverlapON", ((Schwarz)smoother1).EnableOverlapScaling ? 1 : 0, true);
                    //if (iLevel == 0) SetQuery("KcycleSchwarz:OverlapScale", ((Schwarz)smoother1).Overlap, true);

                    levelSolver = new kcycle() {
                        PreSmoother = solve1,
                        PostSmoother = solve1,
                    };

                    if (iLevel > 0) {
                        ((kcycle)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    } else {
                        SetLinItCallback(levelSolver, true);
                        ((kcycle)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    }

                    /*
                    ((kcycle)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " + iter + " : " + renorm);
                        };
                    */
                    // Extended Multigrid Analysis
                    //((OrthonormalizationMultigrid)levelSolver).IterationCallback += MultigridAnalysis;                    

                }
                SolverChain.Add(levelSolver);

                if (iLevel > 0) {

                    ((kcycle)(SolverChain[iLevel - 1])).CoarserLevelSolver = levelSolver;

                }

                if (useDirect) {
                    Console.WriteLine("INFO: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }
            }

            return SolverChain[0];
        }

        ISolverSmootherTemplate KcycleMultiSchwarz_4Rheology() {

            // my tests show that the ideal block size may be around 10'000
            int DirectKickIn = m_lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();

            Console.WriteLine("experimental MG configuration for rheology");

            int MPIsize = this.MPIsize;

            for (int iLevel = 0; iLevel < m_lc.NoOfMultigridLevels; iLevel++) {
                MaxUsedMGLevel = iLevel;
                int SysSize = _LocalDOF[iLevel].MPISum();
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                //useDirect |= (SysSize < DirectKickIn);
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (iLevel == 0) {
                    useDirect = false;
                    NoOfBlocks = 8 / MPIsize;
                } else {
                    break;
                    //useDirect = true;
                }

                if (useDirect)
                    Console.WriteLine("   KcycleMultiSchwarz: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("   KcycleMultiSchwarz: lv {0}, no of blocks {1} : ", iLevel, NoOfBlocks);

                ISolverSmootherTemplate levelSolver;
                if (useDirect) {
                    levelSolver = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                        TestSolution = false
                    };
                } else {

                   
                    var smoother1 = new Schwarz() {
                        FixedNoOfIterations = 1,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsOnCurrentProcess = NoOfBlocks
                        },
                        //m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                        //    Depth = 1
                        //},
                        Overlap = 2, 
                        EnableOverlapScaling = false,
                        UsePMGinBlocks = false
                    };

                   
                    var smoother2 = smoother1;
                    /*
                    var smoother2 = new Schwarz() {
                        FixedNoOfIterations = 1,
                        CoarseSolver = null,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy {
                            NoOfPartsPerProcess = NoOfBlocks * 2
                        },
                        Overlap = 2,
                        EnableOverlapScaling = false,
                        UsePMGinBlocks = false
                    };
                    */

                    var CoarseSolver = new LevelPmg() {
                        UseHiOrderSmoothing = true,
                        OrderOfCoarseSystem = m_lc.pMaxOfCoarseSolver,
                    };

                    levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother2,
                        CoarserLevelSolver = CoarseSolver
                    };

                    if (iLevel > 0) {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    }

                    /*
                    ((OrthonormalizationMultigrid)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " + iter + " : " + renorm);
                        };
                    */

                }
                SolverChain.Add(levelSolver);

                if (iLevel > 0) {

                    ((OrthonormalizationMultigrid)(SolverChain[iLevel - 1])).CoarserLevelSolver = levelSolver;

                }

                if (useDirect) {
                    Console.WriteLine("Kswz: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, SysSize, DirectKickIn);
                    break;
                }

                //Current = Current.CoarserLevel;
            }



            return SolverChain[0];
        }
        #endregion

        #region check 'n' clear methods
        /// <summary>
        /// clears overgiven selfmade solvers
        /// </summary>
        public void Clear() {
            this.m_linsolver = null;
            this.m_nonlinsolver = null;
            this.m_precond = null;
        }

        private void Check_NonLinearSolver(NonlinearSolver NLSolver) {
            if (NLSolver is Newton newtonsolver)
                Check_Newton(newtonsolver);
            if (NLSolver is FixpointIterator picardsolver)
                Check_Picard(picardsolver);
        }

        private void Check_Newton(Newton NewtonSolver) {
            bool check = true;
            check = m_nc.ConvergenceCriterion == NewtonSolver.ConvCrit &&
            m_nc.MaxSolverIterations == NewtonSolver.MaxIter &&
            m_nc.MinSolverIterations == NewtonSolver.MinIter;
            //m_nc.PrecondSolver.Equals(NewtonSolver.Precond);
            Debug.Assert(check);
        }

        private void Check_Picard(FixpointIterator FPSolver) {
            bool check = true;
            check = m_nc.ConvergenceCriterion == FPSolver.ConvCrit &&
            m_nc.MaxSolverIterations == FPSolver.MaxIter &&
            m_nc.MinSolverIterations == FPSolver.MinIter &&
            //m_nc.PrecondSolver.Equals(FPSolver.Precond) &&
            m_nc.UnderRelax == FPSolver.UnderRelax;
            Debug.Assert(check);
        }


        /// <summary>
        /// Checks overgiven selfmade linear solver
        /// </summary>
        /// <returns></returns>
        private void Check_linsolver(ISolverSmootherTemplate solver) {
            int DirectKickIn = m_lc.TargetBlockSize;

            if (MaxUsedMGLevel >= 2) { //indicates that a mutligrid solver was build ...
                if (DirectKickIn < _LocalDOF.Last()) {
                    Console.WriteLine("WARNING: target blocksize ({0}) < smallest blocksize of MG ({1}).", DirectKickIn, _LocalDOF.Last());
                    Console.WriteLine("\tYour Choice of LinearSolverConfig.TargetBlockSize will not influence the MG depth. ");
                    Console.WriteLine("\t Choose LinearSolverConfig.TargetBlockSize > {0} or Raise LinearSolverConfig.NoOfMultigridLevels", _LocalDOF.Last());
                    DirectKickIn = _LocalDOF.Last();
                }
            }

            switch (m_lc.SolverCode) {
                //case LinearSolverCode.exp_Kcycle_schwarz:

                //    Schwarz kcycleSchwarz = null;
                //    if (solver.GetType() == typeof(DirectSolver))
                //        break; //we meet PARDISO here, because no MG descend
                //    try {
                //        kcycleSchwarz = (Schwarz)((OrthonormalizationMultigrid)solver).PreSmoother;
                //    } catch (Exception e) {
                //        throw new ApplicationException("someone messed up kcycle settings");
                //    }

                //    CompareAttributes("FixedNoOfIterations", 1, kcycleSchwarz.FixedNoOfIterations);
                //    Debug.Assert(kcycleSchwarz.CoarseSolver==null);
                //    CompareAttributes("m_BlockingStrategy", typeof(Schwarz.METISBlockingStrategy), kcycleSchwarz.m_BlockingStrategy.GetType());
                //    CompareAttributes("Overlap", 1, kcycleSchwarz.Overlap);
                //    CompareAttributes("EnableOverlapScaling", true, kcycleSchwarz.EnableOverlapScaling);
                //    CompareAttributes("UsePMGinBlocks", true, kcycleSchwarz.UsePMGinBlocks);
                //    break;
                case LinearSolverCode.exp_gmres_levelpmg:
                    LevelPmg TGP = null;
                    try {
                        TGP = (LevelPmg)((SoftGMRES)solver).Precond;
                    } catch (Exception) {
                        throw new ApplicationException("levelpmg setting is messed up");
                    }

                    CompareAttributes("UseHiOrderSmoothing", true, TGP.UseHiOrderSmoothing);
                    CompareAttributes("CoarseLowOrder", 1, TGP.OrderOfCoarseSystem);
                    CompareAttributes("UseDiagonalPmg", true, TGP.UseDiagonalPmg);
                    break;
                case LinearSolverCode.classic_pardiso:
                    DirectSolver sparsesolver = null;
                    try {
                        sparsesolver = (DirectSolver)solver;
                    } catch (Exception) {
                        throw new ApplicationException("someone messed up classic pardiso settings");
                    }

                    if (sparsesolver.WhichSolver != DirectSolver._whichSolver.PARDISO)
                        throw new ApplicationException("someone messed up classic pardiso settings");
                    CompareAttributes("SolverVersion", Parallelism.OMP.ToString(), sparsesolver.SolverVersion.ToString());
                    break;
            }

        }

        private void CompareAttributes(string property, string defaultatt, string checkatt) {
            if (defaultatt != checkatt)
                Console.WriteLine("WARNING : {0} is {1}, default is {2}", property, checkatt, defaultatt);
        }

        private void CompareAttributes(string property, int defaultatt, int checkatt) {
            if (defaultatt != checkatt)
                Console.WriteLine("WARNING : {0} is {1}, default is {2}", property, checkatt, defaultatt);
        }

        private void CompareAttributes(string property, Type defaultatt, Type checkatt) {
            if (defaultatt != checkatt)
                Console.WriteLine("WARNING : {0} is {1}, default is {2}", property, checkatt, defaultatt);
        }

        private void CompareAttributes(string property, bool defaultatt, bool checkatt) {
            if (defaultatt != checkatt)
                Console.WriteLine("WARNING : {0} is {1}, default is {2}", property, checkatt, defaultatt);
        }

        #endregion

    }
}