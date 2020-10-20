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

        private int m_MaxMGLevel = -1;

        /// <summary>
        /// Maximal used Mg Level is saved here.
        /// This will be saved in queries
        /// </summary>
        private int MaxMGLevel {
            get { return m_MaxMGLevel; }
            set { m_MaxMGLevel = Math.Max(value, m_MaxMGLevel); }
        }

        private bool SetQuery(string desciption, double value, bool setonce) {
            if (m_qh != null) {
                m_qh.ValueQuery(desciption, value, !setonce);
            }
            return m_qh != null;
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
            if (m_linsolver != null) {
                m_lc.SolverCode = LinearSolverCode.selfmade;
            }

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

                    if (linsolver.GetType() == typeof(SoftGMRES)) {
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
                            printLambda = nc.printLambda,
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
                            printLambda = nc.printLambda,
                            Globalization = nc.Globalization,
                            ApproxJac = Newton.ApproxInvJacobianOptions.ExternalSolver,
                            Precond = linsolver,
                            ConvCrit = nc.ConvergenceCriterion,
                            constant_newton_it = nc.constantNewtonIterations
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
            int[] LocalDOF = GetLocalDOF(MultigridBasis, MultigridOperatorConfig);
            int NoOfBlocks = (int)Math.Max(1, Math.Round(LocalDOF[0] / (double)lc.TargetBlockSize));
            //int SpaceDim = MultigridSequence[0].SpatialDimension;
            int MaxMGDepth = MultigridBasis.Count();

            switch (lc.SolverCode) {

                //no preconditioner used ...
                case LinearSolverCode.classic_mumps:
                case LinearSolverCode.classic_pardiso:
                case LinearSolverCode.classic_cg:
                case LinearSolverCode.exp_softgmres:
                case LinearSolverCode.exp_Kcycle_schwarz:
                case LinearSolverCode.exp_decomposedMG_OrthoScheme:
                case LinearSolverCode.exp_Kcycle_schwarz_4Rheology:
                case LinearSolverCode.exp_AS:
                case LinearSolverCode.exp_AS_MG:
                case LinearSolverCode.automatic:
                    precond[0] = null;
                    break;

                case LinearSolverCode.exp_gmres_Schur:
                    precond[0] = new SchurPrecond() {
                        SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                    };
                    break;

                case LinearSolverCode.exp_gmres_Simple:
                    precond[0] = new SchurPrecond() {
                        SchurOpt = SchurPrecond.SchurOptions.SIMPLE
                    };
                    break;


                case LinearSolverCode.exp_gmres_AS_MG:
                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = lc.NoOfMultigridLevels-1,
                        },
                        CoarseSolver = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS,    //PARDISO
                            LinConfig = lc
                        },

                        Overlap = 1
                    };
                    break;


                case LinearSolverCode.exp_gmres_localPrec:
                    precond[0] = new LocalizedOperatorPrec() {
                        m_dt = lc.exp_localPrec_Min_dt,
                        m_muA = lc.exp_localPrec_muA,
                    };
                    break;
                case LinearSolverCode.exp_gmres_AS:

                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks,
                        },
                        CoarseSolver = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS,
                            LinConfig = lc
                        },
                        Overlap = 1
                    };
                    break;

                case LinearSolverCode.exp_softpcg_schwarz:

                    precond[0] = new Schwarz() {
                        FixedNoOfIterations = 1,
                        CoarseSolver = null,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy {
                            NoOfPartsPerProcess = NoOfBlocks
                        },
                        Overlap = 1
                    };
                    break;

                case LinearSolverCode.exp_softpcg_schwarz_directcoarse:

                    precond[0] = new Schwarz() {
                        CoarseSolver = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS
                        },
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks
                        },
                        Overlap = 1,
                    };
                    break;

                case LinearSolverCode.exp_gmres_levelpmg:
                    precond[0] = new LevelPmg() { UseHiOrderSmoothing = true, CoarseLowOrder = 1, AssignXdGCellsToLowBlocks = true};
                    SetQuery("XdgCellsToLowBlock", ((LevelPmg)precond[0]).AssignXdGCellsToLowBlocks ? 1 : 0, true);
                    break;

                case LinearSolverCode.exp_gmres_schwarz_pmg:
                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks
                        },
                        CoarseSolver = null,
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

                    precond[0] = My_MG_Precond(lc, LocalDOF, MaxMGDepth, _prechain, _postchain, toppre, toppst);
                    break;

                case LinearSolverCode.exp_OrthoS_pMG:
                    precond = new ISolverSmootherTemplate[]{
                            new LevelPmg() {
                                UseHiOrderSmoothing =true,
                                CoarseLowOrder=1
                            },

                            new Schwarz() {
                                FixedNoOfIterations = 1,
                                CoarseSolver = null,
                                Overlap=1,
                                m_BlockingStrategy = new Schwarz.METISBlockingStrategy()          {
                                    NoOfPartsPerProcess = NoOfBlocks
                                },
                                UsePMGinBlocks=false,
                                EnableOverlapScaling=false,
                                CoarseLowOrder=1,
                            },
                    };
                    break;

                case LinearSolverCode.selfmade:
                    Console.WriteLine("INFO: Selfmade Preconditioner is used!");
                    precond[0] = m_precond;
                    break;




                default:
                    throw new NotImplementedException("Preconditioner not available");
            }


            foreach (var pre in precond) {
                if (pre == null)
                    continue;
                if(lc.verbose)
                    GenerateInfoMessageAtSetup(pre,NoOfBlocks);
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
            int[] LocalDOF = GetLocalDOF(MultigridBasis, MultigridOperatorConfig);
            int NoCellsLoc = MultigridBasis.First()[0].AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int NoOfBlocks = (int)Math.Max(1, Math.Round(LocalDOF[0] / (double)lc.TargetBlockSize));
            int SpaceDim = MultigridBasis.First()[0].AggGrid.ParentGrid.SpatialDimension;
            int MaxMGDepth = MultigridBasis.Count();

            switch (lc.SolverCode) {
                case LinearSolverCode.automatic:
                    if (m_nc != null) {
                       templinearSolve = AutomaticSolver(lc, LocalDOF, SpaceDim, NoCellsLoc,
                           precond);
                    }
                    break;

                case LinearSolverCode.classic_mumps:
                    templinearSolve = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.MUMPS,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverCode.classic_pardiso:
                    templinearSolve = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                        SolverVersion = Parallelism.OMP,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverCode.exp_AS:

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks,
                        },
                        Overlap = 1,
                        CoarseSolver = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.PARDISO,
                            LinConfig = lc
                        }
                    };
                    break;

                case LinearSolverCode.exp_AS_MG:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = lc.NoOfMultigridLevels - 1
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
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
                case LinearSolverCode.exp_gmres_localPrec:
                case LinearSolverCode.exp_gmres_Schur:
                case LinearSolverCode.exp_gmres_Simple:
                    templinearSolve = new SoftGMRES() {
                        //m_Tolerance = lc.ConvergenceCriterion,
                        //m_MaxIterations = lc.MaxSolverIterations,
                        MaxKrylovDim = lc.MaxKrylovDim,
                        Precond = precond[0]
                    };
                    break;

                case LinearSolverCode.exp_Kcycle_schwarz:
                    Func<int, int> SblkSizeFunc = delegate (int iLevel) { return m_lc.TargetBlockSize;  };
                    templinearSolve = KcycleMultiSchwarz(MaxMGDepth, LocalDOF, SblkSizeFunc);
                    break;

                case LinearSolverCode.exp_Kcycle_schwarz_4Rheology:
                    templinearSolve = KcycleMultiSchwarz_4Rheology(lc, LocalDOF);
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

                    templinearSolve = MakeOrthoNormMGDecomp(lc,LocalDOF, MaxMGDepth, new SolverSquence() {SolverChain= subsmoother}, new SolverSquence() { SolverChain = topsmoother });
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

                case LinearSolverCode.selfmade:
                    Console.WriteLine("INFO: Selfmade LinearSolver is used!");
                    templinearSolve = m_linsolver;
                    break;

                default:
                    throw new NotImplementedException("Linear solver option not available");
            }
            Debug.Assert(templinearSolve != null);
            SetLinItCallback(templinearSolve,false);

            // set Queries here which apply to all solvers:
            SetQuery("UsedMGDepth", MaxMGLevel + 1, true);

            if (lc.verbose)
                Console.WriteLine("linear solver : {0}", templinearSolve.ToString().Split('.').Last());

            if(templinearSolve is IProgrammableTermination pt) {
                 // delegate to stop linear solver convergence
                bool LinearConvergence(int iter, double r0_l2, double r_l2) {
                    if (iter <= lc.MinSolverIterations)
                        return true;

                    if (iter > lc.MaxSolverIterations)
                        return false;

                    if (r_l2 < lc.ConvergenceCriterion)
                        return false; // terminate

                    return true; // keep running
                }
                pt.TerminationCriterion = LinearConvergence;
            }

            Check_linsolver(templinearSolve);

            return templinearSolve;
        }
        #endregion

        #region auxiliary methods
        private void GenerateInfoMessageAtSetup(ISolverSmootherTemplate solver, int NoOfBlocks) {
            string solvername = solver.ToString().Split('.').Last();
            Console.WriteLine("preconditioner : {0}", solvername);
            List<string> solverinfotxt = new List<string>();
            if (solver.GetType() == typeof(LevelPmg)) {
                solverinfotxt.Add("entries upto p=" + ((LevelPmg)solver).CoarseLowOrder + " are assigned to low order blocks");
            }
            if (solver.GetType() == typeof(Schwarz)) {
                if (((Schwarz)solver).UsePMGinBlocks) {
                    solverinfotxt.Add("pmg used in Schwarz Blocks");
                    solverinfotxt.Add("entries upto p=" + ((Schwarz)solver).CoarseLowOrder + " are assigned to low order blocks");
                }
                solverinfotxt.Add("Additive Schwarz w. direct coarse, No of blocks: " + NoOfBlocks);
            }
            foreach (string line in solverinfotxt)
                Console.WriteLine(solvername + " INFO:\t" + line);
        }

        private int[] GetLocalDOF(IEnumerable<AggregationGridBasis[]> MultigridBasis, ChangeOfBasisConfig[][] MGChangeOfBasis) {

            //This workaround takes into account, the wierd structure of the ChangeOfBasis-thing
            //prohibits out of bounds exception
            Func<int, int, int[]> getDGs = delegate (int iLevel, int iVar) {
                int tLevel = iLevel < MGChangeOfBasis.Length ? iLevel : MGChangeOfBasis.Length - 1;
                int tVar = iVar < MGChangeOfBasis[tLevel].Length ? iVar : MGChangeOfBasis[tLevel].Length - 1;
                return MGChangeOfBasis[tLevel][tVar].DegreeS;
            };

            var MGBasis = MultigridBasis.ToArray();
            int MGDepth = MGBasis.Length;
            int[] LocalDOF = new int[MGDepth];

            for (int iLevel = 0; iLevel < MGBasis.Length; iLevel++) {
                LocalDOF[iLevel] = 0;
                int NoOfCells = MGBasis[iLevel][0].AggGrid.iLogicalCells.NoOfLocalUpdatedCells;

                for (int iCell = 0; iCell < NoOfCells; iCell++) {
                    for (int iVar = 0; iVar < MGBasis[iLevel].Length; iVar++) {
                        int pmax = getDGs(iLevel, iVar)[0];
                        try {
                            LocalDOF[iLevel] += MGBasis[iLevel][iVar].GetLength(iCell, pmax);
                        }
                        catch (Exception e) {
                            Console.WriteLine("WARNING: internal error occured during DOF calculation. Using estimate instead, which might not be accurate in case of XDG");
                            return SimpleGetLocalDOF(MultigridBasis, MGChangeOfBasis);
                        }
                    }
                }
            }

            return LocalDOF;
        }

        private int[] SimpleGetLocalDOF(IEnumerable<AggregationGridBasis[]> MultigridBasis, ChangeOfBasisConfig[][] MGChangeOfBasis) {
            int NoOfLevels = MultigridBasis.Count();
            int[] DOFperCell = new int[NoOfLevels];
            int[] LocalDOF = new int[NoOfLevels];
            int counter = 0;
            

            for (int iLevel = 0; iLevel < DOFperCell.Length; iLevel++) {
                counter = iLevel;
                if (iLevel > NoOfLevels - 1)
                    counter = NoOfLevels - 1;
                foreach (var cob in MGChangeOfBasis[counter]) {
                    for (int iVar = 0; iVar < cob.VarIndex.Length; iVar++) {
                        int d = ((AggregationGridBasis)MultigridBasis.First()[iVar]).AggGrid.ParentGrid.SpatialDimension;
                        int p = cob.DegreeS[iVar];
                        switch (d) {
                            case 1:
                                DOFperCell[iLevel] += p + 1 + p + 1;
                                break;
                            case 2:
                                DOFperCell[iLevel] += (p * p + 3 * p + 2) / 2;
                                break;
                            case 3:
                                DOFperCell[iLevel] += (p * p * p + 6 * p * p + 11 * p + 6) / 6;
                                break;
                            default:
                                throw new Exception("wtf?Spacialdim=1,2,3 expected");
                        }
                    }
                }
                LocalDOF[iLevel] = ((AggregationGridBasis)MultigridBasis.First()[0]).AggGrid.iLogicalCells.NoOfLocalUpdatedCells* DOFperCell[iLevel];
            }
            return LocalDOF;
        }

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

        private int GetMPIsize {
            get {
                int MPIsize;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MPIsize);
                return MPIsize;
            }
        }

        private int GetMPIrank {
            get {
                int MPIrank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MPIrank);
                return MPIrank;
            }
        }

        #endregion

        #region callback methods
        private void FirstLineinCallBack() {
            if (m_Iterations == null) {
                string FirstLine = "NLinS-, Precond-, LinS-, total-Iterations : Type, SolverName, InfResi, Multigridlevel";
                m_Iterations = new int[4];
                Console.WriteLine(FirstLine);
            }
        }

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

        private void SetNonLinItCallback(NonlinearSolver nonlinsolver) {

            int _caseselect = 0;
            string _type = "NLinSolver";
            DefaultItCallback = GenerateDefaultCallback<NonlinearSolver>(_type, nonlinsolver, _caseselect);
            if (m_nc.verbose) {
                nonlinsolver.IterationCallback += DefaultItCallback;
            }
            nonlinsolver.IterationCallback += CustomizedCallback;
        }

        private void MultigridCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            int currentMGLevel = Mgop.LevelIndex;

            if (m_MG_Counter - currentMGLevel == +1) {
                double residualNormAf = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm after coarse grid correction
                ProlRes[currentMGLevel] = residualNormAf;
                Console.WriteLine("after Prolongation {0}<-{1}: {2}", currentMGLevel, currentMGLevel + 1, ProlRes[currentMGLevel] - ProlRes[currentMGLevel + 1]);
            }

            if (m_MG_Counter - currentMGLevel == 0 && currentMGLevel == MaxMGLevel) {
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
            MaxMGLevel = currentMGLevel;
        }
        #endregion

        #region special configs
        /// <summary>
        /// Automatic choice of linear solver depending on problem size, immersed boundary, polynomial degree, etc. In addition the nonlinearsolver config is considered as well.
        /// </summary>
        private ISolverSmootherTemplate AutomaticSolver( LinearSolverConfig lc, int[] LDOF, int Dim, int NoCellsLoc, ISolverSmootherTemplate[] PreCond) {

            NonLinearSolverConfig nc = null;
            if (m_nc != null) {
                nc = m_nc;
            } else {
                throw new ArgumentException("No NonlinearSolver specified");
            }

            var D = Dim;

            //int pV = Control.FieldOptions["VelocityX"].Degree;
            int pV = LDOF[0];
            int pP = pV - 1; //Control.FieldOptions["Pressure"].Degree;


            // Detecting variables for solver determination 

            var cellsLoc = NoCellsLoc;
            var cellsGlo = NoCellsLoc.MPISum();

            ISolverSmootherTemplate tempsolve = null;

                switch (D) {
                    case 1:
                        tempsolve = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS,
                            LinConfig = lc
                        };
                        break;
                    case 2:
                        tempsolve = new DirectSolver() {
                            WhichSolver = DirectSolver._whichSolver.MUMPS,
                            LinConfig = lc
                        };
                        break;
                    case 3:
                        var dofsPerCell3D = LDOF[0] / NoCellsLoc;
                        var dofsLoc = dofsPerCell3D * cellsLoc;
                        var dofsGlo = dofsPerCell3D * cellsGlo;

                        if (dofsGlo > 10000) {

                            if (lc.NoOfMultigridLevels < 2)
                                throw new ApplicationException("At least 2 Multigridlevels are required");

                            tempsolve = new SoftGMRES() {
                                MaxKrylovDim = lc.MaxKrylovDim,
                                //m_Tolerance = lc.ConvergenceCriterion,
                                TerminationCriterion = ((iter, r0_l2, r_l2) => r_l2 > lc.ConvergenceCriterion && iter < lc.MaxSolverIterations),
                                Precond = new Schwarz() {
                                    m_BlockingStrategy = new Schwarz.SimpleBlocking() {
                                        NoOfPartsPerProcess = (int)Math.Ceiling(dofsLoc / 6500.0),
                                    },
                                    Overlap = 1,
                                    CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                                },
                            };
                        } else {
                            tempsolve = new DirectSolver() {
                                WhichSolver = DirectSolver._whichSolver.MUMPS,
                                LinConfig = lc
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
        /// Determines a solver sequence depending on MGlevels
        /// </summary>
        private ISolverSmootherTemplate DetermineMGSquence(int MGlevels, LinearSolverConfig lc) {
            ISolverSmootherTemplate solver;
            if (MGlevels > 0) {
                solver = new ClassicMultigrid() { CoarserLevelSolver = DetermineMGSquence(MGlevels - 1, lc) };
            } else {
                solver = new DirectSolver() {
                    WhichSolver = DirectSolver._whichSolver.MUMPS,
                    LinConfig = lc
                };
            }
            return solver;
        }

        /// <summary>
        /// experimental. Is connected to Decomposed MG OrthoScheme. Can be deleted if not used anymore ...
        /// </summary>
        private ISolverSmootherTemplate BareMGSquence(int MGlevels, ISolverSmootherTemplate coarseSolver, ISolverSmootherTemplate smoother=null)
        {
            ISolverSmootherTemplate solver;
            if (MGlevels > 0)
            {
                solver = new ClassicMultigrid() {
                    CoarserLevelSolver = BareMGSquence(MGlevels - 1, coarseSolver, smoother),
                    PreSmoother= smoother.CloneAs(),
                    PostSmoother= smoother.CloneAs(),
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
                    CoarserLevelSolver = BareMGSquence(MGlevels - 1, coarseSolver, smoother),
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

        private ISolverSmootherTemplate My_MG_Precond(LinearSolverConfig _lc, int[] _LocalDOF, int MGLength, ISolverSmootherTemplate[] prechain, ISolverSmootherTemplate[] postchain, ISolverSmootherTemplate toplevelpre, ISolverSmootherTemplate toplevelpst)
        {

            int MGDepth = Math.Min(MGLength, _lc.NoOfMultigridLevels);

            bool isLinPrecond = true;
            

            int DirectKickIn = _lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            ISolverSmootherTemplate[] MultigridChain = new ISolverSmootherTemplate[MGDepth];
            for (int iLevel = 0; iLevel < MGDepth; iLevel++) {
                MaxMGLevel = iLevel;
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
                MaxMGLevel = iLevel;
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
                            CoarseSolver = null,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsPerProcess = NoOfBlocks
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


        private ISolverSmootherTemplate MakeOrthoNormMGDecomp(LinearSolverConfig _lc, int[] _LocalDOF, int MSLength, ISolverSmootherTemplate subsmootherchain, ISolverSmootherTemplate toplevelsmootherchain)
        {

            int MGDepth = Math.Min(MSLength, m_lc.NoOfMultigridLevels);

            List<ISolverSmootherTemplate> MG_list = new List<ISolverSmootherTemplate>();

            int DirectKickIn = _lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            MaxMGLevel = GetMaxMGLevel(DirectKickIn, MGDepth, _LocalDOF);
            SetLinItCallback(subsmootherchain, true);

            foreach(var solversmoother in ((SolverSquence)subsmootherchain).SolverChain)
                SetLinItCallback(solversmoother, true);

            foreach (var toplevelsmoother in ((SolverSquence)toplevelsmootherchain).SolverChain)
                SetLinItCallback(toplevelsmoother, true);

            for (
                
                int iDepth = MaxMGLevel; iDepth >= 0; iDepth--)
            {
                ISolverSmootherTemplate solvertoinject;

                if (iDepth == 0) {
                    solvertoinject = toplevelsmootherchain.CloneAs();
                }
                else if (iDepth == MaxMGLevel) {
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
                MaxKrylovDim = _lc.MaxKrylovDim,
                MaxIter = 30,
                Tolerance = _lc.ConvergenceCriterion
            };
            SetLinItCallback(orthosolve, true);

            return orthosolve;
        }



        /// <summary>
        /// I want to ride my k-cycle ...
        /// </summary>
        ISolverSmootherTemplate KcycleMultiSchwarz(int MaxMGDepth, int[] _LocalDOF, Func<int,int> SchwarzblockSize) {

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
                MaxMGLevel = iLevel;
                double SizeFraction = (double)_LocalDOF[iLevel] / (double)SchwarzblockSize(iLevel);
                int SysSize = _LocalDOF[iLevel].MPISum();
                if (SizeFraction < 1 && iLevel == 0)
                    Console.WriteLine("WARNING: Schwarzblock size ({0}) exceeds local system size ({1}); \n resetting local number of Schwarzblocks to 1.", SchwarzblockSize, _LocalDOF[iLevel]);
                int LocalNoOfSchwarzBlocks = Math.Max(1, (int)Math.Ceiling(SizeFraction));
                int TotalNoOfSchwarzBlocks = LocalNoOfSchwarzBlocks.MPISum();
                //SetQuery("LocalSblocks at Lvl" + iLevel, LocalNoOfSchwarzBlocks, true);
                SetQuery("GlobalSblocks at Lvl" + iLevel, TotalNoOfSchwarzBlocks, true);
                SetQuery("SblockSize at Lvl"+ iLevel, SchwarzblockSize(iLevel), true);

                bool useDirect = false;
                // It has to be ensured, that directKickin takes place on all ranks at same level
                // therefore only global criterions have to be used here !!!
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == m_lc.NoOfMultigridLevels - 1;
                useDirect |= TotalNoOfSchwarzBlocks < GetMPIsize;

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
                } else {

                    var smoother1 = new Schwarz() {
                        FixedNoOfIterations = 1,
                        CoarseSolver = null,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = LocalNoOfSchwarzBlocks
                        },
                        Overlap = 1, // overlap seems to help; more overlap seems to help more
                        EnableOverlapScaling = true,
                        UsePMGinBlocks = true,
                        CoarseSolveOfCutcells = true,
                    };

                    if (iLevel == 0) SetQuery("KcycleSchwarz:XdgCellsToLowBlock", ((Schwarz)smoother1).CoarseSolveOfCutcells ? 1 : 0, true);
                    if (iLevel == 0) SetQuery("KcycleSchwarz:OverlapON", ((Schwarz)smoother1).EnableOverlapScaling ? 1 : 0, true);
                    if (iLevel == 0) SetQuery("KcycleSchwarz:OverlapScale", ((Schwarz)smoother1).Overlap, true);

                    levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother1,
                        m_omega = 1,
                    };

                    if (iLevel > 0) {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    }

                    ((OrthonormalizationMultigrid)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " + iter + " : " + renorm);
                        };

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


        ISolverSmootherTemplate KcycleMultiSchwarz_4Rheology(LinearSolverConfig _lc, int[] _LocalDOF) {

            // my tests show that the ideal block size may be around 10'000
            int DirectKickIn = _lc.TargetBlockSize;
            if (DirectKickIn < _LocalDOF.Last()) {
                Console.WriteLine("WARNING: target blocksize: {0} < smallest blocksize of MG: {1}; Resetting target blocksize to: {1}", DirectKickIn, _LocalDOF.Last());
                DirectKickIn = _LocalDOF.Last();
            }

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();

            Console.WriteLine("experimental MG configuration for rheology");

            int MPIsize = GetMPIsize;

            for (int iLevel = 0; iLevel < _lc.NoOfMultigridLevels; iLevel++) {
                MaxMGLevel = iLevel;
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
                        CoarseSolver = null,
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks
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
                        CoarseLowOrder = 1
                    };

                    levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother2,
                        CoarserLevelSolver = CoarseSolver
                    };

                    if (iLevel > 0) {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    }


                    ((OrthonormalizationMultigrid)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " + iter + " : " + renorm);
                        };


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
            switch (m_lc.SolverCode) {
                case LinearSolverCode.exp_Kcycle_schwarz:

                    Schwarz kcycleSchwarz = null;
                    if (solver.GetType() == typeof(DirectSolver))
                        break; //we meet PARDISO here, because no MG descend
                    try {
                        kcycleSchwarz = (Schwarz)((OrthonormalizationMultigrid)solver).PreSmoother;
                    } catch (Exception e) {
                        throw new ApplicationException("someone messed up kcycle settings");
                    }

                    CompareAttributes("FixedNoOfIterations", 1, kcycleSchwarz.FixedNoOfIterations);
                    Debug.Assert(kcycleSchwarz.CoarseSolver==null);
                    CompareAttributes("m_BlockingStrategy", typeof(Schwarz.METISBlockingStrategy), kcycleSchwarz.m_BlockingStrategy.GetType());
                    CompareAttributes("Overlap", 1, kcycleSchwarz.Overlap);
                    CompareAttributes("EnableOverlapScaling", true, kcycleSchwarz.EnableOverlapScaling);
                    CompareAttributes("UsePMGinBlocks", true, kcycleSchwarz.UsePMGinBlocks);
                    break;
                case LinearSolverCode.exp_gmres_levelpmg:
                    LevelPmg TGP = null;
                    try {
                        TGP = (LevelPmg)((SoftGMRES)solver).Precond;
                    } catch (Exception e) {
                        throw new ApplicationException("levelpmg setting is messed up");
                    }

                    CompareAttributes("UseHiOrderSmoothing", true, TGP.UseHiOrderSmoothing);
                    CompareAttributes("CoarseLowOrder", 1, TGP.CoarseLowOrder);
                    CompareAttributes("UseDiagonalPmg", true, TGP.UseDiagonalPmg);
                    break;
                case LinearSolverCode.classic_pardiso:
                    DirectSolver sparsesolver = null;
                    try {
                        sparsesolver = (DirectSolver)solver;
                    } catch (Exception e) {
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