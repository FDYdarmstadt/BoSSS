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

namespace BoSSS.Solution {

    public class SolverFactory {

        /// <summary>
        /// This <see cref="SolverFactory"/> enables creation of linear and nonlinear solver objects. The configuration <see cref="LinearSolverConfig"/> and <see cref="NonLinearSolverConfig"/> can be set in Control-Files (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        public SolverFactory(NonLinearSolverConfig nc, LinearSolverConfig lc) {
            lc_check = lc.CloneAs();
            nc_check = nc.CloneAs();
            m_lc = lc;
            m_nc = nc;
        }

        /// <summary>
        /// This will return <code>linear</code> and <code>nonlinear</code> solver objects, which are configured according to <see cref="LinearSolverConfig"/> and <see cref="NonLinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="nonlinSolver"></param>
        /// <param name="linsolver"></param>
        /// <param name="ts_AssembleMatrixCallback"></param>
        /// <param name="ts_MultigridBasis"></param>
        /// <param name="ts_MultigridOperatorConfig"></param>
        /// <param name="ts_SessionPath"></param>
        /// <param name="ts_MGS"></param>
        public void GenerateNonLin(out NonlinearSolver nonlinSolver, out ISolverSmootherTemplate linsolver, OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis,  MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig, string ts_SessionPath, AggregationGridData[] ts_MGS) {

            if (m_nonlinsolver != null && (m_linsolver == null || m_precond == null))
                throw new NotImplementedException("an uncomplete nonlinear solver is overgiven.");

            if (m_nonlinsolver != null)
                m_nc.SolverCode = NonLinearSolverCode.selfmade;
            if (m_linsolver != null)
                m_lc.SolverCode = LinearSolverCode.selfmade;

            linsolver = null;
            nonlinSolver = null;
            ISolverSmootherTemplate precondsolver = null;

            linsolver = GenerateLinear_body(m_lc, ts_MGS, ts_MultigridOperatorConfig);
            
            var pretmp = GeneratePrecond(m_lc, ts_MGS, ts_MultigridOperatorConfig);//only interesting for GMRES-Newton, look into _body for more info
            if (pretmp != null) {
                precondsolver = pretmp[0];
            }
            Debug.Assert(linsolver != null);

            nonlinSolver = GenerateNonLin_body(ts_AssembleMatrixCallback, ts_MultigridBasis, m_nc, m_lc, linsolver, precondsolver, ts_MultigridOperatorConfig, ts_SessionPath);

            Debug.Assert(nonlinSolver != null);
            return;
        }

        /// <summary>
        /// This one is the method-body of <see cref="GenerateNonLin"/> and shall not be called from the outside. The parameters are mainly handed over to the NonLinearSolver object, which lives in <see cref="AdvancedSolvers.NonlinearSolver"/>.
        /// </summary>
        /// <param name="ts_AssembleMatrixCallback"></param>
        /// <param name="ts_MultigridBasis"></param>
        /// <param name="nc"></param>
        /// <param name="lc"></param>
        /// <param name="LinSolver"></param>
        /// <param name="PrecondSolver"></param>
        /// <param name="MultigridOperatorConfig"></param>
        /// <param name="SessionPath"></param>
        /// <returns></returns>
        private NonlinearSolver GenerateNonLin_body(OperatorEvalOrLin ts_AssembleMatrixCallback, IEnumerable<AggregationGridBasis[]> ts_MultigridBasis, NonLinearSolverConfig nc, LinearSolverConfig lc, ISolverSmootherTemplate LinSolver, ISolverSmootherTemplate PrecondSolver, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig, string SessionPath) {

            //Timestepper.Config_MultigridOperator;
            //Timestepper.SessionPath;

            // +++++++++++++++++++++++++++++++++++++++++++++
            // the nonlinear solvers:
            // +++++++++++++++++++++++++++++++++++++++++++++

            NonlinearSolver nonlinSolver;

            switch (nc.SolverCode) {
                case NonLinearSolverCode.Picard:

                    nonlinSolver = new FixpointIterator(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                        MaxIter = nc.MaxSolverIterations,
                        MinIter = nc.MinSolverIterations,
                        m_LinearSolver = LinSolver,
                        m_SessionPath = SessionPath, //is needed for Debug purposes, output of inter-timesteps
                        ConvCrit = nc.ConvergenceCriterion,
                        UnderRelax = nc.UnderRelax,
                    };
                    break;

                //Besides NonLinearSolverConfig Newton needs also LinearSolverConfig
                //Newton uses MUMPS as linearsolver by default
                case NonLinearSolverCode.Newton:

                    if (LinSolver.GetType() == typeof(SoftGMRES)) {
                        nonlinSolver = new Newton(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                            maxKrylovDim = lc.MaxKrylovDim,
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            ApproxJac = Newton.ApproxInvJacobianOptions.GMRES,
                            Precond = PrecondSolver,
                            GMRESConvCrit = lc.ConvergenceCriterion,
                            ConvCrit = nc.ConvergenceCriterion,
                            m_SessionPath = SessionPath,
                            UsePresRefPoint = nc.UsePresRefPoint,
                        };
                    } else {
                        nonlinSolver = new Newton(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                            maxKrylovDim = lc.MaxKrylovDim,
                            MaxIter = nc.MaxSolverIterations,
                            MinIter = nc.MinSolverIterations,
                            UsePresRefPoint = nc.UsePresRefPoint,
                            ApproxJac = Newton.ApproxInvJacobianOptions.DirectSolver,
                            linsolver = LinSolver,
                            Precond = PrecondSolver,
                            GMRESConvCrit = lc.ConvergenceCriterion,
                            ConvCrit = nc.ConvergenceCriterion,
                            m_SessionPath = SessionPath,
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
                        m_LinearSolver = LinSolver,
                        m_SessionPath = SessionPath, //is needed for Debug purposes, output of inter-timesteps
                        ConvCrit = nc.ConvergenceCriterion,
                        UnderRelax = nc.UnderRelax,
                    };
                    SetNonLinItCallback(myFixPoint);
                    var myNewton = new Newton(
                        ts_AssembleMatrixCallback,
                        ts_MultigridBasis,
                        MultigridOperatorConfig) {
                        maxKrylovDim = lc.MaxKrylovDim,
                        MaxIter = nc.MaxSolverIterations,
                        MinIter = nc.MinSolverIterations,
                        ApproxJac = Newton.ApproxInvJacobianOptions.DirectSolver,
                        linsolver = LinSolver,
                        Precond = PrecondSolver,
                        GMRESConvCrit = lc.ConvergenceCriterion,
                        ConvCrit = nc.ConvergenceCriterion,
                        m_SessionPath = SessionPath,
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
                    Console.WriteLine("Selfmade Nonlinear Solver used!");
                    nonlinSolver = m_nonlinsolver;
                    break;
                default:
                    throw new NotImplementedException();
            }

            Check_NonLinearSolver(nonlinSolver);

            SetNonLinItCallback(nonlinSolver);
#if DEBUG
            Console.WriteLine("nonlinear Solver: {0}", nc.SolverCode.ToString());
#endif
            Check_nonlinsolver();
            return nonlinSolver;
        }

        /// <summary>
        /// This will return a linear solver object, which is configured according to <see cref="LinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>). 
        /// </summary>
        /// <param name="linSolve">the linear solver object</param>
        /// <param name="ts_MultigridSequence"></param>
        /// <param name="ts_MultigridOperatorConfig"></param>
        public void GenerateLinear(out ISolverSmootherTemplate linSolve, AggregationGridData[] ts_MultigridSequence, MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig) {
            if (m_linsolver != null) {
                m_lc.SolverCode = LinearSolverCode.selfmade;
            }
            linSolve = GenerateLinear_body(m_lc, ts_MultigridSequence, ts_MultigridOperatorConfig);
            
            
            Debug.Assert(linSolve != null);
        }

        /// <summary>
        /// This will return a linear solver object, which is configured according to <see cref="LinearSolverConfig"/>, which can be adjusted from Controlfile (defined in <see cref="AppControl"/>).
        /// </summary>
        /// <param name="linSolve">the linear solver object</param>
        /// <param name="ts_MultigridSequence"></param>
        /// <param name="ts_MultigridOperatorConfig"></param>
        /// <param name="IterationCallbacks">insert callback specified by user here</param>
        public void GenerateLinear(out ISolverSmootherTemplate linSolve, AggregationGridData[] ts_MultigridSequence, MultigridOperator.ChangeOfBasisConfig[][] ts_MultigridOperatorConfig, List<Action<int, double[], double[], MultigridOperator>> IterationCallbacks) {
            IterationCallbacks.ForEach(ICB => this.CustomizedCallback += ICB);
            GenerateLinear(out linSolve, ts_MultigridSequence, ts_MultigridOperatorConfig);
        }


        private int[] GetLocalDOF(AggregationGridData[] MgSeq, MultigridOperator.ChangeOfBasisConfig[][] MgOpConfig) {

            //Calculate number of local DOF for every Multigridlevel
            int[] DOFperCell = new int[MgSeq.Length];
            int[] LocalDOF = new int[MgSeq.Length];
            int counter = 0;
            for (int iLevel = 0; iLevel < DOFperCell.Length; iLevel++) {
                counter = iLevel;
                if (iLevel > MgOpConfig.Length - 1)
                    counter = MgOpConfig.Length - 1;

                int d = MgSeq[iLevel].SpatialDimension;
                foreach (var variable in MgOpConfig[counter]) {
                    int p = variable.DegreeS[0];
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
                //CHANGE THIS ASAP: This is not right in case of XDG ... because blocks have truncated size at this code level
                LocalDOF[iLevel] = MgSeq[iLevel].CellPartitioning.LocalLength * DOFperCell[iLevel];
            }
            return LocalDOF;
        }

        private ISolverSmootherTemplate[] GeneratePrecond(LinearSolverConfig lc, AggregationGridData[] MultigridSequence, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig) {

            var precond = new ISolverSmootherTemplate[1];

            int[] LocalDOF = GetLocalDOF(MultigridSequence, MultigridOperatorConfig);

            ////these values are acquired for some solvers
            //int NoCellsLoc = MultigridSequence[0].CellPartitioning.LocalLength;
            //int NoCellsGlob = MultigridSequence[0].CellPartitioning.TotalLength;
            int NoOfBlocks = (int)Math.Max(1, Math.Round(LocalDOF[0] / (double)lc.TargetBlockSize));
            //int SpaceDim = MultigridSequence[0].SpatialDimension;
            int MultigridSeqLength = MultigridSequence.Length;

            switch (lc.SolverCode) {
                
                //no preconditioner used ...
                case LinearSolverCode.classic_mumps:
                case LinearSolverCode.classic_pardiso:
                case LinearSolverCode.classic_cg:
                case LinearSolverCode.exp_schwarz_directcoarse_overlap:
                case LinearSolverCode.exp_schwarz_directcoarse:
                case LinearSolverCode.exp_schwarz_MG_directcoarse:
                case LinearSolverCode.exp_schwarz_MG_directcoarse_overlap:
                case LinearSolverCode.exp_softgmres:
                case LinearSolverCode.exp_ILU:
                case LinearSolverCode.exp_Schur:
                case LinearSolverCode.exp_Simple:
                case LinearSolverCode.exp_AS_1000:
                case LinearSolverCode.exp_AS_5000:
                case LinearSolverCode.exp_AS_10000:
                case LinearSolverCode.exp_AS_MG:
                case LinearSolverCode.exp_localPrec:
                case LinearSolverCode.exp_direct_lapack:
                case LinearSolverCode.exp_Kcycle_schwarz:
                case LinearSolverCode.exp_decomposedMG_OrthoScheme:
                    precond = null;
                    break;
                case LinearSolverCode.automatic:
                    precond[0] = AutomaticPrecond(lc);
                    break;
                case LinearSolverCode.exp_softgmres_schwarz_Kcycle_directcoarse_overlap:
                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = lc.NoOfMultigridLevels - 1
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverCode.exp_softgmres_schwarz_directcoarse_overlap:
                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverCode.exp_multigrid:
                    throw new NotImplementedException("advertising space");

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
                        CoarseSolver = new SparseSolver() {
                            WhichSolver = SparseSolver._whichSolver.MUMPS
                        },
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks
                        },
                        Overlap = 1,
                    };
                    break;

                case LinearSolverCode.exp_gmres_levelpmg:
                    precond[0] = new LevelPmg() { UseHiOrderSmoothing = true, CoarseLowOrder = 1 };
                    break;

                case LinearSolverCode.exp_gmres_schwarz_pmg:
                    precond[0] = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks
                        },
                        CoarseSolver = null,
                        Overlap = 1,
                        EnableOverlapScaling = false,
                        UsePMGinBlocks = false,

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

                    precond[0] = My_MG_Precond(lc, LocalDOF, MultigridSeqLength, _prechain, _postchain, toppre, toppst);
                    break;

                case LinearSolverCode.exp_softpcg_schwarz_mg:
                    precond[0] = new Schwarz() {
                        FixedNoOfIterations = 1,
                        CoarseSolver = new SparseSolver() {
                            WhichSolver = SparseSolver._whichSolver.MUMPS,
                            LinConfig = lc
                        },
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = NoOfBlocks
                        },
                        Overlap = 1,
                    };
                    break;

                case LinearSolverCode.exp_softpcg_mg:
                    precond[0] = SpecialMultilevelSchwarz(lc, LocalDOF, MultigridSeqLength, MultigridOperatorConfig);
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
                                pLow=1,
                            },
                    };
                    break;
                case LinearSolverCode.selfmade:
                    Console.WriteLine("Selfmade Preconditioner is used!");
                    precond[0] = m_precond;
                    break;


                default:
                    throw new NotImplementedException("Preconditioner not available");
            }

            if (precond != null) {
                //special stdout stuff for Schwarz Preconds
                if (precond[0] is Schwarz)
                    Console.WriteLine("Additive Schwarz w. direct coarse, No of blocks: " + NoOfBlocks.MPISum());
                if (precond.Length > 1) {
                    foreach (var pre in precond) {
                        SetLinItCallback(pre, true);
                    }
                } else {
                    SetLinItCallback(precond[0], true);
                }
#if DEBUG
                if (precond.Length > 1) {
                    foreach (var pre in precond) {
                        Console.WriteLine("preconditioner: {0}", pre.GetType());
                    }
                } else {
                    Console.WriteLine("preconditioner: {0}", precond[0].GetType());
                }
#endif
                Check_precond();
            }

            return precond;
        }

        /// <summary>
        /// This one is the method-body of <see cref="GenerateLinear"/> and shall not be called from the outside. Some Solver aquire additional information, thus the timestepper is overgiven as well.
        /// </summary>
        private ISolverSmootherTemplate GenerateLinear_body(LinearSolverConfig lc, AggregationGridData[] MultigridSequence, MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig) {

            // +++++++++++++++++++++++++++++++++++++++++++++
            // the linear solvers:
            // +++++++++++++++++++++++++++++++++++++++++++++

            //if (lc == null)
            //    throw new ArgumentNullException();

            ISolverSmootherTemplate templinearSolve = null;



            //Calculate number of local DOF for every Multigridlevel
            int[] LocalDOF = GetLocalDOF(MultigridSequence, MultigridOperatorConfig);
            //int counter = 0;
            //for (int iLevel = 0; iLevel < DOFperCell.Length; iLevel++) {
            //    counter = iLevel;
            //    if (iLevel > MultigridOperatorConfig.Length - 1)
            //        counter = MultigridOperatorConfig.Length - 1;

            //    int d = MultigridSequence[iLevel].SpatialDimension;
            //    foreach (var variable in MultigridOperatorConfig[counter]) {
            //        int p = variable.Degree;
            //        switch (d) {
            //            case 1:
            //                DOFperCell[iLevel] += p + 1 + p + 1;
            //                break;
            //            case 2:
            //                DOFperCell[iLevel] += (p * p + 3 * p + 2) / 2;
            //                break;
            //            case 3:
            //                DOFperCell[iLevel] += (p * p * p + 6 * p * p + 11 * p + 6) / 6;
            //                break;
            //            default:
            //                throw new Exception("wtf?Spacialdim=1,2,3 expected");
            //        }
            //    }
            //    //CHANGE THIS ASAP: This is not right in case of XDG ... because blocks have truncated size at this code level
            //    LocalDOF[iLevel] = MultigridSequence[iLevel].CellPartitioning.LocalLength * DOFperCell[iLevel];
            //}

            //these values are acquired for some solvers
            int NoCellsLoc = MultigridSequence[0].CellPartitioning.LocalLength;
            int NoCellsGlob = MultigridSequence[0].CellPartitioning.TotalLength;
            int NoOfBlocks = (int)Math.Max(1, Math.Round(LocalDOF[0] / (double)lc.TargetBlockSize));
            int SpaceDim = MultigridSequence[0].SpatialDimension;
            int MultigridSeqLength = MultigridSequence.Length;
            ISolverSmootherTemplate[] precond = GeneratePrecond(lc, MultigridSequence, MultigridOperatorConfig);


            switch (lc.SolverCode) {
                case LinearSolverCode.automatic:
                    if (m_nc != null) {
                       templinearSolve = AutomaticSolver(lc, LocalDOF, SpaceDim, NoCellsLoc, NoCellsGlob,
                           precond);
                    }
                    break;

                case LinearSolverCode.classic_mumps:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.MUMPS,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverCode.classic_pardiso:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.PARDISO,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverCode.exp_schwarz_directcoarse_overlap:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 1,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverCode.exp_schwarz_directcoarse:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                            NoOfPartsPerProcess = 1,
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverCode.exp_schwarz_MG_directcoarse:

                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            Depth = lc.NoOfMultigridLevels - 1
                        },
                        Overlap = 0,
                        CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
                    };
                    break;

                case LinearSolverCode.exp_schwarz_MG_directcoarse_overlap:

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

                case LinearSolverCode.exp_softgmres:

                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = lc.MaxKrylovDim
                    };
                    break;

                case LinearSolverCode.exp_softgmres_schwarz_Kcycle_directcoarse_overlap:

                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = lc.MaxKrylovDim,
                        //m_Tolerance = lc.ConvergenceCriterion,
                        Precond = precond[0],
                    };
                    break;

                case LinearSolverCode.exp_softgmres_schwarz_directcoarse_overlap:
                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");

                    templinearSolve = new SoftGMRES() {
                        MaxKrylovDim = lc.MaxKrylovDim,
                        //m_Tolerance = lc.ConvergenceCriterion,
                        Precond = precond[0],
                    };
                    break;

                case LinearSolverCode.exp_multigrid:
                    if (lc.NoOfMultigridLevels < 2)
                        throw new ApplicationException("At least 2 Multigridlevels are required");
                    templinearSolve = new ILU() { };
                    break;

                case LinearSolverCode.exp_ILU:
                    templinearSolve = new ILU() { };
                    break;

                case LinearSolverCode.exp_Schur:
                    templinearSolve = new SchurPrecond() {
                        SchurOpt = SchurPrecond.SchurOptions.decoupledApprox
                    };
                    break;

                case LinearSolverCode.exp_Simple:
                    templinearSolve = new SchurPrecond() {
                        SchurOpt = SchurPrecond.SchurOptions.SIMPLE
                    };
                    break;

                case LinearSolverCode.exp_AS_1000:
                    if (MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF 
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 76,
                                NoOfPartsPerProcess = 213, // Warum 76
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    } else  //2D --> 75088DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 213,
                                NoOfPartsPerProcess = 213,
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    }
                    break;

                case LinearSolverCode.exp_AS_5000:
                    if (MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 43,
                                NoOfPartsPerProcess = 43,
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    } else  //2D --> 75088DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 16,
                                NoOfPartsPerProcess = 43,
                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    }

                    break;

                case LinearSolverCode.exp_AS_10000:
                    if (MultigridSequence[0].SpatialDimension == 3)   //3D --> 212940DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 22,
                                NoOfPartsPerProcess = 22,

                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    } else  //2D --> 75088DoF
                    {
                        templinearSolve = new Schwarz() {
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                //noofparts = 8,
                                NoOfPartsPerProcess = 22, //

                            },
                            CoarseSolver = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                                LinConfig = lc
                            },
                            Overlap = 1
                        };
                    }

                    break;

                case LinearSolverCode.exp_AS_MG:
                    templinearSolve = new Schwarz() {
                        m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                            //depth = asdepth,
                            Depth = 2,
                        },
                        CoarseSolver = new SparseSolver() {
                            WhichSolver = SparseSolver._whichSolver.MUMPS,    //PARDISO
                            LinConfig = lc
                        },

                        Overlap = 1
                    };
                    break;


                case LinearSolverCode.exp_localPrec:
                    templinearSolve = new LocalizedOperatorPrec() {
                        m_dt = lc.exp_localPrec_Min_dt,
                        m_muA = lc.exp_localPrec_muA,
                    };
                    break;

                case LinearSolverCode.classic_cg:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.CG,
                        LinConfig = lc
                    };
                    break;

                case LinearSolverCode.exp_softpcg_mg:
                    templinearSolve = new SoftPCG() {
                        //m_MaxIterations = _lc.MaxSolverIterations,
                        //m_Tolerance = _lc.ConvergenceCriterion,
                        TerminationCriterion = ((iter, r0_l2, r_l2) => r_l2 > lc.ConvergenceCriterion && iter < lc.MaxSolverIterations),
                        Precond = precond[0]
                    };
                    break;

                case LinearSolverCode.exp_softpcg_schwarz:

                    templinearSolve = new SoftPCG() {
                        //m_MaxIterations = lc.MaxSolverIterations,
                        //m_Tolerance = lc.ConvergenceCriterion,
                        Precond = precond[0]
                    };
                    break;

                case LinearSolverCode.exp_direct_lapack:
                    templinearSolve = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.Lapack
                    };
                    break;

                case LinearSolverCode.exp_softpcg_schwarz_directcoarse:

                    templinearSolve = new SoftPCG() {
                        //m_MaxIterations = lc.MaxSolverIterations,
                        //m_Tolerance = lc.ConvergenceCriterion,
                        Precond = precond[0],
                    };
                    break;

                case LinearSolverCode.exp_Kcycle_schwarz:
                    templinearSolve = KcycleMultiSchwarz(lc, LocalDOF);
                    break;

                case LinearSolverCode.exp_Kcycle_schwarz_4Rheology:
                    templinearSolve = KcycleMultiSchwarz_4Rheology(lc, LocalDOF);
                    break;

                case LinearSolverCode.exp_gmres_levelpmg:
                    templinearSolve = new SoftGMRES() {
                        //m_Tolerance = lc.ConvergenceCriterion,
                        //m_MaxIterations = lc.MaxSolverIterations,
                        MaxKrylovDim = lc.MaxKrylovDim,
                        Precond = precond[0]
                    };
                    break;

                case LinearSolverCode.exp_gmres_schwarz_pmg:
                    templinearSolve = new SoftGMRES() {
                        //m_Tolerance = lc.ConvergenceCriterion,
                        //m_MaxIterations = lc.MaxSolverIterations,
                        Precond = precond[0]
                    };

                    break;

                case LinearSolverCode.exp_softpcg_jacobi_mg:

                    templinearSolve = new SoftPCG() {
                        //m_MaxIterations = lc.MaxSolverIterations,
                        //m_Tolerance = lc.ConvergenceCriterion,
                        Precond = precond[0],
                    };
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

                    templinearSolve = MakeOrthoNormMGDecomp(lc,LocalDOF, MultigridSeqLength, new SolverSquence() {SolverChain= subsmoother}, new SolverSquence() { SolverChain = topsmoother });
                    break;

                case LinearSolverCode.exp_softpcg_schwarz_mg:
          
                    templinearSolve = new SoftPCG() {
                        //m_MaxIterations = lc.MaxSolverIterations,
                        //m_Tolerance = lc.ConvergenceCriterion,
                        Precond = precond[0],
                    };
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
                    Console.WriteLine("Selfmade LinearSolver is used!");
                    templinearSolve = m_linsolver;
                    break;

                default:
                    throw new NotImplementedException("Linear solver option not available");
            }
            Debug.Assert(templinearSolve != null);
            SetLinItCallback(templinearSolve,false);
#if DEBUG
            Console.WriteLine("linear Solver: {0}", lc.SolverCode.ToString());
#endif

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



            return templinearSolve;
        }

        private Action<int, double[], double[], MultigridOperator> CustomizedCallback;
        private Action<int, double[], double[], MultigridOperator> DefaultItCallback;

        private void SetLinItCallback(ISolverSmootherTemplate solverwithoutcallback, bool IsLinPrecond) {

            ISolverWithCallback _solverwithcallback = (ISolverWithCallback)solverwithoutcallback;

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

        private double m_ResOfPreviousSolver=0;

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
            DefaultItCallback = GenerateDefaultCallback<NonlinearSolver>(_type, nonlinsolver,_caseselect);
            if (m_nc.verbose) {
                nonlinsolver.IterationCallback += DefaultItCallback;
            }
            nonlinsolver.IterationCallback += CustomizedCallback;
        }

        private void FirstLineinCallBack() {
            if (m_Iterations == null) {
                string FirstLine = "NLinS-, Precond-, LinS-, total-Iterations : Type, SolverName, InfResi, Multigridlevel";
                m_Iterations = new int[4];
                Console.WriteLine(FirstLine);
            }
        }

        private int[] m_Iterations;

        public int[] GetIterationcounter {
            get {
                return m_Iterations;
            }
        }

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

        private LinearSolverConfig lc_check;
        private NonLinearSolverConfig nc_check;

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
        /// Automatic choice of linear solver depending on problem size, immersed boundary, polynomial degree, etc. In addition the nonlinearsolver config is considered as well.
        /// </summary>
        private ISolverSmootherTemplate AutomaticSolver( LinearSolverConfig lc, int[] LDOF, int Dim, int NoCellsLoc, int NoCellsGlob, ISolverSmootherTemplate[] PreCond) {

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
            var cellsGlo = NoCellsGlob;

            ISolverSmootherTemplate tempsolve = null;

            //var size = Timestepper.MultigridSequence[0].CellPartitioning.MpiSize;

            //// !!!!!!!!!!!UNTERSCHEIDUNG OB PICARD ODER NEWTON!!!!!!!!!!!!
            //if (nc.SolverCode == NonLinearSolverCode.NewtonGMRES) {

            //    // Spatial Dimension
            //    switch (D) {
            //        case 1:
            //            break;
            //            throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
            //        //break;

            //        case 2:
            //            throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
            //        //break;

            //        case 3:
            //            //var dofsPerCell3D = (3 * (pV * pV * pV + 6 * pV * pV + 11 * pV + 6) / 6 + 1 * (pP * pP * pP + 6 * pP * pP + 11 * pP + 6) / 6);
            //            int dofsPerCell3D = LDOF[0] / NoCellsLoc;
            //            var dofsLoc = dofsPerCell3D * cellsLoc;
            //            int dofsGlo = dofsPerCell3D * cellsGlo;

            //            var PPP = (int)Math.Ceiling(dofsLoc / 6500.0);

            //            Console.WriteLine("Analysing the problem yields " + PPP + " parts per process.");

            //            if (dofsGlo > 10000) {

            //                if (lc.NoOfMultigridLevels < 2)
            //                    throw new ApplicationException("At least 2 Multigridlevels are required");

            //                tempsolve = new Schwarz() {
            //                    m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
            //                        NoOfPartsPerProcess = PPP,
            //                    },
            //                    Overlap = 1,
            //                    CoarseSolver = DetermineMGSquence(lc.NoOfMultigridLevels - 2, lc)
            //                };
            //            } else {
            //                tempsolve = new SparseSolver() {
            //                    WhichSolver = SparseSolver._whichSolver.MUMPS,
            //                    LinConfig = lc
            //                };
            //            }
            //            break;

            //        default:
            //            throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
            //    }
            //} else {
                // Spatial Dimension
                switch (D) {
                    case 1:
                        break;
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                    //break;

                    case 2:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                    //break;

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
                            tempsolve = new SparseSolver() {
                                WhichSolver = SparseSolver._whichSolver.MUMPS,
                                LinConfig = lc
                            };
                        }
                        break;

                    default:
                        throw new NotImplementedException("Currently not implemented for " + D + " Dimensions");
                //}

            }

            return tempsolve;
            //Timestepper.

            // Wenn Gesamtproblem in 2D < 100000 DoFs -> Direct Solver
            // Wenn Gesamtproblem in 3D < 10000 DoFs -> Direct Solver 

            // Block Solve 3D ca. 6000 DoFs per Process -> Adjust Blocks per Process
            // Coarse Solve ca. 5000 bis 10000 DoFs. -> Adjust Multigrid Levels

        }

        private ISolverSmootherTemplate AutomaticPrecond(LinearSolverConfig lc) {
            ISolverSmootherTemplate precond=null;
            return precond;
        }

            /// <summary>
            /// Determines a solver sequence depending on MGlevels
            /// </summary>
            private ISolverSmootherTemplate DetermineMGSquence(int MGlevels, LinearSolverConfig lc) {
            ISolverSmootherTemplate solver;
            if (MGlevels > 0) {
                solver = new ClassicMultigrid() { CoarserLevelSolver = DetermineMGSquence(MGlevels - 1, lc) };
            } else {
                solver = new SparseSolver() {
                    WhichSolver = SparseSolver._whichSolver.MUMPS,
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

        /// <summary>
        /// clears overgiven selfmade solvers
        /// </summary>
        public void Clear() {
            this.m_linsolver = null;
            this.m_nonlinsolver = null;
            this.m_precond = null;
        }

        private void Check_NonLinearSolver(NonlinearSolver NLSolver)
        {
            if (NLSolver is Newton newtonsolver)
                Check_Newton(newtonsolver);
            if (NLSolver is FixpointIterator picardsolver)
                Check_Picard(picardsolver);
        }

        private void Check_Newton(Newton NewtonSolver)
        {
            bool check=true;
            check = m_nc.ConvergenceCriterion == NewtonSolver.ConvCrit &&
            m_nc.MaxSolverIterations == NewtonSolver.MaxIter &&
            m_nc.MinSolverIterations == NewtonSolver.MinIter;
            //m_nc.PrecondSolver.Equals(NewtonSolver.Precond);
            Debug.Assert(check);
        }

        private void Check_Picard(FixpointIterator FPSolver)
        {
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
        private void Check_linsolver() {
            if (!m_lc.Equals(lc_check))
                Console.WriteLine("WARNING: LinearSolverConfig was changed since constructor!");
        }

        /// <summary>
        /// Checks overgiven selfmade proconditioner solver
        /// </summary>
        /// <returns></returns>
        private void Check_precond() {
            if (!m_lc.Equals(lc_check))
                Console.WriteLine("WARNING: LinearSolverConfig was changed since constructor!");
        }

        /// <summary>
        /// Checks overgiven selfmade nonlinear solver
        /// </summary>
        /// <returns></returns>
        private void Check_nonlinsolver() {
            if(!m_nc.Equals(nc_check))
                Console.WriteLine("WARNING: NonLinearSolverConfig was changed since constructor!");
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
                    case LinearSolverCode.exp_softgmres_schwarz_directcoarse:
                    case LinearSolverCode.exp_softgmres_schwarz_directcoarse_overlap:
                    case LinearSolverCode.exp_softgmres_schwarz_Kcycle_directcoarse_overlap:
                    case LinearSolverCode.exp_softgmres:
                        LinIsGmres = true;
                        break;
                    default:
                        break;
                }

                return m_nc.SolverCode == NonLinearSolverCode.Newton && LinIsGmres;
                    }
        }

        private int m_MG_Counter = 0;
        private double[] ProlRes = new double[10];
        private double[] RestRes = new double[10];
        //private double[] ProlRes = {1,2,3,4,5};
        //private double[] RestRes = {10,9,8,7,6};

        private int m_MG_Depth=-1;
        private int MG_Depth{
            get{ return m_MG_Depth; }
            set{ m_MG_Depth = Math.Max(value, m_MG_Depth); }
            }

        private void MultigridCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop)
        {
            int currentMGLevel = Mgop.LevelIndex;

            if (m_MG_Counter - currentMGLevel == +1)
            {
                double residualNormAf = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm after coarse grid correction
                ProlRes[currentMGLevel] = residualNormAf;
                Console.WriteLine("after Prolongation {0}<-{1}: {2}", currentMGLevel, currentMGLevel+1 , ProlRes[currentMGLevel] - ProlRes[currentMGLevel+1]);
            }

            if(m_MG_Counter - currentMGLevel ==0 && currentMGLevel== MG_Depth)
            {
                double residualNormAf = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm after coarse grid correction
                ProlRes[currentMGLevel] = residualNormAf;
                Console.WriteLine("after Prolongation {0}<-{1}: {2}", currentMGLevel, currentMGLevel+1, ProlRes[currentMGLevel]- RestRes[currentMGLevel]);
            }

            if (m_MG_Counter - currentMGLevel == 0 && currentMGLevel == 0)
            {
                double residualNormB4 = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm before coarse grid correction
                RestRes[currentMGLevel] = residualNormB4;
                Console.WriteLine("before Restriction {0}->{1}: {2}", currentMGLevel, currentMGLevel + 1, RestRes[currentMGLevel] - ProlRes[currentMGLevel]);
            }

            if (m_MG_Counter - currentMGLevel == -1)
            {
                double residualNormB4 = currentRes.L2Norm().MPISum(Mgop.OperatorMatrix.MPI_Comm); // residual norm before coarse grid correction
                RestRes[currentMGLevel] = residualNormB4;
                Console.WriteLine("before Restriction {0}->{1}: {2}", currentMGLevel, currentMGLevel + 1, RestRes[currentMGLevel] - RestRes[currentMGLevel-1]);
            }
            m_MG_Counter = currentMGLevel;
            MG_Depth = currentMGLevel;
        }

        private ISolverSmootherTemplate My_MG_Precond(LinearSolverConfig _lc, int[] _LocalDOF, int MSLength, ISolverSmootherTemplate[] prechain, ISolverSmootherTemplate[] postchain, ISolverSmootherTemplate toplevelpre, ISolverSmootherTemplate toplevelpst) {

            bool isLinPrecond = true;

            int DirectKickIn = _lc.TargetBlockSize;

            ISolverSmootherTemplate[] MultigridChain = new ISolverSmootherTemplate[MSLength];
            for (int iLevel = 0; iLevel < MSLength; iLevel++) {
                int SysSize = _LocalDOF[iLevel].MPISum();
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == MSLength - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    MultigridChain[iLevel] = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.MUMPS,
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

            //MultigridOperator Current = op;
            ISolverSmootherTemplate[] MultigridChain = new ISolverSmootherTemplate[MSLength];
            for (int iLevel = 0; iLevel < MSLength; iLevel++) {
                int SysSize = _LocalDOF[iLevel].MPISum();
                //int SysSize = Current.Mapping.TotalLength;
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == MSLength - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    MultigridChain[iLevel] = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.MUMPS,
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

            List<ISolverSmootherTemplate> MG_list = new List<ISolverSmootherTemplate>();

            int DirectKickIn = _lc.TargetBlockSize;
            int MaxMGDepth = GetMGdepth(DirectKickIn, MSLength, _LocalDOF);
            SetLinItCallback(subsmootherchain, true);

            foreach(var solversmoother in ((SolverSquence)subsmootherchain).SolverChain)
                SetLinItCallback(solversmoother, true);

            foreach (var toplevelsmoother in ((SolverSquence)toplevelsmootherchain).SolverChain)
                SetLinItCallback(toplevelsmoother, true);

            for (
                
                int iDepth = MaxMGDepth; iDepth >= 0; iDepth--)
            {
                ISolverSmootherTemplate solvertoinject;

                if (iDepth == 0) {
                    solvertoinject = toplevelsmootherchain.CloneAs();
                }
                else if (iDepth == MaxMGDepth) {
                    solvertoinject = new SparseSolver()
                    {
                        WhichSolver = SparseSolver._whichSolver.MUMPS,
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

        private int GetMGdepth(int DirectKickIn, int MSLength, int[] LocalDOF) {
            int MGdepth = -1;
            for (int iLevel = 0; iLevel < MSLength; iLevel++) {
                int SysSize = LocalDOF[iLevel].MPISum();
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                useDirect |= iLevel == MSLength - 1;
                useDirect |= NoOfBlocks.MPISum() <= 1;

                if (useDirect) {
                    MGdepth = iLevel;
                }
            }
            return MGdepth;
        }

        /// <summary>
        /// 
        /// </summary>
        ISolverSmootherTemplate KcycleMultiSchwarz(LinearSolverConfig _lc, int[] _LocalDOF) {

            // my tests show that the ideal block size may be around 10'000
            int DirectKickIn = _lc.TargetBlockSize;
            
            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();
            

            for (int iLevel = 0; iLevel < _lc.NoOfMultigridLevels; iLevel++) {
                int SysSize = _LocalDOF[iLevel].MPISum();
                int NoOfBlocks = (int)Math.Ceiling(((double)SysSize) / ((double)DirectKickIn));

                
                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                //useDirect |= iLevel == _lc.NoOfMultigridLevels - 1; // 
                useDirect |= NoOfBlocks.MPISum() <= 1;

                NoOfBlocks = Math.Max(1, NoOfBlocks); // ensure at least 1 block per MPI core, if direct solver isn't used
                
                if (useDirect)
                    Console.WriteLine("   KcycleMultiSchwarz: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("   KcycleMultiSchwarz: lv {0}, no of blocks {1} : ", iLevel, NoOfBlocks);

                ISolverSmootherTemplate levelSolver;
                if (useDirect) {
                    levelSolver = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.PARDISO,
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
                        Overlap = 1, // overlap seems to help; more overlap seems to help more
                        EnableOverlapScaling = true,
                        UsePMGinBlocks = false
                    };
                    

                    levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother1,
                    };

                    if(iLevel > 0) {
                        ((OrthonormalizationMultigrid)levelSolver).TerminationCriterion = (i, r0, r) => i <= 1;
                    }


                    ((OrthonormalizationMultigrid)levelSolver).IterationCallback =
                        delegate (int iter, double[] X, double[] Res, MultigridOperator op) {
                            double renorm = Res.MPI_L2Norm();
                            Console.WriteLine("      OrthoMg " +iter  + " : " + renorm);
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

            //SolverChain.Add(new DynamicMultigrid());

            /*
            var S = new SoftGMRES() {
                m_Tolerance = _lc.ConvergenceCriterion,
                Precond = new LevelPmg()
            };

            SolverChain.Add(S);
            */


            return SolverChain[0];
        }


        ISolverSmootherTemplate KcycleMultiSchwarz_4Rheology(LinearSolverConfig _lc, int[] _LocalDOF) {

            // my tests show that the ideal block size may be around 10'000
            int DirectKickIn = _lc.TargetBlockSize;

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();


            Console.WriteLine("experimental MG configuration for rheology");

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPIsize);

            for (int iLevel = 0; iLevel < _lc.NoOfMultigridLevels; iLevel++) {
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
                    useDirect = true;
                }

                if (useDirect)
                    Console.WriteLine("   KcycleMultiSchwarz: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("   KcycleMultiSchwarz: lv {0}, no of blocks {1} : ", iLevel, NoOfBlocks);

                ISolverSmootherTemplate levelSolver;
                if (useDirect) {
                    levelSolver = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.PARDISO,
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
                        UseHiOrderSmoothing = false,
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

    }
}