using ilPSP;
using ilPSP.LinSolvers.HYPRE;
using ilPSP.Tracing;
using log4net.Core;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// 
    /// </summary>
    public enum SchwarzImplementation {

        /// <summary>
        /// default option
        /// </summary>
        Auto = 0,

        /// <summary>
        /// force to use the <see cref="SchwarzForCoarseMesh"/> implementation
        /// </summary>
        CoarseMesh = 2,


        /// <summary>
        /// force to use the <see cref="Schwarz"/> implementation
        /// </summary>
        PerProcess = 2
    }

    /// <summary>
    /// Dynamic configuration for the orthonormalization multigrid.
    /// As smoothers, additive Schwarz (<see cref="Schwarz"/>) is used on all mesh levels, except the coarsest one.
    /// There, a direct solver is used.
    /// The coarsest level is determined dynamically, depending on the problem and the mesh size using <see cref="TargetBlockSize"/>.
    /// </summary>
    [Serializable]
    public class OrthoMGSchwarzConfig : IterativeSolverConfig {

        /// <summary>
        /// selection of Schwarz solver;
        /// should not be changed by the user; the main reason for this member is, to force the use of both code (<see cref="Schwarz"/>, <see cref="SchwarzForCoarseMesh"/>) paths in tests.
        /// </summary>
        [DataMember]
        public SchwarzImplementation SchwarzImplementation = SchwarzImplementation.Auto;

        
        /// <summary>
        /// If any blocking is used (Schwarz, block Jacobi), a target for the block size.
        /// Tests show that the ideal block size may be around 10000, but this may depend on computer, DG polynomial order, etc.
        /// 
        /// This also determines the actual number of multigrid levels used in dependence of the problem size;
        /// As soon as the number of DOF's on a certain multigrid level fall below this threshold, a direct 
        /// solver is used and no further multigrid levels are allocated.
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.ExclusiveLowerBound(99.0)]
        public int TargetBlockSize = 10000;


        /// <summary>
        /// Determines maximal DG order within coarse system of a p-Multigrid;
        /// a negative number triggers an automatic, dynamic selection;
        /// </summary>
        [DataMember]
        public int pMaxOfCoarseSolver = -11;

        int Get_pMaxOfCoarseSolver(int maxDeg) {
            if (pMaxOfCoarseSolver >= 0)
                return pMaxOfCoarseSolver;
            else
                return (int)Math.Ceiling(maxDeg*0.5);
        }


        /// <summary>
        /// Sets the **maximum** number of Multigrid levels to be used.
        /// The numbers of levels which are actually used is probably much less, and determined via <see cref="TargetBlockSize"/>.
        /// Multigrid approach is used to get a preconditioner for Krylov solvers, e.g. GMRES.
        /// </summary>
        [DataMember]
        public int NoOfMultigridLevels = 1000000;

        /// <summary>
        /// Use a p-two-grid approach in Additive Schwarz (<see cref="Schwarz.Config.UsePMGinBlocks"/>).
        /// 
        /// This results in a lesser number of Schwarz blocks and more computational cells per Schwarz block.
        /// It seems to accelerate convergence for certain problems (mostly single-phase).
        /// The order of the low-order system can be adjusted using <see cref="pMaxOfCoarseSolver"/>.
        /// </summary>
        [DataMember]
        public bool UsepTG = false;

        /// <summary>
        /// System size at which the multigrid recursion is terminated;
        /// This should be significantly above the break-even point from which on iterative solvers are faster than direct ones, 
        /// since the coarse solver is called very often and the direct solver therefore can use its advantage of saving the factorization.
        /// </summary>
        [DataMember]
        public int CoarseKickIn = 90000;

        /// <summary>
        /// - if true, a p-two-grid method (<see cref="PTGconfig"/>) is used for the coarsest level
        /// - if false, a direct solver (<see cref="DirectSolver"/>) is used at the coarsest level
        /// </summary>
        [DataMember]
        public bool CoarseUsepTG = false;


        /// <summary>
        /// 
        /// </summary>
        public override string Name => "OrthonormalizationMultigrid with Additive Schwarz Smoother";

        /// <summary>
        /// 
        /// </summary>
        public override string Shortname => "OrthoMG w Add Swz";

        /// <summary>
        /// 
        /// </summary>
        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            Func<int, int> SblkSizeFunc = delegate (int iLevel) {
                return TargetBlockSize; 
            };
            var instance = KcycleMultiSchwarz(level, SblkSizeFunc);
            instance.Init(level);
            return instance;
        }
        
        double ComputeInbalance(MultigridOperator level) {
            long LocalLength = level.Mapping.LocalLength;
            long MaxLocalLength = LocalLength.MPIMax();
            long LocalInbalance = MaxLocalLength - LocalLength; ;
            if (LocalInbalance < 0)
                throw new ArithmeticException("Internal error in MPI maximum");
            long TotalInbalance = LocalInbalance.MPISum();

            return ((double)TotalInbalance) / ((double)level.Mapping.TotalLength);
        }

        int[] GetLoOrderDegrees(MultigridOperator level) {
            int[] Degrees = level.Degrees;
            int DegDecrease = Math.Max(Degrees.Max() - this.Get_pMaxOfCoarseSolver(level.Degrees.Max()), 0);
            int[] loDegrees = Degrees.Select(deg => Math.Max(deg - DegDecrease, 0)).ToArray();
            return loDegrees;
        }

        /// <summary>
        /// Number of DG basis polynomials for different spatial dimensions
        /// </summary>
        static int Np(int D, int p) {
            switch (D) {
                case 1: return p + 1;
                case 2: return ((p + 1)*(p + 1) + p + 1)/2;
                case 3: return (p*p*p + 6*p*p + 11*p + 6)/6;
                default: throw new ArgumentOutOfRangeException("unsupported spatial dimension: " + D);
            }
        }

        int GetLowOrderLength(MultigridOperator level) {
            int D = level.GridData.SpatialDimension;

            // dof for some un-cut cell:
            int uncutCellDof = level.Mapping.DgDegree.Sum(deg => Np(D, deg));
            int uncutLowpCellDof = GetLoOrderDegrees(level).Sum(deg => Np(D, deg));

            if (level.Mapping.LocalLength % uncutCellDof != 0)
                throw new ArgumentException($"Internal error: LocalLength = {level.Mapping.LocalLength}, uncut cell dof = {uncutCellDof}");

            return uncutLowpCellDof * level.Mapping.LocalLength / uncutCellDof;
        }


        const double INBALANCE_THRESHOLD = 0.3;

        const int PROCESSLOCAL_SCHWARZBLOCK_MINIMUM = 2;


        /// <summary>
        /// - if true, the <see cref="SchwarzForCoarseMesh"/> smoother should be used
        /// - if false, the <see cref="Schwarz"/> smoother should be used
        /// </summary>
        (ISolverSmootherTemplate swz, int LevelLocalBlocks, int LevelGlobalBlocks) SelectSchwarzSmoother(MultigridOperator level, int FinerLevelLocalBlocks, int FinerLevelGlobalBlocks) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                int MPIsize = level.Mapping.MpiSize;
                int GlobalNoOfBlocks = Math.Max(2, (int)Math.Round((double)level.Mapping.TotalLength / (double)this.TargetBlockSize)); // for the coarse-mesh-Schwarz
                tr.Info($"GlobalNoOfBlocks = {GlobalNoOfBlocks}");

                int LocalNoOfBlocks; // for the process-local blocking mesh
                {
                    if (this.UsepTG) {
                        var lowDOF = this.GetLowOrderLength(level);
                        tr.Info($"low-order (pMaxOfCoarseSolver = {this.Get_pMaxOfCoarseSolver(level.Degrees.Max())}) DOF: {lowDOF}");
                        LocalNoOfBlocks = lowDOF/TargetBlockSize;
                    } else {
                        LocalNoOfBlocks = level.Mapping.LocalLength/TargetBlockSize;
                    }

                    LocalNoOfBlocks = Math.Max(LocalNoOfBlocks, 1);
                    LocalNoOfBlocks = LocalNoOfBlocks.MPIMax();
                }

                double inbal = ComputeInbalance(level);
                tr.Info("DOF MPI inbalance is " + inbal);

                if (inbal <= INBALANCE_THRESHOLD
                    && LocalNoOfBlocks >= (this.UsepTG ? 1 : PROCESSLOCAL_SCHWARZBLOCK_MINIMUM)
                    && FinerLevelLocalBlocks >= 0
                    && (SchwarzImplementation == SchwarzImplementation.Auto || SchwarzImplementation == SchwarzImplementation.PerProcess)) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // * DOF's seem sufficiently distributed across MPI cores;
                    // * sufficient number of blocks per process
                    // * upper/finer level also uses process-local Schwarz blocks
                    // sufficiently balanced to use Schwarz with process-local blocking strategy
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    bool Reduce = (FinerLevelLocalBlocks > LocalNoOfBlocks).MPIOr();

                    if (Reduce) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Number of blocks on this level is sufficiently lower then on upper/finer level
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        GlobalNoOfBlocks = LocalNoOfBlocks.MPISum();
                        tr.Info($"per-process blocking Schwarz, number of blocks blocking is " + (LocalNoOfBlocks.MPIAllGather().ToConcatString("[", "-", "]")) + ", sum = " + GlobalNoOfBlocks);
                        var r = new Schwarz() {
                            FixedNoOfIterations = 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsOnCurrentProcess = LocalNoOfBlocks
                            },
                            //ActivateCachingOfBlockMatrix = SmootherCaching,

                        };
                        r.config.EnableOverlapScaling = true;
                        r.config.Overlap = 1; // overlap seems to help; more overlap seems to help more
                        r.config.UsePMGinBlocks = this.UsepTG;
                        //*/

                        return (r, LocalNoOfBlocks, GlobalNoOfBlocks);
                    } else {
                        tr.Info("Failing to reduce **local** number of blocks " + (LocalNoOfBlocks.MPIAllGather().ToConcatString("[", "-", "]")) +  " wrt. previous level " + (LocalNoOfBlocks.MPIAllGather().ToConcatString("[", "-", "]")) + ".");
                    }
                }

                if(SchwarzImplementation == SchwarzImplementation.Auto || SchwarzImplementation == SchwarzImplementation.CoarseMesh) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // For some reason, we failed using the per-process-blocking Schwarz;
                    // so, try with the coarse-mesh-Schwarz
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    if (GlobalNoOfBlocks <= 1) {
                        // should use somthing else
                        return (null, -1, -11);
                    } else if (GlobalNoOfBlocks <= MPIsize) {
                        // ok with that
                        // no up-rounding; some processors are allowed to sleep
                    } else {
                        GlobalNoOfBlocks = ((int)Math.Round((double)GlobalNoOfBlocks / MPIsize))*MPIsize;
                    }



                    if (GlobalNoOfBlocks < FinerLevelGlobalBlocks) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Number of blocks on this level is sufficiently lower then on upper/finer level
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        
                        tr.Info($"using coarse-mesh-Schwarz, GlobalNoOfBlocks = {GlobalNoOfBlocks}");

                        var r = new SchwarzForCoarseMesh();
                        r.config.NoOfBlocks = GlobalNoOfBlocks;
                        r.config.EnableOverlapScaling = true;
                        return (r, -1, GlobalNoOfBlocks);
                    } else  {
                        tr.Info("Failing to reduce **global** number of blocks (" + GlobalNoOfBlocks + ") wrt. previous level (" + FinerLevelGlobalBlocks + ").");
                    }
                }

                tr.Info($"unable to create Schwarz smoother at level {level.LevelIndex}");
                return (null, -1, -1);
            }
        }

        ISolverSmootherTemplate GetAlternateSmoother() {
            var altSmooth3 = new SoftGMRES();
            altSmooth3.MaxKrylovDim = 80;
            altSmooth3.TerminationCriterion = (iIter, r0, ri) => {
                return (iIter <= 82 && ri > 1e-7*r0, ri <= 1e-7*r0);
            };
            altSmooth3.Precond = new BlockJacobi();

            return altSmooth3;
        }




        ISolverSmootherTemplate KcycleMultiSchwarz(MultigridOperator op, Func<int, int> SchwarzblockSize) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                int MSLength = op.NoOfLevels;

                int mpiSz = op.Mapping.MpiSize;
                int maxDG = op.DGpolynomialDegreeHierarchy.First().Max();
                var SolverChain = new List<ISolverSmootherTemplate>();

                int FinerLevelGlobalBlocks = int.MaxValue;
                int FinerLevelLocalBlocks = FinerLevelGlobalBlocks / mpiSz;


                for (MultigridOperator op_lv = op; op_lv != null; op_lv = op_lv.CoarserLevel) {
                    int iLevel = op_lv.LevelIndex;

                    
                    tr.Info($"KcycleMultiSchwarz: lv {iLevel}, DOFs: " + op_lv.Mapping.MpiSize.ForLoop(rnk => op_lv.Mapping.GetLocalLength(rnk)).ToConcatString("[", ",", "]") + ", sum = " + op_lv.Mapping.TotalLength);
                    tr.Info($"KcycleMultiSchwarz: lv {iLevel}, Cells " + op_lv.Mapping.MpiSize.ForLoop(rnk => op_lv.Mapping.GetLocalNoOfBlocks(rnk)).ToConcatString("[", ",", "]") + ", sum = " + op_lv.Mapping.TotalNoOfBlocks);

                    bool TerminateMultigrid = false;
                    if(iLevel >= this.NoOfMultigridLevels) {
                        tr.Info($"KcycleMultiSchwarz: maximum number of multigrid levels ({this.NoOfMultigridLevels}) reached.");
                        TerminateMultigrid = true;
                    }

                    if(this.CoarseUsepTG && op_lv.Mapping.TotalLength > this.CoarseKickIn) {
                        long CoarseLength = ((long)this.GetLowOrderLength(op_lv)).MPISum();

                        if (CoarseLength <= this.CoarseKickIn) {
                            tr.Info($"KcycleMultiSchwarz: lv {iLevel}, using p-two-grid (pTG), low-degree is {this.Get_pMaxOfCoarseSolver(op_lv.Degrees.Max())}, low-order DOF {CoarseLength} <= kick-in-value ({this.CoarseKickIn}). ");
                            TerminateMultigrid = true;
                        }


                    } else {
                        if(op_lv.Mapping.TotalLength <= this.CoarseKickIn) {
                            tr.Info($"KcycleMultiSchwarz: lv {iLevel}, using direct solver, DOF {op_lv.Mapping.TotalLength} <= kick-in-value ({this.CoarseKickIn}). ");
                            TerminateMultigrid = true;
                        }
                    }

                    // instantiate solver at level
                    // ===========================

                    ISolverSmootherTemplate levelSolver;
                    if (TerminateMultigrid) {

                        // ++++++++++++++++++++++++++++++++++++++
                        // terminate multigrid recursion and ...
                        // ++++++++++++++++++++++++++++++++++++++

                        if (this.CoarseUsepTG && op_lv.Mapping.TotalLength > this.CoarseKickIn) {
                            // + + + + + + + + + + + + + + + + + + +
                            // ... use p-two-grid as coarse solver
                            // + + + + + + + + + + + + + + + + + + +


                            var coarseConfig = new PTGconfig() {
                                pMaxOfCoarseSolver = this.Get_pMaxOfCoarseSolver(op_lv.Degrees.Max()),
                            };

                            var _levelSolver2 = coarseConfig.CreateInstance(op_lv);
                            (_levelSolver2 as SoftGMRES).TerminationCriterion = delegate (int i, double r0, double r) {
                                var ret = (i <= 1 || r > r0 * 0.01, true);
                                return ret;
                            };
                            levelSolver = _levelSolver2;

                        } else {
                            // + + + + + + + + + + + + + + + + +
                            // ... use direct as coarse solver
                            // + + + + + + + + + + + + + + + + +

                            var _levelSolver3 = new DirectSolver();
                            _levelSolver3.config.WhichSolver = DirectSolver._whichSolver.PARDISO;
                            _levelSolver3.config.UseDoublePrecision = (iLevel == 0); // only use double precision id Direct solver is top solver
                            _levelSolver3.config.TestSolution = false;
                            _levelSolver3.ActivateCaching = (__1, __2) => true;
                            levelSolver = _levelSolver3;
                        }

                    } else {

                        // +++++++++++++++++++++++++++++++++++
                        // do some further multigrid recursion
                        // +++++++++++++++++++++++++++++++++++



                        ISolverSmootherTemplate smoother1;


                        int locBlk, glbBlk;
                        (smoother1, locBlk, glbBlk) = this.SelectSchwarzSmoother(op_lv, FinerLevelLocalBlocks, FinerLevelGlobalBlocks);
                        FinerLevelGlobalBlocks = glbBlk;
                        FinerLevelLocalBlocks = locBlk;
                        /*
                        var smoother1 = new SchwarzForCoarseMesh();
                        smoother1.config.NoOfBlocks = NoBlocks[iLevel].MPISum();
                        smoother1.config.EnableOverlapScaling = true;
                        //*/

                        SoftGMRES altSmooth3;
                        {   altSmooth3 = new SoftGMRES();
                            altSmooth3.MaxKrylovDim = 80;
                            altSmooth3.TerminationCriterion = (iIter, r0, ri) => {
                                return (iIter <= 82 && ri > 1e-7*r0, ri <= 1e-7*r0);
                            };
                            altSmooth3.Precond = new BlockJacobi();
                        }

                        if(smoother1 == null) {
                            // failed to init Schwarz smoother 
                            tr.Info($"KcycleMultiSchwarz: lv {iLevel}, no Schwarz smoother available; using GMRES+BlockJacobi. ");
                            smoother1 = altSmooth3;
                            altSmooth3 = null;
                        }


                        var _levelSolver4 = new OrthonormalizationMultigrid() {
                            PreSmoother = smoother1,
                            PostSmoother = smoother1
                        };
                        _levelSolver4.config.m_omega = 1; // v-cycle
                        //_levelSolver4.config.m_omega = 2; // w-cycle


                        if(altSmooth3 != null) {
                            _levelSolver4.AdditionalPostSmoothers = new[] { altSmooth3 }; 
                        }


                        if (iLevel == 0) {
                            _levelSolver4.config.NoOfPostSmootherSweeps = 20;
                        } else {
                            if (glbBlk > 0)
                                _levelSolver4.config.NoOfPostSmootherSweeps = (int)Math.Max(2, Math.Max(maxDG, Math.Round(Math.Log10(glbBlk) * 3.0) + maxDG - 2));
                            else
                                _levelSolver4.config.NoOfPostSmootherSweeps = 2;
                        }
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, NoOfPostSmootherSweeps = {_levelSolver4.config.NoOfPostSmootherSweeps}");

                        levelSolver = _levelSolver4;
                    }

                    // Set termination criterion for respective level
                    // ==============================================

                    if (levelSolver is IProgrammableTermination _levelSolver) { 
                        if (iLevel == 0) {
                            (_levelSolver).TerminationCriterion = this.DefaultTermination;
                            //if(maxDG > 3)
                            //    _levelSolver.config.NoOfPostSmootherSweeps = 6;
                        } else if (iLevel == 1) {

                          
                            (_levelSolver).TerminationCriterion = delegate (int i, double r0, double r) {
                                var ret = (i <= 1 || r > r0 * 0.01, true);
                                return ret;
                            };

                            /*
                            

                            //smoother1.config.Overlap = 0;

                            _levelSolver.AdditionalPostSmoothers = new ISolverSmootherTemplate[] {
                              altSmooth3
                                };
                            */

                        } else {
                            (_levelSolver).TerminationCriterion = (i, r0, r) => (i <= 1, true);
                        }            
                    }
                    SolverChain.Add(levelSolver);


                    // Plug multigrid together: 
                    // ========================
                    if (iLevel > 0) {
                        if (SolverChain[iLevel - 1] is OrthonormalizationMultigrid omg_fine)
                            omg_fine.CoarserLevelSolver = levelSolver;
                        else if (SolverChain[iLevel - 1] is GenericRestriction rest)
                            rest.CoarserLevelSolver = levelSolver;
                    }

                    if (TerminateMultigrid) {
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, terminating multigrid recursion.");
                        break;
                    }

                    if(op_lv.CoarserLevel == null) {
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, terminating multigrid recursion, no further level available.");

                    }
                }

                tr.Info($"KcycleMultiSchwarz: depth of multigrid configuration: {SolverChain.Count}.");
                return SolverChain[0];
            }
        }

                 /*
        ISolverSmootherTemplate KcycleMultiSchwarz(MultigridOperator op, Func<int, int> SchwarzblockSize) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                int MSLength = op.NoOfLevels;

                int mpiSz = op.Mapping.MpiSize;

                var SolverChain = new List<ISolverSmootherTemplate>();

                int maxDG = op.DGpolynomialDegreeHierarchy.First().Max();
                int minDG = Math.Max(op.DGpolynomialDegreeHierarchy.Last().Max(), pMaxOfCoarseSolver);
                minDG = Math.Min(minDG, maxDG); // make sure that minDG <= MaxDG;

                int pCoarsest = UsepTG ? minDG : -1;
                int[] NoBlocks = NoOfSchwarzBlocks(op, pCoarsest, SchwarzblockSize);

                Console.WriteLine("REM: block size reduced");
                for (int i = 0; i < NoBlocks.Length; i++)
                    NoBlocks[i] = Math.Max(1, NoBlocks[i]/3);

                int[] GlobalNoBlocks = NoBlocks.MPISum(op.Mapping.MPI_Comm);

                if (NoBlocks.Any(no => no <= 0))
                    throw new ApplicationException("Error in algorithm: No Of blocks cannot be 0.");
                var OneBlockPerNode = NoBlocks.Select(nb => nb <= 1).MPIAnd(op.Mapping.MPI_Comm);


                for (MultigridOperator op_lv = op; op_lv != null; op_lv = op_lv.CoarserLevel) {
                    int iLevel = op_lv.LevelIndex;

                    bool useDirect = false;
                    useDirect |= op_lv.Mapping.TotalLength <= TargetBlockSize;
                    useDirect |= GlobalNoBlocks[iLevel] <= 1;
                    //useDirect |= (double)PrevSize / (double)SysSize < 1.5 && SysSize < 50000; // degenerated MG-Agglomeration, because too few candidates
                    useDirect |= op_lv.CoarserLevel == null;

                    bool terminateMGRecursion = false;
                    if (useDirect == false && iLevel > 0) {
                        if (OneBlockPerNode[iLevel] && OneBlockPerNode[iLevel - 1]) {
                            // $"Only one block per node on level {iLevel}, and should be skipped, but cannot be skipped be skipped, but not supported yet."
                            //throw new NotSupportedException($"Only one block per node on level {iLevel}, and should be skipped, but cannot be skipped be skipped, but not supported yet.");
                            terminateMGRecursion = true;
                        }
                    }
                    Debug.Assert(useDirect.MPIEquals(op.Mapping.MPI_Comm));
                    Debug.Assert(terminateMGRecursion.MPIEquals(op.Mapping.MPI_Comm));

                    Console.WriteLine("Rem: testcode aktiv !!!!!!!! (hardcoded direct solver @ lv 2) ");
                    terminateMGRecursion = op_lv.LevelIndex >= 1;
                    useDirect = op_lv.LevelIndex >= 2;


                    //if (terminateMGRecursion) {
                    //    // skipping of levels makes convergence bad.
                    //    useDirect = true;
                    //    terminateMGRecursion = false;
                    //}

                    tr.Info($"KcycleMultiSchwarz: lv {iLevel}, DOFs: " + op_lv.Mapping.MpiSize.ForLoop(rnk => op_lv.Mapping.GetLocalLength(rnk)).ToConcatString("[", ",", "]") + ", sum = " + op_lv.Mapping.TotalLength);
                    tr.Info($"KcycleMultiSchwarz: lv {iLevel}, Cells " + op_lv.Mapping.MpiSize.ForLoop(rnk => op_lv.Mapping.GetLocalNoOfBlocks(rnk)).ToConcatString("[", ",", "]") + ", sum = " + op_lv.Mapping.TotalNoOfBlocks);

                    if (useDirect)
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, Direct solver ");
                    else if (terminateMGRecursion) {
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, GMRES with Level PMG");
                    } else {
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, Schwarz, no of blocks : " + (new int[] { NoBlocks[iLevel] }).MPIAllGather().ToConcatString("[", ",", "]") + " sum = " + GlobalNoBlocks[iLevel]);
                    }


                    ISolverSmootherTemplate levelSolver;

                    if (useDirect) {

                        var _levelSolver = new DirectSolver();
                        levelSolver = _levelSolver;
                        _levelSolver.config.WhichSolver = DirectSolver._whichSolver.PARDISO;
                        _levelSolver.config.TestSolution = false;
                        _levelSolver.ActivateCaching = (__1, __2) => true;

                        // Not really recommend ILU for coarse: for Botti2DStokes, 3x int time, 2x in iterations
                        //
                        //var _levelSolver = new CellILU();
                        //levelSolver = _levelSolver;
                        //_levelSolver.ILU_level = 0;

                    } else if (terminateMGRecursion) {

                        var coarseConfig = new PTGconfig() {
                            pMaxOfCoarseSolver = this.pMaxOfCoarseSolver,
                        };

                        var _levelSolver2 = coarseConfig.CreateInstance(op_lv);
                        (_levelSolver2 as SoftGMRES).TerminationCriterion = delegate (int i, double r0, double r) {
                            var ret = (i <= 1 || r > r0 * 0.01, true);
                            return ret;
                        };
                        levelSolver = _levelSolver2;



                    } else {


                        var smoother1 = new Schwarz() {
                            FixedNoOfIterations = 1, //iLevel == 1 ? 2 : 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsOnCurrentProcess = NoBlocks[iLevel],
                            },
                            //ActivateCachingOfBlockMatrix = SmootherCaching,

                        };
                        smoother1.config.EnableOverlapScaling = true;
                        smoother1.config.Overlap = 1; // overlap seems to help; more overlap seems to help more
                        smoother1.config.UsePMGinBlocks = this.UsepTG;


                        //var smoother1 = new SchwarzForCoarseMesh();
                        //smoother1.config.NoOfBlocks = NoBlocks[iLevel].MPISum();
                        //smoother1.config.EnableOverlapScaling = true;
                        //

                        var _levelSolver = new OrthonormalizationMultigrid() {
                            PreSmoother = smoother1,
                            PostSmoother = smoother1
                        };
                        _levelSolver.config.m_omega = 1; // v-cycle


                        _levelSolver.config.NoOfPostSmootherSweeps = (int)Math.Max(2, Math.Max(maxDG, Math.Round(Math.Log10(GlobalNoBlocks[iLevel]) * 3.0) + maxDG - 2));
                        if (iLevel == 0)
                            _levelSolver.config.NoOfPostSmootherSweeps = 20;
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, NoOfPostSmootherSweeps = {_levelSolver.config.NoOfPostSmootherSweeps}");


                        if (iLevel == 0) {
                            (_levelSolver).TerminationCriterion = this.DefaultTermination;
                            //if(maxDG > 3)
                            //    _levelSolver.config.NoOfPostSmootherSweeps = 6;
                        } else if (iLevel == 1) {

                            _levelSolver.config.m_omega = 2;

                            (_levelSolver).TerminationCriterion = delegate (int i, double r0, double r) {
                                var ret = (i <= 1 || r > r0 * 0.01, true);
                                return ret;
                            };

                            var altSmooth3 = new SoftGMRES();
                            altSmooth3.MaxKrylovDim = 80;
                            altSmooth3.TerminationCriterion = (iIter, r0, ri) => {
                                return (iIter <= 82 && ri > 1e-7*r0, ri <= 1e-7*r0);
                            };
                            altSmooth3.Precond = new BlockJacobi();


                            _levelSolver.AdditionalPostSmoothers = new ISolverSmootherTemplate[] {
                                  altSmooth3
                                    };

                        } else {
                            (_levelSolver).TerminationCriterion = (i, r0, r) => (i <= 1, true);
                        }
                        levelSolver = _levelSolver;


                        // Extended Multigrid Analysis
                        //((OrthonormalizationMultigrid)levelSolver).IterationCallback += MultigridAnalysis;                    
                    }
                    SolverChain.Add(levelSolver);

                    if (iLevel > 0) {
                        if (SolverChain[iLevel - 1] is OrthonormalizationMultigrid omg_fine)
                            omg_fine.CoarserLevelSolver = levelSolver;
                        else if (SolverChain[iLevel - 1] is GenericRestriction rest)
                            rest.CoarserLevelSolver = levelSolver;
                    }

                    if (useDirect || terminateMGRecursion)
                        break;
                }

                return SolverChain[0];
            }
        } */

    }


  
}
