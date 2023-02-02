using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Dynamic configuration for the orthonormalization multigrid.
    /// As smoothers, additive Schwarz (<see cref="Schwarz"/>) is used on all mesh levels, except the coarsest one.
    /// There, a direct solver is used.
    /// The coarsest level is determined dynamically, depending on the problem and the mesh size using <see cref="TargetBlockSize"/>.
    /// </summary>
    [Serializable]
    public class OrthoMGSchwarzConfig : IterativeSolverConfig {
        
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
        /// Determines maximal DG order within coarse system of a p-Multigrid. 
        /// </summary>
        [DataMember]
        public int pMaxOfCoarseSolver = 1;

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


        int GetLocalDOF(MultigridOperator op, int pOfLowOrderSystem) {
            if(pOfLowOrderSystem <= -1) {
                return op.Mapping.LocalLength;
            } else {

                
                var map = op.Mapping;
                int D = map.AggGrid.SpatialDimension;
                int acc = 0;
                int NoOfVars = map.NoOfVariables;
                int J = map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                var bS = map.AggBasis;
                for(int j = 0; j < J; j++) {
                    for(int iVar = 0; iVar < NoOfVars; iVar++) {
                        int p = (iVar == D ? pOfLowOrderSystem - 1 : pOfLowOrderSystem);
                        acc += bS[iVar].GetLength(j, p);
                    }
                }
                return acc;
            }
        }

        /// <summary>
        /// Number of Blocks at this level. Minimum 1 per core
        /// </summary>
        private int NoOfBlocksAtLevel(MultigridOperator op, int MGLevel = 0, int pCoarsest = -1, Func<int, int> targetblocksize = null) {
            if(targetblocksize == null)
                targetblocksize = (int iLevel) => TargetBlockSize;
            int SchwarzblockSize = targetblocksize(MGLevel);
            int LocalDOF4directSolver = GetLocalDOF(op, pCoarsest);
            double SizeFraction = (double)LocalDOF4directSolver / (double)SchwarzblockSize;
            //SizeFraction = SizeFraction.MPIMax(); [Toprak] the bug cannot be reproduced, therefore, it is reverted to the previous version.
            //when the number of blocks are not the same on the processors, it causes deadlock in MPI processes
            //Hard solution would be checking all MPI processes written in SubBlocks and modify them to allow different numbered block partitioning 

            //if(SizeFraction < 1) {
            //    Console.WriteLine($"WARNING: local system size ({LocalDOF4directSolver}) < Schwarz-Block size ({SchwarzblockSize});");
            //    Console.WriteLine($"resetting local number of Schwarz-Blocks to 1.");
            //}
            int LocalNoOfSchwarzBlocks = Math.Max(1, (int)Math.Floor(SizeFraction));
            return LocalNoOfSchwarzBlocks;
        }


        /// <summary>
        /// No Of Local Schwarz Blocks (process local) for all MG levels
        /// </summary>
        private int[] NoOfSchwarzBlocks(MultigridOperator op, int pCoarsest = -1, Func<int, int> targetblocksize = null) {
            int MSLength = op.NoOfLevels;
            var NoOfBlocks = MSLength.ForLoop(level => -1);
            for(int iLevel = 0; iLevel < MSLength; iLevel++) {
                int LocalNoOfSchwarzBlocks = NoOfBlocksAtLevel(op.GetLevel(iLevel), iLevel, pCoarsest, targetblocksize);
                NoOfBlocks[iLevel] = LocalNoOfSchwarzBlocks;
            }
            return NoOfBlocks;
        }
        
        ISolverSmootherTemplate KcycleMultiSchwarz(MultigridOperator op, Func<int,int> SchwarzblockSize) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                int MSLength = op.NoOfLevels;

                int mpiSz = op.Mapping.MpiSize;

                var SolverChain = new List<ISolverSmootherTemplate>();

                int maxDG = op.DGpolynomialDegreeHierarchy.First().Max();
                int minDG = Math.Max(op.DGpolynomialDegreeHierarchy.Last().Max(), pMaxOfCoarseSolver);
                minDG = Math.Min(minDG, maxDG); // make sure that minDG <= MaxDG;

                int pCoarsest = UsepTG ? minDG : -1;
                int[] NoBlocks = NoOfSchwarzBlocks(op, pCoarsest, SchwarzblockSize);
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

                    bool skipLevel = false;
                    if(useDirect == false && iLevel > 0) {
                        if (OneBlockPerNode[iLevel] && OneBlockPerNode[iLevel - 1]) {
                            // $"Only one block per node on level {iLevel}, and should be skipped, but cannot be skipped be skipped, but not supported yet."
                            //throw new NotSupportedException($"Only one block per node on level {iLevel}, and should be skipped, but cannot be skipped be skipped, but not supported yet.");
                            skipLevel = true;
                        }
                    }
                    Debug.Assert(useDirect.MPIEquals(op.Mapping.MPI_Comm));
                    Debug.Assert(skipLevel.MPIEquals(op.Mapping.MPI_Comm));

                    //useDirect = op_lv.LevelIndex >= 1;


                    if (skipLevel) {
                        useDirect = true;
                        skipLevel = false;
                    }

                    if (useDirect)
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, L = {op_lv.Mapping.TotalLength}, Direct solver ");
                    else if (skipLevel) {
                        tr.Info($"Only one block per node on levels {iLevel} and {iLevel - 1}, so {iLevel} will be skipped.");
                    } else {
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, L = {op_lv.Mapping.TotalLength}, no of blocks total: {GlobalNoBlocks[iLevel]}");                       
                    }


                    ISolverSmootherTemplate levelSolver;

                    if(useDirect) {
                        
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
                        

                    } else if (skipLevel) {
                        var _levelSolver = new GenericRestriction();
                        levelSolver = _levelSolver;
                    } else {

                        var smoother1 = new Schwarz() {
                            FixedNoOfIterations = 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsOnCurrentProcess = NoBlocks[iLevel],
                            },
                            //ActivateCachingOfBlockMatrix = SmootherCaching,
                            
                        };
                        smoother1.config.EnableOverlapScaling = true;
                        smoother1.config.Overlap = 1; // overlap seems to help; more overlap seems to help more
                        smoother1.config.UsePMGinBlocks = UsepTG;

                        

                        var _levelSolver = new OrthonormalizationMultigrid() {
                            PreSmoother = smoother1,
                            PostSmoother = smoother1,
                           
                        };
                        _levelSolver.config.m_omega = 1; // v-cycle
                        /*
                        if (iLevel == 0) {
                            tr.Info("Using additional ILU smoother on level " + iLevel);
                            var otherSmoother = new CellILU() {

                            };
                            _levelSolver.config.NoOfPostSmootherSweeps = 4;
                            _levelSolver.AdditionalPostSmoothers = new ISolverSmootherTemplate[] { otherSmoother };
                        }
                        */
                        //if (iLevel == 0) {
                        //    tr.Info("Using additional post smoother sweeps on level " + iLevel);
                        //    _levelSolver.config.NoOfPostSmootherSweeps = 8;
                        //}

                        _levelSolver.config.NoOfPostSmootherSweeps = (int)Math.Max(2, Math.Max(maxDG, Math.Round(Math.Log10(GlobalNoBlocks[iLevel]) * 3.0) + maxDG - 2));
                        tr.Info($"KcycleMultiSchwarz: lv {iLevel}, NoOfPostSmootherSweeps = {_levelSolver.config.NoOfPostSmootherSweeps}");


                        if (iLevel == 0) {
                            (_levelSolver).TerminationCriterion = this.DefaultTermination;
                            //if(maxDG > 3)
                            //    _levelSolver.config.NoOfPostSmootherSweeps = 6;
                        } else if(iLevel == 1) {
                            (_levelSolver).TerminationCriterion = delegate (int i, double r0, double r) {
                                var ret = (i <= 1 || r > r0 * 0.01, true);
                                //if(ret.Item1 == false)
                                //    Console.WriteLine($"level 1: {i} {r / r0:0.##e-00} {ret}");
                                return ret;

                            };

                        } else {
                            (_levelSolver).TerminationCriterion = (i, r0, r) => (i <= 1, true);
                        }
                        levelSolver = _levelSolver;

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

                    if(iLevel > 0) {
                        if(SolverChain[iLevel - 1] is OrthonormalizationMultigrid omg_fine)
                            omg_fine.CoarserLevelSolver = levelSolver;
                        else if(SolverChain[iLevel - 1] is GenericRestriction rest)
                            rest.CoarserLevelSolver = levelSolver;
                    }

                    if (useDirect)
                        break;
                }

                return SolverChain[0];
            }
        }

    }


  
}
