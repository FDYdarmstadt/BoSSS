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
        /// Use a p-two-grid approach in Additive Schwarz (<see cref="Schwarz.UsePMGinBlocks"/>).
        /// 
        /// This results in a lesser number of Schwarz blocks and more computational cells per Schwarz block.
        /// It seems to accelerate convergence for certain problems (mostly single-phase).
        /// The order of the low-order system can be adjusted using <see cref="pMaxOfCoarseSolver"/>.
        /// </summary>
        [DataMember]
        bool UsepTG = false; 

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
            Func<int, int> SblkSizeFunc = delegate (int iLevel) { return TargetBlockSize; };
            var instance = KcycleMultiSchwarz(level, SblkSizeFunc);
            instance.Init(level);
            return instance;
        }



        /// <summary>
        /// Determine max MG level from target blocksize
        /// </summary>
        private int GetMGDepth(MultigridOperator op) {
            int MSLength = op.NoOfLevels;
            int DirectKickIn = this.TargetBlockSize;
            int MGDepth = Math.Min(MSLength, NoOfMultigridLevels);
            //var PrevSize = GlobalDOF[0]; // system size on finer level
            long PrevSize = op.FinestLevel.Mapping.TotalLength;

            int iLevel = 0;
            for (var level = op; level != null && iLevel < MGDepth; iLevel++) {
                long SysSize = op.Mapping.TotalLength;
                //Console.WriteLine("DOF on L{0}: {1}", iLevel, SysSize);
                bool useDirect = false;
                useDirect |= (SysSize < DirectKickIn);
                //useDirect |= (double)PrevSize / (double)SysSize < 1.5 && SysSize < 50000; // degenerated MG-Agglomeration, because too few candidates
                useDirect |= iLevel == NoOfMultigridLevels - 1;
                useDirect |= iLevel == MGDepth - 1; // otherwise error, due to missing coarse solver
                Debug.Assert(useDirect.MPIEquals());

                if (useDirect) {
                    break;
                }

                PrevSize = SysSize;
                level = op.CoarserLevel;
            }
            //Console.WriteLine("INFO: using {0} levels, lowest level DOF is {1}, target size is {2}.", iLevel + 1, GlobalDOF[iLevel], DirectKickIn);
            //MaxUsedMGLevel = iLevel; // remember this for saving in query
            return iLevel + 1;
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

            if(SizeFraction < 1) {
                Console.WriteLine($"WARNING: local system size ({LocalDOF4directSolver}) < Schwarz-Block size ({SchwarzblockSize});");
                Console.WriteLine($"resetting local number of Schwarz-Blocks to 1.");
            }
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
                int LocalNoOfSchwarzBlocks = NoOfBlocksAtLevel(op, iLevel, pCoarsest, targetblocksize);
                int TotalNoOfSchwarzBlocks = LocalNoOfSchwarzBlocks.MPISum();
                //SetQuery("GlobalSblocks at Lvl" + iLevel, TotalNoOfSchwarzBlocks, true);
                NoOfBlocks[iLevel] = LocalNoOfSchwarzBlocks;
            }
            return NoOfBlocks;
        }

        private int getMaxDG (MultigridOperator op, int iLevel, int iVar) {
            //This workaround takes into account, the wierd structure of the ChangeOfBasis-thing
            //prohibits out of bounds exception

            var MGChangeOfBasis = op.GetLevel(iLevel).Config;

            int tVar = iVar < MGChangeOfBasis.Length ? iVar : MGChangeOfBasis.Length - 1;
            return MGChangeOfBasis[tVar].DegreeS[0];
        }

        
        ISolverSmootherTemplate KcycleMultiSchwarz(MultigridOperator op, Func<int,int> SchwarzblockSize) {
            using(var tr = new FuncTrace()) {
                int MSLength = op.NoOfLevels;


                var SolverChain = new List<ISolverSmootherTemplate>();

                int maxDG = getMaxDG(op, 0, 0);
                
                int pCoarsest = pMaxOfCoarseSolver < maxDG && UsepTG ? pMaxOfCoarseSolver : -1;
                int[] NoBlocks = NoOfSchwarzBlocks(op, pCoarsest, SchwarzblockSize);
                int[] GlobalNoBlocks = NoBlocks.Select(No => No.MPISum()).ToArray();
                int MGDepth = GetMGDepth(op);
                int LastLevel = MGDepth - 1;

                Func<int, int, int, bool> SmootherCaching = delegate (int Iter, int MgLevel, int iBlock) {
                    //return Iter >= ((MaxMGLevel - MgLevel) * 3) * (1);
                    return true;
                };
                Func<int, int, bool> CoarseCaching = delegate (int Iter, int MgLevel) {
                    //return Iter >= MaxMGDepth+1;
                    return true;
                };

                for(int iLevel = 0; iLevel < MGDepth; iLevel++) {

                    bool useDirect = iLevel == LastLevel;


                    if(useDirect)
                        Console.WriteLine("KcycleMultiSchwarz: lv {0}, Direct solver ", iLevel);
                    else
                        Console.WriteLine("KcycleMultiSchwarz: lv {0}, no of blocks {1} : ", iLevel, GlobalNoBlocks[iLevel]);

                    ISolverSmootherTemplate levelSolver;

                    if(useDirect) {
                        var _levelSolver = new DirectSolver();
                        levelSolver = _levelSolver;
                        _levelSolver.config.WhichSolver = DirectSolver._whichSolver.PARDISO;
                        _levelSolver.config.TestSolution = false;
                        _levelSolver.ActivateCaching = CoarseCaching;


                    } else {

                        var smoother1 = new Schwarz() {
                            FixedNoOfIterations = 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsOnCurrentProcess = NoBlocks[iLevel],
                            },
                            ActivateCachingOfBlockMatrix = SmootherCaching,
                            Overlap = 1, // overlap seems to help; more overlap seems to help more
                            EnableOverlapScaling = true,
                            UsePMGinBlocks = UsepTG
                        };

                        var _levelSolver = new OrthonormalizationMultigrid() {
                            PreSmoother = smoother1,
                            PostSmoother = smoother1
                        };
                        _levelSolver.config.m_omega = 1;

                        if(iLevel > 0) {
                            (_levelSolver).TerminationCriterion = (i, r0, r) => (i <= 1, true);
                        } else {
                            (_levelSolver).TerminationCriterion = this.DefaultTermination;
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
                        ((OrthonormalizationMultigrid)(SolverChain[iLevel - 1])).CoarserLevelSolver = levelSolver;
                    }
                }

                return SolverChain[0];
            }
        }

    }
}
