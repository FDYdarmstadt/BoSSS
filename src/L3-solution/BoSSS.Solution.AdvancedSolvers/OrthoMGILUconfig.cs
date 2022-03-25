using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Dynamic configuration for the orthonormalization multigrid.
    /// As smoothers, <see cref="CellILU"/> is used on all mesh levels, except the coarsest one.
    /// There, a direct solver is used.
    /// The coarsest level is determined dynamically, depending on the problem and the mesh size using <see cref="TargetBlockSize"/>.
    /// </summary>
    public class OrthoMGILUconfig : IterativeSolverConfig {
        
        /// <summary>
        /// 
        /// </summary>
        public override string Name => "OrthonormalizationMultigrid with CellILU Smoother";


        /// <summary>
        /// 
        /// </summary>
        public override string Shortname => "OrthoMG w ILU";

        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            var instance = KcycleMultiILU(level);
            instance.Init(level);
            return instance;
        }

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


        ///// <summary>
        ///// Sets the **maximum** number of Multigrid levels to be used.
        ///// The numbers of levels which are actually used is probably much less, and determined via <see cref="TargetBlockSize"/>.
        ///// Multigrid approach is used to get a preconditioner for Krylov solvers, e.g. GMRES.
        ///// </summary>
        //[DataMember]
        //public int NoOfMultigridLevels = 1000000;

        ISolverSmootherTemplate KcycleMultiILU(MultigridOperator level) {

            //MultigridOperator Current = op;
            var SolverChain = new List<ISolverSmootherTemplate>();
            
            for (MultigridOperator mgop = level; mgop != null; mgop = mgop.CoarserLevel) {
                int iLevel = mgop.LevelIndex;

                bool useDirect = false;
                useDirect |= mgop.CoarserLevel == null;
                useDirect |= (mgop.Mapping.TotalLength < this.TargetBlockSize);
                useDirect = useDirect.MPIOr();

                if (useDirect)
                    Console.WriteLine("KcycleMultiILU: lv {0}, Direct solver ", iLevel);
                else
                    Console.WriteLine("KcycleMultiILU: lv {0}, ", iLevel);

                ISolverSmootherTemplate levelSolver;
                if (useDirect) {
                    var _levelSolver = new DirectSolver();
                    _levelSolver.config.WhichSolver = DirectSolver._whichSolver.PARDISO;
                    _levelSolver.config.TestSolution = false;
                    levelSolver = _levelSolver;

                } else {
                    //var smoother1 = new HypreILU() {
                    //    LocalPreconditioning = true
                    //};

                    var smoother1 = new CellILU() {
                        ILU_level = iLevel == 0 ? 1 : 0 // Use ILU(1) on fine mesh and ILU(0) on all others, because:
                        //                                 - ILU(2) already seems to be to expensive and not very beneficial for high DG degrees
                        //                                 - on coarse mesh levels, where cells have more and more neighbors, we can only afford ILU(0)
                    };


                    var _levelSolver = new OrthonormalizationMultigrid() {
                        PreSmoother = smoother1,
                        PostSmoother = smoother1,
                    };
                    _levelSolver.config.m_omega = 1;
                    levelSolver = _levelSolver;

                    if (iLevel > 0) {
                        _levelSolver.TerminationCriterion = (i, r0, r) => (i <= 1, true);
                    } else {
                        _levelSolver.TerminationCriterion = this.DefaultTermination;
                    }

                }

                if(iLevel > 0) {
                    ((OrthonormalizationMultigrid)(SolverChain[iLevel - 1])).CoarserLevelSolver = levelSolver;

                }
                SolverChain.Add(levelSolver);

                if(useDirect)
                    break;
               
            }

            return SolverChain[0];
        }
    }
}
