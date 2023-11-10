using ilPSP;
using ilPSP.Tracing;
using log4net.Core;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// GMRES with p-two-grid preconditioner (<see cref="LevelPmg"/>).
    /// 
    /// Should work well for intermediate-size-systems, 
    /// where the embedded low-order system on the finest mesh is sufficiently small for the 
    /// direct solver.
    /// </summary>
    [Serializable]
    public class PTGconfig : IterativeSolverConfig {

        /// <summary>
        /// Determines maximal DG order within coarse system of the p-Multigrid. 
        /// </summary>
        public int pMaxOfCoarseSolver = 1;

        /// <summary>
        /// If true, the high order system is solved cell-by-cell (i.e. a Block-Jacobi)
        /// </summary>
        public bool UseDiagonalPmg = true;

        /// <summary>
        /// If true, cell-local solvers will be used to approximate a solution to high-order modes
        /// (this is equivalent to a Block-Jacobi approach)
        /// </summary>
        public bool UseHiOrderSmoothing = true;

        /// <summary>
        /// If true blocks/cells containing more than one species are completely assigned to low order block solver.
        /// This hopefully is better than the default approach
        /// </summary>
        public bool FullSolveOfCutcells = true;

        /// <summary>
        /// number of iterations between restarts
        /// </summary>
        public int MaxKrylovDim = 80;

        /// <summary>
        /// 
        /// </summary>
        public override string Name => "GMRES with p-two-grid preconditioner";

        /// <summary>
        /// 
        /// </summary>
        public override string Shortname => "GMRES w p2G";



        /// <summary>
        /// 
        /// </summary>
        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {

            //var precond = CreatePrecond_alt(level);
            var precond = CreatePrecond();

            var templinearSolve = new SoftGMRES() {
                MaxKrylovDim = MaxKrylovDim,
                Precond = precond,
                TerminationCriterion = this.DefaultTermination
            };


            templinearSolve.Init(level);
            return templinearSolve;
        }


        private ISolverSmootherTemplate CreatePrecond() {
            var precond = new LevelPmg();
            precond.config.UseHiOrderSmoothing = this.UseHiOrderSmoothing;
            precond.config.OrderOfCoarseSystem = this.pMaxOfCoarseSolver;
            precond.config.FullSolveOfCutcells = this.FullSolveOfCutcells;
            precond.config.UseDiagonalPmg = this.UseDiagonalPmg;
            return precond;
        }

        int[] GetLoOrderDegrees(MultigridOperator level) {
            int[] Degrees = level.Degrees;
            int DegDecrease = Math.Max(Degrees.Max() - this.pMaxOfCoarseSolver, 0);
            int[] loDegrees = Degrees.Select(deg => Math.Max(deg - DegDecrease, 0)).ToArray();
            return loDegrees;
        }

        private ISolverSmootherTemplate CreatePrecond_alt(MultigridOperator level) {

            var lowSolver = new DirectSolver();
            lowSolver.config.WhichSolver = DirectSolver._whichSolver.PARDISO;
            lowSolver.config.UseDoublePrecision = false;
            lowSolver.config.TestSolution = false;
            lowSolver.ActivateCaching = (int NoIter, int MgLevel) => true;

            var coarseSolver = new GridAndDegRestriction() {
                RestrictedDeg = GetLoOrderDegrees(level),
                LowerPSolver = lowSolver
            };

            /*
            var postSmoother = new BlockJacobi() {
                NoOfIterations = 1,
            };
            //*/

            var postSmoother = new Schwarz() {
                FixedNoOfIterations = 1,
                m_BlockingStrategy = new Schwarz.MultigridBlocks() { 
                    Depth = 0
                }
            };
            postSmoother.config.EnableOverlapScaling = true;
            postSmoother.config.Overlap = 1;
            postSmoother.config.UsePMGinBlocks = false;
            //*/
            //Debugger.Launch();

            var precond = new ClassicMultigrid();
            precond.PreSmoother = coarseSolver;
            precond.PostSmoother = postSmoother;
            precond.config.m_omega = 1;
            precond.TerminationCriterion = (iter, R0_l2, R_l2) => (iter <= 1, false);



            return precond;
        }

    }

 
}
