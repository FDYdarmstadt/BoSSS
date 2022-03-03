using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {
    
    /// <summary>
    /// GMRES with p-two-grid preconditioner (<see cref="LevelPmg"/>).
    /// 
    /// Should work well for intermediate-size-systems, where the embedded low-order system on the finest mesh is sufficiently small for the 
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
        /// Hack, for the treatment of incompressible flows:
        /// - false (default): the D-th variable, where D is the spatial dimension (2 or 3), is assumed to be the pressure; the order is one lower than for velocity
        /// - true: no special treatment of individual variables
        /// </summary>
        public bool EqualOrder = false;

        /// <summary>
        /// If true, cell-local solvers will be used to approximate a solution to high-order modes
        /// </summary>
        public bool UseHiOrderSmoothing = true;

        ///// <summary>
        ///// DG degree at low order blocks. This degree is the border, which divides into low order and high order blocks
        ///// </summary>
        //public int OrderOfCoarseSystem = 1;

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
            var precond = new LevelPmg();
            precond.config.UseHiOrderSmoothing = this.UseHiOrderSmoothing;
            precond.config.OrderOfCoarseSystem = this.pMaxOfCoarseSolver;
            precond.config.FullSolveOfCutcells = this.FullSolveOfCutcells;
            precond.config.UseDiagonalPmg = this.UseDiagonalPmg;
            precond.config.EqualOrder = this.EqualOrder;

            var templinearSolve = new SoftGMRES() {
                MaxKrylovDim = MaxKrylovDim,
                Precond = precond,
                TerminationCriterion = this.DefaultTermination
            };


            templinearSolve.Init(level);
            return templinearSolve;
        }

    }


    public class PTGconfig2 : PTGconfig {


        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            var precond = new LevelPmg();
            precond.config.UseHiOrderSmoothing = false;
            precond.config.OrderOfCoarseSystem = this.pMaxOfCoarseSolver;
            precond.config.FullSolveOfCutcells = this.FullSolveOfCutcells;
            precond.config.UseDiagonalPmg = this.UseDiagonalPmg;
            precond.config.EqualOrder = this.EqualOrder;
            

            var bj = new BlockJacobi() {
                NoOfIterations = 1,
                
            };

            /*
            var smoother = new SoftGMRES() {
                MaxKrylovDim = 10,
                TerminationCriterion = (int iter, double r0, double ri) => (iter <= 10, true),
                Precond = bj
            };
            */

            var loPsol = new DirectSolver() {

            };

            var coarseSolver = new PRestriction() {
                LowerPSolver = loPsol
            };


            var templinearSolve = new OrthonormalizationMultigrid() {
                //CoarserLevelSolver = loPsol,
                PreSmoother = coarseSolver,
                TerminationCriterion = this.DefaultTermination,
                PostSmoother = bj,
                //CoarserLevelSolver = loPsol
            };




            templinearSolve.Init(level);
            return templinearSolve;
        }
    }
}
