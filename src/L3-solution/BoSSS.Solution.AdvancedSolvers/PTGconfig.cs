using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
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

        ///// <summary>
        ///// Hack, for the treatment of incompressible flows:
        ///// - false (default): the D-th variable, where D is the spatial dimension (2 or 3), is assumed to be the pressure; the order is one lower than for velocity
        ///// - true: no special treatment of individual variables
        ///// </summary>
        //public bool EqualOrder = false;

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

            var templinearSolve = new SoftGMRES() {
                MaxKrylovDim = MaxKrylovDim,
                Precond = precond,
                TerminationCriterion = this.DefaultTermination
            };


            templinearSolve.Init(level);
            return templinearSolve;
        }

    }

    /// <summary>
    /// Orthonormalization (<see cref="OrthonormalizationMultigrid"/>) scheme with full p-multigrid (i.e. multigrid over DG polynomial degree) preconditioning.
    /// In comparison to <see cref="PTGconfig"/>, where only two p-levels are used, this uses all available p-levels.
    /// </summary>
    public class PmgConfig : IterativeSolverConfig {


        /// <summary>
        /// 
        /// </summary>
        public override string Name => "Orthonormalization with p-multi-grid preconditioner";

        /// <summary>
        /// 
        /// </summary>
        public override string Shortname => "Ortho w pmG";
  
        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            //var precond = new LevelPmg();
            //precond.config.UseHiOrderSmoothing = false;
            //precond.config.OrderOfCoarseSystem = this.pMaxOfCoarseSolver;
            //precond.config.FullSolveOfCutcells = this.FullSolveOfCutcells;
            //precond.config.UseDiagonalPmg = this.UseDiagonalPmg;
            //precond.config.EqualOrder = this.EqualOrder;

            
           

            int[][] DegreeHierarchy = level.DGpolynomialDegreeHierarchy;
            
            ISubsystemSolver CreateLevelRecursive(int iLevel) {

                //level.Mapping.

                //long globalDOF;
                //int localDOF;


               

                OrthonormalizationMultigrid createOMG() {

                    Console.WriteLine("OrthoMG solver on level " + iLevel + ", restricting to Degree " + DegreeHierarchy[iLevel + 1].ToConcatString("[", ",", "]"));


                    var coarseSolver = new PRestriction() {
                        RestrictedDeg = DegreeHierarchy[iLevel + 1].CloneAs(),
                        LowerPSolver = CreateLevelRecursive(iLevel + 1)
                    };

                    ISolverSmootherTemplate post;

                    //var bj = new BlockJacobi() {
                    //    NoOfIterations = 1,

                    //};

                    //GMRES bringt nix als PC
                    //var post = new SoftGMRES() {
                    //    //Precond = bj
                    //};
                    //post.TerminationCriterion = (int iter, double R0_l2, double R_l2) => (iter <= 1, true);

                    if(iLevel == 0) {
                        post = new CellILU() {
                            ILU_level = 0
                        };
                    } else {
                        post = new BlockJacobi() {
                            NoOfIterations = 1,
                        };
                    }


                    var omg = new OrthonormalizationMultigrid() {
                        //CoarserLevelSolver = loPsol,
                        PreSmoother = coarseSolver,
                        PostSmoother = post,
                        //CoarserLevelSolver = loPsol
                    };
                    omg.config.NoOfPostSmootherSweeps = 10;
                    return omg;
                }


                if(iLevel >= DegreeHierarchy.Length - 1) {

                    Console.WriteLine("direct solver on level " + iLevel);

                    // ++++++++++++
                    // lowest level 
                    // ++++++++++++
                    var s = new DirectSolver();
                    s.config.TestSolution = false;
                    s.ActivateCaching = (a, b) => true;
                    return s;
                } else if(iLevel >= 2) {
                    // +++++++++++++++++++++++++++
                    // intermediate level / coarse
                    // +++++++++++++++++++++++++++

                    Console.WriteLine("GMRES solver on level " + iLevel);

                    var gmRes = new SoftGMRES();
                    gmRes.TerminationCriterion = (int iter, double R0_l2, double R_l2) => (iter <= 4, true);
                    gmRes.Precond = CreateLevelRecursive(iLevel + 1);
                    return gmRes;

                } else if(iLevel >= 1) { 
                    // +++++++++++++++++++++++++++
                    // intermediate level / fine
                    // +++++++++++++++++++++++++++


                    var omg = createOMG();


                    (bool bNotTerminate, bool bSuccess) SubLevelTermination(int iter, double R0_l2, double R_l2) {
                        double reduction = iLevel >= 2 ? 1e-2 : 1e-2;
                        int limit = iLevel >= 2 ? 4 : 100;

                        if(iter <= 1)
                            return (true, false); // keep running

                        if(R_l2 < R0_l2 * reduction) {
                            //Console.WriteLine("succ " + iter + " iterations");
                            return (false, true); // success
                        }


                        if(iter > limit) {
                            //Console.WriteLine("running out of iterations before sufficient residual reduction");
                            return (false, false); // fail
                        }


                        return (true, false); // keep running
                    }

                    omg.TerminationCriterion = SubLevelTermination;
                    return omg;

                } else {
                    // ++++++++++++
                    // finest level
                    // ++++++++++++
                    var omg = createOMG();
                    omg.TerminationCriterion = this.DefaultTermination;
                    return omg;
                }

                    
            }


            var topLevel = CreateLevelRecursive(0);
            topLevel.Init(level);
            return topLevel;            
        }
    
        //*/
    }
}
