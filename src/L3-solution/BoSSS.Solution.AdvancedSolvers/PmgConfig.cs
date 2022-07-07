using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {
    /// <summary>
    /// Orthonormalization (<see cref="OrthonormalizationMultigrid"/>) scheme with full p-multigrid (i.e. multigrid over DG polynomial degree) preconditioning.
    /// In comparison to <see cref="PTGconfig"/>, where only two p-levels are used, this uses all available p-levels.
    /// </summary>
    [Serializable]
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
            return CreateInstanceImpl(level, level.DGpolynomialDegreeHierarchy);
        }

        /// <summary>
        /// - true: <see cref="CellILU"/> is used as smoother on finest p-level
        /// - false: <see cref="BlockJacobi"/> is used as smoother on finest p-level
        /// </summary>
        public bool UseILU = true;

        internal ISubsystemSolver CreateInstanceImpl(IOperatorMappingPair level, int[][] DegreeHierarchy) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                //var precond = new LevelPmg();
                //precond.config.UseHiOrderSmoothing = false;
                //precond.config.OrderOfCoarseSystem = this.pMaxOfCoarseSolver;
                //precond.config.FullSolveOfCutcells = this.FullSolveOfCutcells;
                //precond.config.UseDiagonalPmg = this.UseDiagonalPmg;
                //precond.config.EqualOrder = this.EqualOrder;

                ISubsystemSolver CreateLevelRecursive(int iLevel) {

                    //level.Mapping.

                    //long globalDOF;
                    //int localDOF;




                    OrthonormalizationMultigrid createOMG() {

                        tr.Info("OrthoMG solver on level " + iLevel + ", restricting to Degree " + DegreeHierarchy[iLevel + 2].ToConcatString("[", ",", "]"));


                        var coarseSolver = new GridAndDegRestriction() {
                            RestrictedDeg = DegreeHierarchy[iLevel + 2].CloneAs(),
                            LowerPSolver = CreateLevelRecursive(iLevel + 1)
                        };

                        ISolverSmootherTemplate post, pre;

                        //var bj = new BlockJacobi() {
                        //    NoOfIterations = 1,

                        //};

                        //GMRES bringt nix als PC
                        //var post = new SoftGMRES() {
                        //    //Precond = bj
                        //};
                        //post.TerminationCriterion = (int iter, double R0_l2, double R_l2) => (iter <= 1, true);

                        if (iLevel == 0 && UseILU) {
                            //pre = new BlockJacobi() {
                            //    NoOfIterations = 1,
                            //};
                            pre = null;


                            post = new CellILU() {
                                ILU_level = 0
                            };

                            ((CellILU)post).id = "Lv0";

                            //post = new SparseILU() {
                            //    UsedLibrary = SparseILU.Library.Intel_MKL
                            //};

                        } else {
                            post = new BlockJacobi() {
                                NoOfIterations = 1,
                            };
                            pre = null;
                        }


                        var omg = new OrthonormalizationMultigrid() {
                            CoarserLevelSolver = coarseSolver,
                            PreSmoother = pre,
                            PostSmoother = post,
                        };
                        omg.config.NoOfPostSmootherSweeps = 10;
                        omg.config.m_omega = 1;
                        omg.config.CoarseOnLovwerLevel = false;
                        return omg;
                    }

                    Console.WriteLine("###########################   testcode aktiv ##################################");
                    if (iLevel >= 1/* DegreeHierarchy.Length - 1*/) {

                        tr.Info("direct solver on level " + iLevel);

                        // ++++++++++++
                        // lowest level 
                        // ++++++++++++
                        var s = new DirectSolver();
                        s.config.TestSolution = false;
                        s.ActivateCaching = (a, b) => true;
                        return s;
                    } else if (iLevel >= 2) {
                        // +++++++++++++++++++++++++++
                        // intermediate level / coarse
                        // +++++++++++++++++++++++++++

                        tr.Info("GMRES solver on level " + iLevel);

                        var gmRes = new SoftGMRES();
                        gmRes.TerminationCriterion = (int iter, double R0_l2, double R_l2) => (iter <= 4, true);
                        gmRes.Precond = CreateLevelRecursive(iLevel + 1);
                        return gmRes;

                    } else if (iLevel >= 1) {
                        // +++++++++++++++++++++++++++
                        // intermediate level / fine
                        // +++++++++++++++++++++++++++


                        var omg = createOMG();


                        (bool bNotTerminate, bool bSuccess) SubLevelTermination(int iter, double R0_l2, double R_l2) {
                            double reduction = iLevel >= 2 ? 1e-2 : 1e-2;
                            int limit = iLevel >= 2 ? 4 : 100;

                            if (iter <= 1)
                                return (true, false); // keep running

                            if (R_l2 < R0_l2 * reduction) {
                                //Console.WriteLine("succ " + iter + " iterations");
                                return (false, true); // success
                            }


                            if (iter > limit) {
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
        }

        //*/
    }
}
