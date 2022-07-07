using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
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
            return CreateInstanceImpl__Kummer(level, level.DGpolynomialDegreeHierarchy);
        }

        /// <summary>
        /// - true: <see cref="CellILU"/> is used as smoother on finest p-level
        /// - false: <see cref="BlockJacobi"/> is used as smoother on finest p-level
        /// </summary>
        public bool UseILU = true;


        static private long GetDOF(IOperatorMappingPair level, int[] DgDegrees) {
            int NoOfVar = level.DgMapping.NoOfVariables;
            if (DgDegrees.Length != NoOfVar)
                throw new ArgumentException();

            if (DgDegrees.ListEquals(level.DgMapping.DgDegree))
                return level.DgMapping.TotalLength;

            int GetNp(int p) {
                int SpacDim = level.DgMapping.SpatialDimension;
                switch (SpacDim) {
                    case 1: return p + 1;
                    case 2: return (p * p + 3 * p + 2) / 2;
                    case 3: return (p * p * p + 6 * p * p + 11 * p + 6) / 6;
                    default: throw new ArgumentException("unknown spatial dimension: " + SpacDim);
                }
            }


            long sum = 0;
            int J = level.DgMapping.NoOfLocalUpdatedCells;
            for(int j = 0; j < J; j++) {
                int NoOfSpc = level.DgMapping.GetNoOfSpecies(j);

                for (int iVar = 0; iVar < NoOfVar; iVar++)
                    sum += GetNp(DgDegrees[iVar]) * NoOfSpc;
            }
            return sum;
        }



        internal ISubsystemSolver CreateInstanceImpl__Kummer(IOperatorMappingPair level, int[][] DegreeHierarchy) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                //var precond = new LevelPmg();
                //precond.config.UseHiOrderSmoothing = false;
                //precond.config.OrderOfCoarseSystem = this.pMaxOfCoarseSolver;
                //precond.config.FullSolveOfCutcells = this.FullSolveOfCutcells;
                //precond.config.UseDiagonalPmg = this.UseDiagonalPmg;
                //precond.config.EqualOrder = this.EqualOrder;

                var MPIcomm = level.DgMapping.MPI_Comm;


                long DOF_top = level.DgMapping.TotalLength;

                tr.Info("UseILU = " + UseILU);


                ISubsystemSolver CreateLevelRecursive(int iLevel) {

                    int skip = 1;
                    long len = GetDOF(level, DegreeHierarchy[iLevel]);

                    OrthonormalizationMultigrid createOMG() {

                        tr.Info("OrthoMG solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);


                        var coarseSolver = new GridAndDegRestriction() {
                            RestrictedDeg = DegreeHierarchy[iLevel + skip].CloneAs(),
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

                        if (iLevel <= 2 && UseILU) {
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

                    //Console.WriteLine("###########################   testcode aktiv ##################################");
                    if (iLevel >= 2/*DegreeHierarchy.Length - 1*/) {

                        tr.Info("direct solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);

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

                        tr.Info("GMRES solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);

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
        
        
        internal ISubsystemSolver CreateInstanceImpl__BottiDiPietro(IOperatorMappingPair level, int[][] DegreeHierarchy) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;

                var MPIcomm = level.DgMapping.MPI_Comm;

                long DOF_top = level.DgMapping.TotalLength;

                tr.Info("UseILU = " + UseILU);


                ISubsystemSolver CreateLevelRecursive(int iLevel) {

                    int skip = 1;
                    long len = GetDOF(level, DegreeHierarchy[iLevel]);

                    ClassicMultigrid createMG() {

                        tr.Info("OrthoMG solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);


                        var coarseSolver = new GridAndDegRestriction() {
                            RestrictedDeg = DegreeHierarchy[iLevel + skip].CloneAs(),
                            LowerPSolver = CreateLevelRecursive(iLevel + 1)
                        };

                        ISolverSmootherTemplate post, pre;

                        

                        //GMRES bringt nix als PC
                        var smoother = new SoftGMRES() {
                        //    //Precond = bj
                        };
                        smoother.TerminationCriterion = (int iter, double R0_l2, double R_l2) => (iter <= 1, true);
                        if (UseILU) {
                            var _post_pc = new CellILU() {
                                ILU_level = 0
                            };
                            smoother.Precond = _post_pc;
                        } else {
                            var bj = new BlockJacobi() {
                                NoOfIterations = 1,
                            };
                            smoother.Precond = bj;
                        }

                        pre = smoother;
                        post = smoother;


                        var mg = new ClassicMultigrid() {
                            CoarserLevelSolver = coarseSolver,
                            PreSmoother = pre,
                            PostSmoother = post,
                        };
                        mg.config.m_omega = 1;
                        mg.config.CoarseOnLovwerLevel = false;
                        
                        return mg;
                    }

                    //Console.WriteLine("###########################   testcode aktiv ##################################");
                    if (iLevel >= 2/*DegreeHierarchy.Length - 1*/) {

                        tr.Info("direct solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);

                        // ++++++++++++
                        // lowest level 
                        // ++++++++++++
                        var s = new DirectSolver();
                        s.config.TestSolution = false;
                        s.ActivateCaching = (a, b) => true;
                        return s;
                    } else {
                        // +++++++++++++++++++++++++++
                        // intermediate level / fine
                        // +++++++++++++++++++++++++++

                        tr.Info("V-cycle on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);

                        var omg = createMG();


                        return omg;

                    
                    }


                }


                var topLevel_precond = CreateLevelRecursive(0);

                var topLevel = new OrthonormalizationMultigrid() {
                    PreSmoother = topLevel_precond
                };
                                
                topLevel.Init(level);
                return topLevel;
            }
        }


        //*/
    }
}
