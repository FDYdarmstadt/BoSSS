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
            var r = CreateInstanceImpl__Kummer(level, level.DGpolynomialDegreeHierarchy);
            //var r = CreateInstanceImpl__BottiDiPietro(level, level.DGpolynomialDegreeHierarchy);

            if (r is IProgrammableTermination pt) {
                pt.TerminationCriterion =  this.DefaultTermination;

            }
            return r;
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


        /// <summary>
        /// Our own approach (see also <see cref="CreateInstanceImpl__BottiDiPietro"/>):
        /// p-Multigrid using <see cref="OrthonormalizationMultigrid"/> as a multigrid and <see cref="CellILU"/> or <see cref="BlockJacobi"/> as smoother
        /// </summary>
        internal ISubsystemSolver CreateInstanceImpl__Kummer(IOperatorMappingPair level, int[][] DegreeHierarchy) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                
                var MPIcomm = level.DgMapping.MPI_Comm;

                //weirdo idea:
                //DegreeHierarchy = new int[][] { new int[] { 5, 5, 4 }, new int[] { 4, 4, 4 }, new int[] { 3, 3, 4 } };


                long DOF_top = level.DgMapping.TotalLength;

                tr.Info("UseILU = " + UseILU);


                ISubsystemSolver CreateLevelRecursive(int iLevel) {

                    int skip = 1;
                    long len = GetDOF(level, DegreeHierarchy[iLevel]);

                    OrthonormalizationMultigrid createOMG() {

                       

                        var coarseSolver = new GridAndDegRestriction() {
                            RestrictedDeg = DegreeHierarchy[iLevel + skip].CloneAs(),
                            LowerPSolver = CreateLevelRecursive(iLevel + 1)
                        };

                        ISolverSmootherTemplate post, pre;

                        ISolverSmootherTemplate[] aps = null;

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
                            
                            if (iLevel == 0) {
                                aps = new ISolverSmootherTemplate[] {
                                    //new VelocityPredictionSolver(),
                                    //new PressureCorrectionSolver()
                                };
                            }
                        } else {
                            post = new BlockJacobi() {
                                NoOfIterations = 1,
                            };
                            pre = null;
                        }

                        tr.Info("OrthoMG (pre = " + (pre?.GetType()?.Name ?? "NULL") + " post = " + (post?.GetType()?.Name ?? "NULL") +  ") solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);



                        var omg = new OrthonormalizationMultigrid() {
                            CoarserLevelSolver = coarseSolver,
                            PreSmoother = pre,
                            PostSmoother = post,
                            AdditionalPostSmoothers = aps
                        };
                        omg.config.NoOfPostSmootherSweeps = 10;
                        omg.config.m_omega = 1;
                        omg.config.CoarseOnLovwerLevel = false;
                        return omg;
                    }

                    //Console.WriteLine("###########################   testcode aktiv ##################################");
                    if ((iLevel >= 2) || (DegreeHierarchy.Length == iLevel + 1)) {
                        // ++++++++++++
                        // lowest level 
                        // ++++++++++++
                        tr.Info("direct solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);

                        var s = new DirectSolver();
                        s.config.TestSolution = false;
                        s.ActivateCaching = (a, b) => true;
                        s.config.UseDoublePrecision = false;
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
                            double reduction = iLevel >= 2 ? 1e-1 : 1e-1;
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

        /// <summary>
        /// Approach according to (Botti, Di Pietro, 2021: p‑Multilevel Preconditioners for HHO Discretizations of the Stokes Equations with Static Condensation):
        /// p-Multigrid using <see cref="ClassicMultigrid"/> and <see cref="SoftGMRES"/> with <see cref="CellILU"/> or <see cref="BlockJacobi"/> as smoother.
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// - the Botti and DiPietro approach uses only ILU, so <see cref="UseILU"/> must be set to true to comply with their findings;
        /// - or own approach (<see cref="CreateInstanceImpl__Kummer"/>) out-performs this configuration;
        ///   E.g. in their 2D test-case, polynomial order 5, 576 cells on one processor:
        ///   - Botti and Di Pietro: 32 iter, 24 sec
        ///   - Kummer: 7 iter, 3.3 sec
        /// </remarks>
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

                        tr.Info("ClassicMG solver on level " + iLevel + ", Degrees " + DegreeHierarchy[iLevel].ToConcatString("[", ",", "]") + ", dofs = " + len);


                        var coarseSolver = new GridAndDegRestriction() {
                            RestrictedDeg = DegreeHierarchy[iLevel + skip].CloneAs(),
                            LowerPSolver = CreateLevelRecursive(iLevel + 1)
                        };

                        ISolverSmootherTemplate post, pre;



                        // GMRES as recommended in Botti/DiPietro paper
                        const int FewIter = 2;
                        var smoother = new SoftGMRES() {
                        };
                        smoother.MaxKrylovDim = FewIter + 2;
                        //var smoother = new OrthonormalizationMultigrid() {
                        //};
                        smoother.TerminationCriterion = delegate (int iter, double R0_l2, double R_l2) {
                            //Console.WriteLine("  gmres: iter = " + iter);
                            var RunSucc = (iter <= FewIter, true);
                            return RunSucc;
                        };
                        if (UseILU) {
                            var _post_pc = new CellILU() {
                                ILU_level = 0
                            };
                            smoother.Precond = _post_pc;
                            //smoother.PreSmoother = _post_pc;
                        } else {
                            var bj = new BlockJacobi() {
                                NoOfIterations = 1,
                            };
                            smoother.Precond = bj;
                            //smoother.PreSmoother = bj;
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
