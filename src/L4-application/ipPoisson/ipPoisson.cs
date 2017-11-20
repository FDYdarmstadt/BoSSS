/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Platform;
using ilPSP;
using System.Linq;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Solution.Queries;
using BoSSS.Foundation.Grid.RefElements;
using NUnit.Framework;

namespace ipPoisson {

    class Program : Application<ippControl> {

#pragma warning disable 649
        /// <summary>
        /// dependent variable
        /// </summary>
        [InstantiateFromControlFile("T", "T", IOListOption.ControlFileDetermined)]
        protected SinglePhaseField T;

        /// <summary>
        /// exact solution, to determine L2-Error, see also <see cref="ippControl.ExactSolution_provided"/>.
        /// </summary>
        [InstantiateFromControlFile("Tex", "Tex", IOListOption.Never)]
        protected SinglePhaseField Tex;

        /// <summary>
        /// dependent variable
        /// </summary>
        [InstantiateFromControlFile("RHS", "T", IOListOption.ControlFileDetermined)]
        protected SinglePhaseField RHS;
#pragma warning restore 649

        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            
            _Main(args, false, "", delegate() {
                Program p = new Program();

                Console.WriteLine("ipPoisson: " + ilPSP.Environment.MPIEnv.MPI_Rank + " of " + ilPSP.Environment.MPIEnv.MPI_Size
                    + " on compute node '" + ilPSP.Environment.MPIEnv.Hostname + "';");
                return p;
            });
        }

        /// <summary>
        /// 1D grid that can be blended between linear- and Cosine-spacing
        /// </summary>
        /// <param name="l"></param>
        /// <param name="r"></param>
        /// <param name="a">blending between linear and cosine - spacing: 1.0 = full cosine,
        /// 0.0 = full linear;</param>
        /// <param name="n">number of nodes</param>
        /// <returns></returns>
        private double[] CosLinSpaceing(double l, double r, double a, int n) {
            double[] linnodes = GenericBlas.Linspace(0, Math.PI * 0.5, n);
            double[] linnodes2 = GenericBlas.Linspace(0, 1, n);
            double[] nodes = new double[n];

            for (int i = 0; i < n; i++)
                nodes[i] = linnodes2[i] * (1 - a) + (1.0 - Math.Cos(linnodes[i])) * a;

            

            for (int i = 0; i < n; i++)
                nodes[i] = nodes[i] * (r - l) + l;
            return nodes;

        }



        IMutableMatrixEx LaplaceMtx;
        double[] LaplaceAffine;

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
            using(FuncTrace tr = new FuncTrace()) {

                
                // assemble system, create matrix
                // ------------------------------

                var volQrSch = new CellQuadratureScheme(true, CellMask.GetFullMask(this.GridData));
                var edgQrSch = new EdgeQuadratureScheme(true, EdgeMask.GetFullMask(this.GridData));

                double D = this.GridData.SpatialDimension;
                double penalty_base = (T.Basis.Degree + 1) * (T.Basis.Degree + D) / D;
                double penalty_factor = base.Control.penalty_poisson;

                {
                    // equation assembly
                    // -----------------
                    tr.Info("creating sparse system...");
                    Console.WriteLine("creating sparse system for {0} DOF's ...", T.Mapping.Ntotal);
                    Stopwatch stw = new Stopwatch();
                    stw.Start();

                    SpatialOperator LapaceIp = new SpatialOperator(1, 1,QuadOrderFunc.SumOfMaxDegrees(), "T", "T");
                    var flux = new ipFlux(penalty_base * base.Control.penalty_poisson, this.GridData.Cells.cj, base.Control);
                    LapaceIp.EquationComponents["T"].Add(flux);


                    LapaceIp.Commit();

#if DEBUG
                    var RefLaplaceMtx = new MsrMatrix(T.Mapping);
#endif                    
                    LaplaceMtx = new BlockMsrMatrix(T.Mapping);
                    LaplaceAffine = new double[T.Mapping.LocalLength];

                    LapaceIp.ComputeMatrixEx(T.Mapping, null, T.Mapping,
                                             LaplaceMtx, LaplaceAffine,
                                             volQuadScheme: volQrSch, edgeQuadScheme: edgQrSch);
#if DEBUG
                    LaplaceAffine.ClearEntries();
                    LapaceIp.ComputeMatrixEx(T.Mapping, null, T.Mapping,
                                             RefLaplaceMtx, LaplaceAffine,
                                             volQuadScheme: volQrSch, edgeQuadScheme: edgQrSch);
                    MsrMatrix ErrMtx = RefLaplaceMtx.CloneAs();
                    ErrMtx.Acc(-1.0, LaplaceMtx);
                    double err = ErrMtx.InfNorm();
                    double infNrm = LaplaceMtx.InfNorm();
                    Console.WriteLine("Matrix comparison error: " + err + ", matrix norm is: " + infNrm);
                    Assert.Less(err, infNrm * 1e-10, "MsrMatrix2 comparison failed.");
#endif
                    //int q = LaplaceMtx._GetTotalNoOfNonZeros();
                    //tr.Info("finished: Number of non-zeros: " + q);
                    stw.Stop();
                    Console.WriteLine("done {0} sec.", stw.Elapsed.TotalSeconds);


                    //double condNo = LaplaceMtx.condest(BatchmodeConnector.Flavor.Octave);
                    //Console.WriteLine("condition number: {0:0.####E-00} ",condNo);
                }
            }
        }

        public void WriteSEMMatrices() {
            int kSEM; // nodes per edge

            if (this.Grid.RefElements[0] is Triangle) {
                kSEM = this.T.Basis.Degree + 1;
            } else if (this.Grid.RefElements[0] is Square) {
                switch (this.T.Basis.Degree) {
                    case 2: kSEM = 2; break;
                    case 3: kSEM = 2; break;
                    case 4: kSEM = 3; break;
                    case 5: kSEM = 3; break;
                    case 6: kSEM = 4; break;
                    default: throw new NotSupportedException();
                }
            } else {
                throw new NotImplementedException();
            }
            
            SpecFemBasis SEM_basis = new SpecFemBasis(this.GridData, kSEM);

            SEM_basis.CellNodes[0].SaveToTextFile("NODES_SEM" + kSEM + ".txt");
            SEM_basis.MassMatrix.SaveToTextFileSparse("MASS_SEM" + kSEM + ".txt");

            var Mod2Nod = SEM_basis.GetModal2NodalOperator(this.T.Basis);
            var Nod2Mod = SEM_basis.GetNodal2ModalOperator(this.T.Basis);

            Mod2Nod.SaveToTextFileSparse("MODAL" + this.T.Basis.Degree + "_TO_NODAL" + kSEM + ".txt");
            Nod2Mod.SaveToTextFileSparse("NODAL" + kSEM + "_TO_MODAL" + this.T.Basis.Degree + ".txt");

            this.LaplaceMtx.SaveToTextFileSparse("OPERATOR" + this.T.Basis.Degree + ".txt");


            {
                var TEST = this.T.CloneAs();
                TEST.Clear();
                TEST.ProjectField((_2D)((x, y) => x*y));

                int J = this.GridData.Cells.NoOfCells;
                int N = SEM_basis.NodesPerCell[0];
                MultidimensionalArray TEST_at_NODES = MultidimensionalArray.Create(J, N);

                TEST.Evaluate(0, J, SEM_basis.CellNodes[0], TEST_at_NODES);

                double[] CompareAtNodes = new double[J * N];
                Mod2Nod.SpMVpara(1.0, TEST.CoordinateVector, 0.0, CompareAtNodes);

                double[] gretchen = TEST_at_NODES.ResizeShallow(J * N).To1DArray();

                double fdist = GenericBlas.L2Dist(gretchen, CompareAtNodes);
                Debug.Assert(fdist < 1.0e-8);

                double[] bak = new double[TEST.Mapping.LocalLength];
                Nod2Mod.SpMVpara(1.0, CompareAtNodes, 0.0, bak);

                double hdist = GenericBlas.L2Dist(bak, TEST.CoordinateVector);
                Debug.Assert(hdist < 1.0e-8);
            }


            var Scatter = SEM_basis.GetNodeScatterMatrix();

            Scatter.SaveToTextFileSparse("SCATTER_SEM" + kSEM + ".txt");


            {
                var SEMTEST = new SpecFemField(SEM_basis);
                var csem = SEMTEST.Coordinates;

                Random r = new Random(666);
                for(int i = 0; i < csem.GetLength(0); i++) {
                    csem[i] = r.NextDouble();
                }


                var DG_TEST0 = new SinglePhaseField(SEMTEST.Basis.ContainingDGBasis);
                var DG_TEST = this.T.CloneAs();
                DG_TEST.Clear();
                SEMTEST.AccToDGField(1.0, DG_TEST0);
                DG_TEST.AccLaidBack(1.0, DG_TEST0);

                double[] S2 = new double[Scatter.RowPartitioning.LocalLength];
                double[] S3 = new double[DG_TEST.Mapping.LocalLength];

                Scatter.SpMVpara(1.0, csem.To1DArray(), 0.0, S2);
                Nod2Mod.SpMVpara(1.0, S2, 0.0, S3);

                double gdist = GenericBlas.L2Dist(S3, DG_TEST.CoordinateVector);
                Debug.Assert(gdist < 1.0e-8);
            }
        }


      
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using(new FuncTrace()) {
                //this.WriteSEMMatrices();
                
                base.NoOfTimesteps = -1;
                if(TimestepNo > 1)
                    throw new ApplicationException("steady-state-equation.");

                base.TerminationKey = true;
                                
                // call solver
                // -----------
                double mintime, maxtime;
                bool converged;
                int NoOfIterations;


                if(base.Control.solver_name == null) {
                    ClassicSolve(out mintime, out maxtime, out converged, out NoOfIterations);
                    
                } else {
                    //ExperimentalSolve(out mintime, out maxtime, out converged, out NoOfIterations);
                    throw new NotImplementedException("todo");
                }

                Console.WriteLine("finished; " + NoOfIterations + " iterations.");
                Console.WriteLine("converged? " + converged);
                Console.WriteLine("Timespan: " + mintime + " to " + maxtime + " seconds");


                base.QueryHandler.ValueQuery("minSolRunT", mintime, true);
                base.QueryHandler.ValueQuery("maxSolRunT", maxtime, true);
                base.QueryHandler.ValueQuery("Conv", converged ? 1.0 : 0.0, true);
                base.QueryHandler.ValueQuery("NoIter", NoOfIterations, true);
                base.QueryHandler.ValueQuery("NoOfCells", this.GridData.CellPartitioning.TotalLength);
                base.QueryHandler.ValueQuery("DOFs", T.Mapping.TotalLength, true);
                base.QueryHandler.ValueQuery("BlockSize", T.Basis.Length, true);


                if(base.Control.ExactSolution_provided) {
                    SinglePhaseField ERR;
                    if(Tex.Basis.Degree >= T.Basis.Degree) {
                        ERR = this.Tex.CloneAs();
                        ERR.AccLaidBack(-1.0, T);
                    } else {
                        ERR = this.T.CloneAs();
                        ERR.AccLaidBack(-1.0, Tex);
                    }

                    double L2_ERR = ERR.L2Norm();

                    base.QueryHandler.ValueQuery("SolL2err", L2_ERR, true);

                }


                return 0.0;
            }
        }

        private void ClassicSolve(out double mintime, out double maxtime, out bool Converged, out int NoOfIter) {

            /*
            double rnk, condNo, det;
            bool isPosDef;
            {
                MultidimensionalArray output = MultidimensionalArray.Create(1,5);
                var bmc = new BatchmodeConnector(_Flav:BatchmodeConnector.Flavor.Octave, MatlabExecuteable: "d:\\cygwin64\\bin\\bash.exe");
                //var bmc = new BatchmodeConnector(_Flav: BatchmodeConnector.Flavor.Matlab);
                bmc.PutSparseMatrix(this.LaplaceMtx, "Mtx");
                bmc.Cmd("condNo = condest(Mtx);");
                bmc.Cmd("rnk = rank(Mtx);");
                //bmc.Cmd("rnk = 0;");
                bmc.Cmd("[V,r]=chol(-0.5*(Mtx+Mtx'));");
                //bmc.Cmd("deti = det(Mtx);");
                bmc.Cmd("output = [condNo,rnk,0.0,4.0,r];");
                bmc.GetMatrix(output, "output");
                bmc.Execute(false);
                condNo = output[0, 0];
                rnk = output[0, 1];
                det = output[0, 2];
                Debug.Assert(output[0, 3] == 4.0);
                isPosDef = output[0, 4] == 0;
            }
            Console.WriteLine("Matrix:");
            Console.WriteLine("  size:      " + this.LaplaceMtx.NoOfRows);
            Console.WriteLine("  rank:      " + rnk);
            Console.WriteLine("  det:       " + det);
            Console.WriteLine("  cond. No.: {0:0.####E-00}", condNo);
            Console.WriteLine("  pos. def.: " + isPosDef);
            */

            
            // create sparse solver
            // --------------------
            ISparseSolver ipSolver;
            {
                ipSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

                //ipSolver = new ilPSP.LinSolvers.monkey.CG();
                //((ilPSP.LinSolvers.monkey.CG)ipSolver).MaxIterations = 50000;
                //((ilPSP.LinSolvers.monkey.CG)ipSolver).Tolerance = 1.0e-10;
                                
                ipSolver.DefineMatrix(LaplaceMtx);
            }

            // call solver
            // -----------

            mintime = double.MaxValue;
            maxtime = double.MinValue;
            SolverResult solRes = default(SolverResult);
            int NoOfSolverRuns = base.Control.NoOfSolverRuns;
            for(int i = 0; i < NoOfSolverRuns; i++) {
                T.Clear();

                Console.WriteLine("RUN " + i + ": solving system...");

                var RHSvec = RHS.CoordinateVector.ToArray();
                //RHSvec.SaveToTextFile("DG" + this.T.Basis.Degree + "_RHS.txt");
                BLAS.daxpy(RHSvec.Length, -1.0, this.LaplaceAffine, 1, RHSvec, 1);

                T.Clear();
                T.InitRandom();
                solRes = ipSolver.Solve(T.CoordinateVector, RHSvec);
                mintime = Math.Min(solRes.RunTime.TotalSeconds, mintime);
                maxtime = Math.Max(solRes.RunTime.TotalSeconds, maxtime);
                System.GC.Collect();

                //T.CoordinatesAsVector.SaveToTextFile("DG" + this.T.Basis.Degree + "_SOLUTION.txt");
            }

            Converged = solRes.Converged;
            NoOfIter = solRes.NoOfIterations;

            ipSolver.Dispose();
        }

        List<DGField> MGColoring = new List<DGField>();

        /*
        private void ExperimentalSolve(out double mintime, out double maxtime, out bool Converged, out int NoOfIter) {
            int p = this.T.Basis.Degree;
            var MgSeq = AggregationGrid.CreateSequence(this.GridData,MaxDepth:2);
            AggregationGridBasis[][] AggBasis = MgSeq.Select(aggGrid => new AggregationGridBasis[] { new AggregationGridBasis(this.T.Basis, aggGrid)}).ToArray();

            Console.WriteLine("Setting up multigrid operator...");
            var mgsetup = new Stopwatch();
            mgsetup.Start();
            var MultigridOp = new MultigridOperator(AggBasis, this.T.Mapping, this.LaplaceMtx, null,
                new MultigridOperator.ChangeOfBasisConfig[][] {
                    new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() { VarIndex = new int[] {0}, mode = MultigridOperator.Mode.DiagBlockEquilib, Degree = p }
                    }
                });
            mgsetup.Stop();
            Console.WriteLine("done. (" + mgsetup.Elapsed + ")");

            Stopwatch stw = new Stopwatch();

            ISolverSmootherTemplate solver;
            {
                string solverName = base.Control.solver_name.ToLower();
                switch(solverName) {
                    case "direct":
                        solver = new DirectSolver() {
                            TestSolution = true
                        };
                        break;

                    case "softpcg+schwarz+directcoarse":
                        solver = new SoftPCG() {
                            m_MaxIterations = 50000,
                            m_Tolerance = 1.0e-10,
                            Precond = new Schwarz() {
                                m_MaxIterations = 1,
                                CoarseSolver = new GenericRestriction() {
                                    CoarserLevelSolver = new DirectSolver()
                                },
                                m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                    Depth = 2,
                                },
                                overlap = 2
                            }
                        };
                        break;

                    case "softpcg+schwarz":
                        solver = new SoftPCG() {
                            m_MaxIterations = 50000,
                            m_Tolerance = 1.0e-10,
                            Precond = new Schwarz() {
                                m_MaxIterations = 1,
                                CoarseSolver = null,
                                m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                    Depth = 2,
                                },
                                overlap = 2
                            }
                        };
                        break;

                    case "softpcg":
                        solver = new SoftPCG() {
                            m_MaxIterations = 50000,
                            m_Tolerance = 1.0e-10
                            //Precond = new BlockJacobi() { NoOfIterations = 1, omega = 1 }
                        };
                        break;


                    case "softpcg+direct":
                        solver = new SoftPCG() {
                            m_MaxIterations = 5000,
                            m_Tolerance = 1.0e-10,
                            Precond = new DirectSolver()
                        };
                        break;

                    //case "MultigridPCG":
                    //    solver = SoftPCG.InitMultigridChain(MultigridMap,
                    //        (i, pcg) => {
                    //            if(i == 0) {
                    //                pcg.m_MaxIterations = 50000;
                    //                pcg.m_Tolerance = 1.0e-10;
                    //            } else {
                    //                pcg.m_MaxIterations = 1;
                    //                pcg.m_MinIterations = 1;
                    //                pcg.m_Tolerance = 0;
                    //            }
                    //        },
                    //        () => new DirectSolver());
                    //    break;

                    case "softpcg+mg1":
                        solver = new SoftPCG() {
                            m_MaxIterations = 2000,
                            m_Tolerance = 1.0e-10,
                            Precond = ClassicMultigrid.InitMultigridChain(MultigridOp,
                                i => new BlockJacobi() { NoOfIterations = 5, omega = 0.6 },
                                i => new BlockJacobi() { NoOfIterations = 5, omega = 0.6 },
                                (i, mg) => {
                                    mg.Gamma = 1;
                                    mg.m_MaxIterations = 1;
                                },
                                () => new DirectSolver())
                        };
                        break;
                    case "softpcg+mg2":
                        solver = new SoftPCG() {
                            m_MaxIterations = 2000,
                            m_Tolerance = 1.0e-10,
                            Precond = ClassicMultigrid.InitMultigridChain(MultigridOp,
                                i => new Schwarz() {
                                    m_MaxIterations = 1,
                                    CoarseSolver = null,
                                    m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                        Depth = 2,
                                    },
                                    overlap = 0
                                },
                                i => new Schwarz() {
                                    m_MaxIterations = 1,
                                    CoarseSolver = null,
                                    m_BlockingStrategy = new Schwarz.MultigridBlocks() {
                                        Depth = 2,
                                    },
                                    overlap = 0
                                },
                                (i, mg) => {
                                    mg.Gamma = 1;
                                    mg.m_MaxIterations = 1;
                                },
                                () => new DirectSolver())
                        };
                        break;

         

                    case "softpcg+blockjacobi":

                        solver = new SoftPCG() {
                            m_MaxIterations = 50000,
                            m_Tolerance = 1.0e-10,
                            Precond = new Jacobi() { NoOfIterations = 1, omega = 1 }
                        };
                        break;



                    case "classicmg":
                        solver = ClassicMultigrid.InitMultigridChain(MultigridOp,
                            i => new BlockJacobi() { NoOfIterations = 20, omega = 0.4 },
                            i => new BlockJacobi() { NoOfIterations = 20, omega = 0.4 },
                            (i, mg) => {
                                mg.m_MaxIterations = i == 0 ? 10000 : 1;
                                mg.m_Tolerance = i == 0 ? 1.0e-10 : 1;
                            },
                            () => new DirectSolver());

                        break;

                    default:
                        throw new ApplicationException("unknown solver: " + solverName);




                }

                if(solver is ISolverWithCallback) {
                    ((ISolverWithCallback)solver).IterationCallback = delegate(int iter, double[] xI, double[] rI, MultigridOperator mgOp) {
                        double l2_RES = rI.L2NormPow2().MPISum().Sqrt();
                        double l2_ERR =  GenericBlas.L2DistPow2(xI, T.CoordinateVector).MPISum().Sqrt();
                        Console.WriteLine("Iter: {0}\tRes: {1:0.##E-00}\tErr: {2:0.##E-00}\tRunt: {3:0.##E-00}", iter, l2_RES, l2_ERR, stw.Elapsed.TotalSeconds);
                        //Tjac.CoordinatesAsVector.SetV(xI);
                        //Residual.CoordinatesAsVector.SetV(rI);
                        //PlotCurrentState(iter, new TimestepNumber(iter), 3);
                    }; 
                }
            }


            solver.Init(MultigridOp);
            //Tjac.Acc(1.0, Tcg);
            mintime = double.MaxValue;
            maxtime = double.MinValue;
            double[] T2 = this.T.CoordinateVector.ToArray();

            for(int irun = 0; irun < base.Control.NoOfSolverRuns; irun++) {
                solver.ResetStat();
                T2.Clear();
                var RHSvec = RHS.CoordinateVector.ToArray();
                BLAS.daxpy(RHSvec.Length, -1.0, this.LaplaceAffine, 1, RHSvec, 1);

                stw.Reset();
                stw.Start();
                // solver.Solve(T2, RHSvec);
                MultigridOp.UseSolver(solver, T2, RHSvec);
                stw.Stop();
                
                mintime = Math.Min(stw.Elapsed.TotalSeconds, mintime);
                maxtime = Math.Max(stw.Elapsed.TotalSeconds, maxtime);

                // Console.WriteLine("Number of iterations: {0}, runtime: {1} sec.", solver.NoOfIterations, stw.Elapsed.TotalSeconds);
            }

            T.CoordinateVector.SetV(T2);

            Converged = solver.Converged;
            NoOfIter = solver.ThisLevelIterations;
        }
        */

        protected override void Bye() {
            object SolL2err;
            if (this.QueryHandler.QueryResults.TryGetValue("SolL2err", out SolL2err)) {
                Console.WriteLine("Value of Query 'SolL2err' " + SolL2err.ToString());
            } else {
                Console.WriteLine("query 'SolL2err' not found.");
            }
        }



        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(new DGField[] { T, Tex, RHS }, "poisson_grid" + timestepNo, phystime, 0);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(ArrayTools.Cat(new DGField[] { T, Tex, RHS }, this.MGColoring), "poisson" + timestepNo, phystime, superSampling);
        }

    }

    /// <summary>
    /// Interior Penalty Flux
    /// </summary>
    class ipFlux : BoSSS.Solution.NSECommon.ipLaplace {

        public ipFlux(double penalty_const, MultidimensionalArray cj, ippControl __ctrl)
            : base(penalty_const, cj, "T") //
        {
            ctrl = __ctrl;
        }


        ippControl ctrl;

        protected override double g_Diri(ref CommonParamsBnd inp) {
            double v = ctrl.g_Diri(inp);
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            return ctrl.g_Neum(inp);
        }
        
        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return ctrl.IsDirichlet(inp);
        }
    }


}
