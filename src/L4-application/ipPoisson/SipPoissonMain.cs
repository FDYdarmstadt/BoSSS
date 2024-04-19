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
using NUnit.Framework;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.Control;
using BoSSS.Solution.Statistic;
using System.IO;
using BoSSS.Platform.Utils.Geom;
using System.Threading;
using ilPSP.LinSolvers.PARDISO;
using System.Collections;
using System.Runtime.InteropServices;
using System.ComponentModel;

namespace BoSSS.Application.SipPoisson {


    static class  BoSSSmkl {
        
    }


    /// <summary>
    /// Benchmark application, solves a Poisson problem using the symmetric interior penalty (SIP) method.
    /// </summary>
    public class SipPoissonMain : Application<SipControl> {


        
        /// <summary>
        /// Main routine
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {



            /*

            Debugger.Launch();
            InitMPI(args);
            //BoSSS.Application.SipPoisson.Tests.TestProgram.TestOperatorConvergence3D(1);

            
            int MpiSz;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MpiSz);
            int Rank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out Rank);

            int Nothreads = ilPSP.Environment.NumThreads;
             

            ilPSP.Environment.StdoutOnlyOnRank0 = false;

            


            //Process proc = Process.GetCurrentProcess();
            var CPUidxs = ilPSP.Utils.CPUAffinity.GetAffinity();
            Console.WriteLine("r" + Rank + $" Windows Affinity to {CPUidxs.Count()} cores: {CPUidxs.ToConcatString("", ", ", ";")}");

            //proc.ProcessorAffinity = (IntPtr) 0xFFFFFFFFFFFF;
            //Console.WriteLine("r" + Rank + $"  Reset Affinity {proc.ProcessorAffinity:X}");
                        


            ilPSP.Environment.StdoutOnlyOnRank0 = true;


           

            //SetOMPAffinity(Nothreads, Rank, MpiSz, false);
            //ilPSP.Environment.InitThreading(false, Nothreads);



            int N = 5000;
            var A = Nothreads.ForLoop(ith => MultidimensionalArray.Create(N, N));
            var B = Nothreads.ForLoop(ith => MultidimensionalArray.Create(N, N));
            var C = Nothreads.ForLoop(ith => MultidimensionalArray.Create(N, N));

            ilPSP.Environment.ParallelFor(0, Nothreads, delegate (int ith) {
                A[ith].Storage.FillRandom();
                B[ith].Storage.FillRandom();
                C[ith].Storage.FillRandom();
            });

            bool ompMkl = args.Contains("ompmkl");
            bool tplMkl = args.Contains("tplmkl");
            bool tpl = args.Contains("tpl");

            if ((ompMkl || tpl || tplMkl) == false)
                throw new Exception("you must turn something on");

            var start = DateTime.Now;
            var timeout = new TimeSpan(hours: 0, minutes: 50, seconds: 1);
            Stopwatch s0 = new Stopwatch();
            int cnt = -1;
            while (true) {
                var duration = DateTime.Now - start;
                if (duration > timeout) {
                    //FinalizeMPI();
                    //return;
                    break;
                }

                {
                    cnt = (int)Math.Round(duration.TotalSeconds / 15);
                    if (cnt < 0)
                        cnt = 0;
                    Console.WriteLine("count = " + cnt);

                    int[] CPUs;
                    if (duration.TotalSeconds < 1*60) {
                        CPUs = new int[] { 4, 5, 6, 7, 8, 9, 10, 11 };
                        //CPUs = 40.ForLoop(i => i + 4);
                    } else if (duration.TotalSeconds < 2*60) {
                        CPUs = new int[] { 8, 9, 10, 11, 12, 13, 14, 15 };
                        //CPUs = 40.ForLoop(i => i + 4 + 64);
                    } else if (duration.TotalSeconds < 3*60) {
                        CPUs = new int[] { 0, 1, 2, 3, 4, 5, 6, 7};
                        //CPUs = ArrayTools.Cat(42.ForLoop(i => i + 2), 42.ForLoop(i => i + 64 + 2));
                    } else {
                        CPUs = 4.ForLoop(i => i*2 + 4);
                        //CPUs = new int[] { 64, 65, 66, 70, 71, 72 };
                    }
                    Console.WriteLine($"Binding to {CPUs.Length} CPUs/cores");
                    ilPSP.MKLservice.BindOMPthreads(CPUs);
                }


                s0.Reset();
                s0.Start();
                var AnrmLoc = new double[A.Length];
                var BnrmLoc = new double[B.Length];
                var CnrmLoc = new double[C.Length];

                ilPSP.Environment.ParallelFor(0, Nothreads, delegate (int iThread) {
                    AnrmLoc[iThread] = A[iThread].L2Norm();
                    BnrmLoc[iThread] = B[iThread].L2Norm();
                    CnrmLoc[iThread] = C[iThread].L2Norm();
                });

                s0.Stop();
                Console.WriteLine("   dnrm time: " + s0.Elapsed.TotalSeconds);

                s0.Reset();
                s0.Start();
                var AnrmG = AnrmLoc.MPISum().Sum();
                var BnrmG = BnrmLoc.MPISum().Sum();
                var CnrmG = CnrmLoc.MPISum().Sum();
                
                Console.WriteLine("A norm: " + AnrmG);
                Console.WriteLine("B norm: " + BnrmG);
                Console.WriteLine("C norm: " + CnrmG);
                s0.Stop();
                Console.WriteLine("   MPIsum time: " + s0.Elapsed.TotalSeconds);

                if (ompMkl && Rank == 0) {
                    Console.WriteLine($"   OMP MKL DGEMM {N}x{N}");
                    s0.Reset();
                    var st2 = DateTime.Now;
                    s0.Start();
                    A[0].GEMM(1.0, B[0], C[0], 0.1);
                    s0.Stop();
                    var en2 = DateTime.Now;
                    Console.WriteLine($"   MKL time: {s0.Elapsed.TotalSeconds:0.##} ( {(en2 - st2).TotalSeconds:0.##}");
                }
                if (tplMkl) {
                    Console.WriteLine($"   TPL MKL: {Nothreads}x Serial MKL-DGEMM {N}x{N}");
                    s0.Reset();
                    var st3 = DateTime.Now;
                    s0.Start(); 
                    ilPSP.Environment.ParallelFor(0, Nothreads, delegate (int iThread) {
                        A[iThread].GEMM(1.0, B[iThread], C[iThread], 1.0);
                        A[iThread].Storage.ScaleV(0.01);


                    });
                    s0.Stop();
                    var en3 = DateTime.Now;
                    Console.WriteLine($"   TPL time: {s0.Elapsed.TotalSeconds:0.##} ( {(en3 - st3).TotalSeconds:0.##}");
                }


                if (tpl) {
                    Console.WriteLine($"   TPL only: {Nothreads}x Serial Naive DGEMM {N}x{N}");
                    s0.Reset();
                    s0.Start();
                    ilPSP.Environment.ParallelFor(0, Nothreads, delegate (int iThread) {
                        //A[iThread].GEMM(1.0, B[iThread], C[iThread], 1.0);
                        //A[iThread].Storage.ScaleV(0.01);

                        var _A = A[iThread];
                        var _B = B[iThread];
                        var _C = C[iThread];


                        for (int iRow = 0; iRow < N; iRow++) {
                            for (int iCol = 0; iCol < N; iCol++) {
                                double Acc = _A[iRow, iCol]*0.1;
                                for (int k = 0; k < N; k++)
                                    Acc += _B[iRow, k]*_C[k, iCol];
                                _A[iRow, iCol] = Acc;
                            }
                        }
                        
                    });
                    s0.Stop();
                    Console.WriteLine("   TPL time: " + s0.Elapsed.TotalSeconds);
                }

            }

            //*/
            Debugger.Launch();
            _Main(args, false, delegate () {
                SipPoissonMain p = new SipPoissonMain();
                
                return p;
            });
        }

       

#pragma warning disable 649
        /// <summary>
        /// dependent variable
        /// </summary>
        [InstantiateFromControlFile("T", "T", IOListOption.ControlFileDetermined)]
        protected SinglePhaseField T;

        /// <summary>
        /// exact solution, to determine L2-Error, see also <see cref="SipControl.ExactSolution_provided"/>.
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
        /// Solver residual, for a DG discretization which is one degree higher than the degree of the solution
        /// </summary>
        private SinglePhaseField ResiualKP1;

        /// <summary>
        /// error of the numerical solution
        /// </summary>
        private SinglePhaseField Error;

        /// <summary>
        /// MPI rank coloring
        /// </summary>
        private SinglePhaseField MPIrank;

        /// <summary>
        /// DG field instantiation
        /// </summary>
        protected override void CreateFields() {
            base.CreateFields();

            ResiualKP1 = new SinglePhaseField(new Basis(this.GridData, T.Basis.Degree + 1), "ResidualKP1");
            base.IOFields.Add(ResiualKP1);

            Error = new SinglePhaseField(new Basis(this.GridData, Math.Max(T.Basis.Degree + 1, Tex.Basis.Degree)), "Error");
            base.m_IOFields.Add(Error);

            // MPI rank coloring
            MPIrank = new SinglePhaseField(new Basis(this.GridData, 0), "MPIRank");
            MPIrank.AccConstant(this.MPIRank);

            // mg coloring
            int iLevel = 0;
            this.MGColoring.Clear();
            foreach (var MgL in this.MultigridSequence) {
                SinglePhaseField c = new SinglePhaseField(new Basis(this.GridData, 0), "MgLevel_" + iLevel);
                Foundation.Grid.Aggregation.CoarseningAlgorithms.ColorDGField(MgL, c);
                this.MGColoring.Add(c);
                base.IOFields.Add(c);
                iLevel++;
            }

            Console.WriteLine("Available multi-grid levels: " + this.MGColoring.Count);
        }

        /*
        unsafe static void my_dgemm(int TRANSA, int TRANSB,
                                        int M, int N, int K,
                                        double ALPHA,
                                        double* A, int LDA,
                                        double* B, int LDB,
                                        double BETA,
                                        double* C, int LDC) {
            for(int m = 0; m < M; m++) {
                for(int n = 0; n < N; n++) {
                    double acc = 0;
                    for(int k = 0; k < K; k++) {
                        acc += A[m * K + k] * B[k * N + n];
                    }
                    C[m * N + n] = BETA * C[m * N + n] + ALPHA * acc;
                }
            }

        }
        */


        /*
#if !DEBUG
        static void MyHandler(object sender, UnhandledExceptionEventArgs args) {
            Exception e = (Exception)args.ExceptionObject;
            Console.WriteLine("MyHandler caught : " + e.Message);
            Console.WriteLine("Runtime terminating: {0}", args.IsTerminating);
            System.Environment.Exit(-1234);
        }
#endif
*/
        /// <summary>
        /// Ensures availability of <see cref="BoSSS.Solution.Statistic.ForeignGridValue"/>
        /// </summary>
        public Type EnsureReference = typeof(ForeignGridValue);


       

        /// <summary>
        /// Sets the multigrid coloring
        /// </summary>
        protected override void SetInitial(double t) {
            /*
#if !DEBUG
            //this will suppress exception prompts
            //Workaround to prevent disturbance while executing batch-client
            if (this.Control.SuppressExceptionPrompt) {
                AppDomain currentDomain = AppDomain.CurrentDomain;
                currentDomain.UnhandledException += new UnhandledExceptionEventHandler(MyHandler);
            }
#endif
            */
            base.SetInitial(t);

            

            //TexactFine = (SinglePhaseField)(GetDatabase().Sessions.First().Timesteps.Last().Fields.Where(fi => fi.Identification == "T"));
        }

        /// <summary>
        /// Spatial operator used by <see cref="UniSolver.Solve"/>
        /// </summary>
        DifferentialOperator LapaceIp;

        /// <summary>
        /// Includes assembly of the matrix.
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            using(FuncTrace tr = new FuncTrace()) {

                // create operator
                // ===============
                BoundaryCondMap<BoundaryType> PoissonBcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "T");
                LapaceIp = new DifferentialOperator(1, 1, QuadOrderFunc.SumOfMaxDegrees(), "T", "T");
                var flux = new ipFlux(base.Control.penalty_poisson, PoissonBcMap);
                LapaceIp.EquationComponents["T"].Add(flux);
                LapaceIp.EquationComponents["T"].Add(new RHSSource(this.RHS));
                LapaceIp.IsLinear = true;
                LapaceIp.Commit();
            }
        }
     

        /// <summary>
        /// control of mesh adaptation
        /// </summary>
        protected override void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            if (this.Control.AdaptiveMeshRefinement && TimestepNo > 1) {

                // compute error against fine solution
                if (Control.ExactSolution_provided) {
                    //Error.Clear();
                    //Error.AccLaidBack(1.0, T);

                    /*
                    var eval = new FieldEvaluation((GridData)(TexactFine.GridDat));

                    void FineEval(MultidimensionalArray input, MultidimensionalArray output) {
                        int L = input.GetLength(0);
                        Debug.Assert(output.GetLength(0) == L);

                        eval.Evaluate(1.0, new DGField[] { TexactFine }, input, 0.0, output.ResizeShallow(L, 1));
                    }

                    Error.ProjectField(-1.0, FineEval);
                    */
                    //Error.AccLaidBack(-1.0, Tex);
                }

                long oldJ = this.GridData.CellPartitioning.TotalLength;

                double LocNormPow2 = this.ResiualKP1.CoordinateVector.L2NormPow2(); // norm of residual on this processor
                double TotNormPow2 = LocNormPow2.MPISum(); //                          norm of residual over all processors
                double MeanNormPow2PerCell = TotNormPow2 / oldJ; //                    mean norm per cell

                double maxSoFar = 0;
                int jMax = -1;
                for (int j = 0; j < oldJ; j++) {
                    double CellNorm = Error.Coordinates.GetRow(j).L2NormPow2();

                    if (CellNorm > maxSoFar) {
                        jMax = j;
                        maxSoFar = CellNorm;
                    }
                }


                int MyLevelIndicator(int j, int CurrentLevel) {
                    double CellNorm = this.ResiualKP1.Coordinates.GetRow(j).L2NormPow2();


                    //if (j == 0)
                    //    CurrentLevel = CurrentLevel + 1;

                    //if (CellNorm > MeanNormPow2PerCell * 1.1)
                    //    return CurrentLevel + 1;
                    //else
                    //    return CurrentLevel;
                    if (j == jMax)
                        return CurrentLevel + 1;
                    else
                        return CurrentLevel;
                }

                GridRefinementController gridRefinementController = new GridRefinementController((GridData)this.GridData,null);
                bool AnyChange = gridRefinementController.ComputeGridChange(MyLevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening);
                int NoOfCellsToRefine = 0;
                int NoOfCellsToCoarsen = 0;
                if (AnyChange) {
                    int[] glb = (new int[] {
                        CellsToRefineList.Count,
                        Coarsening.Sum(L => L.Length),
                        //0, 0
                    }).MPISum();

                    NoOfCellsToRefine = glb[0];
                    NoOfCellsToCoarsen = glb[1];
                }
                //*/


                // Update Grid
                // ===========

                if (AnyChange) {


                    Console.WriteLine("       Refining " + NoOfCellsToRefine + " of " + oldJ + " cells");
                    Console.WriteLine("       Coarsening " + NoOfCellsToCoarsen + " of " + oldJ + " cells");

                    newGrid = ((GridData)(this.GridData)).Adapt(CellsToRefineList, Coarsening, out old2NewGrid);

                } else {

                    newGrid = null;
                    old2NewGrid = null;
                }
            } else {

                newGrid = null;
                old2NewGrid = null;
            }
        }

        /// <summary>
        /// Single run of the solver
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {

                int L = 1000;
                double[] a = new double[L];
                int[] a_acc = L.ForLoop(i => i);
                double[] b = new double[L];
                int[] b_acc = L.ForLoop(i => i);
                a.AccV(1.0, b, acc_index: a_acc, b_index: b_acc);



                if (Control.ExactSolution_provided) {
                    Tex.Clear();
                    Tex.ProjectField(this.Control.InitialValues_Evaluators["Tex"]);

                    RHS.Clear();
                    RHS.ProjectField(this.Control.InitialValues_Evaluators["RHS"]);
                }

                if (Control.AdaptiveMeshRefinement == false) {
                    base.NoOfTimesteps = -1;
                    if (TimestepNo > 1)
                        throw new ApplicationException("steady-state-equation.");
                    base.TerminationKey = true;
                }


                // call solver
                // -----------
                //LastMatrix = this.LapaceIp.GetMatrix(T.Mapping, MgConfig: this.MgConfig);
                //Console.WriteLine("Remember to re-activate solver !!!!!!!");
                this.LapaceIp.Solve(T.Mapping, MgConfig: this.MgConfig, lsc: this.Control.LinearSolver, verbose: true, queryHandler: base.QueryHandler);

                //long J = this.GridData.CellPartitioning.TotalLength;
                //LastMatrix.SaveToTextFileSparse($"LaplaceMtx-J{J}.txt");
                //double condNo = LastMatrix.condest();
                //Console.WriteLine($"Matlab condition number estimate {J} cells: " + condNo);

                /*Stuff for testing OpenMP
                {
                    var mgOp = this.LapaceIp.GetMultigridOperator(T.Mapping, MgConfig: this.MgConfig);
                    var rhs = new double[T.Mapping.LocalLength];
                    rhs.FillRandom();


                    var slvrs = new List<(int num,PARDISOSolver s)>();


                    foreach (int numThreads in new int[] { 8, 2, 4, 8, 16, 4 }) {
                        var _start = DateTime.Now;
                        Console.WriteLine("First solve with " + numThreads + " threads...");
                        var mtx = mgOp.OperatorMatrix;

                        //System.Environment.SetEnvironmentVariable("OMP_NUM_THREADS", numThreads.ToString());
                        //System.Environment.SetEnvironmentVariable("MKL_NUM_THREADS", numThreads.ToString());
                        //System.Environment.SetEnvironmentVariable("MKL_DYNAMIC", "false");

                        ThreadBLAS.OMP_SET_NUM_THREADS(numThreads);

                        var pardiso = new PARDISOSolver();
                        pardiso.Parallelism = Parallelism.OMP;
                        pardiso.CacheFactorization = false;
                        pardiso.DefineMatrix(mtx);
                        pardiso.Solve(new double[rhs.Length], rhs);


                        slvrs.Add((numThreads, pardiso));

                        Console.WriteLine($"done. (took {DateTime.Now - _start})");
                    }


                    var start = DateTime.Now;
                    foreach(var kv in slvrs) {
                        Console.WriteLine($"Now solving with {kv.num} threads:");
                        ThreadBLAS.OMP_SET_NUM_THREADS(kv.num);

                        int cnt = 0;
                        while(DateTime.Now - start < new TimeSpan(hours:0, minutes:0, seconds:59)) {
                            Console.Write($"Do run #{cnt} ... ");
                            kv.s.Solve(new double[rhs.Length], rhs);
                            Console.WriteLine($" done. {(DateTime.Now - start)} passed.");
                        }

                        start = DateTime.Now;
                    }




                    foreach (var s in slvrs)
                        s.s.Dispose();
                }
                */



                if (base.Control.ExactSolution_provided) {
                    Error.Clear();
                    Error.AccLaidBack(1.0, Tex);
                    Error.AccLaidBack(-1.0, T);

                    double L2_ERR = Error.L2Norm();
                    Console.WriteLine("\t\tL2 error on " + this.Grid.NumberOfCells + " cells: " + L2_ERR);
                    last_L2_ERR = L2_ERR;
                    base.QueryHandler.ValueQuery("SolL2err", L2_ERR, true);

                }

                // evaluate residual
                // =================
                {
                    //this.ResiualKP1.Clear();
                    //if (this.Control.InitialValues_Evaluators.ContainsKey("RHS")) {
                    //    this.ResiualKP1.ProjectField(this.Control.InitialValues_Evaluators["RHS"]);
                    //}



                    var ev = this.LapaceIp.GetEvaluatorEx(T.Mapping, new[] { RHS }, ResiualKP1.Mapping);
                    ev.Evaluate(-1.0, 1.0, ResiualKP1.CoordinateVector);
                }

                // return
                // ======

                return 0.0;
            }
        }

        internal double last_L2_ERR;

        List<DGField> MGColoring = new List<DGField>();


        MultigridOperator.Mode m_MgConfig = MultigridOperator.Mode.DiagBlockEquilib;

        MultigridOperator.ChangeOfBasisConfig[][] MgConfig {
            get {
                int p = this.T.Basis.Degree;
                int NoOfLevels = this.MultigridSequence.Length;
                var config = new MultigridOperator.ChangeOfBasisConfig[NoOfLevels][];

                for (int iLevel = 0; iLevel < NoOfLevels; iLevel++) {

                    config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() {
                            VarIndex = new int[] {0},
                            mode = m_MgConfig,
                            DegreeS = new int[] { p }
                            //Degree = Math.Max(1, p - iLevel)
                        }
                    };

                }

                return config;
            }

        }

        /// <summary>
        /// Shutdown function
        /// </summary>
        protected override void Bye() {
            object SolL2err;
            if (this.QueryHandler != null) {
                if (this.QueryHandler.QueryResults.TryGetValue("SolL2err", out SolL2err)) {
                    Console.WriteLine("Value of Query 'SolL2err' " + SolL2err.ToString());
                } else {
                    Console.WriteLine("query 'SolL2err' not found.");
                }
            }
        }

        /// <summary>
        /// Operator stability analysis
        /// </summary>
        override public IDictionary<string,double> OperatorAnalysis(OperatorAnalysisConfig config) {
            using(new FuncTrace()) {
                return this.LapaceIp.OperatorAnalysis(this.T.Mapping, config, this.MgConfig); 
            }
        }

        /// <summary>
        /// default plotting
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0) {
            string caseStr = "";
            if (base.Control.Paramstudy_CaseIdentification != null) {
                var pstudy_case = base.Control.Paramstudy_CaseIdentification.FirstOrDefault(tt => tt.Item1 == "pstudy_case");
                if (pstudy_case != null) {
                    caseStr = "." + pstudy_case.Item2;
                }
            }

            DGField[] Fields = new DGField[] { T, Tex, RHS, ResiualKP1, Error };
            Fields = Fields.Cat(this.MGColoring);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "poisson_MG_coloring" + timestepNo + caseStr, phystime, superSampling);
            //BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, Path.Combine(AnalyseOutputpath, "poisson_MG_coloring" + timestepNo + caseStr), phystime, superSampling);
        }


       
    }

    /// <summary>
    /// Interior Penalty Flux
    /// </summary>
    class ipFlux : BoSSS.Solution.NSECommon.SIPLaplace {

        public ipFlux(double penalty_const, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, "T") //
        {
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap.bndFunction["T"];
        }

        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        protected override double g_Diri(ref CommonParamsBnd inp) {
            double v = m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            double v = m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            return v;
        }

        //public override double Nu(double[] x, double[] p, int jCell) {
        //    return 0.001;
        //}

        //public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        //public override TermActivationFlags InnerEdgeTerms => base.InnerEdgeTerms;

        //public override TermActivationFlags VolTerms => base.VolTerms;


        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType) {
                case BoundaryType.Dirichlet:
                    return true;

                case BoundaryType.Neumann:
                    return false;

                default:
                    throw new NotImplementedException();
            }
        }
    }

    /// <summary>
    /// source term on the RHS
    /// </summary>
    class RHSSource : IVolumeForm, IParameterHandling {

        public RHSSource(DGField rhsSourceField) {
            m_rhsSourceField = rhsSourceField;
        }


        DGField m_rhsSourceField;

        public TermActivationFlags VolTerms => TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new string[0];

        public IList<string> ParameterOrdering => new[] { "RHSsource" };

        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            return new[] { m_rhsSourceField };
        }

        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {
            if(!object.ReferenceEquals(m_rhsSourceField,Parameters[0])) {
                Parameters[0].Clear();
                Parameters[0].Acc(1.0, m_rhsSourceField);
            }
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double rhsVal = cpv.Parameters[0];
            return -rhsVal * V; // actually, the RHS is added on the left, therefore minus!
        }
    }


}
