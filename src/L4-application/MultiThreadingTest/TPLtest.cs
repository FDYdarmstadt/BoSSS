using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiThreadingTest {
    internal class TPLtest {

        static List<(MsrMatrix M, double[] X, double[] B)> GenerateProblem(int n, int count) {
            List<(MsrMatrix M, double[] X, double[] B)> ret = new List<(MsrMatrix M, double[] X, double[] B)>();

            for(int k = 0; k < count; k++) {
                var M = new MsrMatrix(new Partitioning(n, csMPI.Raw._COMM.SELF), new Partitioning(n, csMPI.Raw._COMM.SELF));
                var X = new double[n];
                var B = new double[n];
                for(int i = 0; i < n; i++) {
                    M[i, i] = 4.0;

                    if(i > 0)
                        M[i, i - 1] = -1.0;
                    if(i < n - 1)
                        M[i, i + 1] = -1.0;

                    B[i] = Math.Sin((i + 1.0) * Math.PI / (n + 1)) + 0.1 * Math.Cos((k + 1.0) * i / (n + 1));
                }
                ret.Add((M, X, B));
            }
            return ret;
        }

        static void ParInit(List<(MsrMatrix M, double[] X, double[] B)> problem, PARDISOSolver[] solvers, int i, Parallelism parallelism) {
            using(var tr = new FuncTrace("Init")) {
                var slv_iBlk = new PARDISOSolver() {
                    CacheFactorization = true,
                    UseDoublePrecision = false,
                    Parallelism = parallelism
                };
                slv_iBlk.DefineMatrix(problem[i].M);
                solvers[i] = slv_iBlk;
            }
        }

        static void ParSolve(List<(MsrMatrix M, double[] X, double[] B)> problem, PARDISOSolver[] solvers, int i) {
            Console.WriteLine($"PRank {System.Environment.GetEnvironmentVariable("SLURM_PROCID")} - Task {i} on CPU {MultiThreadingTestMain.sched_getcpu()}");
            using(var tr = new FuncTrace("Solve")) {
                solvers[i].Solve(problem[i].X, problem[i].B);
            }
        }

        public static void ParTestSeq(int n, int count) {
            var Problem = GenerateProblem(n, count);
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                tr.Info("############# PARDISO with Seq in TPL Test #####################");
                 tr.Info("# Foreach");
                PARDISOSolver[] Solvers2 = new PARDISOSolver[count];
                TestThreading4((i) => ParInit(Problem, Solvers2, i, Parallelism.SEQ), 0, count);
                TestThreading4((i) => ParSolve(Problem, Solvers2, i), 0, count);
                CleanCache();

                tr.Info("# BoSSS");
                PARDISOSolver[] Solvers = new PARDISOSolver[count];
                TestThreading((i) => ParInit(Problem, Solvers, i, Parallelism.SEQ), 0, count);
                TestThreading((i) => ParSolve(Problem, Solvers, i), 0, count);
                CleanCache();

                tr.Info("# Special");
                PARDISOSolver[] Solvers3 = new PARDISOSolver[count];
                TestThreading3((i) => ParInit(Problem, Solvers3, i, Parallelism.SEQ), 0, count);
                TestThreading3((i) => ParSolve(Problem, Solvers3, i), 0, count);
                CleanCache();
                tr.Info("############# Finish ####################");
            }
        }

        public static void ParTestOMP(int n, int count) {
            var Problem = GenerateProblem(n, count);
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                tr.Info("############# PARDISO OpenMP Test #####################");
                tr.Info("# OpenMP");
                tr.InfoToConsole = true;
                tr.Info($"Rank {ilPSP.Environment.MPIEnv.MPI_Rank}: ProcessorCount = {System.Environment.ProcessorCount} and NumThreads = {ilPSP.Environment.NumThreads}");
                var sw = Stopwatch.StartNew();
                PARDISOSolver[] Solvers3 = new PARDISOSolver[count];
                for(int i = 0; i < count; i++) {
                    ParInit(Problem, Solvers3, i, Parallelism.OMP);
                }
                sw.Stop();
                tr.Info($"Elapsed: {sw.ElapsedMilliseconds} ms");
                sw.Restart();
                for(int i = 0; i < count; i++) {
                    ParSolve(Problem, Solvers3, i);
                }
                tr.Info($"Elapsed: {sw.ElapsedMilliseconds} ms");
                CleanCache();

                tr.Info("############# Finish ####################");
            }
        }

        static public void TestThreading3(Action<int> task, int i0 = 0, int iE = 10000) {
            task ??= SimpleTestTask;
            using(var tr = new FuncTrace("TestThreading")) {
                tr.InfoToConsole = true;
                tr.Info($"Rank {ilPSP.Environment.MPIEnv.MPI_Rank}: ProcessorCount = {System.Environment.ProcessorCount} and NumThreads = {ilPSP.Environment.NumThreads}");
                var sw = Stopwatch.StartNew();
                ilPSP.Environment.ParallelForWithPartitioner(i0, iE, ilPSP.Environment.NumThreads, (ThreadInfo _, int i) => { task(i); });
                sw.Stop();
                tr.Info($"Elapsed: {sw.ElapsedMilliseconds} ms");
            }
        }

        static public void TestThreading4(Action<int> task, int i0 = 0, int iE = 10000) {
            task ??= SimpleTestTask;
            using(var tr = new FuncTrace("TestThreading")) {
                tr.InfoToConsole = true;
                tr.Info($"Rank {ilPSP.Environment.MPIEnv.MPI_Rank}: ProcessorCount = {System.Environment.ProcessorCount} and NumThreads = {ilPSP.Environment.NumThreads}");

                var sw = Stopwatch.StartNew();

                // Create data range
                var data = Enumerable.Range(i0, iE - i0);

                // Use your ParallelForEach
                ilPSP.Environment.ParallelForEach(data, (ti, i) =>
                {
                    task(i);
                });

                sw.Stop();
                tr.Info($"Elapsed: {sw.ElapsedMilliseconds} ms");
            }
        }

        static public void TestThreading2(Action<int> task, int i0 = 0, int iE = 10000) {
            task ??= SimpleTestTask;
            using(var tr = new FuncTrace("TestThreading")) {
                tr.InfoToConsole = true;
                tr.Info($"Rank {ilPSP.Environment.MPIEnv.MPI_Rank}: ProcessorCount = {System.Environment.ProcessorCount} and NumThreads = {ilPSP.Environment.NumThreads}");
                var sw = Stopwatch.StartNew();
                Parallel.For(i0, iE, new ParallelOptions {
                    MaxDegreeOfParallelism = ilPSP.Environment.NumThreads,
                }, (i) => { task(i); });
                sw.Stop();
                tr.Info($"Elapsed: {sw.ElapsedMilliseconds} ms");
            }
        }

        static public void TestThreading(Action<int> task, int i0 = 0, int iE = 10000) {
            task ??= SimpleTestTask;
            using(var tr = new FuncTrace("TestThreading")) {
                tr.InfoToConsole = true;
                tr.Info($"Rank {ilPSP.Environment.MPIEnv.MPI_Rank}: ProcessorCount = {System.Environment.ProcessorCount} and NumThreads = {ilPSP.Environment.NumThreads}");
                var sw = Stopwatch.StartNew();
                ilPSP.Environment.ParallelFor(i0, iE, (i) => { task(i); });
                sw.Stop();
                tr.Info($"Elapsed: {sw.ElapsedMilliseconds} ms");
            }
        }

        public static void SimpleTestTask(int i) {
            double x = 0;
            for(int j = 0; j < 1000; j++)
                x += Math.Sqrt(i + j);
        }

        static void CleanCache() {
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();
            Thread.Sleep(1000);
        }


        static void SimpleTest() {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                tr.Info("############# Simple Test #####################");
                tr.Info("# Special");
                TestThreading3(SimpleTestTask);
                CleanCache();

                tr.Info("# TPL");
                TestThreading2(SimpleTestTask);
                CleanCache();

                tr.Info("# BoSSS");
                TestThreading(SimpleTestTask);
                CleanCache();
                tr.Info("############# Finish ####################");
            }
        }
    }
}
