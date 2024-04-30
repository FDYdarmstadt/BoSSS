using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading;

namespace ilPSP {

    /// <summary>
    /// 
    /// </summary>
    [Serializable]
    [DataContract]
    public class OnlinePerformanceLog {

        /// <summary>
        /// Collected normalized benchmark results over the application runtime
        /// - Key: benchmark name
        /// - Value: normalized results, i.e. 1.0 represents expected behavior, smaller values denote under-performance, larger values denote over-performance
        /// </summary>
        [DataMember]
        public Dictionary<string, List<double>> BenchResults;


        [DataMember]
        public OMPbindingStrategy OMPbindingStrategy;


        public void WriteStatistics(TextWriter tw) {
            foreach(var br in BenchResults) {
                tw.Write(br.Key);
                tw.Write(": ");
                double[] res = br.Value.ToArray();
                double avg = res.Sum()/res.Count();
                double min = res.Min();
                double max = res.Max();
                tw.Write($"[{min:g5} | {avg:g5} | {max:g5}]");
                tw.Write("  ");
                tw.Write(res.ToConcatString("[", ", ", "]", "g4"));
                tw.WriteLine();
            }
        }

        public override string ToString() {
            using(var tw = new StringWriter()) {
                WriteStatistics(tw);
                return tw.ToString();
            }
        }

    }


    /// <summary>
    /// Online (i.e., executed while BoSSS performs in a production run)
    /// evaluation of micro-benchmarks, to detect any problems with 
    /// computational performance.
    /// </summary>
    public static class OnlinePerformanceMeasurement {

        public const int MaxNumOfBenchmarks = 10;

        static int ExeCount; 

        public static void ExecuteBenchmarks() {
            if (ExeCount >= MaxNumOfBenchmarks)
                return;
            ExeCount++;
                
            using (var tr = new FuncTrace("ExecuteBenchmarks")) {
                //tr.InfoToConsole = true;
                foreach (var b in AllBenchmarks) {
                    if (Log.BenchResults == null)
                        Log.BenchResults = new Dictionary<string, List<double>>();

                    if (!Log.BenchResults.TryGetValue(b.Key, out var results)) {
                        results = new List<double>();
                        Log.BenchResults.Add(b.Key, results);
                    }

                    var result = b.Value.BenchmarkEval(b.Value.Benchmark);
                    tr.Info($"Benchmark '{b.Key}' result: {result}");
                    results.Add(result);
                }

            }
        }

        /// <summary>
        /// results collected during Application runtime
        /// </summary>
        public static OnlinePerformanceLog Log = new OnlinePerformanceLog();


        public static Dictionary<string, (Ibench Benchmark, Func<Ibench, double> BenchmarkEval)> AllBenchmarks;


        /// <summary>
        /// Runtime for the 1024*1024 DGEMM measured on some reference machine
        /// </summary>
        const double DGEMM1024_serialBaselineRuntime = 0.009044966666666666;


        static OnlinePerformanceMeasurement() {
            AllBenchmarks = new Dictionary<string, (Ibench, Func<Ibench, double>)> {
                { "Serial-GEMMbsline", (new GEMMbench(1024, 5), bench => CompareAgainstBaselineValSerial(bench, DGEMM1024_serialBaselineRuntime)) },
                { "OpenMP-GEMMbsline", (new GEMMbench(1024, 5), bench => CompareAgainstBaselineValParallel(bench, DGEMM1024_serialBaselineRuntime/ilPSP.Environment.NumThreads)) },
                { "OpenMP-GEMMaccel", (new GEMMbench(1024, 5), MeasureAcceleration) },
                { "OpenMP-GEMMblock", (new GEMMbench(1024, 5), CheckThreadBlocking) },
                { "TPL-TPLaccel", (new TPLbench(1024, 5), MeasureAcceleration) },
                { "TPL-TPLblock", (new TPLbench(1024, 5), CheckThreadBlocking) }
            };

        }



      
        public interface Ibench {

            string Name { get; }

            (double minTime, double avgTime, double maxTime) DoParallel();
            (double minTime, double avgTime, double maxTime) DoSerial();
        }

        /// <summary>
        /// Benchmarking of OpenMP/BLAS employing double precision, general Matrix-Matrix-Multiplication (DGEMM)
        /// </summary>
        class GEMMbench : Ibench {

            public GEMMbench(int N, int Runs) {
                this.Runs = Runs;

                A = MultidimensionalArray.Create(N, N);
                B = MultidimensionalArray.Create(N, N);
                C = MultidimensionalArray.Create(N, N);

                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        A[i, j] = Math.Sin(i + j);
                        B[i, j] = Math.Cos(i + j);
                        C[i, j] = A[i, j]*B[i, j];

                    }
                }
            }

            public string Name => "GEMM";

            int Runs;

            MultidimensionalArray A, B, C;

            public (double minTime, double avgTime, double maxTime) DoParallel() {
                if (ilPSP.Environment.InParallelSection)
                    throw new ApplicationException("No nested calls allowed");

                return Do();
            }
            public (double minTime, double avgTime, double maxTime) DoSerial() {
                if (ilPSP.Environment.InParallelSection)
                    throw new ApplicationException("No nested calls allowed");
                BLAS.ActivateSEQ();
                var r = Do();
                BLAS.ActivateOMP();
                return r;
            }

            (double minTime, double avgTime, double maxTime) Do() {
                double mintime = double.MaxValue;
                double maxtime = 0.0;

                var RunTimes = new System.Collections.Generic.List<double>();

                for (int i = 0; i < Runs; i++) {
                    var start = DateTime.Now;
                    A.GEMM(1.0, B, C, 0.1);
                    var end = DateTime.Now;

                    double secs = (end - start).TotalSeconds;

                    RunTimes.Add(secs);
                    mintime = Math.Min(mintime, secs);
                    maxtime = Math.Max(maxtime, secs);
                }

                // remove the outliers before computing the average:
                RunTimes.Sort();
                RunTimes.RemoveAt(0);
                RunTimes.RemoveAt(RunTimes.Count - 1);
                double avgTime = RunTimes.Sum() / RunTimes.Count;

                avgTime /= Runs;

                return (mintime, avgTime, maxtime);
            }
        }


        /// <summary>
        /// Benchmarking of Thread Parallel Library: this uses the C# random number generator as benchmark load,
        /// therefore memory bandwidth should not play a significant role and 
        /// therefore the scaling should be quite ideal
        /// </summary>
        class TPLbench : Ibench {

            public TPLbench(int N, int Runs) {
                this.Runs = Runs;

                //A = MultidimensionalArray.Create(N, N);
                N2 = N*N;

                
            }

            public string Name => "TPL";

            int Runs;

            //MultidimensionalArray A;
            int N2;

            public (double minTime, double avgTime, double maxTime) DoParallel() {
                if (ilPSP.Environment.InParallelSection)
                    throw new ApplicationException("No nested calls allowed");
                return Do(true);
            }

            public (double minTime, double avgTime, double maxTime) DoSerial() {
                if (ilPSP.Environment.InParallelSection)
                    throw new ApplicationException("No nested calls allowed");
                return Do(false);
            }


            /// <summary>
            /// Fills an vector with random entries
            /// </summary>
            double BenchOp(bool par) {
                int L = N2;

                double[] globAcc = new double[ilPSP.Environment.NumThreads];
                ilPSP.Environment.ParallelFor(0, L, delegate (int iThread, int i0, int iE) {
                    double locAcc = 0;
                    //Random rnd = new Random();
                    for (int i = i0; i < iE; i++) {
                        locAcc += Math.Sin(i);// + rnd.NextDouble();
                    }
                    globAcc[iThread] = locAcc;
                }, enablePar:par);

                return globAcc.Sum();
            }

            (double minTime, double avgTime, double maxTime) Do(bool par) {
                double mintime = double.MaxValue;
                double maxtime = 0.0;

                var RunTimes = new System.Collections.Generic.List<double>();

                for (int i = 0; i < Runs; i++) {
                    var start = DateTime.Now;
                    BenchOp(par);
                    var end = DateTime.Now;

                    double secs = (end - start).TotalSeconds;

                    RunTimes.Add(secs);
                    mintime = Math.Min(mintime, secs);
                    maxtime = Math.Max(maxtime, secs);
                }

                // remove the outliers before computing the average:
                RunTimes.Sort();
                RunTimes.RemoveAt(0);
                RunTimes.RemoveAt(RunTimes.Count - 1);
                double avgTime = RunTimes.Sum() / RunTimes.Count;

                avgTime /= Runs;

                return (mintime, avgTime, maxtime);
            }
        }



        /// <summary>
        /// Compares the runtime of serial Matrix-Matrix-Multiplication with parallel one.
        /// </summary>
        static double MeasureAcceleration(Ibench gb) {
            using (var tr = new FuncTrace("MeasureAcceleration")) {
                //tr.InfoToConsole = false; 
                if (Environment.OpenMPenabled == false || Environment.NumThreads <= 1) {
                    return -1;
                }

                int NumThreads = Environment.NumThreads;



                
                var SerialTimes = gb.DoSerial();
                tr.Info($"serial benchmark {gb.Name} : {SerialTimes.avgTime} sec");


                //MKLservice.SetNumThreads(NumThreads);
                //Console.WriteLine("OMP Max threads = " + MKLservice.GetMaxThreads() + ", BoSSS numth = " + NumThreads + ", BoSSS max omp = " + MaxNumOpenMPthreads + ", omp numth = " + MKLservice.omp_get_num_threads());
                //Console.WriteLine("MKL_dynamic (s) " + MKLservice.Dynamic);
                //SetOMPbinding();
                ////MKLservice.Dynamic = true;
                //MKLservice.omp_set_num_threads(NumThreads);
                //Console.WriteLine("     ..... omp numth = " + MKLservice.omp_get_num_threads());

                //MKLservice.SetNumThreads(NumThreads);
                //MKLservice.om
                //Console.WriteLine("MKL_dynamic (r) " + MKLservice.Dynamic);

                /*
                for (int i = 0; i < 1000; i++) {
                    Console.Write(i + ".");
                    gb.Do();
                }
                Console.WriteLine();
                */

                var ParallelTimes = gb.DoParallel();
                tr.Info($"parallel benchmark {gb.Name} : {ParallelTimes.avgTime} sec");

                double Accel = SerialTimes.avgTime / ParallelTimes.avgTime; // kleine parallele Laufzeit == gut => große `Accel`
                double RelAccel = Accel / NumThreads;
                tr.Info($"Parallel acceleration of {gb.Name} using {NumThreads}: {Accel}, relative factor {RelAccel}");

                
                //Tracer.Current.UpdateTime();
                //MethodCallRecordExtension.GetMostExpensiveCalls(Console.Error, Tracer.Root);
                //Console.Error.WriteLine(Tracer.Root.TimeExclusive);
                return RelAccel;
            }
        }


        static double CompareAgainstBaselineValSerial(Ibench b, double BaselineRuntime) {
            using (var tr = new FuncTrace("CompareAgainstBaselineValSerial")) {
                //tr.InfoToConsole = true;
                var SerialTimes = b.DoSerial();



                double RelFactor = BaselineRuntime / SerialTimes.avgTime; // große Laufzeit => Ergebnis < 1 => Underperfomer; kleine Laufzeit => Ergebnis > 1 => Overperformer
                tr.Info($"Reference machine comparison (serial) of {b.Name}: relative factor {RelFactor}; abs. runtime: {SerialTimes.avgTime}");

                return RelFactor;
            }
        }
        static double CompareAgainstBaselineValParallel(Ibench b, double BaselineRuntime) {
            using (var tr = new FuncTrace("CompareAgainstBaselineValParallel")) {
                if (Environment.OpenMPenabled == false || Environment.NumThreads <= 1) {
                    return -1;
                }



                var ParallelTimes = b.DoParallel();

                double RelFactor = BaselineRuntime / ParallelTimes.avgTime; // große Laufzeit => Ergebnis < 1 => Underperfomer; kleine Laufzeit => Ergebnis > 1 => Overperformer
                tr.Info($"Reference machine comparison ({Environment.NumThreads} threads) of {b.Name}: relative factor {RelFactor}; abs. runtime: {ParallelTimes.avgTime}");

                return RelFactor;
            }
        }

        /// <summary>
        /// loops through all strategies in <see cref="OMPbindingStrategy"/> and tries to select the fastest one
        /// </summary>
        static public OMPbindingStrategy FindBestOMPstrategy(int[] CPUIndices, out bool OMPisSlow) {
            OMPbindingStrategy[] strats = (OMPbindingStrategy[]) Enum.GetValues(typeof(OMPbindingStrategy));
            double[] performance = new double[strats.Length];

            var bench = new GEMMbench(2048, 5);

            for (int i = 0; i < strats.Length; i++) {
                double measure() {
                    double accel = MeasureAcceleration(bench).MPIMin();
                    double block = CheckThreadBlocking(bench);
                    Console.WriteLine($"strat {strats[i]} [{accel:g4} | {block:g4}]");

                    return Math.Min(accel, block);
                }


                Console.WriteLine("-----------------------   Performance of strategy " +  strats[i] + ": ");
                MKLservice.BindOMPthreads(CPUIndices, strats[i]);
                performance[i] = measure();
                Console.WriteLine("                    =====");
            }

            int IBest = performance.IndexOfMax();
            Console.WriteLine("selecting " + strats[IBest]);
            OMPisSlow = performance[IBest]*CPUIndices.Length <= 0.8;
            if (OMPisSlow) {
                Console.WriteLine("Insufficient OpenMP acceleration // should be disabled!");
            }
            return strats[IBest];
        }




        /// <summary>
        /// We are trying to identify if external libraries from different MPI ranks dead-lock each other;
        /// Therefore:
        /// 1. A reference measurement of a GEMM operation is performed on rank 0
        /// 2. The same operation is then performed on rank [0], [0,1], [0,1,2], ...
        /// 3. The runtime measurements are then compared to the reference measurement
        /// 4. an exception is thrown if the parallel runs take much longer than the reference run
        /// </summary>
        /// <returns>
        /// A "blocking factor": describes how much slower an OpenMP operation becomes if it is performed by multiple threads.
        /// The ideal factor is 1, i.e. no slow-down if multiple threads are active
        /// </returns>
        static double CheckThreadBlocking(Ibench b) {
            using (var tr = new FuncTrace("CheckThreadBlocking")) {
                /* SOME TESTS ON LICHTENBERG: 
                 * 
                 * there seems to be no clear advantage in setting OMP_PLACES;
                 * 
                 * without any OMP_PLACES: -------------------------------------------

                01/25/2024 13:57:37  Running with 4 MPI processes 
                Working path: /work/home/fk69umer/XNSE-25jan24
                Binary path: /work/home/fk69umer/XNSE-25jan24
                Current commit hash: 1ffd4df84aec90b393ee393933b8853b88c85b7a
                User: fk69umer
                Node: mpsc0536 (ranks 0, 1, 2, 3)

                MPI Rank 2: Value for OMP_PLACES: 
                MPI Rank 3: Value for OMP_PLACES: 
                MPI Rank 3: Value for OMP_PROC_BIND: 
                MPI Rank 1: Value for OMP_PLACES: 
                MPI Rank 2: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: Value for OMP_PLACES: 
                MPI Rank 1: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 1: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 2: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 3: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                Ref run: (min|avg|max) : (	7.879E-02 |	8.43E-02 |	1.328E-01)
                Now, doing work on 1 ranks ...
                R0: 1 workers: (min|avg|max) : (	7.878E-02 |	7.9E-02   |	7.924E-02)  --- 		( 1E00     |	9.371E-01 |	5.97E-01)
                Now, doing work on 2 ranks ...
                R0: 2 workers: (min|avg|max) : (	7.579E-02 |	7.853E-02 |	8.01E-02)  --- 		    ( 9.62E-01 |	9.316E-01 |	6.03E-01)
                R1: 2 workers: (min|avg|max) : (	3.419E-02 |	8.897E-02 |	4.743E-01)  --- 		( 4.34E-01 |	1.055E00  |	3.57E00)
                Now, doing work on 3 ranks ...
                R1: 3 workers: (min|avg|max) : (	5.337E-02 |	7.226E-02 |	8.145E-02)  --- 		( 6.77E-01 |	8.572E-01 |	6.14E-01)
                R0: 3 workers: (min|avg|max) : (	7.728E-02 |	8.364E-02 |	9.387E-02)  --- 		( 9.81E-01 |	9.922E-01 |	7.07E-01)
                R2: 3 workers: (min|avg|max) : (	7.862E-02 |	8.826E-02 |	1.322E-01)  --- 		( 9.98E-01 |	1.047E00  |	9.96E-01)
                Now, doing work on 4 ranks ...
                R2: 4 workers: (min|avg|max) : (	7.767E-02 |	8.378E-02 |	9.692E-02)  --- 		( 9.86E-01 |	9.939E-01 |	7.3E-01)
                R0: 4 workers: (min|avg|max) : (	8.063E-02 |	1.092E-01 |	1.83E-01)  --- 		    ( 1.02E00  |	1.296E00  |	1.38E00)
                R1: 4 workers: (min|avg|max) : (	6.777E-02 |	1.142E-01 |	3.435E-01)  --- 		( 8.6E-01  |	1.354E00  |	2.59E00)
                R3: 4 workers: (min|avg|max) : (	3.407E-02 |	1.279E-01 |	6.426E-01)  --- 		( 4.32E-01 |	1.517E00  |	4.84E00)
                            -------------------------------------

                and now, with OMP_PLACES set: -----------------------------------
                01/25/2024 14:05:01  Running with 4 MPI processes 
                Working path: /work/home/fk69umer/XNSE-25jan24
                Binary path: /work/home/fk69umer/XNSE-25jan24
                Current commit hash: 1ffd4df84aec90b393ee393933b8853b88c85b7a
                User: fk69umer
                Node: mpsd0114 (ranks 0, 1, 2, 3)


                MPI Rank 3: Value for OMP_PLACES: 
                MPI Rank 1: Value for OMP_PLACES: 
                MPI Rank 2: Value for OMP_PLACES: 
                MPI Rank 3: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 1: Value for OMP_PROC_BIND: 
                MPI Rank 2: Value for OMP_PROC_BIND: 
                MPI Rank 0: Value for OMP_PLACES: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 2: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 0: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 3: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 1: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                R0, SMP rank 0: setting OMP_PLACES = {0,1,2,3}
                R1, SMP rank 1: setting OMP_PLACES = {4,5,6,7}
                R3, SMP rank 3: setting OMP_PLACES = {12,13,14,15}
                R2, SMP rank 2: setting OMP_PLACES = {8,9,10,11}
                Ref run: (min|avg|max) : (	3.31E-02 |	5.071E-02 |	1.917E-01)
                Now, doing work on 1 ranks ...
                R0: 1 workers: (min|avg|max) : (	3.297E-02 |	3.804E-02 |	5.443E-02)  --- 		( 9.96E-01 |	7.502E-01 |	2.84E-01)
                Now, doing work on 2 ranks ...
                R0: 2 workers: (min|avg|max) : (	3.623E-02 |	4.84E-02 |	6.778E-02)  --- 		( 1.09E00 |	9.545E-01 |	3.54E-01)
                Now, doing work on 3 ranks ...
                R1: 2 workers: (min|avg|max) : (	4.857E-02 |	7.535E-02 |	2.967E-01)  --- 		( 1.47E00 |	1.486E00 |	1.55E00)
                R1: 3 workers: (min|avg|max) : (	5.121E-02 |	5.61E-02 |	7.05E-02)  --- 		( 1.55E00 |	1.106E00 |	3.68E-01)
                R0: 3 workers: (min|avg|max) : (	4.167E-02 |	6.363E-02 |	7.679E-02)  --- 		( 1.26E00 |	1.255E00 |	4E-01)
                R2: 3 workers: (min|avg|max) : (	4.825E-02 |	7.676E-02 |	2.685E-01)  --- 		( 1.46E00 |	1.514E00 |	1.4E00)
                Now, doing work on 4 ranks ...
                R1: 4 workers: (min|avg|max) : (	5.134E-02 |	5.854E-02 |	9.231E-02)  --- 		( 1.55E00 |	1.154E00 |	4.81E-01)
                R2: 4 workers: (min|avg|max) : (	5.128E-02 |	6.41E-02 |	8.576E-02)  --- 		( 1.55E00 |	1.264E00 |	4.47E-01)
                R3: 4 workers: (min|avg|max) : (	5.178E-02 |	8.39E-02 |	3.045E-01)  --- 		( 1.56E00 |	1.655E00 |	1.59E00)
                R0: 4 workers: (min|avg|max) : (	3.749E-02 |	8.64E-02 |	1.526E-01)  --- 		( 1.13E00 |	1.704E00 |	7.96E-01)
                -------------------------------------

                */





                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSz);
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int Rank);


                (double minTime, double avgTime, double maxTime) TimeRef0 = (0, 0, 0);
                if (Rank == 0) {
                    tr.Info("Now, doing reference run on rank 0...");
                    TimeRef0 = b.DoParallel();
                    tr.Info($"Ref run: (min|avg|max) : (\t{TimeRef0.minTime:0.###E-00} |\t{TimeRef0.avgTime:0.###E-00} |\t{TimeRef0.maxTime:0.###E-00})");
                }

                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                var TimeRef = TimeRef0.MPIBroadcast(0);

                double worstAvgFactor = 0;
                for (int ranksToBench = 0; ranksToBench < MpiSz; ranksToBench++) {
                    tr.Info("Now, doing work on " + (ranksToBench + 1) + " ranks ...");

                    (double minTime, double avgTime, double maxTime) TimeX = (BLAS.MachineEps, BLAS.MachineEps, BLAS.MachineEps);
                    if (Rank <= ranksToBench) {
                        TimeX = b.DoParallel();
                    }


                    double minFactor = TimeX.minTime / TimeRef.minTime;
                    double avgFactor = TimeX.avgTime / TimeRef.avgTime; // TimeX klein == gut/schnell => Factor < 1; Muss also noch invertiert werden
                    double maxFactor = TimeX.maxTime / TimeRef.maxTime;
                    worstAvgFactor = Math.Max(worstAvgFactor, avgFactor);

                    if (minFactor.MPIMax() > 10 || maxFactor.MPIMax() > 5 || avgFactor.MPIMax() > 10) {

                        string scaling = $"R{Rank}: {ranksToBench + 1} workers: (min|avg|max) : (\t{TimeX.minTime:0.###E-00} |\t{TimeX.avgTime:0.###E-00} |\t{TimeX.maxTime:0.###E-00})  --- \t\t( {minFactor:0.##E-00} |\t{avgFactor:0.###E-00} |\t{maxFactor:0.##E-00})";
                        tr.Info("Scaling involving " + (ranksToBench + 1) + " ranks: " + scaling);
                        
                        

                        //if(avgFactor > 7)
                        //    throw new ApplicationException("Some very slow processor detected -- maybe some OpenMP locking: " + scaling);
                    }



                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                }

                return 1/worstAvgFactor; // inversion, for return value: small factor = under-performance, large factor = over-performance;
            }
        }

    }
}
