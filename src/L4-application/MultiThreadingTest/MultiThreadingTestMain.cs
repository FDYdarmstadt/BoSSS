using System.Runtime.InteropServices;

namespace MultiThreadingTest {

    public class MultiThreadingTestMain {

        [DllImport("libc")]
        public static extern int sched_getcpu();

        static int Rank => ilPSP.Environment.MPIEnv.MPI_Rank;
        static int NumThreads = ilPSP.Environment.NumThreads;

        static void Main(string[] args) {
            Console.WriteLine("Given input parameters: ");
            foreach(var arg in args) {
                Console.Write(arg);
            }

            BoSSS.Solution.Application.InitMPI(args);
            
            //string logFile = $"stdout_rank{ilPSP.Environment.MPIEnv.MPI_Rank}.txt";
            //Console.SetOut(new StreamWriter(logFile) { AutoFlush = true });

            Console.WriteLine("Starting multithread tests");
            Parallel.For(0, System.Environment.ProcessorCount - 1, new ParallelOptions { MaxDegreeOfParallelism = System.Environment.ProcessorCount - 1 }, i => {
                Console.WriteLine($"PRank {System.Environment.GetEnvironmentVariable("SLURM_PROCID")} - Task {i} on CPU {sched_getcpu()}");
            });
            //TPLtest.SimpleTestTask(1000);

            TPLtest.ParTestSeq(1000, 1000);
            SwitchToTPL(false);
            TPLtest.ParTestOMP(1000, 1000);
            SwitchToTPL(true);
            TPLtest.ParTestSeq(1000, 1000);

            MPI.Wrappers.csMPI.Raw.mpiFinalize();

        }

        public static void SwitchToTPL(bool enableTPL) {
            if(enableTPL) {
                // --- Switch to TPL (C#) mode ---
                ilPSP.Environment.DisableOpenMP();
                Environment.SetEnvironmentVariable("OMP_NUM_THREADS", "1");
                Environment.SetEnvironmentVariable("MKL_NUM_THREADS", "1");
                Environment.SetEnvironmentVariable("MKL_THREADING_LAYER", "SEQUENTIAL");
                Environment.SetEnvironmentVariable("OMP_PROC_BIND", "false");
                Environment.SetEnvironmentVariable("KMP_AFFINITY", "disabled");

                //// Use all logical cores assigned by SLURM for TPL
                //int threads = int.Parse(Environment.GetEnvironmentVariable("SLURM_CPUS_PER_TASK") ?? "1");
                ThreadPool.SetMinThreads(NumThreads, NumThreads);
                ThreadPool.SetMaxThreads(NumThreads, NumThreads);

                Console.WriteLine($"[Rank {Rank}] Switched to TPL mode using {NumThreads} threads.");
            } else {
                // --- Switch to MKL / OpenMP mode ---
                ilPSP.Environment.EnableOpenMP();
                Environment.SetEnvironmentVariable("OMP_NUM_THREADS", NumThreads.ToString());
                Environment.SetEnvironmentVariable("MKL_NUM_THREADS", NumThreads.ToString());
                Environment.SetEnvironmentVariable("MKL_THREADING_LAYER", "INTEL");
                Environment.SetEnvironmentVariable("OMP_PROC_BIND", "true");
                Environment.SetEnvironmentVariable("OMP_PLACES", "cores");
                Environment.SetEnvironmentVariable("KMP_AFFINITY", "granularity=fine,compact,1,0");

                // Limit TPL interference
                ThreadPool.SetMinThreads(1, 1);
                ThreadPool.SetMaxThreads(1, 1);

                Console.WriteLine($"[Rank {Rank}] Switched to MKL/OpenMP mode using {NumThreads} threads.");
            }
        }

    }

}