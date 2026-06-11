using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver {

    [TestFixture]
    public static class XNSE_Solver_MultithreadingTest {

        [Test]
        public static void Perf() {
            int threads = GetThreadCountFromEnvironment();

            Console.WriteLine($"Running XNSE thread performance test with {threads} threads.");
            Console.WriteLine($"OMP_NUM_THREADS      = {Environment.GetEnvironmentVariable("OMP_NUM_THREADS")}");

            string file = $"xnse-thread-{threads}.time";
            EnsureTimingHeader(file);

            foreach(string testcase in GetTestCasesFromEnvironment()) {
                Console.WriteLine($"Running testcase '{testcase}' with {threads} threads.");

                var C = CreateControl(testcase);
                // Warm-up
                RunXNSEControl(C);
                CollectProfilingFiles(testcase, threads, delete: true);


                var timing = RunXNSEControl(C);
                Console.WriteLine(
                        $"XNSE_THREAD_RESULT testcase={testcase}, threads={threads}, phase=, time={timing:g6} s"
                    );
                AppendTiming(
                    file: file,
                    testcase: testcase,
                    phase: "",
                    seconds: timing
                );
                CollectProfilingFiles(testcase, threads);
            }

            Console.WriteLine($"Wrote timing file: {file}");
        }

        [Test]
        public static void Check() {
            var t1 = ReadTimings(1);
            var t2 = ReadTimings(2);
            var t4 = ReadTimings(4);

            string[] testcases = t1
                .Select(r => r.Testcase)
                .Intersect(t2.Select(r => r.Testcase))
                .Intersect(t4.Select(r => r.Testcase))
                .Distinct()
                .OrderBy(s => s)
                .ToArray();

            Assert.That(testcases.Length, Is.GreaterThan(0), "No common testcases found in timing files.");

            foreach(string testcase in testcases) {
                double T1 = GetPhaseTime(t1, testcase, "");
                double T2 = GetPhaseTime(t2, testcase, "");
                double T4 = GetPhaseTime(t4, testcase, "");

                double speedup2 = T1 / T2;
                double speedup4 = T1 / T4;

                Console.WriteLine($"XNSE_THREAD_SCALING_SUMMARY testcase={testcase}");
                Console.WriteLine($"  t1       = {T1:g6} s");
                Console.WriteLine($"  t2       = {T2:g6} s");
                Console.WriteLine($"  t4       = {T4:g6} s");
                Console.WriteLine($"  speedup2 = {speedup2:g6}");
                Console.WriteLine($"  speedup4 = {speedup4:g6}");

                Assert.That(T2, Is.LessThan(0.98 * T1),
                    $"Expected 2 threads to be faster than 1 thread for testcase '{testcase}'. t1={T1:g6}, t2={T2:g6}");

                Assert.That(T4, Is.LessThan(0.95 * T1),
                    $"Expected 4 threads to be at least 5% faster than 1 thread for testcase '{testcase}'. t1={T1:g6}, t4={T4:g6}");

                Assert.That(T4, Is.LessThanOrEqualTo(1.05 * T2),
                    $"Expected 4 threads not to be significantly slower than 2 threads for testcase '{testcase}'. t2={T2:g6}, t4={T4:g6}");
            }
        }

        private static int GetThreadCountFromEnvironment() {
            string omp = Environment.GetEnvironmentVariable("OMP_NUM_THREADS");

            if(string.IsNullOrWhiteSpace(omp))
                Assert.Fail("OMP_NUM_THREADS is not set.");

            if(!int.TryParse(omp, out int threads))
                Assert.Fail($"Unable to parse OMP_NUM_THREADS='{omp}'.");

            if(threads != 1 && threads != 2 && threads != 4)
                Assert.Fail($"Unexpected OMP_NUM_THREADS={threads}. Expected 1, 2, or 4.");

            return threads;
        }

        private static string[] GetTestCasesFromEnvironment() {
            string cases = Environment.GetEnvironmentVariable("XNSE_PERF_CASES");

            if(string.IsNullOrWhiteSpace(cases)) {
                return new[] {
                    "RisingBubble",
                    "RisingBubbleMedium",
                };
            }

            return cases
                .Split(new[] { ',', ';' }, StringSplitOptions.RemoveEmptyEntries)
                .Select(s => s.Trim())
                .Where(s => s.Length > 0)
                .ToArray();
        }



        private static double RunXNSEControl(XNSE_Control control) {
            XNSE solver = null;
            double totalSeconds = 0;
            try {

                solver = new XNSE();
                solver.Init(control);
                Stopwatch total = Stopwatch.StartNew();
                solver.RunSolverMode();
                total.Stop();
                totalSeconds = total.Elapsed.TotalSeconds;
            } catch(Exception e) {
                Assert.Fail($"Exception during XNSE run: {e}");
            } finally {
                solver?.Dispose();
            }
            return totalSeconds;
        }


        private static XNSE_Control CreateControl(string testcase) {
            XNSE_Control C;

            switch(testcase) {
                case "RisingBubble":
                C = XNSE_Solver_MPItest.RisingBubble(p: 2, kelem: 18);
                C.NoOfTimesteps = 4;
                C.ProjectName = "XNSE/MultithreadingScalingTest/Bubble";
                break;

                case "RisingBubbleSmall":
                C = XNSE_Solver_MPItest.RisingBubble(p: 1, kelem: 8);
                C.NoOfTimesteps = 8;
                C.ProjectName = "XNSE/MultithreadingScalingTest/BubbleSmall";
                break;

                case "RisingBubbleMedium":
                C = XNSE_Solver_MPItest.RisingBubble(p: 2, kelem: 16);
                C.NoOfTimesteps = 6;
                C.ProjectName = "XNSE/MultithreadingScalingTest/BubbleMedium";
                break;

                default:
                Assert.Fail($"Unknown XNSE performance testcase '{testcase}'.");
                throw new NotSupportedException();
            }

            C.ImmediatePlotPeriod = -1;
            C.SuperSampling = 0;
            C.savetodb = false;
            C.SessionName = $"XNSE_Perf_{testcase}";

            return C;
        }

        private static void EnsureTimingHeader(string file) {
            if(File.Exists(file))
                return;

            File.WriteAllText(
                file,
                "TimestampUtc,Testcase,Phase,Seconds,OMP_NUM_THREADS" + Environment.NewLine
            );
        }

        private static void AppendTiming(string file, string testcase, string phase, double seconds) {
            string row = string.Join(",",
                DateTime.UtcNow.ToString("o", CultureInfo.InvariantCulture),
                testcase,
                phase,
                seconds.ToString("R", CultureInfo.InvariantCulture),
                Environment.GetEnvironmentVariable("OMP_NUM_THREADS") ?? ""
            );

            File.AppendAllText(file, row + Environment.NewLine);
        }

        private sealed class TimingRow {
            public string Testcase;
            public string Phase;
            public double Seconds;
        }

        private static TimingRow[] ReadTimings(int threads) {
            string file = $"xnse-thread-{threads}.time";

            Assert.That(File.Exists(file), Is.True,
                $"Missing timing file '{file}'.");

            return File.ReadAllLines(file)
                .Skip(1)
                .Where(line => !string.IsNullOrWhiteSpace(line))
                .Select(ParseTimingRow)
                .ToArray();
        }

        private static TimingRow ParseTimingRow(string line) {
            string[] parts = line.Split(',');

            if(parts.Length < 5)
                Assert.Fail($"Invalid timing row: '{line}'.");

            return new TimingRow() {
                Testcase = parts[1],
                Phase = parts[2],
                Seconds = double.Parse(parts[3], CultureInfo.InvariantCulture)
            };
        }

        private static double GetPhaseTime(TimingRow[] rows, string testcase, string phase) {
            double[] values = rows
                .Where(r => r.Testcase == testcase)
                .Where(r => r.Phase == phase)
                .Select(r => r.Seconds)
                .OrderBy(x => x)
                .ToArray();

            Assert.That(values.Length, Is.GreaterThan(0),
                $"No timing data for testcase '{testcase}', phase '{phase}'.");

            return values[values.Length / 2];
        }

        private static void CollectProfilingFiles(string testcase, int threads, bool delete = false) {
            string stamp = DateTime.UtcNow.ToString("yyyyMMdd_HHmmss_fff", CultureInfo.InvariantCulture);

            string targetDir = Path.Combine(
                Directory.GetCurrentDirectory(),
                $"profiling_XNSE_{testcase}_threads{threads}_{stamp}"
            );
            if(!delete)
                Directory.CreateDirectory(targetDir);

            string[] files = Directory
                .GetFiles(Directory.GetCurrentDirectory(), "profiling_bin.*.*.bin")
                .Concat(Directory.GetFiles(Directory.GetCurrentDirectory(), "trace_*.txt"))
                .Distinct()
                .ToArray();

            foreach(string file in files) {
                if(delete) {
                    try {
                        File.Delete(file);
                    } catch(Exception e) {
                        Console.WriteLine($"Could not delete profiling file '{file}': {e.Message}");
                    }
                    continue;
                }

                string dest = Path.Combine(targetDir, Path.GetFileName(file));

                try {
                    File.Copy(file, dest, overwrite: true);
                } catch(Exception e) {
                    Console.WriteLine($"Could not copy profiling file '{file}' to '{dest}': {e.Message}");
                }
            }

            Console.WriteLine($"Collected {files.Length} profiling/trace files in '{targetDir}'.");
        }
    }
}