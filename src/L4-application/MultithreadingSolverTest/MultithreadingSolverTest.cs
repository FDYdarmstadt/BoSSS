using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using BoSSS.Application.BoSSSpad;

namespace BoSSS.Application.MultithreadingSolverTest {

    [TestFixture]
    public static class MultithreadingSolverTest {

        [Test]
        public static void Perf() {
            int threads = GetThreadCountFromEnvironment();

            Console.WriteLine($"Running thread performance test with {threads} threads.");
            Console.WriteLine($"OMP_NUM_THREADS = {Environment.GetEnvironmentVariable("OMP_NUM_THREADS")}");

            string file = $"thread-{threads}.time";
            ResetTimingFile(file);

            foreach(var testcase in GetTestCasesFromEnvironment()) {
                Console.WriteLine($"Running testcase '{testcase.Name}' with {threads} threads.");
                Console.WriteLine($"Measured runs = {testcase.Runs}");

                for(int run = 1; run <= testcase.Runs; run++) {
                    // Remove profiling files from previous runs before starting this run.
                    CollectProfilingFiles(testcase.Name, threads, run, delete: true);

                    double timing = testcase.RunControl(testcase.CreateControl());

                    Console.WriteLine(
                        $"THREAD_RESULT testcase={testcase.Name}, threads={threads}, run={run}, time={timing:g6} s"
                    );

                    AppendTiming(
                        file: file,
                        testcase: testcase.Name,
                        phase: "",
                        run: run,
                        seconds: timing
                    );

                    // Keep one profiling folder per testcase/thread/run combination.
                    CollectProfilingFiles(testcase.Name, threads, run);
                }
            }

            Console.WriteLine($"Wrote timing file: {file}");
        }

        [Test]
        public static void Check() {
            var t1 = ReadTimings(1);
            var t2 = ReadTimings(2);
            var t4 = ReadTimings(4);

            string[] commonTestcases = t1
                .Select(r => r.Testcase)
                .Intersect(t2.Select(r => r.Testcase))
                .Intersect(t4.Select(r => r.Testcase))
                .Distinct()
                .OrderBy(s => s)
                .ToArray();

            Console.WriteLine("COMMON_THREAD_SCALING_TESTCASES:");
            foreach(var testcase in commonTestcases) {
                Console.WriteLine($"  {testcase}");
            }

            Assert.That(
                commonTestcases.Length,
                Is.GreaterThan(0),
                "No common testcases found in timing files."
            );

            CreateExtraPFilesForAllTestcases();

            var selectedTestcases = GetTestCasesFromEnvironment()
                .Where(tc => commonTestcases.Contains(tc.Name))
                .OrderBy(tc => tc.Name)
                .ToArray();

            Console.WriteLine("SELECTED_THREAD_SCALING_TESTCASES:");
            foreach(var testcase in selectedTestcases) {
                Console.WriteLine($"  {testcase.Name}");
            }

            Assert.That(
                selectedTestcases.Length,
                Is.GreaterThan(0),
                "No configured testcases from the environment matched the common timing testcases."
            );

            var failures = new List<string>();

            foreach(var testcase in selectedTestcases) {
                double T1 = GetAveragePhaseTime(t1, testcase.Name, "");
                double T2 = GetAveragePhaseTime(t2, testcase.Name, "");
                double T4 = GetAveragePhaseTime(t4, testcase.Name, "");

                double speedup12 = T1 / T2;
                double speedup14 = T1 / T4;
                double speedup24 = T2 / T4;

                Console.WriteLine($"THREAD_SCALING_SUMMARY testcase={testcase.Name}");
                Console.WriteLine($"  t1 avg       = {T1:g6} s");
                Console.WriteLine($"  t2 avg       = {T2:g6} s");
                Console.WriteLine($"  t4 avg       = {T4:g6} s");
                Console.WriteLine($"  speedup 1->2 = {speedup12:g6}, expected >= {testcase.MinSpeedup12:g6}");
                Console.WriteLine($"  speedup 1->4 = {speedup14:g6}, expected >= {testcase.MinSpeedup14:g6}");
                Console.WriteLine($"  speedup 2->4 = {speedup24:g6}, expected >= {testcase.MinSpeedup24:g6}");

                if(speedup12 < testcase.MinSpeedup12) {
                    failures.Add(
                        $"Expected stronger 1->2 thread speedup for testcase '{testcase.Name}'. " +
                        $"T1={T1:g6}, T2={T2:g6}, speedup={speedup12:g6}, expected >= {testcase.MinSpeedup12:g6}"
                    );
                }

                if(speedup14 < testcase.MinSpeedup14) {
                    failures.Add(
                        $"Expected stronger 1->4 thread speedup for testcase '{testcase.Name}'. " +
                        $"T1={T1:g6}, T4={T4:g6}, speedup={speedup14:g6}, expected >= {testcase.MinSpeedup14:g6}"
                    );
                }

                if(speedup24 < testcase.MinSpeedup24) {
                    failures.Add(
                        $"Expected stronger 2->4 thread speedup for testcase '{testcase.Name}'. " +
                        $"T2={T2:g6}, T4={T4:g6}, speedup={speedup24:g6}, expected >= {testcase.MinSpeedup24:g6}"
                    );
                }
            }

            Assert.That(
                failures,
                Is.Empty,
                "Thread scaling failures:" + Environment.NewLine + string.Join(Environment.NewLine, failures)
            );
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

        private static ThreadScalingTestcase[] GetTestCasesFromEnvironment() {
            string cases = Environment.GetEnvironmentVariable("THREAD_PERF_CASES");

            if(string.IsNullOrWhiteSpace(cases))
                return MultithreadingSolverTestcases.All;

            string[] selectedNames = cases
                .Split(new[] { ',', ';' }, StringSplitOptions.RemoveEmptyEntries)
                .Select(s => s.Trim())
                .Where(s => s.Length > 0)
                .ToArray();

            var selected = MultithreadingSolverTestcases.All
                .Where(tc => selectedNames.Contains(tc.Name))
                .ToArray();

            string[] missing = selectedNames
                .Except(selected.Select(tc => tc.Name))
                .ToArray();

            if(missing.Length > 0)
                Assert.Fail("Unknown thread performance testcase(s): " + string.Join(", ", missing));

            return selected;
        }

        private static void ResetTimingFile(string file) {
            File.WriteAllText(
                file,
                "TimestampUtc,Testcase,Phase,Run,Seconds,OMP_NUM_THREADS" + Environment.NewLine
            );
        }

        private static void AppendTiming(
            string file,
            string testcase,
            string phase,
            int run,
            double seconds
        ) {
            string row = string.Join(",",
                DateTime.UtcNow.ToString("o", CultureInfo.InvariantCulture),
                testcase,
                phase,
                run.ToString(CultureInfo.InvariantCulture),
                seconds.ToString("R", CultureInfo.InvariantCulture),
                Environment.GetEnvironmentVariable("OMP_NUM_THREADS") ?? ""
            );

            File.AppendAllText(file, row + Environment.NewLine);
        }

        private sealed class TimingRow {
            public string Testcase;
            public string Phase;
            public int Run;
            public double Seconds;
        }

        private static TimingRow[] ReadTimings(int threads) {
            string file = $"thread-{threads}.time";

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

            if(parts.Length < 6)
                Assert.Fail($"Invalid timing row: '{line}'.");

            return new TimingRow {
                Testcase = parts[1],
                Phase = parts[2],
                Run = int.Parse(parts[3], CultureInfo.InvariantCulture),
                Seconds = double.Parse(parts[4], CultureInfo.InvariantCulture)
            };
        }

        private static double GetAveragePhaseTime(TimingRow[] rows, string testcase, string phase) {
            double[] values = rows
                .Where(r => r.Testcase == testcase)
                .Where(r => r.Phase == phase)
                .Select(r => r.Seconds)
                .ToArray();

            Assert.That(values.Length, Is.GreaterThan(0),
                $"No timing data for testcase '{testcase}', phase '{phase}'.");

            return values.Average();
        }

        private static void CollectProfilingFiles(string testcase, int threads, int run, bool delete = false) {
            string cwd = Directory.GetCurrentDirectory();

            string targetDir = Path.Combine(
                cwd,
                $"profiling_{testcase}_threads{threads}_run{run}"
            );

            if(!delete)
                Directory.CreateDirectory(targetDir);

            string profilingBin = Directory
                .GetFiles(cwd, "profiling_bin.json")
                .OrderByDescending(File.GetLastWriteTimeUtc)
                .FirstOrDefault();

            var files = new List<string>();

            if(profilingBin != null)
                files.Add(profilingBin);
            else if(!delete)
                Console.WriteLine($"No profiling_bin.json found for testcase '{testcase}', threads={threads}, run={run}.");

            files.AddRange(Directory.GetFiles(cwd, "trace_*.txt"));
            files.AddRange(Directory.GetFiles(cwd, "profiling_summary.*.*.txt"));
            

            files = files
                .Distinct()
                .ToList();

            foreach(string file in files) {
                if(delete) {
                    try {
                        File.Delete(file);
                    } catch(Exception e) {
                        Console.WriteLine($"Could not delete profiling file '{file}': {e.Message}");
                    }
                    continue;
                }

                string destName = Path.GetFileName(file);

                if(file == profilingBin)
                    destName = $"profiling_bin.{testcase}.threads{threads}.run{run}.bin";

                string dest = Path.Combine(targetDir, destName);

                try {
                    File.Copy(file, dest, overwrite: true);
                } catch(Exception e) {
                    Console.WriteLine($"Could not copy profiling file '{file}' to '{dest}': {e.Message}");
                }
            }

            if(!delete) {
                Console.WriteLine(
                    $"Collected {files.Count} profiling/trace files in '{targetDir}'."
                );
            }
        }

        private static void CreateExtraPFilesForAllTestcases(bool exclusive = false) {
            string cwd = Directory.GetCurrentDirectory();

            var folderRegex = new System.Text.RegularExpressions.Regex(
                @"^profiling_(?<testcase>.+)_threads(?<threads>\d+)_run(?<run>\d+)$",
                System.Text.RegularExpressions.RegexOptions.Compiled
            );

            var profilingFolders = Directory
                .GetDirectories(cwd, "profiling_*_threads*_run*")
                .Select(dir => new {
                    Dir = dir,
                    Name = Path.GetFileName(dir),
                    Match = folderRegex.Match(Path.GetFileName(dir))
                })
                .Where(x => x.Match.Success)
                .Select(x => new {
                    x.Dir,
                    Testcase = x.Match.Groups["testcase"].Value,
                    Threads = int.Parse(x.Match.Groups["threads"].Value),
                    Run = int.Parse(x.Match.Groups["run"].Value)
                })
                .GroupBy(x => x.Testcase)
                .OrderBy(g => g.Key);

            foreach(var testcaseGroup in profilingFolders) {
                string testcase = testcaseGroup.Key;
                var datapoints = new List<ExtraPDatapoint>();

                foreach(var entry in testcaseGroup.OrderBy(x => x.Threads).ThenBy(x => x.Run)) {
                    try {
                        var dp = ExtraPAnalysis.LoadDatapointFromProfilingFolder(
                            profilingFolder: entry.Dir,
                            mpi: 1,
                            threads: entry.Threads,
                            p: 1,
                            d: 1,
                            cells: 1,
                            timesteps: 1,
                            sessionName: $"{testcase}_run{entry.Run}"
                        );

                        datapoints.Add(dp);

                        Console.WriteLine(
                            $"Added datapoint: testcase='{testcase}', threads={entry.Threads}, run={entry.Run}"
                        );
                    } catch(Exception e) {
                        Console.WriteLine(
                            $"Could not create datapoint for testcase='{testcase}', folder='{entry.Dir}': {e.Message}"
                        );
                    }
                }

                if(datapoints.Count == 0) {
                    Console.WriteLine($"No datapoints for testcase '{testcase}'. Skipping.");
                    continue;
                }

                string outputFile = Path.Combine(cwd, $"{testcase}.extra-p.txt");

                ExtraPAnalysis.WriteExtraPFile(
                    datapoints,
                    outputFile,
                    exclusive: exclusive
                );

                Console.WriteLine(
                    $"Created Extra-P file '{outputFile}' with {datapoints.Count} datapoints."
                );
            }
        }
    }
}
