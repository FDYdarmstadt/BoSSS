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
            var timingsByThreads = ReadTimingsForAllThreadCounts();
            Assert.That(timingsByThreads.Count, Is.GreaterThanOrEqualTo(2),
                "At least two thread timing files are required for a scaling check.");

            Console.WriteLine("AVAILABLE_THREAD_COUNTS: " +
                string.Join(", ", timingsByThreads.Keys.OrderBy(threads => threads)));

            string[] commonTestcases = timingsByThreads.Values
                .Select(rows => rows.Select(r => r.Testcase).Distinct())
                .Aggregate((common, names) => common.Intersect(names))
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

            var profilingDatapoints = CreateExtraPFilesForAllTestcases();

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
                var scalingChecks = new[] {
                    new {
                        Name = (string)null,
                        Requirements = testcase.Requirements ?? Array.Empty<ThreadScalingRequirement>()
                    }
                }.Concat((testcase.ScalingMethods ?? Array.Empty<ThreadScalingMethod>())
                    .Select(method => new {
                        method.Name,
                        Requirements = method.Requirements ?? Array.Empty<ThreadScalingRequirement>()
                    }));

                foreach(var scalingCheck in scalingChecks) {
                    if(scalingCheck.Requirements.Length == 0)
                        continue;

                    string scalingTarget = scalingCheck.Name == null
                        ? testcase.Name
                        : $"{testcase.Name}, method '{scalingCheck.Name}'";

                    string[] regions = null;
                    if(scalingCheck.Name != null) {
                        Assert.That(profilingDatapoints.ContainsKey(testcase.Name), Is.True,
                            $"No Extra-P profiling data for testcase '{testcase.Name}'.");

                        regions = ResolveProfilingRegions(profilingDatapoints[testcase.Name], scalingCheck.Name);
                    }

                    foreach(ThreadScalingRequirement requirement in scalingCheck.Requirements) {
                        Assert.That(requirement.FromThreads, Is.GreaterThan(0),
                            $"Invalid source thread count for {scalingTarget}.");
                        Assert.That(requirement.ToThreads, Is.GreaterThan(0),
                            $"Invalid target thread count for {scalingTarget}.");
                        Assert.That(requirement.FromThreads, Is.Not.EqualTo(requirement.ToThreads),
                            $"Thread counts must differ for {scalingTarget}.");

                        double fromTime = scalingCheck.Name == null
                            ? GetAveragePhaseTime(GetTimingRows(timingsByThreads, requirement.FromThreads), testcase.Name, "")
                            : GetAverageRegionTime(profilingDatapoints[testcase.Name], requirement.FromThreads, regions);
                        double toTime = scalingCheck.Name == null
                            ? GetAveragePhaseTime(GetTimingRows(timingsByThreads, requirement.ToThreads), testcase.Name, "")
                            : GetAverageRegionTime(profilingDatapoints[testcase.Name], requirement.ToThreads, regions);
                        double speedup = fromTime / toTime;

                        Console.WriteLine($"THREAD_SCALING_SUMMARY testcase={scalingTarget}, threads={requirement.FromThreads}->{requirement.ToThreads}");
                        Console.WriteLine($"  t{requirement.FromThreads} avg = {fromTime:g6} s");
                        Console.WriteLine($"  t{requirement.ToThreads} avg = {toTime:g6} s");
                        PrintSpeedup($"{requirement.FromThreads}->{requirement.ToThreads}", speedup, requirement.MinSpeedup);

                        if(requirement.MinSpeedup > 0 && speedup < requirement.MinSpeedup) {
                            failures.Add(
                                $"Expected stronger {requirement.FromThreads}->{requirement.ToThreads} thread speedup for {scalingTarget}. " +
                                $"T{requirement.FromThreads}={fromTime:g6}, T{requirement.ToThreads}={toTime:g6}, " +
                                $"speedup={speedup:g6}, expected >= {requirement.MinSpeedup:g6}"
                            );
                        }
                    }
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

            if(threads <= 0)
                Assert.Fail($"OMP_NUM_THREADS must be positive, but was {threads}.");

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

        private static void PrintSpeedup(string transition, double speedup, double minimum) {
            if(minimum <= 0) {
                Console.WriteLine($"  speedup {transition} = {speedup:g6}, check disabled");
            } else {
                string result = speedup >= minimum ? "OK" : "FAIL";
                Console.WriteLine($"  speedup {transition} = {speedup:g6}, expected >= {minimum:g6}: {result}");
            }
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

        private static Dictionary<int, TimingRow[]> ReadTimingsForAllThreadCounts() {
            var filenamePattern = new System.Text.RegularExpressions.Regex(
                @"^thread-(?<threads>\d+)\.time$",
                System.Text.RegularExpressions.RegexOptions.Compiled
            );

            var timingFiles = Directory
                .GetFiles(Directory.GetCurrentDirectory(), "thread-*.time")
                .Select(file => new { File = file, Match = filenamePattern.Match(Path.GetFileName(file)) })
                .Where(x => x.Match.Success)
                .Select(x => new {
                    Threads = int.Parse(x.Match.Groups["threads"].Value, CultureInfo.InvariantCulture),
                    x.File
                })
                .ToArray();

            Assert.That(timingFiles.Length, Is.GreaterThan(0), "No thread timing files found.");
            Assert.That(timingFiles.Select(x => x.Threads).Distinct().Count(), Is.EqualTo(timingFiles.Length),
                "Multiple timing files found for the same thread count.");

            return timingFiles.ToDictionary(
                x => x.Threads,
                x => File.ReadAllLines(x.File)
                    .Skip(1)
                    .Where(line => !string.IsNullOrWhiteSpace(line))
                    .Select(ParseTimingRow)
                    .ToArray()
            );
        }

        private static TimingRow[] GetTimingRows(
            IReadOnlyDictionary<int, TimingRow[]> timingsByThreads,
            int threads
        ) {
            Assert.That(timingsByThreads.ContainsKey(threads), Is.True,
                $"Missing timing file for {threads} threads.");

            return timingsByThreads[threads];
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

        private static string[] ResolveProfilingRegions(
            IEnumerable<ExtraPDatapoint> datapoints,
            string methodSelector
        ) {
            string[] matches = datapoints
                .SelectMany(dp => dp.Records?.Keys ?? Enumerable.Empty<string>())
                .Distinct()
                .Where(region =>
                    region.Equals(methodSelector, StringComparison.Ordinal) ||
                    region.EndsWith("->" + methodSelector, StringComparison.Ordinal) ||
                    region.EndsWith("." + methodSelector, StringComparison.Ordinal)
                )
                .ToArray();

            Assert.That(
                matches.Length,
                Is.GreaterThan(0),
                $"Method selector '{methodSelector}' did not match an Extra-P region."
            );

            Console.WriteLine(
                $"THREAD_SCALING_METHOD selector='{methodSelector}', regions={matches.Length}"
            );

            return matches;
        }

        private static double GetAverageRegionTime(
            IEnumerable<ExtraPDatapoint> datapoints,
            int threads,
            IEnumerable<string> regions
        ) {
            string[] selectedRegions = regions.ToArray();
            double[] values = datapoints
                .Where(dp => dp.ThreadSize == threads)
                .Select(dp => {
                    Assert.That(dp.Records, Is.Not.Null,
                        "Profiling records are unavailable.");

                    string[] missingRegions = selectedRegions
                        .Where(region => !dp.Records.ContainsKey(region))
                        .ToArray();
                    Assert.That(missingRegions, Is.Empty,
                        $"Extra-P region(s) missing for {threads} threads in session '{dp.SessionName}': " +
                        string.Join(", ", missingRegions));

                    return selectedRegions.Sum(region =>
                        dp.Records[region].Sum(ExtraPAnalysis.GetInclusiveSeconds)
                    );
                })
                .ToArray();

            Assert.That(values.Length, Is.GreaterThan(0),
                $"No Extra-P profiling data for selected regions with {threads} threads.");

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

        private static Dictionary<string, List<ExtraPDatapoint>> CreateExtraPFilesForAllTestcases(bool exclusive = false) {
            string cwd = Directory.GetCurrentDirectory();
            var result = new Dictionary<string, List<ExtraPDatapoint>>();

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

                result.Add(testcase, datapoints);

                Console.WriteLine(
                    $"Created Extra-P file '{outputFile}' with {datapoints.Count} datapoints."
                );
            }

            return result;
        }
    }
}
