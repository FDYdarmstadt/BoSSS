#nullable enable

using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Text.Json.Serialization;


namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Represents a single datapoint collected from a simulation session for Extra-P analysis.
    /// </summary>
    [DataContract]
    public class ExtraPDatapoint {
        /// <summary>
        /// Name of the session this datapoint was collected from.
        /// </summary>
        [DataMember]
        public string SessionName = "[Empty]";

        /// <summary>
        /// Name of the solver used in this session.
        /// </summary>
        [DataMember] 
        public string SolverName = "[Empty]";

        /// <summary>
        /// Number of MPI processes used.
        /// </summary>
        [DataMember]
        public int MpiSize;

        /// <summary>
        /// Number of threads used.
        /// </summary>
        [DataMember]
        public int ThreadSize;

        /// <summary>
        /// Number of cells in the computational grid. (At the first timestep)
        /// </summary>
        [DataMember]
        public int Cells;

        /// <summary>
        /// Polynomial degree of the basis.
        /// </summary>
        [DataMember]
        public int P;

        /// <summary>
        /// Spatial dimension of the grid.
        /// </summary>
        [DataMember]
        public int D;

        /// <summary>
        /// Number of timesteps in the simulation.
        /// </summary>
        [DataMember]
        public int TimeSteps;

        /// <summary>
        /// Total runtime in seconds.
        /// </summary>
        [DataMember]
        public double TotalTime = -1.0f;

        /// <summary>
        /// Method call records collected during profiling.
        /// </summary>
        [JsonIgnore]
        public Dictionary<string, List<MethodCallRecord>>? Records = new();

        /// <summary>
        /// Strips profiling records and computes the total time from the root frame.
        /// </summary>
        public void Strip() {
            if(Records == null)
                return;
            Records.TryGetValue("root_frame", out var calls);
            if (calls == null || calls.Count == 0)
                return;
            TotalTime = ExtraPAnalysis.GetInclusiveSeconds(calls.FirstOrDefault());
            Records = null;
        }
    }

    /// <summary>
    /// Provides utility methods for analyzing Extra-P profiling data.
    /// Includes methods to calculate inclusive/exclusive time and generate ExtraPDatapoint lists/files.
    /// </summary>
    public static class ExtraPAnalysis {
        /// <summary>
        /// Returns the inclusive time in seconds for a given method call record.
        /// </summary>
        public static double GetInclusiveSeconds(MethodCallRecord? mcr) {
            if(mcr == null) return 0.0;
            return mcr.TimeSpentInMethod.TotalSeconds;
        }

        /// <summary>
        /// Returns the exclusive time in seconds for a given method call record.
        /// </summary>
        public static double GetExclusiveSeconds(MethodCallRecord? mcr) {
            if(mcr == null) return 0.0;
            return mcr.TimeExclusive.TotalSeconds;
        }

        static bool IsRecursiveDescendant(MethodCallRecord mcr) {
            var parent = mcr.ParrentCall;
            while(parent != null) {
                if(parent == mcr) return false; // paranoia
                if(parent.Name == mcr.Name) return true;
                parent = parent.ParrentCall;
            }
            return false;
        }

        static string CallPath(MethodCallRecord call) {
            var s = new Stack<string>();
            for(var n = call; n != null; n = n.ParrentCall)
                s.Push(n.Name);
            return string.Join("->", s);
        }

        static Dictionary<string, List<MethodCallRecord>> TakeFirstProfile(OnlineProfiling profile) {
            var map = new Dictionary<string, List<MethodCallRecord>>();

            foreach(var call in profile.RootCall.Flatten()) {
                var key = CallPath(call);
                (map.TryGetValue(key, out var l) ? l : map[key] = new()).Add(call);
            }

            return map;
        }

        static HashSet<string> ExpandWithParents(Dictionary<string, List<MethodCallRecord>> dict, HashSet<string> include) {
            var result = new HashSet<string>(include);

            foreach(var calls in dict.Values) {
                foreach(var r in calls) {
                    if(!include.Contains(r.Name))
                        continue;

                    var parent = r.ParrentCall;
                    while(parent != null) {
                        result.Add(parent.Name);
                        parent = parent.ParrentCall;
                    }
                }
            }
            return result;
        }


        static List<KeyValuePair<string, List<MethodCallRecord>>> SelectSignificant(Dictionary<string, List<MethodCallRecord>> dict, IEnumerable<string>? includeMethods = null, int max = 100) {

            HashSet<string>? allowed = null;

            if(includeMethods != null) {
                allowed = ExpandWithParents(
                    dict,
                    new HashSet<string>(includeMethods)
                );
            }

            // Compute times for all methods
            var timesAll = dict
                .Select(kv => new {
                    Key = kv.Key,
                    Exclusive = kv.Value.Sum(r => GetExclusiveSeconds(r)),
                    Inclusive = kv.Value
                        .Where(r => !IsRecursiveDescendant(r))
                        .Sum(r => GetInclusiveSeconds(r))
                })
                .OrderByDescending(x => x.Inclusive)
                .ToList();

            // First take methods from includeMethods
            var selected = timesAll
                .Where(t => allowed != null && allowed.Contains(t.Key))
                .ToList();

            // If not enough, fill with remaining methods
            if(selected.Count < max) {
                var remaining = timesAll
                    .Where(t => !selected.Any(s => s.Key == t.Key))
                    .Take(max - selected.Count);
                selected.AddRange(remaining);
            }

            // Convert to KeyValuePair list
            var sel = new List<KeyValuePair<string, List<MethodCallRecord>>>();
            foreach(var t in selected) {
                sel.Add(new KeyValuePair<string, List<MethodCallRecord>>(t.Key, dict[t.Key]));
            }
            return sel;
        }

        static ExtraPDatapoint CreateDatapoint(ISessionInfo session, OnlineProfiling profiling) {
            var dp = new ExtraPDatapoint();
            try {
                var firstTimestep = session.Timesteps.First();
                var fields = firstTimestep.FieldInitializers;
                IGridInfo grid = firstTimestep.Grid;
                dp.SessionName = session.Name;
                Type? solverType = session.GetControl()?.GetSolverType();
                dp.SolverName = solverType != null ? solverType.Namespace ?? "" + solverType.Name : "[NULL]";
                dp.MpiSize = profiling.MPISize;
                dp.ThreadSize = profiling.NumThreads;
                dp.Cells = grid.NumberOfCells;
                dp.P = 0;
                dp.D = session.Timesteps.First().Grid.SpatialDimension;
                foreach(var field in fields) {
                    //dp.Dofs += session.GetDOF(field.Identification); //this is expensive
                    dp.P = Math.Max(dp.P, field.BasisInfo.Degree);
                }
                dp.TimeSteps = session.Timesteps.Count;
                //dp.SessionRuntime = session.GetApproximateRunTime().Seconds;
            } catch (Exception e) { 
                System.Console.WriteLine($"Could not extract info from session {session.Name}: {e.Message}");
            }
            return dp;
        }


        /// <summary>
        /// Creates list of <see cref="ExtraPDatapoint"/> for a given list of sessions. 
        /// </summary>
        public static List<ExtraPDatapoint> LoadExtraPDatapoints(
            IEnumerable<ISessionInfo> sessions,
            bool strip = false,
            int max_count = 100,
            IEnumerable<string>? mustIncludeMethods = null
            ) {

            var dps = new List<ExtraPDatapoint>();
            var allMethods = new HashSet<string>();

            foreach(var s in sessions) {
                if(!s.SuccessfulTermination)
                    continue;
                var p0 = s.GetProfiling([0]).First();
                var dp = CreateDatapoint(s, p0);
                {
                    var rec = TakeFirstProfile(p0);
                    List<KeyValuePair<string, List<MethodCallRecord>>> sig;

                    if(strip) {
                        // Only select 'root_frame' when strip is true
                        sig = rec.TryGetValue("root_frame", out var rootCalls)
                            ? new List<KeyValuePair<string, List<MethodCallRecord>>> { new("root_frame", rootCalls) }
                            : new List<KeyValuePair<string, List<MethodCallRecord>>>();
                    } else {
                        sig = SelectSignificant(rec, mustIncludeMethods, max_count);
                    }

                    foreach(var kv in sig) {
                        if(dp.Records != null)
                            dp.Records[kv.Key] = kv.Value;
                        if(!strip)
                            allMethods.Add(kv.Key);
                    }
                }

                dps.Add(dp);
                if(strip)
                    dp.Strip();
                System.Console.WriteLine($"Session: {s.Name} added to the datapoints.");
            }
            if(!strip) { 
            // normalize missing regions
                foreach(var dp in dps) {
                    foreach(var m in allMethods) {
                        dp.Records?.TryAdd(m, new());
                    }
                }
            }
            return dps;
        }

        /// <summary>
        /// Writes a list of <see cref="ExtraPDatapoint"/> to an extra p txt file.
        /// </summary>
        public static void WriteExtraPFile(List<ExtraPDatapoint> datapoints, string outputFile, bool exclusive = true) {
            if(datapoints == null || datapoints.Count == 0)
                throw new ArgumentException("No datapoints available.");
            var grouped = datapoints
                .GroupBy(dp => new {
                    dp.MpiSize,
                    dp.ThreadSize,
                    dp.Cells
                })
                .OrderBy(g => g.Key.MpiSize)
                .ThenBy(g => g.Key.ThreadSize)
                .ThenBy(g => g.Key.Cells);

            // Clean output file
            if(File.Exists(outputFile))
                File.Delete(outputFile);

            var utf8NoBom = new UTF8Encoding(false);
            using var sw = new StreamWriter(outputFile, false, utf8NoBom);

            // PARAMETERS
            sw.WriteLine("PARAMETER mpi_procs");
            sw.WriteLine("PARAMETER threads");
            sw.WriteLine("PARAMETER mesh");

            sw.WriteLine("POINTS " + string.Join(" ", grouped.Select(g =>$"({g.Key.MpiSize} {g.Key.ThreadSize} {g.Key.Cells})")));
            sw.WriteLine();

            // METHODS
            if(datapoints[0].Records == null) { 
                Console.WriteLine("Stripped mode: No method records available in datapoints! Writing basic Extra-P-File.");
                sw.WriteLine($"REGION root_frame");
                sw.WriteLine("METRIC time");
                foreach(var group in grouped) {
                    var values = group.Select(dp => dp.TotalTime);
                    sw.WriteLine($"DATA {string.Join(" ", values)}");
                }
            } else {
                var allMethods = datapoints
                    .SelectMany(dp => dp.Records?.Keys ?? Enumerable.Empty<string>())
                    .Distinct()
                    .OrderBy(m => m)
                    .ToList();

                foreach(var method in allMethods) {
                    sw.WriteLine($"REGION {method}");
                    sw.WriteLine("METRIC time");

                    foreach(var group in grouped) {

                        var values = group.Select(dp => {
                            if(dp.Records != null &&
                                dp.Records.TryGetValue(method, out var calls)) {
                                if(exclusive) {
                                    return calls.Sum(c => GetExclusiveSeconds(c));
                                } else {
                                    // inclusive time without recursive double counting
                                    return calls
                                        .Where(c => !IsRecursiveDescendant(c))
                                        .Sum(c => GetInclusiveSeconds(c));
                                }
                            }
                            return 0.0;
                        });

                        sw.WriteLine($"DATA {string.Join(" ", values)}");
                    }

                    sw.WriteLine();
                }
            }
            Console.WriteLine($"Extra-P file written to {outputFile} containing {datapoints.Count} datapoints");
        }
    }

    /// <summary>
    /// Represents a performance model for simulations, allowing prediction of runtime 
    /// and recommendation of MPI + thread configurations based on training data.
    /// </summary>
    public class ExtraPModel {
        /// <summary>
        /// Scaling factor for the compute time.
        /// </summary>
        public double Ccompute { get; private set; }

        /// <summary>
        /// Fixed setup time for the simulation.
        /// </summary>
        public double Csetup { get; private set; }

        /// <summary>
        /// Exponent for polynomial degree scaling.
        /// </summary>
        public double AlphaP { get; private set; }

        /// <summary>
        /// Exponent for MPI size scaling.
        /// </summary>
        public double AlphaMPI { get; private set; }

        /// <summary>
        /// Exponent for thread count scaling.
        /// </summary>
        public double AlphaThreads { get; private set; }

        /// <summary>
        /// Maximum trusted number of MPI processes for strong scaling.
        /// </summary>
        public int MaxTrustedMPI { get; private set; } = int.MaxValue;

        /// <summary>
        /// Maximum number of threads to consider.
        /// </summary>
        public double MaxThreads { get; private set; } = int.MaxValue;
        
        /// <summary>
        /// Minimum number of cells per MPI process to ensure efficiency.
        /// </summary>
        public double MinCellsPerMPI { get; private set; } = 0.0;

        //public List<ExtraPDatapoint> TrainingData { get; private set; }

        /// <summary>
        /// Constructs an ExtraPModel with all parameters.
        /// </summary>
        public ExtraPModel(double ccompute, double csetup, double alphaP, double alphaMPI, double alphaThreads, int maxTrustedMPI, double maxThreads, double minCellsPerMPI
        ) {
            Ccompute = ccompute;
            Csetup = csetup;
            AlphaP = alphaP;
            AlphaMPI = alphaMPI;
            AlphaThreads = alphaThreads;
            MaxTrustedMPI = maxTrustedMPI;
            MaxThreads = maxThreads;
            MinCellsPerMPI = minCellsPerMPI;
        }

        private ExtraPModel() { }

        /// <summary>
        /// Fits an ExtraPModel from a list of datapoints.
        /// Detects limits, filters data, and fits the smooth scaling model.
        /// </summary>
        /// <param name="data">Training data consisting of ExtraPDatapoints.</param>
        /// <returns>A fitted ExtraPModel instance.</returns>
        public static ExtraPModel Fit(List<ExtraPDatapoint> data) {

            if(data.Count < 5)
                throw new ArgumentException("Not enough datapoints.");

            // --- 1. detect hard limits ---
            int maxMPI = DetectMaxTrustedMPI(data);
            int maxThreads = Math.Max(data.Max(d => d.ThreadSize), 16);
            double minCellsPerMPI = DetectMinCellsPerMPI(data, maxMPI);

            // --- 2. filter data ---
            var filtered = data.Where(dp =>
                dp.MpiSize <= maxMPI &&
                CellsPerMPI(dp) >= minCellsPerMPI
            ).ToList();

            // --- 3. fit smooth model ---
            var model = FitSmoothModel(filtered);

            model.MaxTrustedMPI = maxMPI;
            model.MinCellsPerMPI = minCellsPerMPI;
            model.MaxThreads = maxThreads;
            //model.TrainingData = filtered;

            return model;
        }

        static int DetectMaxTrustedMPI(List<ExtraPDatapoint> data) {

            // group by MPI, check monotonic scaling
            var groups = data.GroupBy(d => d.MpiSize)
                             .OrderBy(g => g.Key)
                             .ToList();

            int lastGood = groups.First().Key;

            for(int i = 1; i < groups.Count; i++) {
                double prevMedian = Median(groups[i - 1].Select(d => d.TotalTime));
                double currMedian = Median(groups[i].Select(d => d.TotalTime));

                // strong-scaling violation -> spike
                if(currMedian > 1.2 * prevMedian)
                    break;

                lastGood = groups[i].Key;
            }

            return lastGood;
        }

        static double DetectMinCellsPerMPI(List<ExtraPDatapoint> data, int maxMPI) {

            var valid = data.Where(d => d.MpiSize <= maxMPI)
                            .Select(d => CellsPerMPI(d))
                            .OrderBy(x => x)
                            .ToList();

            // conservative: lower 0.1 quantile
            int idx = (int)(0.1 * valid.Count);
            return valid[Math.Max(idx, 0)];
        }
        static ExtraPModel FitSmoothModel(List<ExtraPDatapoint> data) {

            double[] alphaPGrid = { 2.0, 2.5, 3.0 };
            double[] alphaMPIGrid = { 0.7, 0.8, 0.9 };
            double[] alphaThreadGrid = { 0.5, 0.6, 0.7 };

            ExtraPModel best = FitLinear(data, alphaPGrid[0], alphaMPIGrid[0], alphaThreadGrid[0]);
            double bestErr = double.PositiveInfinity;

            foreach(var aP in alphaPGrid)
                foreach(var aM in alphaMPIGrid)
                    foreach(var aT in alphaThreadGrid) {

                        var model = FitLinear(data, aP, aM, aT);
                        double err = RelativeRMSError(model, data);

                        if(err < bestErr) {
                            bestErr = err;
                            best = model;
                        }
                    }

            return best;
        }

        static ExtraPModel FitLinear(List<ExtraPDatapoint> data,  double alphaP, double alphaMPI, double alphaThreads) {
            int n = data.Count;
            var A = MultidimensionalArray.Create(n, 2);
            var b = new double[n];

            for(int i = 0; i < n; i++) {
                A[i, 0] = Work(data[i], alphaP, alphaMPI, alphaThreads);
                A[i, 1] = 1.0;
                b[i] = data[i].TotalTime;
            }

            double[] c = NNLS2(A, b);

            return new ExtraPModel {
                Ccompute = c[0],
                Csetup = c[1],
                AlphaP = alphaP,
                AlphaMPI = alphaMPI,
                AlphaThreads = alphaThreads
            };
        }

        /// <summary>
        /// Predicts the runtime for a given datapoint.
        /// </summary>
        /// <param name="dp">The ExtraPDatapoint to predict for.</param>
        /// <param name="ignore_limits">If false, checks MPI size and cells per MPI.</param>
        /// <returns>Predicted runtime.</returns>
        public double Predict(ExtraPDatapoint dp, bool ignore_limits = true) {

            if(!ignore_limits && dp.MpiSize > MaxTrustedMPI)
                throw new InvalidOperationException("MPI outside trusted range.");

            if(!ignore_limits && CellsPerMPI(dp) < MinCellsPerMPI)
                throw new InvalidOperationException("Too little work per core.");

            return Ccompute * Work(dp, AlphaP, AlphaMPI, AlphaThreads) + Csetup;
        }

        /// <summary>
        /// Recommends the MPI and thread configuration for a given problem size.
        /// Uses the efficiency parameter to filter near-optimal configurations.
        /// </summary>
        /// <param name="cells">Number of cells in the simulation.</param>
        /// <param name="p">Polynomial degree.</param>
        /// <param name="timesteps">Number of timesteps.</param>
        /// <param name="efficiency">
        /// Efficiency threshold (0..1). Configurations with predicted time ≤ t_min / efficiency are considered.
        /// </param>
        /// <returns>Tuple of recommended MPI, Threads, and predicted Time.</returns>
        public (int MPI, int Threads, double Time) RecommendConfiguration(int cells, int p, int timesteps, double efficiency = 0.7) {

            var candidates = new List<(int mpi, int th, double t)>();

            for(int mpi = 1; mpi <= MaxTrustedMPI; mpi++)
                for(int th = 1; th <= MaxThreads; th++) {

                    var dp = new ExtraPDatapoint {
                        Cells = cells,
                        P = p,
                        TimeSteps = timesteps,
                        MpiSize = mpi,
                        ThreadSize = th
                    };

                    if(CellsPerMPI(dp) < MinCellsPerMPI)
                        continue;

                    double t = Predict(dp);
                    candidates.Add((mpi, th, t));
                }

            if(!candidates.Any())
                throw new InvalidOperationException("No valid configurations found for the given problem size.");

            double tMin = candidates.Min(c => c.t);

            var efficient = candidates
                .Where(c => c.t <= tMin / efficiency)
                .OrderBy(c => c.t)
                .First();

            return efficient;
        }

        // HELPERS
        static double Work(ExtraPDatapoint dp, double aP, double aMPI, double aT) =>
            dp.Cells * Math.Pow(dp.P, aP) * dp.TimeSteps /
            (Math.Pow(dp.MpiSize, aMPI) * Math.Pow(dp.ThreadSize, aT));

        static double CellsPerMPI(ExtraPDatapoint dp) =>
            dp.Cells / (double)(dp.MpiSize);

        static double RelativeRMSError(ExtraPModel m, List<ExtraPDatapoint> d) =>
            Math.Sqrt(d.Average(dp => {
                double r = (m.Predict(dp) - dp.TotalTime) / dp.TotalTime;
                return r * r;
            }));

        static double Median(IEnumerable<double> x) {
            var a = x.OrderBy(v => v).ToArray();
            return a[a.Length / 2];
        }

        static double[] NNLS2(MultidimensionalArray A, double[] b) {
            int n = b.Length;
            var AtA = MultidimensionalArray.Create(2, 2);
            var Atb = new double[2];

            for(int i = 0; i < n; i++) {
                Atb[0] += A[i, 0] * b[i];
                Atb[1] += b[i];
                AtA[0, 0] += A[i, 0] * A[i, 0];
                AtA[0, 1] += A[i, 0];
                AtA[1, 0] += A[i, 0];
                AtA[1, 1] += 1.0;
            }

            double det = AtA[0, 0] * AtA[1, 1] - AtA[0, 1] * AtA[1, 0];
            return new[] {
                Math.Max(0, ( Atb[0]*AtA[1,1] - Atb[1]*AtA[0,1]) / det),
                Math.Max(0, (-Atb[0]*AtA[1,0] + Atb[1]*AtA[0,0]) / det)
            };
        }
    }

    /// <summary>
    /// Adds some Extensions to a AppControl 
    /// </summary>
    public static class AppControlExtraPExtensions {

        /// <summary>
        /// Predicts the number of cells
        /// </summary>
        public static int? PredictCells<AppControlType>(this AppControlType ctrl, IDatabaseInfo? db = null) where AppControlType : AppControl {
            int? cells = null;
            // Try GridFunc first
            if(ctrl.GridFunc != null) {
                try {
                    var grid = ctrl.GridFunc();
                    cells = grid?.NumberOfCells;
                } catch {
                    // ignore
                }
            }
            if(cells == null && ctrl.GridGuid != Guid.Empty && db == null) {
                System.Console.WriteLine("The grid appears to be stored in the database but you didn't provide a database!");
            }

            if(cells == null && ctrl.GridGuid != Guid.Empty && db != null) {
                try {
                    var gridInfo = db.Grids.Find(ctrl.GridGuid);
                    cells = gridInfo?.NumberOfCells;
                } catch {
                    System.Console.WriteLine("The grid could not be found in the database!");
                }
            }

            return cells;
        }

        /// <summary>
        /// Predicts the number of timesteps the simmulation may take
        /// </summary>
        public static (int, int, int) PredictNoOfTimesteps<AppControlType>(this AppControlType ctrl) where AppControlType : AppControl {
            int avg = -1, min = -1, max = -1;

            if(ctrl.NoOfTimesteps > 0) {
                avg = ctrl.NoOfTimesteps;
            } else {
                // Estimate timesteps from Endtime and dt
                double dt = double.NaN;
                if(ctrl.dtMin > 0 && ctrl.dtMax > 0) {
                    if(ctrl.dtMin == ctrl.dtMax) {
                        dt = ctrl.dtMin;
                    } else {
                        dt = 0.5 * (ctrl.dtMin + ctrl.dtMax);
                    }
                } else if(ctrl.dtFixed > 0 && !double.IsNaN(ctrl.dtFixed)) {
                    dt = ctrl.dtFixed;
                }
                if(dt > 0 && ctrl.Endtime > 0 && !double.IsNaN(ctrl.Endtime)) {
                    avg = (int)Math.Ceiling(ctrl.Endtime / dt);
                }
            }

            return (avg, min, max);
        }



        /// <summary>
        /// Predicts the runtime for a given MPI and thread count using the ExtraPModel for this control.
        /// Uses the polynomial degree, mesh cells, and timesteps from the control object.
        /// Returns null if no model is available or required information is missing.
        /// </summary>
        /// <param name="ctrl">The control object.</param>
        /// <param name="mpi">MPI process count.</param>
        /// <param name="threads">Thread count.</param>
        /// <param name="db">If the grid is stored in the DB, we need this</param>
        /// <returns>Predicted runtime in seconds, or null if no model or required info is available.</returns>
        public static double? PredictRuntime<AppControlType>(this AppControlType ctrl, int mpi, int threads, IDatabaseInfo? db = null) where AppControlType : AppControl {
            var model = ctrl.GetExtraPModel();
            if(model == null) {
                Console.WriteLine($"No ExtraPModel available for this control {ctrl.GetType()}.");
                return null;
            }

            int? cells = ctrl.PredictCells(db);
            if(cells == null) {
                Console.WriteLine("Could not determine number of cells.");
                return null;
            }

            int? p = null;
            if(ctrl.FieldOptions != null && ctrl.FieldOptions.Count > 0) {
                p = ctrl.FieldOptions.Values.Max(opt => opt.Degree);
            }
            if(p == null) {
                Console.WriteLine("Could not determine polynomial degree.");
                return null;
            }

            var (timestepsAvg, timestepsMin, timestepsMax) = ctrl.PredictNoOfTimesteps();
            if(timestepsAvg <= 0) {
                Console.WriteLine("Could not determine number of timesteps.");
                return null;
            }

            Console.WriteLine($"[ExtraP] Predicting runtime for:");
            Console.WriteLine($"  Cells: {cells}");
            Console.WriteLine($"  Polynomial degree (p): {p}");
            Console.WriteLine($"  Timesteps: avg={timestepsAvg}, min={timestepsMin}, max={timestepsMax}");
            Console.WriteLine($"  MPI: {mpi}");
            Console.WriteLine($"  Threads: {threads}");

            var dp = new ExtraPDatapoint {
                Cells = cells.Value,
                P = p.Value,
                TimeSteps = timestepsAvg,
                MpiSize = mpi,
                ThreadSize = threads
            };
            double predicted = model.Predict(dp, ignore_limits: false);
            Console.WriteLine($"[ExtraP] Predicted runtime (avg): {predicted} seconds");

            if(timestepsMin > 0 && timestepsMin != timestepsAvg) {
                var dpMin = new ExtraPDatapoint {
                    Cells = cells.Value,
                    P = p.Value,
                    TimeSteps = timestepsMin,
                    MpiSize = mpi,
                    ThreadSize = threads
                };
                double predictedMin = model.Predict(dpMin, ignore_limits: false);
                Console.WriteLine($"[ExtraP] Predicted runtime (min): {predictedMin} seconds");
            }
            if(timestepsMax > 0 && timestepsMax != timestepsAvg) {
                var dpMax = new ExtraPDatapoint {
                    Cells = cells.Value,
                    P = p.Value,
                    TimeSteps = timestepsMax,
                    MpiSize = mpi,
                    ThreadSize = threads
                };
                double predictedMax = model.Predict(dpMax, ignore_limits: false);
                Console.WriteLine($"[ExtraP] Predicted runtime (max): {predictedMax} seconds");
            }

            return predicted;
        }


        /// <summary>
        /// Returns the best (recommended) MPI/thread configuration and predicted runtime for the current control.
        /// Uses the polynomial degree, mesh cells, and timesteps from the control object.
        /// </summary>
        /// <param name="ctrl">The control object.</param>
        /// <param name="efficiency">Efficiency threshold (0..1), default 0.7.</param>
        /// <param name="db">If the grid is stored in the DB, we need this</param>
        /// <returns>Tuple (MPI, Threads, Time) or null if not possible.</returns>
        public static (int MPI, int Threads, double Time)? PredictBestConfig<AppControlType>(this AppControlType ctrl,double efficiency = 0.7,IDatabaseInfo? db = null) where AppControlType : AppControl {
            var model = ctrl.GetExtraPModel();
            if(model == null) {
                Console.WriteLine($"No ExtraPModel available for this control {ctrl.GetType()}.");
                return null;
            }

            int? cells = ctrl.PredictCells(db);
            if(cells == null) {
                Console.WriteLine("Could not determine number of cells.");
                return null;
            }

            int? p = null;
            if(ctrl.FieldOptions != null && ctrl.FieldOptions.Count > 0) {
                p = ctrl.FieldOptions.Values.Max(opt => opt.Degree);
            }
            if(p == null) {
                Console.WriteLine("Could not determine polynomial degree.");
                return null;
            }

            var (timestepsAvg, _, _) = ctrl.PredictNoOfTimesteps();
            if(timestepsAvg <= 0) {
                Console.WriteLine("Could not determine number of timesteps.");
                return null;
            }

            var best = model.RecommendConfiguration(
                cells: cells.Value,
                p: p.Value,
                timesteps: timestepsAvg,
                efficiency: efficiency
            );

            Console.WriteLine($"[ExtraP] Best configuration (efficiency={efficiency}):");
            Console.WriteLine($"  MPI: {best.MPI}");
            Console.WriteLine($"  Threads: {best.Threads}");
            Console.WriteLine($"  Predicted runtime: {best.Time} seconds");

            return best;
        }


        /// <summary>
        /// Returns an ExtraPModel for the current solver, if available.
        /// </summary>
        public static ExtraPModel? GetExtraPModel<AppControlType>(this AppControlType ctrl) where AppControlType : AppControl {
            return ctrl switch {
                XNSE_Control => new ExtraPModel(
                    ccompute: 0.0046165029629496495,
                    csetup: 27.961561804953913,
                    alphaP: 2.0,
                    alphaMPI: 0.7,
                    alphaThreads: 0.5,
                    maxTrustedMPI: 8,
                    maxThreads: 16,
                    minCellsPerMPI: 100
                ),
                _ => null
            };
        }

    }

}