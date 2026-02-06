#nullable enable

using BoSSS.Foundation.IO;
using BoSSS.Solution;
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
                    var sig = SelectSignificant(rec, mustIncludeMethods, max_count);

                    foreach(var kv in sig) {
                        if(dp.Records != null)
                            dp.Records[kv.Key] = kv.Value;
                        if (!strip)
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

            sw.WriteLine("POINTS " + string.Join(" ", datapoints.Select(dp => $"({dp.MpiSize} {dp.ThreadSize} {dp.Cells})")));
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
}