using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution {



    [Serializable]
    [DataContract]
    public class OnlineProfiling {

        /// <summary>
        /// empty ctor to aid serialization
        /// </summary>
        private OnlineProfiling() {

        }

        public OnlineProfiling(AppControl __ctrl) {
            this.RootCall = Tracer.Root;
            this.AppStartTime = DateTime.Now;
            this.Computer =  ilPSP.Environment.MPIEnv.Hostname;
            this.UserName = System.Environment.UserName; 
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out this.MPIRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out this.MPISize);
            this.NumThreads = ilPSP.Environment.NumThreads;
            this.GitCommitHash = Application.GitCommitHash;
            {
                var _EnvironmentVars = new List<(string Name, string Value)>();
                var ev = System.Environment.GetEnvironmentVariables();
                foreach (string k in ev.Keys) {
                    _EnvironmentVars.Add((k, ev[k] as string));
                }
                this.EnvironmentVars = _EnvironmentVars.ToArray();
            }
            AppDirectory = Directory.GetCurrentDirectory();

            bool serializationWorks = true;
            try {
                string check = __ctrl.Serialize();
                AppControl.Deserialize(check);
            } catch (Exception) {
                serializationWorks = false;
            }
            this.Ctrl = serializationWorks ? __ctrl : null;
            CommandLine = Application.LatestCmdLineArgs?.CloneAs();

            OnlinePerformanceLog = ilPSP.OnlinePerformanceMeasurement.Log;
        }

        [DataMember]
        public ilPSP.OnlinePerformanceLog OnlinePerformanceLog;

        [DataMember]
        public string[] CommandLine;

        [DataMember]
        public (string Name, string Value)[] EnvironmentVars;

        [DataMember]
        public string AppDirectory;

        [DataMember]
        public DateTime AppStartTime;

        [DataMember]
        public DateTime AppEndTime;

        [DataMember]
        public AppControl Ctrl;


        [JsonIgnore]
        public TimeSpan AppRunTime {
            get {
                return (this.AppEndTime - this.AppStartTime);
            }
        }

        [DataMember]
        public string Computer;

        [DataMember]
        public string UserName;

        [DataMember]
        public string GitCommitHash;

        [DataMember]
        public int NumThreads;

        [DataMember]
        public int MPIRank;

        [DataMember]
        public int MPISize;



        /// <summary>
        /// Serializes this object and its childs into a JSON string.
        /// </summary>
        public string Serialize() {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };

            using (var wrt = new StringWriter()) {
                using (JsonWriter writer = new JsonTextWriter(wrt)) {  // Alternative: binary writer: JsonTextWriter
                    formatter.Serialize(writer, this);
                }

                return wrt.ToString();
            }

        }

        /// <summary>
        /// De-Serializes from a JSON string.
        /// </summary>
        public static OnlineProfiling Deserialize(string JsonString) {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };

            using (var rd = new StringReader(JsonString)) {
                OnlineProfiling profiling;
                using (JsonReader reader = new JsonTextReader(rd)) {  // Alternative: binary writer: JsonTextWriter
                    profiling = formatter.Deserialize<OnlineProfiling>(reader);

                    if(profiling.RootCall == null && JsonString.IsNonEmpty()) {
                        // try to load legacy file
                        var mcr = MethodCallRecord.Deserialize(JsonString);
                        profiling.RootCall = mcr;
                    }

                }

                profiling?.RootCall?.FixParrent();

                return profiling;
            }
        }

        /// <summary>
        /// Stack information 
        /// </summary>
        [DataMember]
        public MethodCallRecord RootCall {
            get;
            private set;
        }


        public void UpdateDGInfo(IGrid g, IEnumerable<DGField> m_RegisteredFields) {
            if (g != null) {
                this.CellPartitioning = g.CellPartitioning.MpiSize.ForLoop(rank => g.CellPartitioning.GetLocalLength(rank));
                this.LocalNumberOfCells = g.CellPartitioning.LocalLength;
            }

            var _DgInfo = new List<(string FieldName, int Degree, bool isXdg)>();
            if(m_RegisteredFields != null) {
                foreach (var f in m_RegisteredFields) {
                    _DgInfo.Add((f.Identification, f.Basis.Degree, f.Basis is XDGBasis));
                }
            }
            DGinfo = _DgInfo.ToArray();
        }

        [DataMember]
        public (string FieldName, int Degree, bool isXdg)[] DGinfo;


        [DataMember]
        public int[] CellPartitioning;

        [JsonIgnore] 
        public long NumberOfCells => CellPartitioning.Sum();

        [DataMember]
        public int LocalNumberOfCells;


        private void WriteProfilingHeader(TextWriter stw) {
            stw.WriteLine($"Start Time  : {this.AppStartTime}");
            stw.WriteLine($"End Time    : {this.AppEndTime}");
            stw.WriteLine($"Run Time    : {this.AppRunTime}");
            stw.WriteLine($"Computer    : {this.Computer} ");
            stw.WriteLine($"User Name   : {this.UserName}");
            stw.WriteLine($"MPI rank    : {MPIRank}");
            stw.WriteLine($"MPI size    : {MPISize}");
            stw.WriteLine($"Num threads : {this.NumThreads}");
            stw.WriteLine($"GIT hash    : {this.GitCommitHash}");
            stw.WriteLine($"Command line: {this.CommandLine.ToConcatString("", " ", "")}");

            void TryWrite(string s, Func<object> o) {
                try {
                    stw.WriteLine(s + o());
                } catch (Exception e) {
                    stw.WriteLine(s + $"{e.GetType().Name}, {e.Message}");
                }

            }

            TryWrite("Number of cells (last mesh)      : ", () => this.NumberOfCells);
            TryWrite("Number of local cells (last mesh): ", () => this.LocalNumberOfCells);
            TryWrite("Cell partitioning (last mesh)    : ", () => this.CellPartitioning.ToConcatString("[", ", ", "]"));
            if (DGinfo != null) {
                foreach (var f in DGinfo) {
                    //                                         :
                    TryWrite("    Field: ", () => $"{f.FieldName}, degree {f.Degree}, XDG: {f.isXdg}");
                }
            }
            Console.WriteLine();
        }


        /// <summary>
        /// creates a human-readable performance report from the profiling information stored in <see cref="RootCall"/>.
        /// </summary>
        public void WriteProfilingReport(TextWriter wrt) {
            var R = RootCall;

            WriteProfilingHeader(wrt);

            //if (ilPSP.Environment.OpenMPenabled) {
            //    var acel = ilPSP.Environment.MeasureOpenMPAcceleration();
            //    var blk = ilPSP.Environment.CheckOMPThreading();
            //    wrt.WriteLine($"OpenMP scaling        : absolute: {acel.AbsScaling}, relative: {acel.RelScaling}");
            //    wrt.WriteLine($"Hybrid Parallelization: OpenMP blocking: {blk}");

            //} else {
            //    wrt.WriteLine("No OpenMP acceleration measured, since OpenMP is disabled.");
            //}
            //wrt.WriteLine();

            wrt.WriteLine();
            wrt.WriteLine("Micro-Benchmarking results");
            wrt.WriteLine("=========================================================");
            this.OnlinePerformanceLog.WriteStatistics(wrt);


            wrt.WriteLine();
            wrt.WriteLine("Most expensive calls and blocks (sort by exclusive time):");
            wrt.WriteLine("(sum over all calling parents)");
            wrt.WriteLine("=========================================================");

            MethodCallRecordExtension.GetMostExpensiveCalls(wrt, R);

            wrt.WriteLine();
            wrt.WriteLine("Most expensive calls and blocks (sort by exclusive time):");
            wrt.WriteLine("(distinction by parent call)");
            wrt.WriteLine("=========================================================");

            MethodCallRecordExtension.GetMostExpensiveCallsDetails(wrt, R);

            wrt.WriteLine();
            wrt.WriteLine("Most memory consuming calls and blocks (sort by exclusive allocation size):");
            wrt.WriteLine("(sum over all calling parents)");
            wrt.WriteLine("===========================================================================");

            MethodCallRecordExtension.GetMostMemoryConsumingCalls(wrt, R);

            wrt.WriteLine();
            wrt.WriteLine("Most memory consuming calls and blocks (sort by exclusive allocation size):");
            wrt.WriteLine("(distinction by parent call)");
            wrt.WriteLine("==========================================================================");

            MethodCallRecordExtension.GetMostMemoryConsumingCallsDetails(wrt, R);


            wrt.WriteLine();
            wrt.WriteLine("Environment Variables");
            wrt.WriteLine("=====================");
            foreach(var ev in this.EnvironmentVars) {
                wrt.WriteLine(ev.Name + ": " + ev.Value);
            }


            /*
            wrt.WriteLine();
            wrt.WriteLine("Details on nonlinear operator evaluation:");
            wrt.WriteLine("=========================================");

            var OpEval = R.FindChildren("BoSSS.Foundation.SpatialOperator*Evaluator*Evaluate*");
            if (OpEval.Count() == 0) {
                wrt.WriteLine("not called.");
            } else {
                try {
                    wrt.WriteLine((new CollectionReport(OpEval.ToArray())).ToString());

                    wrt.WriteLine("Blocks:");
                    wrt.WriteLine("-------");

                    List<MethodCallRecord>[] OpEval_Blocks = ((int)2).ForLoop(iii => new List<MethodCallRecord>());
                    Dictionary<string, List<MethodCallRecord>>[] QuadratureExecuteBlocks = ((int)2).ForLoop(iii => new Dictionary<string, List<MethodCallRecord>>());

                    foreach (MethodCallRecord mcr in OpEval) {
                        OpEval_Blocks[0].AddRange(mcr.FindChildren("Edge_Integration_NonLin"));
                        OpEval_Blocks[1].AddRange(mcr.FindChildren("Volume_Integration_NonLin"));
                    }

                    for (int ii = 0; ii < 2; ii++) {
                        var L = OpEval_Blocks[ii];

                        foreach (MethodCallRecord mcr in L) {
                            MethodCallRecord quadCall = mcr.FindChild("*Execute*");

                            foreach (var subBlock in quadCall.Calls.Values) {
                                List<MethodCallRecord> col;
                                if (!QuadratureExecuteBlocks[ii].TryGetValue(subBlock.Name, out col)) {
                                    col = new List<MethodCallRecord>();
                                    QuadratureExecuteBlocks[ii].Add(subBlock.Name, col);
                                }
                                col.Add(subBlock);
                            }
                        }
                    }

                    for (int ii = 0; ii < 2; ii++) {
                        wrt.WriteLine((new CollectionReport(OpEval_Blocks[ii])).ToString());
                        foreach (var col in QuadratureExecuteBlocks[ii].Values) {
                            wrt.Write("  ");
                            wrt.WriteLine((new CollectionReport(col)).ToString());
                        }
                    }


                } catch (Exception e) {
                    wrt.WriteLine(e.GetType().Name + ": " + e.Message);
                    wrt.WriteLine(e.StackTrace);
                }
            }



            wrt.WriteLine();
            wrt.WriteLine("Details on Matrix compilation:");
            wrt.WriteLine("==============================");

            var Matrix = R.FindChildren("BoSSS.Foundation.SpatialOperator*ComputeMatrix*");
            if (Matrix.Count() == 0) {
                wrt.WriteLine("not called.");
            } else {
                try {
                    wrt.WriteLine((new CollectionReport(Matrix.ToArray())).ToString());

                    wrt.WriteLine("Blocks:");
                    wrt.WriteLine("-------");

                    List<MethodCallRecord>[] Matrix_Blocks = ((int)4).ForLoop(iii => new List<MethodCallRecord>());
                    Dictionary<string, List<MethodCallRecord>>[] QuadratureExecuteBlocks = ((int)4).ForLoop(iii => new Dictionary<string, List<MethodCallRecord>>());

                    foreach (MethodCallRecord mcr in Matrix) {
                        Matrix_Blocks[0].AddRange(mcr.FindChildren("Edge_Integration_(legacy)"));
                        Matrix_Blocks[1].AddRange(mcr.FindChildren("Edge_Integration_(new)"));
                        Matrix_Blocks[2].AddRange(mcr.FindChildren("Volume_Integration_(legacy)"));
                        Matrix_Blocks[3].AddRange(mcr.FindChildren("Volume_Integration_(new)"));
                    }

                    for (int ii = 0; ii < 4; ii++) {
                        var L = Matrix_Blocks[ii];

                        foreach (MethodCallRecord mcr in L) {
                            MethodCallRecord quadCall = mcr.FindChild("*Execute*");
                            if(quadCall != null) {
                                foreach(var subBlock in quadCall.Calls.Values) {
                                    List<MethodCallRecord> col;
                                    if(!QuadratureExecuteBlocks[ii].TryGetValue(subBlock.Name, out col)) {
                                        col = new List<MethodCallRecord>();
                                        QuadratureExecuteBlocks[ii].Add(subBlock.Name, col);
                                    }
                                    col.Add(subBlock);
                                }
                            }
                        }
                    }

                    for (int ii = 0; ii < 4; ii++) {
                        if (Matrix_Blocks[ii].Any()) {
                            wrt.WriteLine((new CollectionReport(Matrix_Blocks[ii])).ToString());
                            foreach (var col in QuadratureExecuteBlocks[ii].Values) {
                                wrt.Write("  ");
                                wrt.WriteLine((new CollectionReport(col)).ToString());
                            }
                        }
                    }
                } catch (Exception e) {
                    wrt.WriteLine(e.GetType().Name + ": " + e.Message);
                    wrt.WriteLine(e.StackTrace);
                }
            }


            wrt.WriteLine();
            wrt.WriteLine("Details on XDG Matrix compilation:");
            wrt.WriteLine("==================================");

            var XMatrix = R.FindChildren("BoSSS.Foundation.XDG.XSpatialOperator*ComputeMatrix*");
            if (XMatrix.Count() == 0) {
                wrt.WriteLine("not called.");
            } else {
                try {
                    wrt.WriteLine((new CollectionReport(XMatrix.ToArray())).ToString());

                    wrt.WriteLine("Blocks (coarse):");
                    wrt.WriteLine("----------------");

                    List<MethodCallRecord> XMatrix_agglomeration = new List<MethodCallRecord>();
                    List<MethodCallRecord> XMatrix_surface_integration = new List<MethodCallRecord>();
                    List<MethodCallRecord> XMatrix_bulk_integration = new List<MethodCallRecord>();
                    List<MethodCallRecord> XMatrix_QuadRule_compilation = new List<MethodCallRecord>();

                    foreach (MethodCallRecord mcr in XMatrix) {
                        XMatrix_agglomeration.AddRange(mcr.FindChildren("agglomeration"));
                        XMatrix_surface_integration.AddRange(mcr.FindChildren("surface_integration"));
                        XMatrix_bulk_integration.AddRange(mcr.FindChildren("bulk_integration"));
                        XMatrix_QuadRule_compilation.AddRange(mcr.FindChildren("QuadRule-compilation"));
                    }

                    wrt.WriteLine((new CollectionReport(XMatrix_agglomeration)).ToString());
                    wrt.WriteLine((new CollectionReport(XMatrix_surface_integration)).ToString());
                    wrt.WriteLine((new CollectionReport(XMatrix_bulk_integration)).ToString());
                    wrt.WriteLine((new CollectionReport(XMatrix_QuadRule_compilation)).ToString());

                    wrt.WriteLine("Blocks (fine):");
                    wrt.WriteLine("--------------");

                    List<MethodCallRecord>[] Matrix_Blocks = ((int)5).ForLoop(iii => new List<MethodCallRecord>());
                    Dictionary<string, List<MethodCallRecord>>[] QuadratureExecuteBlocks = ((int)5).ForLoop(iii => new Dictionary<string, List<MethodCallRecord>>());

                    foreach (MethodCallRecord mcr in XMatrix) {
                        Matrix_Blocks[0].AddRange(mcr.FindChildren("Edge_Integration_(legacy)"));
                        Matrix_Blocks[1].AddRange(mcr.FindChildren("Edge_Integration_(new)"));
                        Matrix_Blocks[2].AddRange(mcr.FindChildren("Volume_Integration_(legacy)"));
                        Matrix_Blocks[3].AddRange(mcr.FindChildren("Volume_Integration_(new)"));
                        Matrix_Blocks[4].AddRange(mcr.FindChildren("surface_integration"));
                    }

                    for (int ii = 0; ii < 5; ii++) {
                        var L = Matrix_Blocks[ii];

                        foreach (MethodCallRecord mcr in L) {
                            MethodCallRecord quadCall = mcr.FindChild("*Execute*");

                            foreach (var subBlock in quadCall.Calls.Values) {
                                List<MethodCallRecord> col;
                                if (!QuadratureExecuteBlocks[ii].TryGetValue(subBlock.Name, out col)) {
                                    col = new List<MethodCallRecord>();
                                    QuadratureExecuteBlocks[ii].Add(subBlock.Name, col);
                                }
                                col.Add(subBlock);
                            }
                        }
                    }

                    for (int ii = 0; ii < 5; ii++) {
                        wrt.WriteLine((new CollectionReport(Matrix_Blocks[ii])).ToString());
                        foreach (var col in QuadratureExecuteBlocks[ii].Values) {
                            wrt.Write("  ");
                            wrt.WriteLine((new CollectionReport(col)).ToString());
                        }
                    }

                } catch (Exception e) {
                    wrt.WriteLine(e.GetType().Name + ": " + e.Message);
                    wrt.WriteLine(e.StackTrace);
                }
            }
            */
        }


        public override string ToString() {
            using (var tw = new StringWriter()) {
                WriteProfilingReport(tw);
                return tw.ToString();
            }
        }

    }
}
