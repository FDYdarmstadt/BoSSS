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


using BoSSS.Foundation;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.LoadBalancing;
using BoSSS.Solution.Queries;
using CommandLine;
using CommandLine.Text;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net;
using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Threading;
using System.Xml;

namespace BoSSS.Solution {

    /// <summary>
    /// Configuration of <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector"/>
    /// </summary>
    [DataContract]
    public class MatlabConnectorConfig {
        /// <summary>
        /// <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector.MatlabExecuteable"/>
        /// </summary>
        [DataMember]
        public string MatlabExecuteable;

        /// <summary>
        /// <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector.Flav"/>
        /// </summary>
        [DataMember]
        public string Flav = ilPSP.Connectors.Matlab.BatchmodeConnector.Flavor.Matlab.ToString();


        /// <summary>
        /// Loads configuration from default location in user directory
        /// </summary>
        public static MatlabConnectorConfig LoadDefault(string userDir) {
            string ConfigFile = Path.Combine(userDir, "etc", "MatlabConnectorConfig.json");
            if (!File.Exists(ConfigFile)) {
                var r = new MatlabConnectorConfig();
                r.SaveDefault(userDir);
                return r;
            }
            string str = File.ReadAllText(ConfigFile);

            return Deserialize(str);
        }

        /// <summary>
        /// Saves configuration in default location in user directory
        /// </summary>
        public void SaveDefault(string userDir) {
            string ConfigFile = Path.Combine(userDir, "etc", "MatlabConnectorConfig.json");
            var Str = this.Serialize();
            File.WriteAllText(ConfigFile, Str);
        }


        /// <summary>
        /// JSON deserialization
        /// </summary>
        public static MatlabConnectorConfig Deserialize(string Str) {


            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };


            using (var tr = new StringReader(Str)) {
                //string typeName = tr.ReadLine();
                Type ControlObjectType = typeof(MatlabConnectorConfig); //Type.GetType(typeName);
                using (JsonReader reader = new JsonTextReader(tr)) {
                    var obj = formatter.Deserialize(reader, ControlObjectType);

                    MatlabConnectorConfig ctrl = (MatlabConnectorConfig)obj;
                    return ctrl;
                }

            }
        }

        /// <summary>
        /// JSON serialization
        /// </summary>
        public string Serialize() {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error,
                Formatting = Newtonsoft.Json.Formatting.Indented
            };

            using (var tw = new StringWriter()) {
                //tw.WriteLine(this.GetType().AssemblyQualifiedName);
                using (JsonWriter writer = new JsonTextWriter(tw)) {  // Alternative: binary writer: BsonWriter
                    formatter.Serialize(writer, this);
                }

                string Ret = tw.ToString();
                return Ret;
            }
        }

    }


    /// <summary>
    /// Non-generic version of the <see cref="Application{T}"/>-class for
    /// backward compatibility.
    /// </summary>
    public abstract class Application : Application<Application.EmptyAppControl> {

        /// <summary>
        /// Dummy-implementation of the <see cref="AppControl"/>-class.
        /// </summary>
        public class EmptyAppControl : AppControl {

            /// <summary>
            /// Ctor.
            /// </summary>
            public EmptyAppControl() {
                base.savetodb = false;
            }


            /// <summary>
            /// nix supported.
            /// </summary>
            public override Type GetSolverType() {
                throw new NotSupportedException();
            }

        }
    }

    /// <summary>
    /// Base class for BoSSS applications that helps with the organization of
    /// the general work-flow and offers a simple control file handling. The
    /// standard mode of execution is defined by <see cref="RunSolverMode"/>. 
    /// The <see cref="_Main"/> method
    /// offers a convenient way to start a BoSSS application with minimal
    /// effort.
    /// </summary>
    public abstract class Application<T> : IDisposable, IApplication, IApplication<T>
        where T : AppControl, new() {

        /// <summary>
        /// Logs useful information about the execution.
        /// </summary>
        private static ILog m_Logger = LogManager.GetLogger(typeof(Application<T>));

        ///// <summary>
        ///// Indicates whether a running application must finalize MPI
        ///// </summary>
        //private static bool m_MustFinalizeMPI = false;

        /// <summary>
        /// Set this variable to false if database IO is desired, but no
        /// control file is given; Must be set before <see cref="InitMPI"/>;
        /// </summary>
        protected bool passiveIo = true;

        /// <summary>
        /// User configuration input.
        /// </summary>
        public T Control {
            get;
            private set;
        }

        /// <summary>
        /// See <see cref="Control"/>.
        /// </summary>
        public AppControl ControlBase {
            get {
                return Control;
            }
        }


        /// <summary>
        /// <see cref="QueryHandler"/>
        /// </summary>
        protected QueryHandler m_queryHandler;

        /// <summary>
        /// New 'table' to log query results
        /// </summary>
        private QueryResultTable m_QueryResultTable = new QueryResultTable();

        /// <summary>
        /// New 'table' to log query results (and other things that the application is doing)
        /// </summary>
        public QueryResultTable QueryResultTable {
            get {
                return m_QueryResultTable;
            }
        }

        /// <summary>
        /// The query handler holding all queries to be performed during a run.
        /// </summary>
        public QueryHandler QueryHandler {
            get {
                return m_queryHandler;
            }
        }

        /// <summary>
        /// respective location of native libraries
        /// </summary>
        public static string GetNativeLibraryDir() {
            return BoSSS.Foundation.IO.Utils.GetNativeLibraryDir(m_Logger);
        }




        /// <summary>
        /// Re-loads optional user-configuration file for Matlab connector <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector"/>.
        /// </summary>
        public static void ReadBatchModeConnectorConfig() {
            using (var tr = new FuncTrace()) {
                string userDir = BoSSS.Foundation.IO.Utils.GetBoSSSUserSettingsPath();
                if (userDir.IsEmptyOrWhite())
                    return;

                try {
                    var o = MatlabConnectorConfig.LoadDefault(userDir);

                    ilPSP.Connectors.Matlab.BatchmodeConnector.Flav = (ilPSP.Connectors.Matlab.BatchmodeConnector.Flavor)System.Enum.Parse(typeof(ilPSP.Connectors.Matlab.BatchmodeConnector.Flavor), o.Flav);
                    ilPSP.Connectors.Matlab.BatchmodeConnector.MatlabExecuteable = o.MatlabExecuteable;

                } catch (Exception e) {
                    var errStr = $"{e.GetType().Name} while reading/saving Matlab connector configuration file: {e.Message}";
                    tr.Logger.Error(errStr);
                    Console.Error.WriteLine(errStr);
                }
            }
        }

        /// <summary>
        /// Application startup. Performs bootstrapping of unmanaged resources and initializes MPI.
        /// This method may be called multiple times in an application lifetime -- the MPI init is only performed once.
        /// </summary>
        /// <param name="args">
        /// command line arguments
        /// </param>
        /// <returns>
        /// Whether this call actually initialized MPI
        /// - true, if this routine actually called <see cref="IMPIdriver.Init"/>; then, the call should be 
        ///   an other call to <see cref="FinalizeMPI"/>.
        /// - false, if not.
        /// </returns>
        public static bool InitMPI(string[] args = null) {
            if (args == null)
                args = new string[0];

            m_Logger.Info("In BoSSS MPI init.");

            System.Threading.Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            m_Logger.Info("Bootstrapping.");


            ilPSP.Environment.Bootstrap(
                args,
                GetNativeLibraryDir(),
                out bool _MustFinalizeMPI);

            m_Logger.Info("_MustFinalizeMPI = " + _MustFinalizeMPI);

            int rank, size;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            if (_MustFinalizeMPI) {
                if (rank == 0) {

                    Console.WriteLine(@"      ___           ___           ___           ___           ___     ");
                    Console.WriteLine(@"     /\  \         /\  \         /\  \         /\  \         /\  \    ");
                    Console.WriteLine(@"    /::\  \       /::\  \       /::\  \       /::\  \       /::\  \   ");
                    Console.WriteLine(@"   /:/\:\  \     /:/\:\  \     /:/\ \  \     /:/\ \  \     /:/\ \  \  ");
                    Console.WriteLine(@"  /::\~\:\__\   /:/  \:\  \   _\:\~\ \  \   _\:\~\ \  \   _\:\~\ \  \ ");
                    Console.WriteLine(@" /:/\:\ \:|__| /:/__/ \:\__\ /\ \:\ \ \__\ /\ \:\ \ \__\ /\ \:\ \ \__\");
                    Console.WriteLine(@" \:\~\:\/:/  / \:\  \ /:/  / \:\ \:\ \/__/ \:\ \:\ \/__/ \:\ \:\ \/__/");
                    Console.WriteLine(@"  \:\ \::/  /   \:\  /:/  /   \:\ \:\__\    \:\ \:\__\    \:\ \:\__\  ");
                    Console.WriteLine(@"   \:\/:/  /     \:\/:/  /     \:\/:/  /     \:\/:/  /     \:\/:/  /  ");
                    Console.WriteLine(@"    \::/__/       \::/  /       \::/  /       \::/  /       \::/  /   ");
                    Console.WriteLine(@"     ~~            \/__/         \/__/         \/__/         \/__/    ");
                    Console.WriteLine(@"                                                                      ");

                    Console.Write(DateTime.Now);
                    Console.Write("  ");
                    if (size <= 1)
                        Console.WriteLine("Running with 1 MPI process (single core)");
                    else
                        Console.WriteLine("Running with " + size + " MPI processes ");

                    Console.WriteLine("User: " + System.Environment.UserName);

                    using (var stw = new StringWriter()) {
                        var HostName = ilPSP.Environment.MPIEnv.AllHostNames;
                        int I = HostName.Count;
                        if (I > 1)
                            stw.Write("Nodes: ");
                        else
                            stw.Write("Node: ");

                        int i = 0;
                        foreach (var kv in HostName) {
                            stw.Write(kv.Key);
                            if (kv.Value.Length == 1)
                                stw.Write(" (rank ");
                            else
                                stw.Write(" (ranks ");

                            int J = kv.Value.Length;
                            for (int j = 0; j < J; j++) {
                                stw.Write(kv.Value[j]);
                                if (j < J - 1)
                                    stw.Write(", ");
                            }
                            stw.Write(")");

                            int rest = I - i - 1;
                            if (i >= 1023 && rest >= 1) {
                                stw.Write($" ... and {rest} ,more ... ");
                                break;
                            }

                            if (i < I - 1)
                                stw.Write(", ");

                            i++;

                        }


                        stw.Flush();
                        Console.WriteLine(stw.ToString());
                    }

                }
            }

            m_Logger.Info("Hello from MPI rank " + rank + " of " + size + "!");

            ReadBatchModeConnectorConfig();

            System.Threading.Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            return _MustFinalizeMPI;
        }

        /// <summary>
        /// Recursive collection of all dependencies of some assembly.
        /// </summary>
        public static IEnumerable<Assembly> GetAllAssemblies(Type EntryType) {
            StackTrace stackTrace = new StackTrace();
            List<Assembly> allAssis = new List<Assembly>();

            if (EntryType == null) {
                // Start from entry assembly and _not_
                // - don't use typeof(T).Assembly since T might be located different assembly than the control file
                // - don't use Assembly.GetEntryAssembly() as it is undefined if called by Nunit
                // - exclude mscorlib (does not even appear on Windows, but makes problems using mono)

                for (int i = 0; i < stackTrace.FrameCount; i++) {
                    Type Tuepe = stackTrace.GetFrame(i).GetMethod().DeclaringType;
                    if (Tuepe != null) {
                        Assembly entryAssembly = Tuepe.Assembly;
                        GetAllAssembliesRecursive(entryAssembly, allAssis);
                    }
                }
            } else {
                GetAllAssembliesRecursive(EntryType.Assembly, allAssis);
            }


            return allAssis.Where(a => !a.FullName.StartsWith("mscorlib")).Where(a => a.IsDynamic == false).ToArray();
        }

        /// <summary>
        /// Recursive collection of all dependencies of some assembly.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="assiList">
        /// Output, list where all dependent assemblies are collected.
        /// </param>
        static void GetAllAssembliesRecursive(Assembly a, List<Assembly> assiList) {
            if (a.GlobalAssemblyCache)
                return;


            if (assiList.Contains(a))
                return;



            assiList.Add(a);

            foreach (var _b in a.GetReferencedAssemblies()) {
                Assembly b;
                try {
                    b = Assembly.Load(_b);
                } catch (Exception) {
                    continue;
                }

                if (b.GlobalAssemblyCache)
                    continue;
                GetAllAssembliesRecursive(b, assiList);
            }
        }

        /// <summary>
        /// parses command line arguments, parses control file, runs the
        /// application
        /// </summary>
        /// <param name="args">command line arguments</param>
        /// <param name="noControlFile"></param>
        /// <param name="ApplicationFactory">
        /// A factory that returns the application to be run.
        /// </param>
        public static void _Main(
            string[] args,
            bool noControlFile,
            Func<Application<T>> ApplicationFactory) {

            m_Logger.Info("Entering _Main routine...");
            Tracer.MemoryInstrumentationLevel = MemoryInstrumentationLevel.None;

            int MPIrank = int.MinValue;
            //#if DEBUG
            //            {


            //#else
            //            try {
            //#endif
            CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
            CultureInfo.DefaultThreadCurrentUICulture = CultureInfo.InvariantCulture;

            bool _MustFinalizeMPI = InitMPI(args);
            ReadBatchModeConnectorConfig();

            MPIrank = ilPSP.Environment.MPIEnv.MPI_Rank;

            // lets see if we have environment variables which override command line arguments
            // (environment variables are usually more robust w.r.t. e.g. escape characters)
            args = ArgsFromEnvironmentVars(args);


            // parse arguments
            CommandLineOptions opt = new CommandLineOptions();
            //ICommandLineParser parser = new CommandLine.CommandLineParser(new CommandLineParserSettings(Console.Error));
            
            bool argsParseSuccess;
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var CmdlineParseRes = Parser.Default.ParseArguments<CommandLineOptions>(args);
                opt = CmdlineParseRes.Value;
                argsParseSuccess = CmdlineParseRes.Errors.IsNullOrEmpty();
                //argsParseSuccess = parser.ParseArguments(args, opt);
                argsParseSuccess = argsParseSuccess.MPIBroadcast<bool>(0, csMPI.Raw._COMM.WORLD);
            } else {
                argsParseSuccess = false;
                argsParseSuccess = argsParseSuccess.MPIBroadcast<bool>(0, csMPI.Raw._COMM.WORLD);
            }

            if (!argsParseSuccess) {
                MPI.Wrappers.csMPI.Raw.mpiFinalize();
                _MustFinalizeMPI = false;
                m_Logger.Error("Unable to parse arguments - exiting.");
                System.Environment.Exit(-1);
            }

            if (opt.ControlfilePath != null) {
                opt.ControlfilePath = opt.ControlfilePath.Trim();
            }
            m_Logger.Info("Braodcasting command line options..");
            opt = opt.MPIBroadcast(0, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
            m_Logger.Info("done.");

            // Delete old plots if requested
            if (opt.delPlt) {
                m_Logger.Info("Deletion of old plot files...");
                DeleteOldPlotFiles();
                m_Logger.Info("done.");
            }


            // load control file
            T ctrlV2 = null;
            T[] ctrlV2_ParameterStudy = null;
            if (!noControlFile) {
                m_Logger.Info("Loading control object...");
                LoadControlFile(opt.ControlfilePath, out ctrlV2, out ctrlV2_ParameterStudy);
                m_Logger.Info("done.");
            } else {
                ctrlV2 = new T();
                m_Logger.Info("No control Object.");
            }

            AppEntry(ApplicationFactory, opt, ctrlV2, ctrlV2_ParameterStudy);

            if (_MustFinalizeMPI)
                FinalizeMPI();

            //#if DEBUG
            //            }
            //#else
            //            } catch (Exception e) {
            //                // handle exception logging, for automated processing in WorkflowMgm etc.
            //                // make sure the exception is appended to the stderr in session directory!

            //                Console.Error.WriteLine(e.StackTrace);
            //                Console.Error.WriteLine();
            //                Console.Error.WriteLine();
            //                Console.Error.WriteLine("========================================");
            //                Console.Error.WriteLine("========================================");
            //                Console.Error.WriteLine($"MPI rank {MPIrank}: {e.GetType().Name } :");
            //                Console.Error.WriteLine(e.Message);
            //                Console.Error.WriteLine("========================================");
            //                Console.Error.WriteLine("========================================");
            //                Console.Error.WriteLine();
            //                Console.Error.Flush();
            //                //System.Environment.Exit(-1);
            //            }
            //#endif
        }

        /// <summary>
        /// Appends specially named environment variables (BOSSS_ARG_n) as the n-th command line argument
        /// </summary>
        public static string[] ArgsFromEnvironmentVars(string[] args) {
            {
                List<string> _args = new List<string>(args);
                int ArgCounter = 0;
                while (true) {
                    string ArgOverrideName = "BOSSS_ARG_" + ArgCounter;
                    string ArgValue = System.Environment.GetEnvironmentVariable(ArgOverrideName);
                    if (ArgValue == null)
                        break;

                    System.Environment.SetEnvironmentVariable(ArgOverrideName, null); // delete the envvar
                    // many test internally call the _Main function with arguments;
                    // this would be overridden (and thus not work properly) if we don't delete the variable here and now.

                    if (ArgCounter < _args.Count) {
                        _args[ArgCounter] = ArgValue;
                    } else {
                        _args.Add(ArgValue);
                    }

                    Console.WriteLine("arg #{0} override from environment variable '{1}': {2}", ArgCounter, ArgOverrideName, ArgValue);

                    ArgCounter++;
                }
                args = _args.ToArray();
            }

            return args;
        }

        /// <summary>
        /// Loads the control object from the given file path
        /// </summary>
        /// <param name="ControlFilePath">
        /// This could be either
        ///  - the path to a cs - script, or an object file
        ///  - an instruction, starting with cs:, which refers to an internal function which generates the control object
        /// </param>
        /// <param name="ctrlV2">
        /// a single control object, in the case of a normal solver run
        /// </param>
        /// <param name="ctrlV2_ParameterStudy">
        /// a list of control objects, in the case of a parameter study
        /// </param>
        public static void LoadControlFile(string ControlFilePath, out T ctrlV2, out T[] ctrlV2_ParameterStudy) {
            ctrlV2 = null;
            ctrlV2_ParameterStudy = null;

            if (ControlFilePath.IsNullOrEmpty()) {
                throw new Exception("No control file specified.");
            }

            if (ControlFilePath.ToLower().StartsWith("cs:")) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++
                // C#-instruction provided as argument 
                // usually redicts to some pre-compiled static function
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++


                var StringwithoutPrefix = ControlFilePath.Substring(3);



                ControlObjFromCode(StringwithoutPrefix, out ctrlV2, out ctrlV2_ParameterStudy);

            } else if (ControlFilePath.ToLower().EndsWith(".cs") || ControlFilePath.ToLower().EndsWith(".bws")) {
                // +++++++++
                // C#-script
                // +++++++++
                string ctrlfileContent = "";
                if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {

                    if (ControlFilePath != null) {
                        StreamReader rd = new StreamReader(ControlFilePath);
                        // Trim trailing empty lines since they would
                        // return an empty string instead of a control object
                        ctrlfileContent = rd.ReadToEnd().TrimEnd();
                        rd.Close();
                    }

                }

                ctrlfileContent = ctrlfileContent.MPIBroadcast(0, csMPI.Raw._COMM.WORLD);
                ControlObjFromCode(ctrlfileContent, out ctrlV2, out ctrlV2_ParameterStudy);
            } else if (ControlFilePath.ToLower().EndsWith(".obj")) {
                // +++++++++++++++++++++
                // control object
                // +++++++++++++++++++++

            
                string JSON = File.ReadAllText(ControlFilePath);
                object controlObj = AppControl.Deserialize(JSON);
                ctrlV2 = controlObj as T;


                if (ctrlV2 == null) {
                    throw new ApplicationException(string.Format(
                        "Invalid control instruction: unable to cast the last result of the control file/cs-script of type {0} to type {1}",
                        controlObj.GetType().FullName,
                        typeof(T).FullName));
                }

            } else {
                throw new ArgumentException("unable to interpret: " + ControlFilePath);
            }
        }

        /// <summary>
        /// On process rank 0, deletes all plot files in the current directory
        /// </summary>
        public static void DeleteOldPlotFiles() {
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                Console.Write("rm");
                foreach (var pltFile in dir.GetFiles("*.plt").Concat(dir.GetFiles("*.curve"))) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }
        }
        /// <summary>
        /// On process rank 0, deletes all txt and csv files in the current directory
        /// </summary>
        public static void DeleteOldTextFiles() {
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                Console.Write("rm");
                foreach (var pltFile in dir.GetFiles("*.txt").Concat(dir.GetFiles("*.csv"))) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }
        }

        /// <summary>
        /// Loads a control object, resp. a series of control objects (in the case of a parameter study)
        /// form a C#-script.
        /// </summary>
        /// <param name="ctrlfileContent">the script.</param>
        /// <param name="ctrl">output, for a singe control object.</param>
        /// <param name="ctrl_ParameterStudy">output, for a series of control objects.</param>
        static void ControlObjFromCode(string ctrlfileContent, out T ctrl, out T[] ctrl_ParameterStudy) {
            //object controlObj = ControlObjFromCode(ctrlfileContent, typeof(T));
            AppControl.FromCode(ctrlfileContent, typeof(T), out AppControl _ctrl, out AppControl[] _ctrl_ParamStudy);

            Debug.Assert((_ctrl == null) != (_ctrl_ParamStudy == null));

            if (_ctrl != null) {
                ctrl = (T)_ctrl;

                ctrl_ParameterStudy = null;
            } else if (_ctrl_ParamStudy != null) {
                //ctrl_ParameterStudy = ((IEnumerable<T>)controlObj).ToArray();
                //foreach (var c in ctrl_ParameterStudy)
                //    c.ControlFileText = ctrlfileContent;

                ctrl_ParameterStudy = new T[_ctrl_ParamStudy.Length];
                for (int i = 0; i < ctrl_ParameterStudy.Length; i++) {
                    ctrl_ParameterStudy[i] = (T)(_ctrl_ParamStudy[i]);
                }

                ctrl = null;
            } else {
                //throw new ApplicationException(string.Format(
                //"Invalid control instruction: unable to cast the last result of the control file/cs-script of type {0} to type {1} or IEnumerable<{1}>",
                //controlObj.GetType().FullName,
                //typeof(T).FullName));
                throw new ApplicationException();
            }
        }


        private static void AppEntry(
            Func<Application<T>> ApplicationFactory,
            CommandLineOptions opt, T ctrlV2, T[] ctrlV2_ParameterStudy) {

            if (opt == null)
                throw new ArgumentNullException();

            // run application
            // ===============
            if (ctrlV2_ParameterStudy != null) {

                ParameterStudyModeV2(opt, ctrlV2_ParameterStudy, ApplicationFactory);
            } else if (ctrlV2 != null) {

                // control file overrides from command line
                // ----------------------------------------

                // load control file, parse args
                if (opt.ProjectName != null)
                    ctrlV2.ProjectName = opt.ProjectName;
                if (opt.SessionName != null && ctrlV2.SessionName.IsEmptyOrWhite())
                    ctrlV2.SessionName = opt.SessionName;

                if (opt.ImmediatePlotPeriod != null) {
                    ctrlV2.ImmediatePlotPeriod = opt.ImmediatePlotPeriod.Value;
                }
                if (opt.SuperSampling != null) {
                    ctrlV2.SuperSampling = opt.SuperSampling.Value;
                }
                if (opt.fullTracing) {
                    ctrlV2.TracingNamespaces = "*";
                }

                // ad-hoc added tags (added by command-line option)
                if (opt != null && (opt.TagsToAdd != null && opt.TagsToAdd.Length > 0)) {
                    string[] AdHocTags;
                    AdHocTags = opt.TagsToAdd.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    ctrlV2.Tags.AddRange(AdHocTags);
                }


                // run app
                // -------

                using (Application<T> app = ApplicationFactory()) {
                    m_Logger.Info("Running application...");                    
                    app.Init(ctrlV2);
                    app.RunSolverMode();
                    m_Logger.Info("Application finished.");
                }
            } else {
                // no control file 


                throw new ArgumentException();
            }
        }

        /// <summary>
        /// the very end of any BoSSS application.
        /// </summary>
        public static void FinalizeMPI() {
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// Initializes the environment of the application
        /// </summary>
        /// <param name="control">
        /// control object
        /// </param>
        public virtual void Init(AppControl control) {

            if (control != null) {
                this.Control = (T)control;
            } else {
                this.Control = new T();
            }

            ReadBatchModeConnectorConfig();

            // set . as decimal separator:
            // ===========================
            CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
            CultureInfo.DefaultThreadCurrentUICulture = CultureInfo.InvariantCulture;


            // set a few switches
            // ==================
            if (this.Control != null) {
                this.Control.Verify();



                if (this.Control.NoOfTimesteps >= 0) {
                    this.NoOfTimesteps = this.Control.NoOfTimesteps;
                }

                if (this.Control.Endtime >= 0) {
                    this.EndTime = this.Control.Endtime;
                }

                if (this.Control.saveperiod > 0) {
                    this.SavePeriod = this.Control.saveperiod;
                }

                if (this.Control.rollingSaves == true) {
                    this.RollingSave = this.Control.rollingSaves;
                }
            }
        }

        /// <summary>
        /// Size of the MPI world communicator
        /// </summary>
        public int MPISize {
            get {
                int size;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                return size;
            }
        }


        /// <summary>
        /// rank of this process within the MPI world communicator
        /// </summary>
        public int MPIRank {
            get {
                int rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                return rank;
            }
        }

        /// <summary>
        /// Prints a text summary of the control object to <paramref name="txt"/>.
        /// </summary>
        static void PrintObject(TextWriter txt, object obj, string spaces = "  ") {
            if (obj == null)
                txt.WriteLine(" null;");

            if (spaces.Length > 20)
                return;

            string desc = obj.ToString();
            desc = desc.Replace(txt.NewLine, txt.NewLine + "//" + spaces + "     ");
            txt.WriteLine(desc);

            Type objT = obj.GetType();
            if (objT.IsPrimitive || objT == typeof(string)) {
                return;
            }

            if (objT.IsSubclassOf(typeof(System.Delegate))) {
                return;
            }


            if (obj is System.Collections.IEnumerable) {
                System.Collections.IEnumerable objEnu = (System.Collections.IEnumerable)obj;
                int cnt = 0;
                foreach (var objE in objEnu) {
                    txt.Write("//" + spaces + "   ");
                    txt.Write("[{0}]: ", cnt);
                    cnt++;
                    PrintObject(txt, objE, spaces + "      ");
                }
                return;
            }

            BindingFlags biFlags = BindingFlags.FlattenHierarchy | BindingFlags.Instance | BindingFlags.Public | BindingFlags.GetProperty;

            PropertyInfo[] PIs = objT.GetProperties(biFlags);
            foreach (var pi in PIs) {
                object piVal = pi.GetValue(obj, biFlags, null, null, null);
                txt.Write("//");
                txt.Write(spaces);
                txt.Write(pi.Name);
                txt.Write(": ");

                if (objT.IsSubclassOf(typeof(AppControl)) && pi.Name == "ControlFileText") {
                    txt.WriteLine("...");
                    //txt.Write("//");
                    //txt.Write(spaces);
                } else {
                    if (piVal != null) {
                        PrintObject(txt, piVal, spaces + "  ");
                    } else {
                        txt.WriteLine(" null;");
                    }
                }
            }

            FieldInfo[] FIs = objT.GetFields(biFlags);
            foreach (var fi in FIs) {
                object fiVal = fi.GetValue(obj);
                txt.Write("//");
                txt.Write(spaces);
                txt.Write(fi.Name);
                txt.Write(": ");
                if (fiVal != null) {
                    PrintObject(txt, fiVal, spaces + "  ");
                } else {
                    txt.WriteLine(" null;");
                }
            }
        }


        static private bool TestSerializibility(object o) {
            try {
                using (var ms = new MemoryStream()) {
                    ms.Position = 0;


                    JsonSerializer formi = new JsonSerializer() {
                        NullValueHandling = NullValueHandling.Ignore,
                        TypeNameHandling = TypeNameHandling.Auto,
                        ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                        ReferenceLoopHandling = ReferenceLoopHandling.Ignore

                    };

                    JsonWriter jWrt = new JsonTextWriter(new StreamWriter(ms));
                    formi.Serialize(jWrt, o);
                    jWrt.Flush();

                    ms.Position = 0;
                    JsonReader jRed = new JsonTextReader(new StreamReader(ms));
                    object backAgain = formi.Deserialize(jRed, o.GetType());

                    return (o.Equals(backAgain));
                }
            } catch (Exception) {
                return false;
            }
        }


        /// <summary>
        /// Generates key/value pairs from control objects to identify sessions.
        /// </summary>
        public static void FindKeys(IDictionary<string, object> Keys, AppControl ctrl) {
            foreach (var fldOpt in ctrl.FieldOptions) {
                string KeyName = "DGdegree:" + fldOpt.Key;
                int FldDeg = fldOpt.Value.Degree;

                if (!Keys.ContainsKey(KeyName))
                    Keys.Add(KeyName, FldDeg);
            }

            foreach (var bc in ctrl.BoundaryValues) {
                string EdgeTagName = bc.Key;
                AppControl.BoundaryValueCollection bcCol = bc.Value;
                string bcType = bcCol.type;
                string KeyName = "Bndtype:" + EdgeTagName;

                if (bcType == null)
                    bcType = "auto:" + EdgeTagName;

                if (!Keys.ContainsKey(KeyName))
                    Keys.Add(KeyName, bcType);
            }

            if (ctrl.Paramstudy_CaseIdentification != null) {
                foreach (var kv in ctrl.Paramstudy_CaseIdentification) {
                    string name = kv.Item1;
                    object value = kv.Item2;

                    name = "id:" + name;
                    if (!Keys.ContainsKey(name) && TestSerializibility(value))
                        Keys.Add(name, value);
                }
            }

            FindKeysRecursive(Keys, ctrl, 0, "");

            foreach (string k in Keys.Keys.ToArray()) {
                if (Keys[k] == null)
                    Keys.Remove(k);
            }
        }

        /// <summary>
        /// Generates key/value pairs to identify sessions.
        /// </summary>
        static void FindKeysRecursive(IDictionary<string, object> Keys, object obj, int RecursionDepth, string KeyName = "") {
            if (obj == null)
                return;

            if (RecursionDepth > 20)
                return;

            Type objT = obj.GetType();
            if ((objT.IsPrimitive || objT.IsEnum || objT == typeof(string))) {
                if (!Keys.ContainsKey(KeyName) && TestSerializibility(obj))
                    Keys.Add(KeyName, obj);
                return;
            }

            if (objT.IsSubclassOf(typeof(System.Delegate))) {
                // unable to log delegates
                return;
            }


            if (obj is System.Collections.IEnumerable) {
                System.Collections.IEnumerable objEnu = (System.Collections.IEnumerable)obj;
                int cnt = 0;
                foreach (var objE in objEnu) {
                    //txt.Write("//" + spaces + "   ");
                    //txt.Write("[{0}]: ", cnt);
                    //cnt++;
                    //PrintObject(txt, objE, spaces + "      ");
                    string KeyNameIdx = KeyName + "[" + cnt + "]";
                    FindKeysRecursive(Keys, objE, RecursionDepth + 1, KeyNameIdx);
                    cnt++;
                }
                return;
            }

            BindingFlags biFlags = BindingFlags.FlattenHierarchy | BindingFlags.Instance | BindingFlags.Public | BindingFlags.GetProperty;

            MemberInfo[] PIs = objT.GetProperties(biFlags);
            MemberInfo[] FIs = objT.GetFields(biFlags);

            foreach (MemberInfo mi in ArrayTools.Cat(PIs, FIs)) {
                object Val;
                if (mi is PropertyInfo) {
                    PropertyInfo pi = ((PropertyInfo)mi);
                    if (!pi.CanRead)
                        continue;
                    if (pi.GetIndexParameters() != null && pi.GetIndexParameters().Length > 0)
                        // no suport for indexed properties.
                        continue;

                    //pi.GetIndexParameters
                    try {
                        Val = pi.GetValue(obj, biFlags, null, null, null);
                    } catch (TargetInvocationException) {
                        Val = null;
                    }
                } else if (mi is FieldInfo) {
                    Val = ((FieldInfo)mi).GetValue(obj);
                } else {
                    continue;
                }
                if (objT.IsSubclassOf(typeof(AppControl)) &&
                    (mi.Name == "ControlFileText"
                    || mi.Name == "DbPath"
                    || mi.Name == "Paramstudy_CaseIdentification"
                    || mi.Name == "ProjectDescription"
                    || mi.Name == "Tags"
                    || mi.Name == "FieldOptions"
                    || mi.Name == "InitialValues"
                    || mi.Name == "BoundaryValues"
                    || mi.Name == "InitialValues_Evaluators"
                    || mi.Name == "m_Grid")) {


                    // these guys are filtered...

                } else {
                    if (Val != null) {
                        string txt;
                        if (KeyName.Length <= 0)
                            txt = mi.Name;
                        else
                            txt = KeyName + "." + mi.Name;
                        FindKeysRecursive(Keys, Val, RecursionDepth + 1, txt);

                        //PrintObject(txt, piVal, spaces + "  ");
                    } else {
                        //txt.WriteLine(" null;");
                    }
                }
            }
        }


        /// <summary>
        /// Sets up the environment for the current run of the application.
        /// Depending on whether IO is active, this includes the following
        /// steps
        /// <list type="number">
        ///     <item>Opening a database <see cref="GetDatabase"/></item>
        ///     <item>Loading a grid via <see cref="CreateOrLoadGrid"/></item>
        ///     <item>Initializing <see cref="GridData"/></item>
        ///     <item>Creating fields via <see cref="CreateFields"/></item>
        ///     <item>
        ///     Loading the query handler via <see cref="QueryHandlerFactory"/>
        ///     </item>
        /// </list>
        /// </summary>
        virtual protected void SetUpEnvironment() {
            // init database
            // =============

            m_Logger.Info("Try to open database...");
            m_Database = GetDatabase();
            if (m_Database == null)
                m_Logger.Info("Database is NULL.");
            else
                m_Logger.Info("Opened database: " + m_Database.ToString());

            if (this.Control != null) {
                this.passiveIo = !this.Control.savetodb;
                if (this.DatabaseDriver.FsDriver is NullFileSystemDriver)
                    this.passiveIo = true;

                if (this.Control.TracingNamespaces == null) {
                    Tracer.NamespacesToLog = new string[0];
                } else {
                    var NamespacesToLog = this.Control.TracingNamespaces.Split(new char[] { ',', ' ', '\n', '\t', '\r' }, StringSplitOptions.RemoveEmptyEntries);

                    if (NamespacesToLog.Any(s => s.Equals("*"))) {
                        Tracer.NamespacesToLog = new string[] { "" };
                        m_Logger.Info("Logging ALL namespaces / full tracing");
                    } else {
                        Tracer.NamespacesToLog = NamespacesToLog;
                        m_Logger.Info("Logging " + NamespacesToLog.Length + " namespaces.");
                        foreach (var n in NamespacesToLog) {
                            m_Logger.Info("Logging namespace: " + n);
                        }
                    }
                }

                Tracer.MemoryInstrumentationLevel = this.Control.MemoryInstrumentationLevel;
            } else {
                this.passiveIo = true;
            }
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            if (!passiveIo) {
                CurrentSessionInfo = m_Database.Controller.DBDriver.CreateNewSession(m_Database);

                if (this.Control != null) {
                    CurrentSessionInfo.Name = String.IsNullOrWhiteSpace(this.Control.SessionName) ? "empty-session-name" : this.Control.SessionName;
                    CurrentSessionInfo.ProjectName = String.IsNullOrWhiteSpace(this.Control.ProjectName) ? "empty-project-name" : this.Control.ProjectName;

                    CurrentSessionInfo.Description = this.Control.ProjectDescription;

                    if (this.Control.ProjectName != null)
                        CurrentSessionInfo.KeysAndQueries.Add(PROJECTNAME_KEY, this.Control.ProjectName);
                    if (this.Control.SessionName != null)
                        CurrentSessionInfo.KeysAndQueries.Add(SESSIONNAME_KEY, this.Control.SessionName);
                }

                //System.Environment.CurrentDirectory
                //this.GetType().Assembly.Location

                //TODO: tags, projectname, ...
            } else {
                CurrentSessionInfo = new SessionInfo(Guid.Empty, this.m_Database);
            }
            InitCurrentSessionInfo();
            {
                var tags = CurrentSessionInfo.Tags.ToList();
                tags.Add(SessionInfo.NOT_TERMINATED_TAG);
                CurrentSessionInfo.Tags = tags;
            }


            this.DatabaseDriver.InitTraceFile(CurrentSessionInfo);
            using (var ht = new FuncTrace()) {
                // create or load grid
                //====================

                Grid = CreateOrLoadGrid();
                if (Grid == null) {
                    throw new ApplicationException("No grid loaded through CreateOrLoadGrid");
                }
                // access the GridData object here, to enforce its creation:
                ht.Info("loaded grid with " + this.GridData.CellPartitioning.TotalLength + " cells");

               

                bool DoDbLogging = !passiveIo
                    && this.Control != null
                    && MPIRank == 0
                    && DatabaseDriver != null
                    && (!this.CurrentSessionInfo.ID.Equals(Guid.Empty));
                if (DoDbLogging && this.Control != null) {
                    //TextWriter tw = DatabaseDriver.FsDriver.GetNewLog("Control", this.CurrentSessionInfo.ID);
                    //if (this.Control.ControlFileText != null)
                    //    tw.WriteLine(this.Control.ControlFileText);
                    //else
                    //    tw.WriteLine("// (no string representation of control object available.)");
                    //tw.WriteLine();
                    //tw.WriteLine("///////////////////////////////////////////");
                    //tw.WriteLine("// SUMMARY");
                    //tw.WriteLine("///////////////////////////////////////////");
                    //tw.Write("//");
                    //try {
                    //    PrintObject(tw, this.Control, "  ");
                    //} catch (Exception exc) {
                    //    tw.Flush();
                    //    tw.WriteLine("//" + exc.GetType().FullName + ": '" + exc.Message);
                    //}
                    //tw.Close();

                    if (this.Control.GeneratedFromCode) {
                        using (var tw = DatabaseDriver.FsDriver.GetNewLog("Control-script", this.CurrentSessionInfo.ID)) {
                            tw.WriteLine("//" + this.Control.GetType().AssemblyQualifiedName);
                            tw.Write(this.Control.ControlFileText);
                            tw.Close();
                        }
                    } else {
                        using (var tw = DatabaseDriver.FsDriver.GetNewLog("Control-obj", this.CurrentSessionInfo.ID)) {
                            tw.Write(this.Control.Serialize());
                            tw.Close();
                        }
                    }



                    Dictionary<string, object> KV = new Dictionary<string, object>();
                    FindKeys(KV, this.Control);

                    foreach (var kv in KV) {
                        try {
                            if (!this.CurrentSessionInfo.KeysAndQueries.ContainsKey(kv.Key))
                                this.CurrentSessionInfo.KeysAndQueries.Add(kv.Key, kv.Value);
                        } catch (Exception e) {
                            Console.WriteLine("Exception " + e.GetType().Name + " during logging of key/value [" + kv.Key + "," + kv.Value + "].");
                        }

                    }
                }

                // Basic session init
                if (DatabaseDriver != null) {
                    string DBpath = this.CurrentSessionInfo.Database.Path;
                    if (DBpath.Length <= 0)
                        DBpath = "EMPTY";
                    if (DBpath == null)
                        DBpath = "NULL";
                    if (DatabaseDriver.MyRank == 0) {
                        Console.WriteLine("Session ID: {0}, DB path: '{1}'.", this.CurrentSessionInfo.ID.ToString(), DBpath);
                    }
                } else {
                    Console.WriteLine("IO deactivated.");
                }

                // kernel setup
                //====================
               {
                    Grid.Redistribute(DatabaseDriver, Control.GridPartType, Control.GridPartOptions);
                    if (!passiveIo && !DatabaseDriver.GridExists(Grid.ID)) {

                        DatabaseDriver.SaveGrid(this.Grid, this.m_Database);
                        //DatabaseDriver.SaveGridIfUnique(ref _grid, out GridReplaced, this.m_Database);
                    }


                    if (this.Control == null || this.Control.NoOfMultigridLevels > 0) {
                        this.MultigridSequence = CoarseningAlgorithms.CreateSequence(this.GridData, MaxDepth: (this.Control != null ? this.Control.NoOfMultigridLevels : 1));
                    } else {
                        this.MultigridSequence = new AggregationGridData[0];
                    }
                }
                
               
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                if (this.CurrentSessionInfo != null) {
                    if (!this.CurrentSessionInfo.KeysAndQueries.ContainsKey("Grid:NoOfCells"))
                        this.CurrentSessionInfo.KeysAndQueries.Add("Grid:NoOfCells", Grid.NumberOfCells);
                    if (!this.CurrentSessionInfo.KeysAndQueries.ContainsKey("Grid:SpatialDimension"))
                        this.CurrentSessionInfo.KeysAndQueries.Add("Grid:SpatialDimension", GridData.SpatialDimension);

                    try //ToDo
                    {
                        if (!this.CurrentSessionInfo.KeysAndQueries.ContainsKey("Grid:hMax"))
                            this.CurrentSessionInfo.KeysAndQueries.Add("Grid:hMax", ((GridData)GridData).Cells.h_maxGlobal);
                        if (!this.CurrentSessionInfo.KeysAndQueries.ContainsKey("Grid:hMin"))
                            this.CurrentSessionInfo.KeysAndQueries.Add("Grid:hMin", ((GridData)GridData).Cells.h_minGlobal);
                    } catch (InvalidCastException e) {
                        Console.WriteLine("Error: Could not log everything.\n {0}", e);
                    }
                }

                // initialize solver residual logger
                m_ResLogger = new ResidualLogger(this.MPIRank, this.DatabaseDriver, this.CurrentSessionInfo.ID);

                // Make sure everything that is loaded from disk uses this grid
                // data object (if it corresponds to the same grid)
                m_Database.Controller.AddGridInitializationContext(GridData);

                // create fields
                // =============
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                if (this.Control != null) {
                    InitFromAttributes.CreateFieldsAuto(
                        this, GridData, this.Control.FieldOptions, this.Control.CutCellQuadratureType, this.m_IOFields, this.m_RegisteredFields);
                }
                CreateTracker();
                using(new BlockTrace("CreateFieldsBlock", ht)) {
                    CreateFields(); // full user control
                }



                // load queries from control file
                //========================
                m_queryHandler = QueryHandlerFactory();

                if (this.Control != null && this.Control.Queries != null) {
                    foreach (var queryIdPair in this.Control.Queries) {
                        m_queryHandler.AddQuery(queryIdPair.Key, queryIdPair.Value);
                    }
                }
                //this.QueryHandler.ValueQuery("UsedNoOfMultigridLevels", this.MultigridSequence.Length, true);

                // logging
                // ==================================
                m_PostprocessingModules.Clear();
                if (this.Control != null && this.Control.PostprocessingModules != null) {
                    m_PostprocessingModules.AddRange(this.Control.PostprocessingModules);
                }

                // save session information
                // ========================
                if (DatabaseDriver.FsDriver != null
                    && !this.CurrentSessionInfo.ID.Equals(Guid.Empty)) {
                    this.CurrentSessionInfo.Save();
                }
            }
        }

        List<InSituPostProcessingModule> m_PostprocessingModules = new List<InSituPostProcessingModule>();

        /// <summary>
        /// <see cref="Control.AppControl.PostprocessingModules"/>
        /// </summary>
        public IList<InSituPostProcessingModule> PostprocessingModules {
            get {
                return m_PostprocessingModules;
            }
        }



        /// <summary>
        /// Information about the currently active session.
        /// </summary>
        public SessionInfo CurrentSessionInfo {
            get;
            private set;
        }

        /// <summary>
        /// Initializes <see cref="CurrentSessionInfo"/>
        /// </summary>
        private void InitCurrentSessionInfo() {
            if (MPIRank == 0) {
                // set info that comes from the user
                if (this.Control != null) {
                    CurrentSessionInfo.Description = this.Control.ProjectDescription;
                    CurrentSessionInfo.Tags = this.Control.Tags.ToArray();
                }

                // set master git commit
                //CurrentSessionInfo.MasterGitCommit = Properties.Resources.MasterGitCommit;
                CurrentSessionInfo.MasterGitCommit = ((AssemblyInformationalVersionAttribute)
                  (Assembly.GetAssembly(typeof(BoSSS.Solution.Application))
                  .GetCustomAttributes(typeof(AssemblyInformationalVersionAttribute), false)[0]))
                  .InformationalVersion;

                // set deploy directory path
                string path = typeof(BoSSS.Solution.Application).Assembly.Location;
                string outpath = Path.GetDirectoryName(path); // skip "\BoSSS.Solution.dll" at end of path

                CurrentSessionInfo.DeployPath = outpath;

                // set computeNode - names:
                CurrentSessionInfo.ComputeNodeNames.Clear();
                CurrentSessionInfo.ComputeNodeNames.AddRange(ilPSP.Environment.MPIEnv.HostnameForRank);

                // save
                this.CurrentSessionInfo.Save();
            }
        }

        /// <summary>
        /// Returns a "lightweight" object that stores information about the
        /// current database. Override this method to create a different type
        /// of database; Especially, if you want to use IO but don't use the
        /// control file, you must override this method to return a valid
        /// <see cref="IDatabaseInfo"/>;
        /// </summary>
        /// <returns></returns>
        protected virtual IDatabaseInfo GetDatabase() {
            if (this.Control == null || (this.Control.DbPath.IsNullOrEmpty() && (this.Control.AlternateDbPaths == null || this.Control.AlternateDbPaths.Length <= 0))) {
                return NullDatabaseInfo.Instance;
            } else {
                return DatabaseInfo.Open(this.Control.DbPath, this.Control.AlternateDbPaths);
            }
        }


        /// <summary>
        /// Optional interface tracker if XDG is used.
        /// </summary>
        public LevelSetTracker LsTrk {
            get;
            protected set;
        }

        /// <summary>
        /// Extended grid information.
        /// </summary>
        public IGridData GridData {
            get {
                return Grid.iGridData;
            }
        }

        /// <summary>
        /// Multigrid levels, sorted from fine to coarse, i.e. the 0-th entry contains the finest grid.
        /// </summary>
        public AggregationGridData[] MultigridSequence {
            get;
            private set;
        }

        /// <summary>
        /// BoSSS grid.
        /// </summary>
        public IGrid Grid {
            get;
            private set;
        }

        ///// <summary>
        ///// Provisional alternative grid;
        ///// </summary>
        //protected AggregationGrid AggGrid;


        /// <summary>
        /// <see cref="DatabaseDriver"/>
        /// </summary>
        private IDatabaseInfo m_Database;

        /// <summary>
        /// interface to the database driver
        /// </summary>
        public IDatabaseDriver DatabaseDriver {
            get {
                if (m_Database == null) {
                    return null;
                } else {
                    return m_Database.Controller.DBDriver;
                }
            }
        }

        /// <summary>
        /// use this method to either create a grid or to load a grid from
        /// the IO system (See <see cref="IDatabaseDriver"/>);
        /// </summary>
        protected virtual IGrid CreateOrLoadGrid() {
            using (var ht = new FuncTrace()) {

                


                if (this.Control != null) {
                    if (this.Control.GridFunc != null && this.Control.GridGuid != Guid.Empty)
                        throw new ApplicationException("Control object error: 'AppControl.GridFunc' and 'AppControl.GridGuid' are exclusive, cannot be unequal null at the same time.");
                    if (this.Control.GridFunc != null && this.Control.RestartInfo != null)
                        throw new ApplicationException("Control object error: 'AppControl.GridFunc' and 'AppControl.RestartInfo' are exclusive, cannot be unequal null at the same time.");
                    if (this.Control.GridFunc == null && this.Control.GridGuid == Guid.Empty && this.Control.RestartInfo == null)
                        throw new ApplicationException("Control object error: No grid specified -- either 'AppControl.GridFunc' or 'AppControl.GridGuid' or 'AppControl.RestartInfo' must be specified.");

                    if (this.Control.RestartInfo != null) {
                        // +++++++++++++++++++++++++++++++++++++++++++
                        // restart - use grid guid from last time-step
                        // +++++++++++++++++++++++++++++++++++++++++++

                        ISessionInfo session = m_Database.Controller.GetSessionInfo(this.Control.RestartInfo.Item1);
                        TimestepNumber timestep = this.Control.RestartInfo.Item2;
                        ITimestepInfo tsi_toLoad;
                        if (timestep == null || timestep.MajorNumber < 0) {
                            tsi_toLoad = session.Timesteps.OrderBy(tsi => tsi.PhysicalTime).Last();
                        } else {
                            tsi_toLoad = session.Timesteps.Single(t => t.TimeStepNumber.Equals(timestep));
                        }

                        if (this.Control.GridGuid != null && this.Control.GridGuid != default(Guid) && this.Control.GridGuid != Guid.Empty) {
                            if (this.Control.GridGuid != tsi_toLoad.GridID)
                                throw new ArgumentException($"Grid Guid mismatch for restart: 'Control.GridGuid' is set to {this.Control.GridGuid}, but grid for restart-timestep has id {tsi_toLoad.GridID}.");
                        }


                        var _Grid = DatabaseDriver.LoadGrid(tsi_toLoad.GridID, m_Database);

                        if (_Grid is GridCommons) {
                            GridCommons __Grid = (GridCommons)_Grid;
                            foreach (string oldBndy in this.Control.BoundaryValueChanges.Keys) {
                                int bndyInd = __Grid.EdgeTagNames.Values.FirstIndexWhere(bndyVal => bndyVal.Equals(oldBndy, StringComparison.InvariantCultureIgnoreCase));
                                if (bndyInd > -1) {
                                    __Grid.EdgeTagNames[__Grid.EdgeTagNames.Keys.ElementAt(bndyInd)] = this.Control.BoundaryValueChanges[oldBndy];
                                } else {
                                    throw new ArgumentException("Boundary " + oldBndy + " is not found in EdgeTagNames of the loaded Grid");
                                }
                            }
                        }

                        ht.LogMemoryStat();
                        return _Grid;
                    } else if (this.Control.GridGuid != null && this.Control.GridGuid != default(Guid)) {
                        // ++++++++++++++++++++
                        // load grid by GridID
                        // ++++++++++++++++++++
                        var _Grid = DatabaseDriver.LoadGrid(this.Control.GridGuid, m_Database);
                        ht.LogMemoryStat();
                        return _Grid;

                    } else if (this.Control.GridFunc != null) {
                        // ++++++++++++++++++++++++++++
                        // use grid-generating function
                        // ++++++++++++++++++++++++++++
                        var g = this.Control.GridFunc();
                        return g;
                    } else {
                        throw new ApplicationException("Unable to create grid from control object. 'AppControl.GridFunc' and 'AppControl.GridGuid' and 'AppControl.RestartInfo' are all null.");
                    }
                }

                throw new ApplicationException("No control object was given -- cannot load grid, you must override this method.");
            }
        }


        /// <summary>
        /// Adds some DG field to <see cref="m_RegisteredFields"/> and, optionally, to <see cref="m_IOFields"/>.
        /// </summary>
        public void RegisterField(DGField f, IOListOption ioOpt = IOListOption.ControlFileDetermined) {
            foreach (var _f in this.m_RegisteredFields)
                if (object.ReferenceEquals(_f, f))
                    return;

            this.m_RegisteredFields.Add(f);

            IDictionary<string, FieldOpts> FieldOptions = null;
            if (this.Control != null) {
                FieldOptions = this.Control.FieldOptions;
            }

            //FieldOpts fopts;
            //bool isSpec = FieldOptions.TryGetValue(f.Identification, out fopts);
            FieldOpts fopts = FieldOptions.Where(kv => kv.Key.WildcardMatch(f.Identification)).SingleOrDefault().Value;

            if (fopts != null) {
                if (ioOpt == IOListOption.Always && fopts.SaveToDB == FieldOpts.SaveToDBOpt.FALSE)
                    throw new ApplicationException("IO for field '" + f.Identification + "' cannot be turned OFF, i.e. 'SaveToDB==false' is illegal.");
                if (ioOpt == IOListOption.Never && fopts.SaveToDB == FieldOpts.SaveToDBOpt.TRUE)
                    throw new ApplicationException("IO for field '" + f.Identification + "' cannot be turned ON, i.e. 'SaveToDB==true' is illegal");

                if (ioOpt == IOListOption.Always || fopts.SaveToDB == FieldOpts.SaveToDBOpt.TRUE)
                    IOFields.Add(f);

            } else {
                if (ioOpt == IOListOption.Always)
                    IOFields.Add(f);
            }

        }

        /// <summary>
        /// Adds some DG fields to <see cref="m_RegisteredFields"/> and, optionally, to <see cref="m_IOFields"/>.
        /// </summary>
        public void RegisterField(IEnumerable<DGField> Lf, IOListOption ioOpt = IOListOption.ControlFileDetermined) {
            foreach (DGField f in Lf)
                RegisterField(f, ioOpt);
        }


        /// <summary>
        /// all DG fields that were decorated by an <see cref="InstantiateFromControlFileAttribute"/>.
        /// </summary>
        protected List<DGField> m_RegisteredFields = new List<DGField>();

        /// <summary>
        /// List of DG Fields which were somehow made known to the app
        /// </summary>
        public ICollection<DGField> RegisteredFields {
            get {
                return m_RegisteredFields.AsReadOnly();
            }
        }

        /// <summary>
        /// see <see cref="IOFields"/>
        /// </summary>        
        protected List<DGField> m_IOFields = new List<DGField>();

        /// <summary>
        /// All fields, for which IO should be performed by
        /// <see cref="SaveToDatabase"/> and <see cref="RestartFromDatabase"/>
        /// should be put into this collection; Also, the
        /// <see cref="SetInitial"/>-method takes into account only fields
        /// listed in this collection; that can happen during the
        /// <see cref="CreateFields"/>-call;
        /// </summary>
        public ICollection<DGField> IOFields {
            get {
                return m_IOFields;
            }
        }

        /// <summary>
        /// creates the <see cref="TimestepInfo"/> object that will be saved during <see cref="SaveToDatabase"/>;
        /// override this method to add user-data to a time-step
        /// </summary>
        /// <param name="t">
        /// time value which will be associated with the field
        /// </param>
        /// <param name="timestepno">time-step number</param>
        protected virtual TimestepInfo GetCurrentTimestepInfo(TimestepNumber timestepno, double t) {

            // store the MPI rank
            var _fields = m_IOFields.ToArray();

            {
                string ID_MPIrank = "MPIrank";
                {
                    int no = 2;
                    while (IOFields.Any(f => f.Identification == ID_MPIrank)) {
                        ID_MPIrank = "MPIrank_" + no;
                        no++;
                    }
                }

                SinglePhaseField MPIrnk;
                {
                    MPIrnk = new SinglePhaseField(new Basis(this.GridData, 0), ID_MPIrank);
                    int rnk = MPIRank;
                    MPIrnk.AccConstant(this.MPIRank);
                }

               
                MPIrnk.AddToArray(ref _fields);
            }

            if(this.IOFields.Any(f => f is XDGField) && LsTrk != null) {

                string ID_CutCells = "CutCells";
                {
                    int no = 2;
                    while (IOFields.Any(f => f.Identification == ID_CutCells)) {
                        ID_CutCells = "CutCells_" + no;
                        no++;
                    }
                }

                SinglePhaseField CutCells;
                {
                    CutCells = new SinglePhaseField(new Basis(this.GridData, 0), ID_CutCells);
                    CutCells.AccConstant(1.0, LsTrk.Regions.GetCutCellMask());
                }

                CutCells.AddToArray(ref _fields);
            }

            //
            return new TimestepInfo(t, this.CurrentSessionInfo, timestepno, _fields);
        }

        /// <summary>
        /// If data logging is turned on, saves all fields in
        /// <see cref="m_IOFields"/> to the database 
        /// </summary>
        /// <param name="t">
        /// time value which will be associated with the field
        /// </param>
        /// <param name="timestepno">time-step number</param>
        protected virtual ITimestepInfo SaveToDatabase(TimestepNumber timestepno, double t) {
            using (var ht = new FuncTrace()) {

                if (DatabaseDriver.FsDriver == null)
                    return null;
                if (this.CurrentSessionInfo.ID.Equals(Guid.Empty))
                    return null;

                TimestepInfo tsi = GetCurrentTimestepInfo(timestepno, t);
                //Exception e = null;
                try {
                    this.DatabaseDriver.SaveTimestep(tsi);
                } catch (Exception ee) {
                    Console.Error.WriteLine(ee.GetType().Name + " on rank " + this.MPIRank + " saving time-step " + timestepno + ": " + ee.Message);
                    Console.Error.WriteLine(ee.StackTrace);
                    //tsi = null;
                    //e = ee;

                    if (ContinueOnIOError) {
                        Console.WriteLine("Ignoring IO error: " + DateTime.Now);

                    } else {
                        throw ee;
                    }

                    tsi = null;
                }

                // e.ExceptionBcast();
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                return tsi;
            }
        }

        /// <summary>
        /// Calculation is not stopped if an I/O exception is thrown in <see cref="SaveToDatabase(TimestepNumber, double)"/>,
        /// see also <see cref="AppControl.ContinueOnIoError"/>.
        /// </summary>
        protected virtual bool ContinueOnIOError {
            get {
                if (this.Control != null)
                    return this.Control.ContinueOnIoError;
                else
                    return true;
            }
        }


        /// <summary>
        /// Loads all fields in <see cref="m_IOFields"/> from the database
        /// using the given <see cref="AppControl.RestartInfo"/> (<see cref="Control"/>)
        /// </summary>
        /// <param name="time">
        /// On exit, contains the physical time represented by the time-step
        /// </param>
        /// <returns>
        /// Returns the actual time-step number of the loaded time-step
        /// </returns>
        protected virtual TimestepNumber RestartFromDatabase(out double time) {
            using (var tr = new FuncTrace()) {

                // obtain session timesteps:
                var sessionToLoad = this.Control.RestartInfo.Item1;
                ISessionInfo session = m_Database.Controller.GetSessionInfo(sessionToLoad);

                var all_ts = session.Timesteps;

                // find timestep to load
                Guid tsi_toLoad_ID;
                tsi_toLoad_ID = GetRestartTimestepID();
                ITimestepInfo tsi_toLoad = all_ts.Single(t => t.ID.Equals(tsi_toLoad_ID));

                time = tsi_toLoad.PhysicalTime;

                bool MustResetTime = false;
                if(this.Control.TimesteppingMode == AppControl._TimesteppingMode.Transient) {
                    double dtMin = this.Control.dtMin;
                    if(dtMin > 0) {
                        if(time + dtMin == time) { // time is so advanced, that the timestep becomes invisible
                            MustResetTime = true;
                        }
                    }
                }

                if (session.KeysAndQueries.TryGetValue("TimesteppingMode", out object mode) || MustResetTime) {
                    if (Convert.ToInt32(mode) == (int)AppControl._TimesteppingMode.Steady) {
                        Console.WriteLine("Restarting from steady-state, resetting time ...");
                        time = 0.0; // Former simulation is steady-state, this should be restarted with time = 0.0
                    }
                }

                if (tsi_toLoad is BoSSS.Foundation.IO.TimestepProxy tp) {
                    var tsiI = tp.GetInternal() as TimestepInfo;
                    if (tsiI != null) {
                        OnRestartTimestepInfo(tsiI);
                    }
                }

                DatabaseDriver.LoadFieldData(tsi_toLoad, ((GridData)(this.GridData)), this.IOFields);

                // return
                return tsi_toLoad.TimeStepNumber;
            }
        }

        /// <summary>
        /// called on a restart, after the <see cref="TimestepInfo"/> is loaded from database.
        /// </summary>
        protected virtual void OnRestartTimestepInfo(TimestepInfo tsi) {

        }

        /// <summary>
        /// Returns, in the case of a restart (<see cref="AppControl.RestartInfo"/>), the ID of the time-step to restart.
        /// </summary>
        protected Guid GetRestartTimestepID() {
            var sessionToLoad = this.Control.RestartInfo.Item1;
            var timestep = this.Control.RestartInfo.Item2;

            ISessionInfo session = m_Database.Controller.GetSessionInfo(sessionToLoad);
            var all_ts = session.Timesteps;

            Guid tsi_toLoad_ID;
            if (timestep == null || timestep.MajorNumber < 0) {
                tsi_toLoad_ID = all_ts.OrderBy(tsi => tsi.PhysicalTime).Last().ID;
            } else {
                tsi_toLoad_ID = all_ts.Single(t => t.TimeStepNumber.Equals(timestep)).ID;
            }

            return tsi_toLoad_ID;
        }

        /// <summary>
        /// override this method to create the level-set tracker
        /// </summary>
        protected virtual void CreateTracker() {
        }

        /// <summary>
        /// override this method to create fields, if the decoration with
        /// <see cref="InstantiateFromControlFileAttribute"/>-attributes is not
        /// sufficient.
        /// </summary>
        protected virtual void CreateFields() {
        }

        /// <summary>
        /// override this method to create Equations
        /// and Solvers, time-steppers, ...
        /// </summary>
        /// <param name="L">
        /// If restarted after dynamic load balancing, the respective data.
        /// </param>
        protected abstract void CreateEquationsAndSolvers(GridUpdateDataVaultBase L);


        /// <summary>
        /// sets initial values as defined in the control file. Override this
        /// method to set initial values for the fields;
        /// </summary>
        protected virtual void SetInitial(double time) {
            using (var tr = new FuncTrace()) {

                this.QueryResultTable.UpdateKey("Timestep", ((int)0));

                if (this.Control == null) {
                    tr.Info("no control - file loaded; don't do anything;");
                    return;
                }

                var relevantFields = m_RegisteredFields.Union(
                    m_IOFields, ReferenceComparer.Instance).ToArray();

                foreach (var f in relevantFields) {
                    if (f is XDGField) {
                        var xdgf = (XDGField)f;
                        if (!object.ReferenceEquals(xdgf.Basis.Tracker, this.LsTrk))
                            throw new ApplicationException("XDG is defined against unknown level-set tracker!");
                    }
                }


                // pass 1: single phase fields
                // ===========================

                var Pass2_Evaluators = new Dictionary<string, ScalarFunction>();
                foreach (var val in this.Control.InitialValues_EvaluatorsVec) {
                    string DesiredFieldName = val.Key;
                    //ScalarFunction Function = Utils.NonVectorizedScalarFunction.Vectorize(val.Value);
                    ScalarFunction Function = val.Value.SetTime(time);

                    bool found = false;
                    foreach (DGField f in relevantFields) {
                        if (f.Identification.Equals(DesiredFieldName)) {
                            tr.Info("projecting field \"" + f.Identification + "\"");
                            f.ProjectField(Function);
                            found = true;
                            break;
                        } else {

                            // now, the XDG hack:
                            var NameAndSpc = DesiredFieldName.Split(new string[] { "#" }, StringSplitOptions.RemoveEmptyEntries);
                            if (NameAndSpc.Length == 2 && f.Identification.Equals(NameAndSpc[0])) {
                                tr.Info("projecting XDG-field \"" + f.Identification + "\"");
                                string spc = NameAndSpc[1];
                                var xdgf = (XDGField)f;
                                //var SpeciesOnlyField = xdgf.GetSpeciesShadowField(spc);
                                //SpeciesOnlyField.ProjectField(Function);

                                Pass2_Evaluators.Add(val.Key, val.Value.SetTime(time));

                                found = true;
                                break;
                            }
                        }
                    }

                    if (!found) {
                        Console.WriteLine("Warning: " +
                            "initial value specified for a field named \"" + DesiredFieldName +
                            "\", but no field with that identification exists in context.");
                    }
                }

                if (LsTrk != null) {
                    LsTrk.UpdateTracker(time);
                    LsTrk.UpdateTracker(time); // doppeltes Update h�lt besser; 
                }

                // pass 2: XDG fields (after tracker update)
                // =========================================
                if (Pass2_Evaluators.Count > 0) {


                    foreach (var val in Pass2_Evaluators) {
                        string DesiredFieldName = val.Key;
                        //ScalarFunction Function = Utils.NonVectorizedScalarFunction.Vectorize(val.Value);
                        ScalarFunction Function = val.Value;

                        bool found = false;
                        foreach (DGField f in relevantFields) {
                            if (f.Identification.Equals(DesiredFieldName)) {
                                throw new ApplicationException();
                            } else {

                                // now, the XDG hack:
                                var NameAndSpc = DesiredFieldName.Split(new string[] { "#" }, StringSplitOptions.RemoveEmptyEntries);
                                if (NameAndSpc.Length == 2 && f.Identification.Equals(NameAndSpc[0])) {
                                    tr.Info("projecting XDG-field \"" + f.Identification + "\"");
                                    string spc = NameAndSpc[1];
                                    var xdgf = (XDGField)f;
                                    var SpeciesOnlyField = xdgf.GetSpeciesShadowField(spc);
                                    SpeciesOnlyField.ProjectField(Function);
                                    found = true;
                                    break;
                                }
                            }
                        }

                        if (!found) {
                            throw new ApplicationException(
                                "initial value specified for a field named \"" + DesiredFieldName +
                                "\", but no field with that identification exists in context.");
                        }
                    }
                }
            }
        }


        /// <summary>
        /// number of time-steps to be performed; <see cref="RunSolverMode"/>
        /// terminates if the number of time-steps exceeds this number; At startup, initialized equal to  <see cref="AppControl.NoOfTimesteps"/>.
        /// </summary>
        /// <remarks>
        /// This also holds if the configured <see cref="EndTime"/> has not
        /// been reached yet.
        /// </remarks>
        protected int NoOfTimesteps = 1;

        /// <summary>
        /// physical end time at which the simulation should be stopped:
        /// <see cref="RunSolverMode"/> terminates if the physical simulation
        /// time exceeds this value;
        /// </summary>
        /// <remarks>
        /// This also holds if the configured <see cref="NoOfTimesteps"/> has
        /// not been performed yet.
        /// </remarks>
        protected double EndTime = double.MaxValue;

        /// <summary>
        /// if this value is set to n, a restart file is written every n time-steps, <see cref="AppControl.saveperiod"/>
        /// </summary>
        protected int SavePeriod = 1;

        /// <summary>
        /// <see cref="AppControl.rollingSaves"/>
        /// </summary>
        protected bool RollingSave = false;

        /// <summary>
        /// Number of Consecutive timesteps which are saved -- this is intended to be used by BDF or Adams-Bashforth time integrators which require multiple time steps
        /// (e.g. 3 to save time-step 98, 99, 100 for a save-period of 100;)
        /// </summary>
        protected virtual int BurstSave {
            get {
                return Math.Max(1, this.Control.BurstSave);
            }
        }


        /// <summary>
        /// Implement this method by performing a single time-step of the
        /// solution algorithm.
        /// </summary>
        /// <param name="phystime">
        /// physical time prior to the call of this time step
        /// </param>
        /// <param name="TimestepNo">
        /// Time-step number 
        /// </param>
        /// <param name="dt">
        /// A desired time-step size; if greater than 0.0, then exactly this
        /// time should be simulated and passed back as an return value; (This
        /// is used for temporal convergence studies).
        /// </param>
        /// <returns>
        /// The size of the time-step; If <paramref name="dt"/> is greater than
        /// 0, the return value should be equal to this;
        /// </returns>
        protected abstract double RunSolverOneStep(int TimestepNo, double phystime, double dt);

        /// <summary>
        /// Override this method to implement plotting (by e.g. the use of the 
        /// Tecplot interface) of the current state.
        /// </summary>
        /// <param name="physTime">
        /// physical time
        /// </param>
        /// <param name="timestepNo">
        /// Time-step number, can be used by implementer to define file name
        /// for plots
        /// </param>
        /// <param name="superSampling">
        /// Super sampling option (recursive subdivisions of individual grid
        /// cells)
        /// </param>
        protected abstract void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0);

        /// <summary>
        /// For terminating the <see cref="RunSolverMode"/> before reaching the
        /// <see cref="EndTime"/> or <see cref="NoOfTimesteps"/>
        /// </summary>
        protected Boolean TerminationKey = false;

        /// <summary>
        /// If the current simulation has been restarted, <see cref="TimeStepNoRestart"/>
        /// is set by the method <see cref="LoadRestart(out double, out TimestepNumber)"/>.
        /// </summary>
        protected TimestepNumber TimeStepNoRestart = null;

        /// <summary>
        /// 
        /// </summary>
        public virtual void RunSolverMode() {
#if RELEASE
            try {
#endif
                this.RunSolverModeInternal();
#if RELEASE
            } catch (Exception e) {
                SolverExceptionLogger.SaveException(e, this);
                throw;
            }
#endif
        }

        /// <summary>
        /// Runs the application in the "solver"-mode. This method makes
        /// multiple calls to <see cref="RunSolverOneStep"/> method. The
        /// termination behavior is determined by the variables
        /// <see cref="EndTime"/> and <see cref="NoOfTimesteps"/>. To modify
        /// the termination behavior of this method, one can overwrite
        /// <see cref="EndTime"/> and <see cref="NoOfTimesteps"/>, e.g. from
        /// within <see cref="RunSolverOneStep"/>
        /// </summary>
        /// <remarks>
        /// Important methods are called in the following order:
        /// <list type="number">
        ///     <item><see cref="SetUpEnvironment"/></item>
        ///     <item><see cref="SetInitial"/> or <see cref="LoadRestart"/></item>
        ///     <item><see cref="CreateEquationsAndSolvers"/></item>
        ///     <item>
        ///         Possibly multiple times:
        ///         <list type="number">
        ///             <item><see cref="RunSolverOneStep"/></item>
        ///             <item><see cref="SaveToDatabase"/></item>
        ///             <item><see cref="PlotCurrentState"/></item>
        ///         </list>
        ///     </item>
        ///     <item><see cref="SaveToDatabase"/></item>
        ///     <item><see cref="Queries.QueryHandler.EvaluateQueries"/></item>
        /// </list>
        /// </remarks>
        private void RunSolverModeInternal() {

            // =========================================
            // loading grid, initializing database, etc:
            // =========================================
            SetUpEnvironment(); // remark: tracer is not avail before setup

            using (var tr = new FuncTrace()) {
                //tr.InfoToConsole = true;
                var rollingSavesTsi = new List<Tuple<int, ITimestepInfo>>();

                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                double physTime = 0.0;
                TimestepNumber i0 = 0;
                if (this.Control == null) {
                    SetInitial(0.0); // default behavior if no control file is present
                } else {
                    if (this.Control != null) {
                        if (!this.Control.InitialValues_Evaluators.IsNullOrEmpty() && this.Control.RestartInfo != null) {
                            //throw new ApplicationException("Invalid state in control object: the specification of initial values ('AppControl.InitialValues') and restart info ('AppControl.RestartInfo') is exclusive: "
                            //    + " both cannot be unequal null at the same time.");
                            Console.WriteLine("Warning: InitialValues set, while restarting a simulation.");
                        }

                        if (this.Control.RestartInfo != null) {
                            LoadRestart(out physTime, out i0);
                            TimeStepNoRestart = i0;
                        } else {
                            SetInitial(0.0);
                        }
                    }
                }
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                m_queryHandler.QueryResults.Clear();

                if (this.Control.RestartInfo != null) {
                    CreateEquationsAndSolvers(null);
                    tr.LogMemoryStat();
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                    if(LsTrk != null) {
                        if(LsTrk.Regions.Time != physTime)
                            LsTrk.UpdateTracker(physTime);
                        LsTrk.PushStacks();
                    }
                }

                // =========================================================
                // Adaptive-Mesh-Refinement and/or load balancing on startup
                // =========================================================


                // load balancing solo
                tr.Info("DynamicLoadBalancing_RedistributeAtStartup = " + this.Control.DynamicLoadBalancing_RedistributeAtStartup);
                tr.Info("AdaptiveMeshRefinement = " + this.Control.AdaptiveMeshRefinement);
                if (this.Control.DynamicLoadBalancing_RedistributeAtStartup && !this.Control.AdaptiveMeshRefinement) {
                    PlotAndSave(physTime, i0, rollingSavesTsi);
                    MpiRedistributeAndMeshAdaptOnInit(i0.MajorNumber, physTime);
                    PlotAndSave(physTime, i0, rollingSavesTsi);
                }
                          

                // load balancing and adaptive mesh refinement
                if (this.Control.AdaptiveMeshRefinement) {
                    
                    // unprocessed initial value IO
                    if (this.Control != null && this.Control.ImmediatePlotPeriod > 0)
                        PlotCurrentState(physTime, new TimestepNumber(i0.Numbers.Cat(0)), this.Control.SuperSampling);

                    var ts0amr = SaveToDatabase(new TimestepNumber(i0.Numbers.Cat(0)), physTime); // save the initial value
                    if (this.RollingSave)
                        rollingSavesTsi.Add(Tuple.Create(0, ts0amr));


                    bool initialRedist = false;
                    for (int s = 1; s <= this.Control.AMR_startUpSweeps; s++) {
                        initialRedist |= this.MpiRedistributeAndMeshAdaptOnInit(i0.MajorNumber, physTime);

                        if (initialRedist == true) {

                            if (this.Control.ImmediatePlotPeriod > 0)
                                PlotCurrentState(physTime, new TimestepNumber(i0.Numbers.Cat(s)), this.Control.SuperSampling);

                            ts0amr = SaveToDatabase(new TimestepNumber(i0.Numbers.Cat(s)), physTime); // save the AMR'ed initial value
                            if (this.RollingSave)
                                rollingSavesTsi[0] = Tuple.Create(0, ts0amr);

                        }
                    }
                }

                // ================================================================================
                // sometimes, the operators depend on parameters,
                // therefore 'CreateEquationsAndSolvers()' has to be called after ' SetInitial()',
                // resp. 'LoadRestart(..)'!!!
                // ================================================================================

                if (this.Control.RestartInfo == null) {
                    CreateEquationsAndSolvers(null);
                }
                {
                    tr.LogMemoryStat();
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                    if (LsTrk != null)
                        LsTrk.PushStacks();
                }
                // ========================================================================
                // initial value IO:
                // (note: in some apps, the initial values might be tweaked in the 
                // 'CreateEquationsAndSolvers(...)' method; but here we should have the 
                // "true" initial value)
                // ========================================================================


                if (this.Control != null && this.Control.ImmediatePlotPeriod > 0)
                    PlotCurrentState(physTime, i0, this.Control.SuperSampling);

                var ts0 = SaveToDatabase(i0, physTime); // save the initial value
                if (this.RollingSave)
                    rollingSavesTsi.Add(Tuple.Create(0, ts0));

                // =========================================
                // Adaptive-Mesh-Refinement on startup
                // =========================================

                //bool initialRedist = false;
                //for (int s = 0; s < this.Control.AMR_startUpSweeps; s++) {
                //    initialRedist |= this.MpiRedistributeAndMeshAdapt(i0.MajorNumber, physTime);

                //    if (this.Control.ImmediatePlotPeriod > 0 && initialRedist == true)
                //        PlotCurrentState(physTime, new TimestepNumber(i0.Numbers.Cat(s)), this.Control.SuperSampling);

                //    if (initialRedist == true) {
                //        ts0 = SaveToDatabase(new TimestepNumber(i0.Numbers.Cat(s)), physTime); // save the AMR'ed initial value
                //        if (this.RollingSave)
                //            rollingSavesTsi[0] = Tuple.Create(0, ts0);
                //    }
                //}

                // =================================================================================
                // Main/outmost time-stepping loop
                // (in steady-state: only one iteration)
                // =================================================================================
                {


                    // setup of logging
                    foreach (var l in PostprocessingModules) {
                        l.Setup(this);
                        l.DriverTimestepPostProcessing(i0.MajorNumber, physTime);
                    }


                    bool RunLoop(int i) {
                        return (i <= i0.MajorNumber + (long)NoOfTimesteps) && EndTime - physTime > 1.0E-10 && !TerminationKey;
                    }

                    for (int i = i0.MajorNumber + 1; RunLoop(i); i++) {
                        tr.Info("performing timestep " + i + ", physical time = " + physTime);
                        this.MpiRedistributeAndMeshAdapt(i, physTime);
                        this.QueryResultTable.UpdateKey("Timestep", ((int)i));
                        // Call the solver    vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                        double dt = RunSolverOneStep(i, physTime, -1);
                        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        tr.Info("simulated time: " + dt + " timeunits.");
                        tr.LogMemoryStat();
                        physTime += dt;

                        if(LsTrk != null) {
                            if(LsTrk.Regions.Time != physTime) {
                                // correct the level-set tracker time if some solver did not correctly updated it.
                                LsTrk.UpdateTracker(physTime);
                            }
                        }


                        foreach (var l in PostprocessingModules) {
                            l.DriverTimestepPostProcessing(i, physTime);
                        }

                        ITimestepInfo tsi = null;

                        if (this.BurstSave < 1) {
                            throw new NotSupportedException("misconfiguration of burst save variable.");
                        }
                        

                        for (int sb = 0; sb < this.BurstSave; sb++) {
                            if ((i + sb) % SavePeriod == 0 || (!RunLoop(i + 1) && sb == 0)) {
                                tsi = SaveToDatabase(i, physTime);
                                this.ProfilingLog();
                                break;
                            }
                        }
                        
                        if (this.RollingSave) {
                            if (tsi == null) {
                                tsi = SaveToDatabase(i, physTime);
                            }
                            rollingSavesTsi.Add(Tuple.Create(i, tsi));

                            while (rollingSavesTsi.Count > this.BurstSave) { // delete overdue rolling timesteps...
                                var top_i_tsi = rollingSavesTsi[0];

                                rollingSavesTsi.RemoveAt(0);

                                if ((top_i_tsi.Item1 != 0) && (top_i_tsi.Item1 % SavePeriod != 0)) { // ...only if they should not be saved anyway
                                    if (DatabaseDriver.FsDriver != null &&
                                        !this.CurrentSessionInfo.ID.Equals(Guid.Empty)) {
                                        if (MPIRank == 0) {
                                            this.CurrentSessionInfo.RemoveTimestep(top_i_tsi.Item2.ID);
                                            ((DatabaseController)this.m_Database.Controller).DeleteTimestep(top_i_tsi.Item2, false);
                                        }
                                    }
                                }
                            }
                        }

                        if (this.Control != null && this.Control.ImmediatePlotPeriod > 0 && i % this.Control.ImmediatePlotPeriod == 0)
                            PlotCurrentState(physTime, i, this.Control.SuperSampling);
                    }


                    // =================================================================================
                    // Evaluate queries and write log file 
                    // (either to session directory or current directory)
                    // =================================================================================
                    m_queryHandler.EvaluateQueries(this.m_RegisteredFields.Union(m_IOFields), physTime);
                    foreach (var kv in m_queryHandler.QueryResults) {
                        QueryResultTable.LogValue(kv.Key, kv.Value);
                        if (!this.CurrentSessionInfo.KeysAndQueries.ContainsKey(kv.Key))
                            this.CurrentSessionInfo.KeysAndQueries.Add(kv.Key, kv.Value);
                    }

                    if (MPIRank == 0 && m_queryHandler.QueryResults.Count > 0) {
                        TextWriter queryLogFile_Txt;

                        if (Control != null && Control.savetodb) {
                            queryLogFile_Txt = DatabaseDriver.FsDriver.GetNewLog("queryResults", this.CurrentSessionInfo.ID);
                        } else {
                            try {
                                queryLogFile_Txt = new StreamWriter("queryResults.txt");
                            } catch (Exception e) {
                                // in a parameter study, 
                                //     - when running in different processes
                                //     - but simultaneously
                                //     - without database
                                //  two or more processes may try to access queryResults.txt 
                                //  => Exception
                                // this is such a rare case, that I don't implement a smarter solution
                                // (in the parameter study case, the query results will be in the ParameterStudy file anyway

                                Console.WriteLine("WARNING: not writing queryResults.txt file due to exception: {0} \n '{1}'",
                                    e.GetType().Name, e.Message);
                                queryLogFile_Txt = null;
                            }
                        }

                        if (queryLogFile_Txt != null) {
                            QueryHandler.LogQueryResults(m_queryHandler.QueryResults, queryLogFile_Txt);
                            queryLogFile_Txt.Close();
                        }
                    }

                    CorrectlyTerminated = true;

                }
            }
        }

        private void PlotAndSave(double TS, TimestepNumber TSno, List<Tuple<int, ITimestepInfo>> rollingSavesSammeldingens) {
            if (this.Control != null && this.Control.ImmediatePlotPeriod > 0)
                PlotCurrentState(TS, new TimestepNumber(TSno.Numbers.Cat(0)), this.Control.SuperSampling);

            var ts0amr = SaveToDatabase(new TimestepNumber(TSno.Numbers.Cat(0)), TS); // save the initial value
            if (this.RollingSave)
                rollingSavesSammeldingens.Add(Tuple.Create(0, ts0amr));
        }

        /// <summary>
        /// Main routine for dynamic load balancing and adaptive mesh refinement.
        /// </summary>
        /// <returns>
        /// - true if the mesh is actually changed or re-distributed
        /// - false if not
        /// </returns>
        virtual protected bool MpiRedistributeAndMeshAdapt(int TimeStepNo, double physTime, int[] fixedPartition = null, Permutation fixedPermutation = null) {

            DoMeshAdaption(TimeStepNo, physTime);
            DoLoadbalancing(TimeStepNo, physTime, fixedPartition, fixedPermutation);

            //this.QueryHandler.ValueQuery("UsedNoOfMultigridLevels", this.MultigridSequence.Length, true); 
            //PlotCurrentState(physTime, new TimestepNumber(new int[] { TimeStepNo, 12 }), 2);
            return true;
        }

        /// <summary>
        /// routine combining load balancing and adaptive mesh refinement before Timestepper-Init and Create CreateEquationsAndSolvers.
        /// This is a temporary and ugly state, a developer with more insight, may unite these methods
        /// </summary>
        /// <param name="TimeStepNo"></param>
        /// <param name="physTime"></param>
        /// <param name="fixedPartition"></param>
        /// <param name="fixedPermutation"></param>
        /// <returns></returns>
        virtual protected bool MpiRedistributeAndMeshAdaptOnInit(int TimeStepNo, double physTime, int[] fixedPartition = null, Permutation fixedPermutation = null) {
            double tmp = this.Control.DynamicLoadBalancing_ImbalanceThreshold;
            bool IsRestarted = this.Control.RestartInfo != null;
            bool IsInit = !IsRestarted;

            this.Control.DynamicLoadBalancing_ImbalanceThreshold = 0.0; // ensures that there is a redistribution at startup, idependant of threshold
            DoMeshAdaption(TimeStepNo, physTime, IsInit);
            DoLoadbalancing(TimeStepNo, physTime, fixedPartition, fixedPermutation, IsInit);
            this.Control.DynamicLoadBalancing_ImbalanceThreshold = tmp;
            //this.QueryHandler.ValueQuery("UsedNoOfMultigridLevels", this.MultigridSequence.Length, true);
            //PlotCurrentState(physTime, new TimestepNumber(new int[] { TimeStepNo, 11 }), 2);

            return true;
        }

        private bool DoLoadbalancing(int TimeStepNo, double physTime, int[] fixedPartition = null, Permutation fixedPermutation = null, bool IsInit=false) {
            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // no mesh adaptation, but (maybe) grid redistribution (load balancing)
            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            using (new FuncTrace()) {
                int[] NewPartition = fixedPartition ?? ComputeNewCellDistribution(TimeStepNo, physTime);
                if (NewPartition == null)
                    return false; // immediate quit, because there is nothing to do

                int JupOld = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int NoOfRedistCells = CheckPartition(NewPartition, JupOld);

                if(NoOfRedistCells <= 0) {
                    return false;
                } else {
#if DEBUG
                    Console.WriteLine("Re-distribution of " + NoOfRedistCells + " cells.");
#endif
                }

                // backup old data
                // ===============
                GridData oldGridData = ((GridData)(this.GridData));
                Permutation tau;
                GridUpdateDataVault_LoadBal loadbal = new GridUpdateDataVault_LoadBal(oldGridData, this.LsTrk);

                if(IsInit)
                    BackupDataOnInit(oldGridData, this.LsTrk, loadbal, out tau);
                    //BackupData(oldGridData, this.LsTrk, loadbal, out tau);
                else
                    BackupData(oldGridData, this.LsTrk, physTime, loadbal, out tau);

                // create new grid
                // ===============
                GridData newGridData;
                {
                    this.MultigridSequence = null;

                    this.Grid.RedistributeGrid(NewPartition);
                    newGridData = (GridData)this.Grid.iGridData;
                    oldGridData.Invalidate();
                    if (this.LsTrk != null) {
                        this.LsTrk.Invalidate();
                    }

                    if (this.Control == null || this.Control.NoOfMultigridLevels > 0)
                        this.MultigridSequence = CoarseningAlgorithms.CreateSequence(this.GridData, MaxDepth: (this.Control != null ? this.Control.NoOfMultigridLevels : 1));
                    else
                        this.MultigridSequence = new AggregationGridData[0];

                    //Console.WriteLine("P {0}: new grid: {1} cells.", MPIRank, newGridData.iLogicalCells.NoOfLocalUpdatedCells);
                }

                // compute redistribution permutation
                // ==================================

                Permutation Resorting;
                {
                    // sigma is the GlobalID-permutation of the **new** grid
                    Permutation sigma = fixedPermutation ?? newGridData.CurrentGlobalIdPermutation;

                    // compute resorting permutation
                    Permutation invSigma = sigma.Invert();
                    Resorting = invSigma * tau;
                    tau = null;
                    invSigma = null;
                }
                //Console.WriteLine("P {0}: Resorting: {1} entries.", MPIRank, Resorting.LocalLength);
                //ilPSP.Environment.StdoutOnlyOnRank0 = true;
                //Debug.Assert(Resorting.LocalLength == newGridData.iLogicalCells.NoOfLocalUpdatedCells);


                // sent data around the world
                // ==========================
                int newJ = newGridData.CellPartitioning.LocalLength;

                //int[] newTrackerData = null;
                //if(oldTrackerData != null) {
                //    newTrackerData = new int[newJ];
                //    Resorting.ApplyToVector(oldTrackerData, newTrackerData, newGridData.CellPartitioning);
                //    oldTrackerData = null;
                //}

                loadbal.Resort(Resorting, newGridData);

                // re-init simulation
                // ==================

                // release old DG fields
                this.m_RegisteredFields.Clear();
                this.m_IOFields.Clear();
                this.LsTrk = null;

                // re-create fields
                if (this.Control != null) {
                    InitFromAttributes.CreateFieldsAuto(
                        this, GridData, this.Control.FieldOptions, this.Control.CutCellQuadratureType, this.m_IOFields, this.m_RegisteredFields);
                }
                CreateFields(); // full user control   
                PostRestart(physTime, TimeStepNo);


                // re-set Level-Set tracker
                int trackerVersion = loadbal.SetNewTracker(this.LsTrk);
                //if(this.LsTrk != null) {
                //    Debug.Assert(object.ReferenceEquals(this.LsTrk.GridDat, this.GridData));
                //    Debug.Assert(this.LsTrk.Regions.Version == trackerVersion);
                //    foreach(var f in m_RegisteredFields) {
                //        if(f is XDGField) {
                //            ((XDGField)f).Override_TrackerVersionCnt(trackerVersion);
                //        }
                //    }
                //}

                //// skip this for init
                ReCreateEquationAndSolvers(IsInit, loadbal, physTime);
                return true;
            }
        }


        private bool DoMeshAdaption(int TimeStepNo, double physTime, bool IsInit = false) {
            using (var tr = new FuncTrace()) {

                bool plotAdaption = false;

                this.AdaptMesh(TimeStepNo, out var newGrid, out var old2newGridCorr);
                if (newGrid == null)
                    return false;

                using (new BlockTrace("process mesh Adaption", tr)) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // mesh adaptation
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    if (plotAdaption)
                        PlotCurrentState(physTime, new TimestepNumber(new int[] { TimeStepNo, 10 }), 2);

                    // backup old data
                    // ===============

                    GridData oldGridData = (GridData)this.GridData;
                    GridCommons oldGrid = oldGridData.Grid;
                    Guid oldGridId = oldGrid.ID;
                    Permutation tau;
                    GridUpdateDataVault_Adapt remshDat = new GridUpdateDataVault_Adapt(oldGridData, this.LsTrk);
                    
                    if(IsInit)
                        BackupDataOnInit(oldGridData, this.LsTrk, remshDat, out tau);
                        //BackupData(oldGridData, this.LsTrk, remshDat, out tau);
                    else
                        BackupData(oldGridData, this.LsTrk, physTime, remshDat, out tau);
                    
                    // save new grid to database
                    // ==========================

                    if (!passiveIo) {

                        if (newGrid.ID == null || newGrid.ID.Equals(Guid.Empty))
                            throw new ApplicationException();
                        if (newGrid.ID.Equals(oldGridId))
                            throw new ApplicationException();
                        if (DatabaseDriver.GridExists(newGrid.ID))
                            throw new ApplicationException();

                        DatabaseDriver.SaveGrid(newGrid, this.m_Database);
                    }

                    // create new grid
                    // ===============
                    GridData newGridData;
                    {
                        this.MultigridSequence = null;

                        this.Grid = newGrid;
                        newGridData = (GridData)this.Grid.iGridData;
                        oldGridData.Invalidate();
                        if (this.LsTrk != null) {
                            this.LsTrk.Invalidate();
                        }
                        LsTrk = null;
                        oldGridData = null;

                        if (this.Control == null || this.Control.NoOfMultigridLevels > 0)
                            this.MultigridSequence = CoarseningAlgorithms.CreateSequence(this.GridData,
                                MaxDepth: (this.Control != null ? this.Control.NoOfMultigridLevels : 1));
                        else
                            this.MultigridSequence = new AggregationGridData[0];

                        //Console.WriteLine("P {0}: new grid: {1} cells.", MPIRank, newGridData.iLogicalCells.NoOfLocalUpdatedCells);
                    }

                    // compute redistribution permutation
                    // ==================================

                    old2newGridCorr.ComputeDataRedist(newGridData);

                    int newJ = newGridData.Cells.NoOfLocalUpdatedCells;
                    //int[][] TargMappingIdx = old2newGridCorr.GetTargetMappingIndex(newGridData.CellPartitioning);

                    // sent data around the world
                    // ==========================

                    remshDat.Resort(old2newGridCorr, newGridData);

                    // re-init simulation
                    // ==================

                    // release old DG fields
                    this.m_RegisteredFields.Clear();
                    this.m_IOFields.Clear();

                    // re-set Level-Set tracker
                    this.CreateTracker();
                    int trackerVersion = remshDat.SetNewTracker(this.LsTrk);

                    // re-create fields
                    if (this.Control != null) {
                        InitFromAttributes.CreateFieldsAuto(
                            this, GridData, this.Control.FieldOptions, this.Control.CutCellQuadratureType, this.m_IOFields, this.m_RegisteredFields);
                    }
                    CreateFields(); // full user control   
                                    //PostRestart(physTime, TimeStepNo);
                    if(this.LsTrk != null) {
                        if(this.LsTrk.Regions.Time != physTime)
                            this.LsTrk.UpdateTracker(physTime);
                    }

                    if (plotAdaption)
                        PlotCurrentState(physTime, new TimestepNumber(new int[] { TimeStepNo, 11 }), 2);
                    
                    if (IsInit && this.Control.RestartInfo != null)
                        PostRestart(physTime, TimeStepNo);

                    ReCreateEquationAndSolvers(IsInit, remshDat, physTime);
                }
                return true;
            }
        }


        private void ReCreateEquationAndSolvers(bool IsInit, GridUpdateDataVaultBase GDataVault, double physTime) {
            if (IsInit) {
                if (this.Control.RestartInfo == null)
                    SetInitial(physTime);
            } else {
                //set dg coordinates
                foreach (var f in m_RegisteredFields) {
                    if (f is XDGField) {
                        XDGBasis xb = ((XDGField)f).Basis;
                        if (!object.ReferenceEquals(xb.Tracker, this.LsTrk))
                            throw new ApplicationException();
                    }
                    GDataVault.RestoreDGField(f);
                }

                // re-create solvers, etc.
                CreateEquationsAndSolvers(GDataVault);
            }
        }

        private void BackupData(GridData oldGridData, LevelSetTracker oldLsTrk, double physTime,
            GridUpdateDataVaultBase loadbal, out Permutation tau) {
            if(oldLsTrk != null && !object.ReferenceEquals(oldGridData, oldLsTrk.GridDat))
                throw new ApplicationException();
          
            // id's of the fields which we are going to rescue
            string[] FieldIds = m_RegisteredFields.Select(f => f.Identification).ToArray();
            
            // tau   is the GlobalID-permutation of the **old** grid
            tau = oldGridData.CurrentGlobalIdPermutation.CloneAs();

            // backup level-set tracker 
            if (this.LsTrk != null) {
                loadbal.BackupTracker(physTime);
            }

            // backup DG Fields
            foreach(var f in this.m_RegisteredFields) {
                if(f is XDGField) {
                    XDGBasis xb = ((XDGField)f).Basis;
                    if(!object.ReferenceEquals(xb.Tracker, oldLsTrk))
                        throw new ApplicationException();
                }
                if(!object.ReferenceEquals(f.Basis.GridDat, oldGridData))
                    throw new ApplicationException();


                loadbal.BackupField(f);
            }

            // backup user data
            this.DataBackupBeforeBalancing(loadbal);
        }


        private void BackupDataOnInit(GridData oldGridData, LevelSetTracker oldLsTrk, GridUpdateDataVaultBase loadbal, out Permutation tau) {

            tau = oldGridData.CurrentGlobalIdPermutation.CloneAs();

            // backup level-set tracker 
            if (this.LsTrk != null) {
                loadbal.BackupTracker(0.0);
            }
        }


        private static int CheckPartition(int[] NewPartition, int JupOld) {
            int mpiRank;// = oldGridData.CellPartitioning.MpiRank;
            int mpiSize;// = oldGridData.CellPartitioning.MpiSize;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out mpiRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out mpiSize);

            int locNoOfRedistCells = 0;

            // check Partitioning
            if (NewPartition.Length != JupOld)
                throw new ApplicationException("Illegal length of return value of 'ComputeNewCellDistribution'.");

            int[] locNoOfCellsPerProcessor = new int[mpiSize + 1], NoOfCellsPerProcessor;
            for (int j = 0; j < JupOld; j++) {
                if (NewPartition[j] < 0)
                    throw new ApplicationException("Illegal MPI rank '" + NewPartition[j] + "' from 'ComputeNewCellDistribution(...)'.");
                if (NewPartition[j] >= mpiSize)
                    throw new ApplicationException("Illegal MPI rank '" + NewPartition[j] + "' from 'ComputeNewCellDistribution(...)' - larger or equal than MPI size, which is '" + mpiSize + "'.");
                if (NewPartition[j] != mpiRank)
                    locNoOfRedistCells++;

                locNoOfCellsPerProcessor[NewPartition[j]]++;
            }
            locNoOfCellsPerProcessor[mpiSize] = locNoOfRedistCells;

            NoOfCellsPerProcessor = locNoOfCellsPerProcessor.MPISum();
            int NoOfRedistCells = NoOfCellsPerProcessor[mpiSize];

            //for (int proc = 0; proc < mpiSize; proc++)
            //    Console.WriteLine("Proc {0}: {1} cells.", proc, NoOfCellsPerProcessor[proc]);

            for (int proc = 0; proc < mpiSize; proc++) {
                if (NoOfCellsPerProcessor[proc] <= 0)
                    throw new ApplicationException("Zero cells on one processor are nor allowed.");
            }

            return NoOfRedistCells;
        }

        /// <summary>
        /// Additional backup (e.g. internal states of time integrators before grid redistribution)
        /// during dynamic load balancing.
        /// May also be used to invalidate internal states related to the old <see cref="GridData"/> or <see cref="LsTrk"/> objects.
        /// </summary>
        public virtual void DataBackupBeforeBalancing(GridUpdateDataVaultBase L) {

        }


        LoadBalancer m_Balancer;

        /// <summary>
        /// Adaptation of the current mesh (<see cref="Grid"/>).
        /// </summary>
        protected virtual void AdaptMesh(int TimestepNo, out GridCommons newGrid, out GridCorrelation old2NewGrid) {
            newGrid = null;
            old2NewGrid = null;
        }

        /// <summary>
        /// Default implementation (which does nothing) for the computation of the grid partitioning
        /// during runtime (dynamic load balancing).
        /// </summary>
        /// <param name="TimeStepNo"></param>
        /// <param name="physTime"></param>
        /// <returns>
        /// - if null, no re-partitioning of the calculation is done.
        /// - Otherwise, an array which defines, for each cell, the MPI rank of the desired owner processor.
        /// In the case of one MPI process, this method should always return <c>null</c>.
        /// </returns>
        protected virtual int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                if (Control == null
                    || !Control.DynamicLoadBalancing_On
                    || (     TimeStepNo % Control.DynamicLoadBalancing_Period != 0 
                        && !(Control.DynamicLoadBalancing_RedistributeAtStartup && TimeStepNo == TimeStepNoRestart))  // Variant for single partioning at restart
                    || (Control.DynamicLoadBalancing_Period < 0 && !Control.DynamicLoadBalancing_RedistributeAtStartup)
                    || MPISize <= 1) {
                    return null;
                }

                
                if (m_Balancer == null) {
                    m_Balancer = new LoadBalancer(Control.DynamicLoadBalancing_CellCostEstimators, this);

                    
                }

                return m_Balancer.GetNewPartitioning(
                    this,
                    TimeStepNo,
                    Control.GridPartType,
                    Control.GridPartOptions,
                    Control != null ? Control.DynamicLoadBalancing_ImbalanceThreshold : 0.12,
                    Control != null ? Control.DynamicLoadBalancing_Period : 5,
                    redistributeAtStartup: Control.DynamicLoadBalancing_RedistributeAtStartup,
                    TimestepNoRestart: TimeStepNoRestart);
            }
           
        }

        /// <summary>
        /// The name of a specific simulation should be logged in the <see cref="ISessionInfo.KeysAndQueries"/> under this key,
        /// see also <see cref="AppControl.SessionName"/>
        /// </summary>
        public const string SESSIONNAME_KEY = "SessionName";


        /// <summary>
        /// The name of a specific project should be logged in the <see cref="ISessionInfo.KeysAndQueries"/> under this key,
        /// see also <see cref="AppControl.ProjectName"/>
        /// </summary>
        public const string PROJECTNAME_KEY = "ProjectName";

        /// <summary>
        /// Called before application finishes (internal Bye)
        /// </summary>
        void ByeInt() {
            // remove the 'NotTerminated' tag from the session info
            // =====================================================

            // code extra-cautious, since exceptions in Dispose() are, especially in Mono,
            // sometimes not correctly reported and may cause unexplainable segfaults.
            var app = this;
            if (app == null)
                return;

            var csi = app.CurrentSessionInfo;
            IEnumerable<string> tags = csi != null ? csi.Tags : null;
            bool contains_not_terminated = tags != null ? tags.Contains(SessionInfo.NOT_TERMINATED_TAG) : false;
            if (csi != null && tags != null && app.CorrectlyTerminated && contains_not_terminated) {

                Console.WriteLine("Removing tag: " + SessionInfo.NOT_TERMINATED_TAG);
                IList<string> sessTags = tags.ToList();
                sessTags.Remove(SessionInfo.NOT_TERMINATED_TAG);
                this.CurrentSessionInfo.Tags = sessTags;

                if (Tracer.MemtraceFile != null) {
                    try {
                        Tracer.MemtraceFile.Flush();
                        Tracer.MemtraceFile.Close();
                        Tracer.MemtraceFile.Dispose();
                    } catch (IOException) {

                    } finally {
                        Tracer.MemtraceFile = null;
                    }
                }

            }
            if (m_ResLogger != null) {
                m_ResLogger.Close();
                m_ResLogger = null;
            }


        }

        /// <summary>
        /// Called before application finishes. Override this method in order
        /// to implement some finalization-tasks.
        /// </summary>
        protected virtual void Bye() {
        }

        /// <summary>
        /// writes the profiling report 
        /// </summary>
        protected virtual void ProfilingLog() {
            var R = Tracer.Root;

            if (this.DatabaseDriver != null && this.CurrentSessionInfo != null) {
                try {
                    using (Stream stream = this.DatabaseDriver.GetNewLogStream(this.CurrentSessionInfo, "profiling_bin")) {
                        var str = R.Serialize();
                        using (StreamWriter stw = new StreamWriter(stream)) {
                            stw.Write(str);
                            stw.Flush();
                        }

                    }
                } catch (Exception e) {
                    Console.Error.WriteLine(e.GetType().Name + " during writing of profiling_bin: " + e.Message);
                }

                try {
                    using (Stream stream = this.DatabaseDriver.GetNewLogStream(this.CurrentSessionInfo, "profiling_summary")) {
                        using (StreamWriter stw = new StreamWriter(stream)) {
                            WriteProfilingHeader(stw);
                            WriteProfilingReport(stw, R);
                            stw.Flush();
                            stream.Flush();
                            stw.Close();
                        }
                    }
                } catch (Exception e) {
                    Console.Error.WriteLine(e.GetType().Name + " during writing of profiling_summary: " + e.Message);
                }

            }
        }

        private void WriteProfilingHeader(StreamWriter stw ) {
            stw.WriteLine($"Date       : {DateTime.Now}");
            stw.WriteLine($"Computer   : {ilPSP.Environment.MPIEnv.Hostname} ");
            stw.WriteLine($"User Name  : {System.Environment.UserName}");
            stw.WriteLine($"MPI rank   : {this.MPIRank}");
            stw.WriteLine($"MPI size   : {this.MPISize}");

            void TryWrite(string s, Func<object> o) {
                try {
                    stw.WriteLine(s + o());
                } catch(Exception e) {
                    stw.WriteLine(s + $"{e.GetType().Name}, {e.Message}");
                }

            }

            TryWrite("Number of cells (last mesh)      : ", () => this.Grid.NumberOfCells);
            TryWrite("Number of local cells (last mesh): ", () => this.Grid.CellPartitioning.LocalLength);
            if(m_RegisteredFields != null) {
                foreach(var f in m_RegisteredFields) {
                    //                                         :
                    TryWrite("    Field: ", () => $"{f.Identification}, degree {f.Basis.Degree}, XDG: {f.Basis is XDGBasis}");
                }
            }
            Console.WriteLine();
        }


        /// <summary>
        /// Runs this application in parameter study mode which means that it
        /// casts a sequence of runs (using different sessions!) with different
        /// sets of parameters (as specified in the control file).
        /// </summary>
        /// <param name="cases">
        /// enumeration of all parameter study cases
        /// </param>
        /// <param name="opt">%</param>
        /// <param name="ApplicationFactory">%</param>
        static void ParameterStudyModeV2(CommandLineOptions opt, IEnumerable<T> cases, Func<Application<T>> ApplicationFactory) {

            TextWriter log = null;
            QueryResultTable nlog = new QueryResultTable();
            Stream nlog_stream = null;

            if (opt.PstudyCase >= 0
                && opt.PstudyCase < 0 || opt.PstudyCase >= cases.Count()) {
                throw new IndexOutOfRangeException(string.Format("Argument 'pstudy_case' out of range: expected to be in {0} to {1} (both including).", 0, cases.Count() - 1));
            }

            int failCount = 0;
            int numberOfRuns = cases.Count();

            for (int iPstudy = 0; iPstudy < numberOfRuns; iPstudy++) { // loop over cases
                // Run only one case in this process if requested
                if (opt.PstudyCase > 0 && opt.PstudyCase != iPstudy) {
                    continue;
                }

                // select control
                var _control = cases.ElementAt(iPstudy);

                // control file overrides from command line
                if (opt.ProjectName != null)
                    _control.ProjectName = opt.ProjectName;
                if (opt.SessionName != null)
                    _control.SessionName = opt.SessionName;

                if (opt.ImmediatePlotPeriod != null) {
                    _control.ImmediatePlotPeriod = opt.ImmediatePlotPeriod.Value;
                }
                if (opt.SuperSampling != null) {
                    _control.SuperSampling = opt.SuperSampling.Value;
                }

                // ad-hoc added tags (added by command-line option)
                if (opt != null && (opt.TagsToAdd != null && opt.TagsToAdd.Length > 0)) {
                    string[] AdHocTags;
                    AdHocTags = opt.TagsToAdd.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    _control.Tags.AddRange(AdHocTags);
                }

                // add the case index to the 'Paramstudy_CaseIdentification'
                {
                    List<Tuple<string, object>> idl;
                    if (_control.Paramstudy_CaseIdentification != null)
                        idl = _control.Paramstudy_CaseIdentification.ToList();
                    else
                        idl = new List<Tuple<string, object>>();
                    Tuple<string, object> caseId = new Tuple<string, object>("pstudy_case", iPstudy);
                    if (!idl.Contains(caseId, ((Func<Tuple<string, object>, Tuple<string, object>, bool>)((a, b) => a.Item1.Equals(b.Item1))).ToEqualityComparer())) {
                        idl.Add(caseId);
                        _control.Paramstudy_CaseIdentification.AddRange(idl.ToArray());
                    }
                }


                Stopwatch watch = new Stopwatch();
                watch.Start();

                // If everything went fine, this should start the application

                // in solver mode
                Console.WriteLine("Parameter study: Starting run " + (iPstudy) + " of " + numberOfRuns + " (zero-based-index).");
                nlog.CurrentKeyHistory.Clear();
                nlog.UpdateKey("pstudy_case", (iPstudy));
                foreach (var blabla in _control.Paramstudy_CaseIdentification) {
                    nlog.UpdateKey(blabla.Item1, blabla.Item2);
                }


                using (Application<T> app = ApplicationFactory()) {
                    long afterInit = 0;
                    app.m_QueryResultTable = nlog;
                    bool CorrectlyTerminated = false;
#if DEBUG
                    {
#else
                    try {
#endif
                        app.Init(_control);
                        afterInit = watch.ElapsedTicks;
                        app.RunSolverMode();
                        CorrectlyTerminated = true;
                        nlog.LogValue("pstudy_case_successful", true);
                        nlog.LogValue("GrdRes:NumberOfCells", app.Grid.NumberOfCells);
                        if (app.GridData is GridData) {
                            nlog.LogValue("GrdRes:h_min", ((GridData)(app.GridData)).Cells.h_minGlobal);
                            nlog.LogValue("GrdRes:h_max", ((GridData)(app.GridData)).Cells.h_maxGlobal);
                        } else {
                            Console.WriteLine("Warning: unable to obtain grid resolution");
                        }
#if DEBUG
                    }
#else
                    } catch (Exception e) {
                        nlog.LogValue("pstudy_case_successful", false);
                        if (_control.Paramstudy_ContinueOnError) {
                            Console.WriteLine("WARNING: Run" + (iPstudy) + "failed with message '{0}'", e.Message);
                            failCount++;
                        } else {
                            throw;
                        }
                    }
#endif
                    watch.Stop();


                    // Create log during first iteration (happens at this point
                    // because $m_IOFields would undefined before)
                    if (log == null) {
                        log = InitParameterStudyLog(app.m_IOFields, app,
                            out nlog_stream,
                            opt.PstudyCase > 0 ? iPstudy : -1);
                    }

                    double InitTime = (double)afterInit / (double)Stopwatch.Frequency;
                    double solverTime = (double)(watch.ElapsedTicks - afterInit) / (double)Stopwatch.Frequency;
                    nlog.LogValue("InitTime(sec)", InitTime);
                    nlog.LogValue("solverTime(sec)", solverTime);
                    nlog.LogValue("SessionGuid", app.DatabaseDriver.FsDriver == null ? Guid.Empty : app.CurrentSessionInfo.ID);
                    if (app.QueryHandler != null) {
                        foreach (var kv in app.QueryHandler.QueryResults) { // Assume queries have already been evaluated
                            nlog.LogValue(kv.Key, kv.Value);
                        }
                    }

                    // Log only exists on rank 0
                    if (log != null) {

                        // Feed the log
                        if (app.DatabaseDriver.FsDriver == null) {
                            WriteToLog(log, Guid.Empty.ToString(), 40);
                        } else {
                            WriteToLog(log, app.CurrentSessionInfo.ID.ToString(), 40);
                        }
                        WriteToLog(log, InitTime);
                        WriteToLog(log, solverTime);

                        foreach (var blabla in _control.Paramstudy_CaseIdentification) {
                            WriteToLog(log, blabla.Item2 is double ? ((double)blabla.Item2) : double.NaN);
                        }

                        // Assume queries have already been evaluated
                        if (app.QueryHandler != null) {
                            foreach (var kv in app.QueryHandler.QueryResults) {
                                WriteToLog(log, (double)kv.Value);
                            }
                        }
                        log.WriteLine();
                        log.Flush();

                        if (nlog_stream != null) {
                            nlog_stream.Position = 0; // we can only write the whole table, so we have to overwrite the so-far-written stuff
                            var w = new StreamWriter(nlog_stream);
                            nlog.WriteToStream(w);
                            w.Flush();
                        }
                    }

                    // log to session directory
                    if (app.MPIRank == 0 && !app.CurrentSessionInfo.ID.Equals(Guid.Empty)) {
                        var nlog_stream_session = app.DatabaseDriver.GetNewLogStream(
                            app.CurrentSessionInfo, "ParameterStudy.case-" + iPstudy);
                        try {
                            var w = new StreamWriter(nlog_stream_session);
                            nlog.WriteToStream(w, RowFilter: nlog.CurrentKeyHistory);
                            w.Flush();
                        }
                        finally {
                            nlog_stream_session.Flush();
                            nlog_stream_session.Close();
                        }
                    }

                    // finalize
                    Console.WriteLine("Parameter study run " + iPstudy + " successful: " + CorrectlyTerminated);
#if DEBUG
                    {
#else
                    try {
#endif

                        app.ByeInt();
                        app.Bye();
                        app.ProfilingLog();
#if DEBUG
                    }
#else
                    } catch (Exception e) {
                        nlog.LogValue("pstudy_case_successful", false);
                        if (_control.Paramstudy_ContinueOnError) {
                            Console.WriteLine("WARNING: Run" + (iPstudy) + "failed with message '{0}'", e.Message);
                            failCount++;
                        } else {
                            throw;
                        }
                    }
#endif
                }

                System.GC.Collect();
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }

            string failedRunsWarning = string.Format(
                "WARNING: {0} run(s) failed during execution. See program output for more information.",
                failCount);
            if (failCount > 0) {
                Console.WriteLine(failedRunsWarning);
            }

            if (log != null) {
                if (failCount > 0) {
                    log.WriteLine(failedRunsWarning);
                }
                log.Close();
            }
        }

        bool CorrectlyTerminated = false;

        static private StreamWriter InitParameterStudyLog(
            ICollection<DGField> ioFields, Application<T> app,
            out Stream BinIOStream,
            int counter) {
            // Only do something on rank 0
            if (ilPSP.Environment.MPIEnv.MPI_Rank != 0) {
                BinIOStream = null;
                return null;
            }

            // Prepare log file using specified path
            StreamWriter log = null;
            {
                string basePath;
                if (app.Control != null)
                    basePath = app.Control.logFileDirectory;
                else
                    throw new ApplicationException();
                string baseName = "ParameterStudy";
                if (counter >= 0)
                    baseName += ".case-" + counter.ToString();
                string nextName = baseName;
                int i = 1;

                BinIOStream = null;
                while (log == null) {
                    string pathToFile = Path.Combine(basePath, nextName + ".txt");
                    if (!(File.Exists(pathToFile))) {
                        try {
                            var fs = new FileStream(pathToFile, FileMode.CreateNew, FileAccess.Write);
                            log = new StreamWriter(fs);
                            BinIOStream = new FileStream(Path.Combine(basePath, nextName + ".table"), FileMode.CreateNew, FileAccess.Write);
                        } catch (Exception) {
                            log = null;
                        }
                    }

                    // generate next filename
                    nextName = baseName + "." + i;
                    i++;
                }
            }

            // Write header
            WriteToLog(log, "Session", 40);
            WriteToLog(log, "Setup time (s)");
            WriteToLog(log, "Solution time (s)");

            if (app.Control != null) {
                foreach (string id in app.Control.Paramstudy_CaseIdentification.Select(t => t.Item1)) {
                    WriteToLog(log, id);
                }
            } else {
                throw new ApplicationException("should not occur.");
            }

            //foreach (string id in app.m_PerformanceFigure.Keys) {
            //    WriteToLog(log, id);
            //}

            //foreach (Field field in ioFields) {
            //    WriteToLog(log, field.Identification + "_L2");
            //    WriteToLog(log, field.Identification + "_Linf");
            //}

            // Process exact solutions
            //Dictionary<string, AppControl._Base> exactSolutions =
            //    app.m_control.confParameterStudy.ExactSolutions;

            //Dictionary<string, Guid> referenceSolutions =
            //    app.m_control.confParameterStudy.ReferenceSolutions;
            //foreach (Field field in ioFields) {
            //    if (exactSolutions.ContainsKey(field.Identification)
            //        || referenceSolutions.ContainsKey(field.Identification)) {
            //        WriteToLog(log, field.Identification + "_error_L2");
            //        WriteToLog(log, field.Identification + "_error_Linf");
            //        WriteToLog(log, field.Identification + "_relError_L2");
            //        WriteToLog(log, field.Identification + "_relError_Linf");
            //    }
            //}
            log.WriteLine();
            log.Flush();

            return log;
        }

        static private void WriteToLog(TextWriter log, string value, int padding) {
            if (log != null) {
                log.Write(value.PadRight(padding) + "\t");
            }
        }

        static private void WriteToLog(TextWriter log, string value) {
            WriteToLog(log, value, 17);
        }

        static private void WriteToLog(TextWriter log, int value) {
            WriteToLog(log, value.ToString());
        }

        static private void WriteToLog(TextWriter log, double value) {
            string valueString = value.ToString(
                "E",
                System.Globalization.CultureInfo.InvariantCulture);
            WriteToLog(log, valueString);
        }

        /// <summary>
        /// loads a field from the database
        /// </summary>
        protected virtual void LoadField(ITimestepInfo tsi, string fieldName, string newFieldName = null) {
            using (new ilPSP.Tracing.FuncTrace()) {
                DGField field = DatabaseDriver.LoadFields(tsi, (GridData)GridData, new[] { fieldName }).Single();
                field.Identification = newFieldName ?? fieldName;
                m_IOFields.Add(field);
            }
        }

        /// <summary>
        /// performs restart as specified by <see cref="Control"/>;
        /// </summary>
        /// <param name="Time">
        /// on exit, the physical time associated with the field state
        /// </param>
        /// <param name="TimestepNo">
        /// on exit, the physical time associated with the field state
        /// </param>
        protected virtual void LoadRestart(out double Time, out TimestepNumber TimestepNo) {
            Time = 0;
            TimestepNo = 0;

            if (MPIRank == 0) {
                Console.WriteLine("----------------------------------------");
                Console.WriteLine("Loading Restart...");
            }


            if (this.Control != null && this.Control.RestartInfo != null) {
                //if (!this.Control.InitialValues_Evaluators.IsNullOrEmpty())
                //    throw new ApplicationException("control object error: initial values ('AppControl.InitialValues') and restart info ('AppControl.RestartInfo') cannot be specified at the same time.");

                TimestepNo = RestartFromDatabase(out Time);
                this.CurrentSessionInfo.RestartedFrom = this.Control.RestartInfo.Item1;
                this.CurrentSessionInfo.Save();

            } else {
                throw new ApplicationException("restart was not specified");
                //throw new NotImplementedException("unknown restart type.");
            }

            PostRestart(Time, TimestepNo);

            if(this.LsTrk != null) {
                if(this.LsTrk.Regions.Time != Time)
                    this.LsTrk.UpdateTracker(Time);
            }

            System.GC.Collect();
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

            if (MPIRank == 0) {
                Console.WriteLine("Finished loading restart.");
                Console.WriteLine("----------------------------------------");
            }

        }

        /// <summary>
        /// override this method to implement any user-specific tasks which
        /// should be carried out after a restart file has been loaded (e.g.,
        /// setting the correct time for a time-stepper)
        /// </summary>
        public virtual void PostRestart(double time, TimestepNumber timestep) {
        }

        /// <summary>
        /// Override this method in order to construct a
        /// <see cref="QueryHandler"/> specific to your application. This might
        /// be necessary if you want to user application-specific queries in
        /// the control file (see <see cref="QueryHandler"/>)
        /// </summary>
        /// <returns>
        /// A new instance of <see cref="QueryHandler"/>
        /// </returns>
        protected virtual QueryHandler QueryHandlerFactory() {
            return new QueryHandler(this);
        }

        /// <summary>
        /// true after <see cref="Dispose"/> has been called
        /// </summary>
        public bool IsDisposed {
            get;
            private set;
        }

        /// <summary>
        /// disposes database driver (<see cref="IFileSystemDriver"/>-instance) if necessary and calls MPI_Finalize();
        /// </summary>
        public virtual void Dispose() {
            if (!IsDisposed) {
                try {
                    ByeInt();
                    Bye();
                    ProfilingLog();

                    foreach( var l in PostprocessingModules) {
                        l.Dispose();
                    }

                    if (this.CurrentSessionInfo != null)
                        this.CurrentSessionInfo.Dispose();
                    if (DatabaseDriver != null) {
                        DatabaseDriver.Dispose();
                    }

                    if (m_Database != null) {
                        DatabaseInfo.Close(m_Database);
                    }

                    Console.Out.Flush();
                    Console.Error.Flush();

                } catch (Exception e) {
                    Console.WriteLine(e.GetType().Name + " in Dispose() " + e.Message);
                    Console.WriteLine(e.StackTrace);
                }
                IsDisposed = true;
            }
        }

        /// <summary>
        /// returns the size of the fixed timestep (<see cref="AppControl.dtFixed"/>); if a variable timestep is set, this method throws an exception.
        /// </summary>
        virtual public double GetTimestep() {
            if (this.Control != null) {
                return this.Control.dtFixed;
            } else {
                throw new NotSupportedException(
                    "Cannot get time-step from control object, since there is no control object specified.");
            }
        }

        /// <summary>
        /// The solver residual logger for this application.
        /// </summary>
        public ResidualLogger ResLogger {
            get {
                return m_ResLogger;
            }
        }

        ResidualLogger m_ResLogger;

        /// <summary>
        /// creates a human-readable performance report from the profiling information stored in <see cref="Tracer.Root"/>.
        /// </summary>
        public static void WriteProfilingReport(TextWriter wrt, MethodCallRecord Root) {
            var R = Root;

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


        /// <summary>
        /// This method should be overridden to support automatic numerical stability analysis of the PDE's operator
        /// </summary>
        /// <returns>
        /// Pairs of property name and value, e.g. ConditionNumber and the respective value of the operators Jacobian matrix condition number.
        /// </returns>
        virtual public IDictionary<string, double> OperatorAnalysis() {
            throw new NotImplementedException();
        }


        /// <summary>
        /// This method just forces the C# - compiler to integrate
        /// BoSSS.Foundation.Grid.dll into the manifest of BoSSS.Solution.dll;
        /// Otherwise, without this useless method BoSSS.Foundation.Grid.dll
        /// would be loaded delayed (when de-serializing the Grid),
        /// in contrast to other dependencies like BoSSS.Foundation.dll which
        /// are loaded at start-up time. The direct inclusion of
        /// BoSSS.Foundation.Grid.dll into the manifest of BoSSS.Solution.dll  
        /// helps Visual Studio and other tools, e.g. "bcl deploy-at" to find
        /// and copy really all dependencies to the binary directory.
        /// </summary>
        /// <returns>
        /// A useless double array.
        /// </returns>
        /// <remarks>
        /// Must <b>not</b> be private in order to ensure that faking also
        /// works for sub-classes
        /// </remarks>
        protected double[] Fake() {
            return GenericBlas.Linspace(-1, 1, 2);
        }

        /// <summary>
        /// enforce the compiler to integrate Microsoft.CodeAnalysis.dll etc.
        /// </summary>
        public static Type[] DllEnforcer() {
            using(var tr = new FuncTrace()) {
                var types = new Type[] {
                    typeof(Microsoft.CodeAnalysis.Compilation),
                    typeof(Microsoft.CodeAnalysis.CSharp.CSharpCompilation),
                    typeof(Microsoft.CodeAnalysis.Scripting.Script),
                    typeof(Microsoft.CodeAnalysis.CSharp.Scripting.CSharpScript),
                    typeof(System.Configuration.Configuration)
                };

                foreach(var t in types) {
                    tr.Info("Loaded type " + t + " form " + t.Assembly);
                }

                return types;
            }
        }

        ///// <summary>
        ///// 
        ///// </summary>
        //public static Microsoft.CodeAnalysis.Compilation DllEnforcer2() {
        //    return Microsoft.CodeAnalysis.CSharp.CSharpCompilation.C
        //}
    }
}

