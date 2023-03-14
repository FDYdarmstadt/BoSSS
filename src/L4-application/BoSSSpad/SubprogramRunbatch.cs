using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using CommandLine;
using CommandLine.Text;
using ilPSP;

namespace BoSSS.Application.BoSSSpad {
    
    /// <summary>
    /// Routines to deploy and run on some batch queue from command line 
    /// </summary>
    static class SubprogramRunbatch {

        /// <summary>
        /// Variables parsed form command line using <see cref="CommandLine.Parser.ParseArguments{T}(IEnumerable{string})"/>.
        /// </summary>
        [Serializable]
        public sealed class CommandLineOptions {

           
            /*
            /// <summary>
            /// help
            /// </summary>
            [HelpOption(HelpText = "Displays this help screen.")]
            public string GetUsage() {
                var help = new HelpText("BoSSS");
                help.AdditionalNewLineAfterOption = true;
                help.Copyright = new CopyrightInfo("Fachgebiet fuer Stroemungsdynamik, TU Darmstadt", 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021);
                
                help.AddPreOptionsLine("For available queues (aka. HPC servers, clusters, batch processors), see at your ~/.BoSSS/etc/BatchProcessorConfig.json file. ");
                help.AddOptions(this);

                return help;
            }
            */

            /*
            /// <summary>
            /// Immediate plot period: This variable controls immediate
            /// plotting, i.e. plotting during the solver run.<br/>
            /// A positive value indicates that
            /// <see cref="Application{T}.PlotCurrentState"/>"/> will be called every
            /// <see cref="ImmediatePlotPeriod"/>-th time-step;<br/>
            /// A negative value turns immediate plotting off;
            /// </summary>
            [Option("i", "implt", HelpText = "Period (in number of timesteps) for immediate plot output.")]
            public int? ImmediatePlotPeriod;

            /// <summary>
            /// Super sampling: This option controls whether a finer grid
            /// resolution shall be used in the plots created via --implt.
            /// </summary>
            [Option("u", "supersampling", HelpText = "Super sampling (recursive cell divisions) for --implt output")]
            public int? SuperSampling;
            

            /// <summary>
            /// deletion of plot files
            /// </summary>
            [Option("d", "delplt", HelpText = "if set, all *plt files are deleted from the output directory.")]
            public bool delPlt = false;
            */

            /// <summary>
            /// Override for the project name in the control file (<see cref="AppControl.ProjectName"/>), resp. the 
            /// session information (<see cref="ISessionInfo.ProjectName"/>).
            /// </summary>
            [Option('P', "prjnmn", HelpText = "The project name assigned to the job; this is only an aid for post-processing to better filter for jobs in the workflow management.")]
            public string ProjectName { get; set; }


            /// <summary>
            /// Optional name for a computation, override to <see cref="AppControl.SessionName"/>.
            /// </summary>
            [Option('N', "sesnmn", HelpText = "optional name for the compute session.")]
            public string SessionName { get; set; }

            /// <summary>
            /// Tracing of all namespaces on.
            /// </summary>
            [Option('F', "fulltracing", HelpText = "Mainly for debugging purpose, turns tracing of all namespaces on, trace-files will be created.")]
            public bool fullTracing { get; set; }

            /// <summary>
            /// path to control file
            /// </summary>
            [Option('c', "control", HelpText = "path to control file - or  - when starting with the prefix 'cs:', a single line of C#-code which results in a control object.", Required = true)]
            public string ControlfilePath { get; set; }

            /// <summary>
            /// tags which will be added to the session information.
            /// </summary>
            [Option('t', "tags", HelpText = "tags which will be added to the session information when saved in the database, separated by comma (',').")]
            public string TagsToAdd { get; set; }

            /// <summary>
            /// no of MPI procs
            /// </summary>
            [Option('n', "noofMPIprocs", HelpText = "number of MPI processes requested.", Required = true)]
            public int np { get; set; }

            /// <summary>
            /// <see cref="BatchProcessorClient.Name"/>, <see cref="BoSSSshell.ExecutionQueues"/>
            /// </summary>
            [Option('q', "queue", HelpText ="Name or (zero-based) index of the requested batch processor in ~/.BoSSS/etc/BatchProcessorConfig.json", Required = true)]
            public string queue { get; set; }

            /// <summary>
            /// Name of the solver to run.
            /// </summary>
            [Option('s', "solver", HelpText = "Name of the solver to run.", Required = true)]
            public string SolverName { get; set; }
        }
        
        static BatchProcessorClient GetQueue(string id) {
            if(id.IsEmptyOrWhite())
                return BoSSSshell.GetDefaultQueue();

            if(int.TryParse(id, out int idx)) {
                return BoSSSshell.ExecutionQueues[idx];

            } else {
                return BoSSSshell.ExecutionQueues.Where(client => client.Name != null && client.Name.Contains(id)).Single();
            }

        }


        public static int RunBatch(string[] args) {


            //CommandLineOptions opt = new CommandLineOptions();
            //ICommandLineParser parser = new CommandLine.CommandLineParser(new CommandLineParserSettings(Console.Error));
            //bool argsParseSuccess = parser.ParseArguments(args, opt);
            var CmdlineParseRes = Parser.Default.ParseArguments<CommandLineOptions>(args);
            bool argsParseSuccess = !CmdlineParseRes.Errors.Any();
            CommandLineOptions opt = CmdlineParseRes.Value;

            if (!argsParseSuccess) {
                Console.Error.WriteLine("Unable to parse arguments.");
                return -41;
            }

            BatchProcessorClient queue = GetQueue(opt.queue);

            string jobname = opt.SessionName ?? "LousyJob";
            string project = opt.ProjectName ?? "NixProject";


            BoSSSshell.WorkflowMgm.Init(project, queue);

            if(opt.np < 1) {
                throw new ArgumentException("Number of MPI processors must be a positive number (greater or equal than 1).");
            }
          

            string solverName = opt.SolverName;
            Type solver = GetType(solverName); 

            int np = opt.np;
            
            var j = new Job(jobname, solver);

            j.MySetCommandLineArguments("-c", opt.ControlfilePath);

            j.NumberOfMPIProcs = np;


            j.Activate(queue);


            Console.WriteLine("Job Status: " + j.Status);

            return 0;
        }


        static private Type GetType(string partialName) {
            List<Type> foundTypes = new List<Type>();

            var allAssis = new HashSet<Assembly>();
            GetAllAssemblies(typeof(BoSSSpadMain).Assembly, allAssis);

            foreach(var a in allAssis) {
                try {
                    foreach(var t in a.GetTypes()) {
                        //if(t.FullName.StartsWith("BoSSS.Application.XNSE_Solver.XNSE"))
                        //    Console.WriteLine(t.FullName);
                        if(IsApplication(t) && t.Name.EndsWith(partialName, StringComparison.InvariantCultureIgnoreCase)) {
                            Console.WriteLine("Found solver: " + t.Name + "    (" + t.FullName + ")");
                            foundTypes.Add(t);
                        }
                    }
                } catch(Exception) { }
            }

            if(foundTypes.Count <= 0) {
                throw new ArgumentException($"found no solver containing {partialName}");
            }
            if(foundTypes.Count > 1) {
                string AllMatch = "";
                for(int i = 0; i < foundTypes.Count; i++) {
                    AllMatch += foundTypes[i];
                    if(i < foundTypes.Count - 1)
                        AllMatch += ", ";
                }
                
                throw new ArgumentException($"found {foundTypes.Count} matches for {partialName} -- be more specific (matches are: {AllMatch}).");
            }

            return foundTypes[0];

        }


        /// <summary>
        /// true, if <paramref name="t"/> is derived from <see cref="Solution.Application{T}"/> or implements <see cref="Solution.IApplication"/>
        /// </summary>
        static bool IsApplication(Type t) {
            var bbase = typeof(Solution.Application<>);
            var ibase = typeof(Solution.IApplication);

            if(t.Equals(bbase))
                return true;

            if(t.GetInterfaces().Contains(ibase))
                return true;

            if(t.IsGenericType && t.GetGenericTypeDefinition() == bbase)
                return true;


            if(t.BaseType != null)
                return IsApplication(t.BaseType);
            else
                return false;
        }



        /// <summary>
        /// Recursive collection of all dependencies of some assembly.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="assiList">
        /// Output, list where all dependent assemblies are collected.
        /// </param>
        private static void GetAllAssemblies(Assembly a, HashSet<Assembly> assiList) {
            if(a.IsDynamic)
                return;
            if(assiList.Contains(a))
                return;
            assiList.Add(a);

            //string fileName = Path.GetFileName(a.Location);
            //var allMatch = assiList.Where(_a => Path.GetFileName(_a.Location).Equals(fileName)).ToArray();
            //if(allMatch.Length > 1) {
            //    throw new ApplicationException("internal error in assembly collection.");
            //}


            foreach(AssemblyName b in a.GetReferencedAssemblies()) {
                Assembly na;
                try {
                    na = Assembly.Load(b);
                } catch(FileNotFoundException) {
                    //string[] AssiFiles = ArrayTools.Cat(Directory.GetFiles(SearchPath, b.Name + ".dll"), Directory.GetFiles(SearchPath, b.Name + ".exe"));
                    //if(AssiFiles.Length != 1) {
                    //    //throw new FileNotFoundException("Unable to locate assembly '" + b.Name + "'.");
                    //    Console.WriteLine("Skipping: " + b.Name);
                    //    continue;
                    //}
                    //na = Assembly.LoadFile(AssiFiles[0]);
                    continue;
                }

                GetAllAssemblies(na, assiList);
            }
        }
        

    }
}
