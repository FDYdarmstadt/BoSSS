using CommandLine;
using CommandLine.Text;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {

    /// <summary>
    /// Variables parsed form command line using <see cref="CommandLine.Parser"/>.
    /// </summary>
    [Serializable]
    public sealed class CommandLineOptions {

        
        /// <summary>
        /// Case number for parameter study
        /// </summary>
        [Option('p', "pstudy_case", Default = -666, HelpText = "case number for parameter study; if not specified, or smaller than zero, all runs are done in the same process")]
        public int PstudyCase { get; set; }

        /// <summary>
        /// number of threads per MPI process
        /// </summary>
        [Option('T', "num_threads", Default = null, HelpText = "number of threads (OpenMP and C# Task Parallel Library) for each MPI process")]
        public int? NumThreads { get; set; }


        /// <summary>
        /// path to control file
        /// </summary>
        [Option('c', "control", HelpText = "path to control file - or  - when starting with the prefix 'cs:', a single line of C#-code which results in a control object")]
        public string ControlfilePath { get; set; }

        /// <summary>
        /// tags which will be added to the session information.
        /// </summary>
        [Option('t', "tags", HelpText = "tags which will be added to the session information when saved in the database, separated by comma (',').")]
        public string TagsToAdd { get; set; }

        /*
        /// <summary>
        /// help
        /// </summary>
        [CommandLine.Text.HelpText. HelpOption(HelpText = "Displays this help screen.")]
        public string GetUsage() {
            var help = new HelpText("BoSSSpad");
            help.AdditionalNewLineAfterOption = true;
            help.Copyright = new CopyrightInfo("Fachgebiet fuer Stroemungsdynamik, TU Darmstadt", 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021);
            help.AddPreOptionsLine("This is free software. You may redistribute copies of it under the terms of");
            help.AddPreOptionsLine("the Apache License <http://www.apache.org/licenses/LICENSE-2.0>.");
            //help.AddPreOptionsLine("Usage: SampleApp -rMyData.in -wMyData.out --calculate");
            //help.AddPreOptionsLine(string.Format("       SampleApp -rMyData.in -i -j{0} file0.def file1.def", 9.7));
            //help.AddPreOptionsLine("       SampleApp -rMath.xml -wReport.bin -o *;/;+;-");
            //help.AddOptions()

            return help;
        }
        */

        /// <summary>
        /// Immediate plot period: This variable controls immediate
        /// plotting, i.e. plotting during the solver run.<br/>
        /// A positive value indicates that
        /// <see cref="Application{T}.PlotCurrentState"/>"/> will be called every
        /// <see cref="ImmediatePlotPeriod"/>-th time-step;<br/>
        /// A negative value turns immediate plotting off;
        /// </summary>
        [Option('i', "implt", HelpText = "Period (in number of timesteps) for immediate plot output.")]
        public int? ImmediatePlotPeriod { get; set; }

        /// <summary>
        /// Super sampling: This option controls whether a finer grid
        /// resolution shall be used in the plots created via --implt.
        /// </summary>
        [Option('u', "supersampling", HelpText = "Super sampling (recursive cell divisions) for --implt output")]
        public int? SuperSampling { get; set; }

        /// <summary>
        /// deletion of plot files
        /// </summary>
        [Option('d', "delplt", HelpText = "if set, all *plt files are deleted from the output directory.")]
        public bool delPlt { get; set; }

        /// <summary>
        /// Override for the project name in the control file (<see cref="Control.AppControl.ProjectName"/>), resp. the 
        /// session information (<see cref="BoSSS.Foundation.IO.ISessionInfo.ProjectName"/>).
        /// </summary>
        [Option('P', "prjnmn", HelpText = "overrides the project name.")]
        public string ProjectName { get; set; }


        /// <summary>
        /// Optional name for a computation, override to <see cref="Control.AppControl.SessionName"/>.
        /// </summary>
        [Option('n', "sesnmn", HelpText = "optional name for the compute session.")]
        public string SessionName { get; set; }


        /// <summary>
        /// Tracing of all namespaces on.
        /// </summary>
        [Option('F', "fulltracing", HelpText = "Mainly for debugging purpose, turns tracing of all namespaces on, trace-files will be created.")]
        public bool fullTracing { get; set; }
    }
}
