using CommandLine;
using CommandLine.Text;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {
         
        /// <summary>
        /// Variables parsed form command line using <see cref="CommandLine.CommandLineParser"/>.
        /// </summary>
        [Serializable]
        public sealed class CommandLineOptions {

            /// <summary>
            /// Case number for parameter study
            /// </summary>
            [Option("p", "pstudy_case", HelpText = "case number for parameter study; if not specified, or non-positive, all runs are done in the same process")]
            public int PstudyCase = -666;

            /// <summary>
            /// path to control file
            /// </summary>
            [Option("c", "control", HelpText = "path to control file", MutuallyExclusiveSet = "control,session")]
            public string ControlfilePath = null;

            /// <summary>
            /// tags which will be added to the session information.
            /// </summary>
            [Option("t", "tags", HelpText = "tags which will be added to the session information when saved in the database, separated by comma (',').")]
            public string TagsToAdd = null;

            /// <summary>
            /// help
            /// </summary>
            [HelpOption(HelpText = "Displays this help screen.")]
            public string GetUsage() {
                var help = new HelpText("BoSSS");
                help.AdditionalNewLineAfterOption = true;
                help.Copyright = new CopyrightInfo("Fachgebiet fuer Stroemungsdynamik, TU Darmstadt", 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018);
                //help.AddPreOptionsLine("This is free software. You may redistribute copies of it under the terms of");
                //help.AddPreOptionsLine("the MIT License <http://www.opensource.org/licenses/mit-license.php>.");
                //help.AddPreOptionsLine("Usage: SampleApp -rMyData.in -wMyData.out --calculate");
                //help.AddPreOptionsLine(string.Format("       SampleApp -rMyData.in -i -j{0} file0.def file1.def", 9.7));
                //help.AddPreOptionsLine("       SampleApp -rMath.xml -wReport.bin -o *;/;+;-");
                help.AddPreOptionsLine("Bla, Bla, Bla");
                help.AddOptions(this);

                return help;
            }

            /// <summary>
            /// Immediate plot period: This variable controls immediate
            /// plotting, i.e. plotting during the solver run.<br/>
            /// A positive value indicates that
            /// <see cref="Application{T}.PlotCurrentState(double, TimestepNumber, int)"/>"/> will be called every
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

            /// <summary>
            /// Override for the project name in the control file (<see cref="AppControl.ProjectName"/>), resp. the 
            /// session information (<see cref="ISessionInfo.ProjectName"/>).
            /// </summary>
            [Option("P", "prjnmn", HelpText = "overrides the project name.")]
            public string ProjectName = null;


            /// <summary>
            /// Optional name for a computation, override to <see cref="AppControl.SessionName"/>.
            /// </summary>
            [Option("n", "sesnmn", HelpText = "optional name for the compute session.")]
            public string SessionName = null;

        }
}
