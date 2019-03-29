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

using ilPSP;
using BoSSS.Platform;
using BoSSS.Foundation.IO;
using ilPSP.Connectors.Matlab;
using Mono.CSharp;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading;
using System.Diagnostics;
using System.Text.RegularExpressions;
using ilPSP.LinSolvers;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Console version of the BoSSSpad.
    /// </summary>
    static public class ReadEvalPrintLoop {

        private static string[] exitCommands = new[] { "exit", "quit", "q" };
        
        /// <summary>
        /// Reports errors and stuff like that.
        /// </summary>
        public static CompilerContext cmpCont {
            get;
            private set;
        }

        /// <summary>
        /// The current evaluator, which evaluates C# commands.
        /// </summary>
        public static SafeEvaluator eval {
            get;
            private set;
        }
             


        /// <summary>
        /// Creates the evaluator using <see cref="InteractiveShell"/> as
        /// <see cref="Evaluator.InteractiveBaseClass"/>
        /// </summary>
        /// <returns></returns>
        private static Evaluator GetEvaluator() {
            var Settings = new CompilerSettings();
#if DEBUG
            Settings.Optimize = false;
#else
            Settings.Optimize = false;
#endif

            cmpCont = new CompilerContext(
                Settings, new ConsoleReportPrinter());
            Evaluator eval = new Evaluator(cmpCont);
            eval.InteractiveBaseClass = typeof(InteractiveShell);

            return eval;
        }

        /// <summary>
        /// Executes start-up commands
        /// </summary>
        /// <param name="startupCommands">Optional startup commands.</param>
        private static Evaluator Startup(string startupCommands) {
            Evaluator eval = GetEvaluator();


            // Assembly References
            eval.ReferenceAssembly(typeof(System.Data.DataTable).Assembly); // required for session tables
            eval.ReferenceAssembly(typeof(ilPSP.Environment).Assembly);
            eval.ReferenceAssembly(typeof(ilPSP.LinSolvers.SimpleSolversInterface).Assembly);
            eval.ReferenceAssembly(typeof(BatchmodeConnector).Assembly); // Do it this cause connector is not referenced anywhere else, i.e. the assembly will often be missing otherwise
            eval.ReferenceAssembly(typeof(NUnit.Framework.Assert).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.PlotGenerator.PlotApplication).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Platform.Utils.Geom.BoundingBox).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Foundation.Basis).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Foundation.XDG.XDGField).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Foundation.Grid.Classic.Grid1D).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.Application).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.Gnuplot.Gnuplot).Assembly);
			eval.ReferenceAssembly(typeof(BoSSS.Solution.GridImport.Cgns).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.Statistic.CellLocalization).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.Tecplot.Tecplot).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.ASCIIExport.CurveExportDriver).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.AdvancedSolvers.MultigridOperator).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.XNSECommon.CurvatureAlgorithms).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.XdgTimestepping.LevelSetHandling).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Solution.LevelSetTools.ContinuityProjection).Assembly);
            eval.ReferenceAssembly(typeof(BoSSSpad.BoSSSpadMain).Assembly);
            eval.ReferenceAssembly(typeof(Renci.SshNet.SftpClient).Assembly);
            eval.ReferenceAssembly(typeof(MiniBatchProcessor.Client).Assembly);
            eval.ReferenceAssembly(typeof(System.Numerics.Complex).Assembly);
            eval.ReferenceAssembly(typeof(Mono.CSharp.Evaluator).Assembly);
            eval.ReferenceAssembly(typeof(CNS.Program).Assembly);
            eval.ReferenceAssembly(typeof(IBM_Solver.IBM_SolverMain).Assembly);
            eval.ReferenceAssembly(typeof(XNSE_Solver.XNSE_SolverMain).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Application.SipPoisson.SipPoissonMain).Assembly);
            eval.ReferenceAssembly(typeof(Rheology.Rheology).Assembly);
            eval.ReferenceAssembly(typeof(FSI_Solver.FSI_SolverMain).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Foundation.SpecFEM.SpecFemField).Assembly);
            eval.ReferenceAssembly(typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main).Assembly);
          //  eval.ReferenceAssembly(typeof(BoSSS.Application.LowMachCombustionNSE.LowMachCombustionNSE).Assembly);

            //eval.ReferenceAssembly(typeof(FSI_Solver.FSI_SolverMain).Assembly);
            //eval.ReferenceAssembly(typeof(FuelCell.FuelCellMain).Assembly);
            // Helical shit
            eval.ReferenceAssembly(typeof(StokesHelical.HelicalMain).Assembly);
            // eval.ReferenceAssembly(typeof(PosissonScalar3CylinderCoords.PoissonScalar3CCMain).Assembly);
            eval.ReferenceAssembly(typeof(PoissonScalar2CylinderCoords.PoissonScalar2CCMain).Assembly);


            eval.Compile(
                "using System;" + Console.Out.NewLine +
                "using System.Collections.Generic;" + Console.Out.NewLine +
                "using System.Linq;" + Console.Out.NewLine +
                "using ilPSP;" + Console.Out.NewLine +
                "using ilPSP.Utils;" + Console.Out.NewLine +
                "using BoSSS.Platform;" + Console.Out.NewLine +
                "using BoSSS.Foundation;" + Console.Out.NewLine +
                "using BoSSS.Foundation.Grid;" + Console.Out.NewLine +
                "using BoSSS.Foundation.Grid.Classic;" + Console.Out.NewLine +
                "using BoSSS.Foundation.IO;" + Console.Out.NewLine +
                "using BoSSS.Solution;" + Console.Out.NewLine +
                "using BoSSS.Solution.Control; " + Console.Out.NewLine +
                "using BoSSS.Solution.GridImport;" + Console.Out.NewLine +
                "using BoSSS.Solution.Statistic;" + Console.Out.NewLine +
                "using BoSSS.Solution.Utils;" + Console.Out.NewLine +
                "using BoSSS.Solution.Gnuplot;" + Console.Out.NewLine +
                "using BoSSS.Application.BoSSSpad;" + Console.Out.NewLine +
                "using Renci.SshNet;" + Console.Out.NewLine +
                "using Mono.CSharp;"
            );
            

            // Make sure user is asked for password where required
            DatabaseController.PasswordCallback = PromptPassword;

            // Load Databases and DBController. Run console
            using (StringWriter stw = new StringWriter()) {
                ilPSP.Environment.StdOut.WriterS.Add(stw);
                ilPSP.Environment.StdErr.WriterS.Add(stw);
                eval.Run(@"
                    Console.WriteLine("""");
                    Console.WriteLine(""  BoSSSpad C# interpreter"");
                    Console.WriteLine(""  _______________________\n"");

                
                    try {
                        databases = DatabaseController.LoadDatabaseInfosFromXML();
                    
                        string summary = databases.Summary();
                        Console.WriteLine(""Databases loaded:"");
                        Console.WriteLine(summary);
                    } catch (Exception e) {
                        Console.WriteLine();
                        Console.WriteLine(
                            ""{0} occurred with message '{1}' while loading the databases. Type 'LastError' for details."",
                            e.GetType(),
                            e.Message);
                        InteractiveShell.LastError = e;
                    }
                ");

                // Run startup commands
                try {
                    using (StringReader reader = new StringReader(startupCommands)) {
                        string multiline = null;
                        int lineno = 0;
                        for (string line = reader.ReadLine(); line != null; line = reader.ReadLine()) {
                            line = line.TrimEnd();
                            lineno++;

                            if (line.EndsWith("\\")) {
                                if (multiline == null)
                                    multiline = "";

                                multiline += " " + line.Substring(0, line.Length - 1);
                            } else {
                                string completeline;
                                if (multiline == null) {
                                    completeline = line;
                                } else {
                                    completeline = multiline + " " + line;
                                    multiline = null;
                                }

                                try {
                                    eval.Run(completeline);
                                } catch (Exception e) {
                                    throw new AggregateException(
                                        e.GetType().Name + " during the interpretation of DBErc file"
                                        + " line "
                                        + lineno
                                        + "; ",
                                        e);
                                }
                            }
                        }
                    }
                } catch (Exception e) {
                    Console.WriteLine("Running startup commands failed with message '{0}'", e.Message);
                }

                eval.Run(@"Console.WriteLine(""\n Console ready for input. Type 'help' for help."");");

                // Log results of startup
                Console.Out.Flush();
                Console.Error.Flush();
                ilPSP.Environment.StdOut.WriterS.Remove(stw);
                ilPSP.Environment.StdErr.WriterS.Remove(stw);
                Log.Add(new Tuple<string, string>("restart", stw.ToString()));

                return eval;
            }
        }

        /// <summary>
        /// Executes the REPL (Read-Eval-Print-Loop) until terminated.
        /// </summary>
        public static void REPL() {
            EvalPrint("restart", out var dummy1);

            CommandLineReader reader = GetCommandLineReader();

            while (eval != null) {
                Console.WriteLine();
                string line = reader.ReadCommand("> ", "").Trim();

                if (line == null || exitCommands.Contains(line)) {
                    break;
                }

                EvalPrint(line, out var dummy2);
            }
        }

        /// <summary>
        /// Removes all comments from a C#-statement.
        /// </summary>
        public static string RemoveAllComments(string input) {
            // http://stackoverflow.com/questions/3524317/regex-to-strip-line-comments-from-c-sharp/3524689

            var blockComments = @"/\*(.*?)\*/";
            var lineComments = @"//(.*?)\r?\n";
            var strings = @"""((\\[^\n]|[^""\n])*)""";
            var verbatimStrings = @"@(""[^""]*"")+";


            string noComments = Regex.Replace(input,
                                            blockComments + "|" + lineComments + "|" + strings + "|" + verbatimStrings,
                                            me => {
                                                if (me.Value.StartsWith("/*") || me.Value.StartsWith("//"))
                                                    return me.Value.StartsWith("//") ? System.Environment.NewLine : "";
                                                // Keep the literal strings
                                                return me.Value;
                                            },
                                            RegexOptions.Singleline);
            return noComments;
        }

        
        /// <summary>
        /// the 'EP' of <see cref="REPL"/>
        /// </summary>
        public static object EvalPrint(string line, out Assembly AssemblyProduced) {
            
            string lineWithoutWhite = line.TrimStart(new char[] { ' ', '\t', '\r', '\n' }).TrimEnd(new char[] { ' ', '\t', '\r', '\n' });
            if (lineWithoutWhite == "restart" || lineWithoutWhite == "restart;") {
                InteractiveShell.Reset();
                InteractiveShell.LastError = null;
                string runcommands = ParseDBErc();
#if !DEBUG
                try { 
#endif
                    eval = new SafeEvaluator(() => Startup(runcommands));
#if !DEBUG
                } catch (Exception e) {
                    Console.WriteLine(
                        "DBE initialization failed with message '{0}'. Type 'LastError' for details.",
                        e.Message);
                    InteractiveShell.LastError = e;
                    eval = new SafeEvaluator(() => Startup(runcommands));
                }
#endif
                AssemblyProduced = null;
                return null;
            }

            if (eval == null) {
                Console.Error.WriteLine("C# evaluator not initialized: use 'restart'.");
                AssemblyProduced = null;
                return null;
            }

            object result = null;
            bool result_set;
#if !DEBUG
            try {
#endif

            string ans = eval.Evaluator.Evaluate(line, out result, out result_set);
            AssemblyProduced = eval.LatestAssembly;
            
            if (result_set) {
                InteractiveShell.ans = result;
            }
            if (ans != null) {
                string lineWoComment = RemoveAllComments(line);
                if (lineWoComment.IsEmptyOrWhite()) {

                    Console.Error.WriteLine("Incomplete statement - missing closing parentheses?.");
                }
                return null;
            }
            if (cmpCont.Report.Errors > 0 || cmpCont.Report.Warnings > 0) {
                Console.WriteLine("No. of errors|warnings: {0}|{1}.", cmpCont.Report.Errors, cmpCont.Report.Warnings);
            }
            
            if (result_set && result != null) {
                // t.IsAssignableFrom(typeof(char)) => Beware of strings
                // (a.k.a. IEnumerable<char>)
                Type typeArg = result.GetType().GetInterfaces().
                    Where(i => i.GetBaseName() == "IEnumerable").
                    Where(i => i.GetGenericArguments().Length == 1).
                    Select(i => i.GetGenericArguments().First()).
                    SingleOrDefault(t => !t.IsAssignableFrom(typeof(Char)));

                // If type implements IEnumerable<SomeClass>, call
                // IEnumerableExtensions.Summary<SomeClass> (using
                // reflection) instead of ToString in order to get some
                // useful output for the user
                string resultString;
                if (typeArg != null) {
                    var summaryMethod = typeof(BoSSS.Foundation.IO.IEnumerableExtensions).GetMethods().
                        Where(m => m.Name == "Summary").
                        Single(m => m.GetParameters().Length == 1).
                        MakeGenericMethod(typeArg);
                    resultString = summaryMethod.Invoke(result, new object[] { result }).ToString();
                } else {
                    resultString = result.ToString();
                }

                Console.Write(resultString);
                Log.Add(new Tuple<string, string>(line, resultString));
            }
            InteractiveShell.LastError = null;
#if !DEBUG
            } catch (Exception e) {
                Console.WriteLine(String.Format(
                    "{0} occurred: {1}. Type 'LastError' for details.",
                    e.GetType(),
                    e.Message));
                InteractiveShell.LastError = e;
                AssemblyProduced = null;
            }
#endif
            return result;
        }

        private static CommandLineReader GetCommandLineReader() {
            string historyPath = Utils.GetBoSSSUserSettingsPath();
            if (!Directory.Exists(historyPath)) {
                historyPath = null;
            }

            CommandLineReader reader = new CommandLineReader(
                "BoSSSpad", historyPath, 500);
            int timeout = 1000;

            // Get completions from the evaluator knows about
            reader.AutoCompleteEvent += (sender, eventArgs) => {
                string[] completions;
                string prefix;
                bool completed = eval.TryGetCompletions(
                    eventArgs.Text, out completions, out prefix, timeout);

                if (completed && completions != null) {
                    eventArgs.Completions.Add(
                        new CommandLineReader.CompletionList(prefix, completions));
                }
            };

            // Rudimentary handler for completions of InteractiveShell (better
            // than nothing)
            reader.AutoCompleteEvent += (sender, eventArgs) => {
                string text = eventArgs.Text.Trim();
                if (text.IndexOf('.') >= 0) {
                    return;
                }

                string[] completions;
                string prefix;
                bool completed = eval.TryGetCompletions(
                    "InteractiveShell." + text, out completions, out prefix, timeout);

                if (completed && completions != null) {
                    eventArgs.Completions.Add(
                        new CommandLineReader.CompletionList(prefix, completions));
                }
            };

            // Same for InteractiveBase
            reader.AutoCompleteEvent += (sender, eventArgs) => {
                string text = eventArgs.Text.Trim();
                if (text.IndexOf('.') >= 0) {
                    return;
                }

                string[] completions;
                string prefix;
                bool completed = eval.TryGetCompletions(
                    "InteractiveBase." + text, out completions, out prefix, timeout);

                if (completed && completions != null) {
                    eventArgs.Completions.Add(
                        new CommandLineReader.CompletionList(prefix, completions));
                }
            };

            // Standard exit (e.g., user typed quit on the console)
            AppDomain.CurrentDomain.ProcessExit += delegate (object sender, EventArgs e) {
                reader.SaveHistory();
            };

            return reader;
        }

        /// <summary>
        /// Parses the DBErc.cs file which contains customizable startup commands.
        /// </summary>
        /// <returns>The content of the file.</returns>
        private static string ParseDBErc() {
            string dbercPath = Path.Combine(
                Utils.GetBoSSSUserSettingsPath(), "etc", "DBErc.cs");
            string runCommands = "";

            if (File.Exists(dbercPath)) {
                using (StreamReader reader = new StreamReader(dbercPath)) {
                    string content = reader.ReadToEnd();
                    if (!string.IsNullOrWhiteSpace(content)) {
                        runCommands = content;
                    }
                }
            }

            return runCommands;
        }

        /// <summary>
        /// Starts a console prompt for password authentication
        /// </summary>
        /// <param name="context">
        /// The context of the password request (for user information)
        /// </param>
        /// <returns>
        /// The password supplied to the console
        /// </returns>
        private static string PromptPassword(string context) {
            Console.WriteLine();
            Console.WriteLine("Password authentication required for '{0}'", context);

            ConsoleKeyInfo info = Console.ReadKey(true);
            string password = "";
            while (info.Key != ConsoleKey.Enter) {
                if (info.Key == ConsoleKey.Backspace) {
                    if (!string.IsNullOrEmpty(password)) {
                        password = password.Substring(0, password.Length - 1);
                    }
                    Console.Write("\b \b");
                    info = Console.ReadKey(true);
                } else {
                    password += info.KeyChar;
                    Console.Write("*");
                    info = Console.ReadKey(true);
                }
            }

            for (int i = 0; i < password.Length; i++) {
                Console.Write("*");
            }

            return password;
        }

        /// <summary>
        /// A log that keeps track of all entered commands and the
        /// corresponding results. Only used by
        /// <see cref="SaveSessionAsWorksheet"/>.
        /// </summary>
        private static IList<Tuple<string, string>> Log = new List<Tuple<string, string>>();

        /// <summary>
        /// Saves
        /// </summary>
        /// <param name="path"></param>
        public static void SaveSessionAsWorksheet(string path) {
            string extension = ".bws";
            if (!path.EndsWith(extension)) {
                path = path + extension;
            }

            string ResultStartMarker = "**************";
            string ResultEndMarker = "==============";

            using (var fs = File.CreateText(path)) {
                foreach (var t in Log) {
                    string Command = t.Item1;
                    string InterpreterTextOutput = t.Item1;
                    if (Command != null && Command.Length > 0)
                        fs.WriteLine(Command);
                    fs.WriteLine(ResultStartMarker);
                    if (InterpreterTextOutput != null && InterpreterTextOutput.Length > 0)
                        fs.WriteLine(InterpreterTextOutput);
                    fs.WriteLine(ResultEndMarker);
                }
            }

            

        }
    }
}
