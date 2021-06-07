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

using BoSSS.Foundation.IO;
using BoSSS.Platform;
using Microsoft.DotNet.Interactive.Notebook;
using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Threading;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// entry point of the `BoSSSpad` application.
    /// </summary>
    public static class BoSSSpadMain {

        /// <summary>
        /// Modes of operation of BoSSSpad
        /// </summary>
        private enum Modes {

            ///// <summary>
            ///// Classic worksheet mode (with a GUI)
            ///// </summary>
            //Worksheet,

            ///// <summary>
            ///// Interactive console mode (without a GUI)
            ///// </summary>
            //Console,

            ///// <summary>
            ///// Simplified interactive console mode (for embedding into other terminals, experimental)
            ///// </summary>
            //SimpleConsole,

            /// <summary>
            /// Batch execution of .bws files
            /// </summary>
            Batch,

            /// <summary>
            /// Just check the installation
            /// </summary>
            Check,

            /// <summary>
            /// Batch execution of LaTeX files (experimental)
            /// </summary>
            TexBatch,

            /// <summary>
            /// Upgrading some old file format (`bws` or `tex` with BoSSS-macros)
            /// to Jupyter notebooks (`ipnb`)
            /// </summary>
            OldFileUpgrade,

            /// <summary>
            /// Executing a notebook in batch mode
            /// </summary>
            JupyterBatch
        }

   


        /// <summary>
        /// application entry point
        /// </summary>
        //[STAThread]
        public static int Main(string[] args) {
            
            /*
            SshClient_exp ssh = new SshClient_exp("lcluster3.hrz.tu-darmstadt.de", "fk69umer", new PrivateKeyFile("C:\\Users\\flori\\.ssh\\id_rsa"));

            var rr1 = ssh.RunCommand("ls");
            Console.WriteLine(rr1.stdout);
            var rr2 = ssh.RunCommand("ls -l");
            Console.WriteLine(rr2.stdout);
            ssh.Final();
            return 0;
            //*/
            int errCount = 0;

            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            // interpretation of command line options
            // ======================================
            Modes mode;
            {
                bool parseModeSuccesfully = true;
                //if (args.Length == 0) {
                //    // assuming the user wants to run the worksheet mode
                //    mode = Modes.Worksheet;
                //    fileToOpen = null;
                //} else if (args.Length == 1) {
                //    if (args[0].StartsWith("--")) {
                //        parseModeSuccesfully = Enum<Modes>.TryParse(args[0].Substring(2), out mode);
                //        fileToOpen = null;
                //    }
                //    else {
                //        mode = Modes.Worksheet;
                //        fileToOpen = args[0];
                //    }
                //} else 
                if(args.Length <= 0) {
                    PrintUsage();
                    return int.MinValue;
                } else {
                    if(!args[0].StartsWith("--")) {
                        PrintUsage();
                        return int.MinValue;
                    }

                    parseModeSuccesfully = Enum<Modes>.TryParse(args[0].Substring(2), out mode);

                    if(!parseModeSuccesfully) {
                        PrintUsage();
                        return int.MinValue;
                    }
                }

                

               
            }

            // launch the app
            // ==============
            bool IinitializedMPI = BoSSS.Solution.Application.InitMPI();

            switch (mode) {
                //case Modes.Worksheet:
                //throw new NotSupportedException("GUI has been removed; use Jupyter notebook!");
                //var ws = new Worksheet(fileToOpen);
                //ws.Shown += Worksheet.OnShown; // Workaround for wrong word-wrap on start-up of the application
                //System.Windows.Forms.Application.Run(ws);

                //ws.m_ExecutorOfCommandQueue_RegularTermination = false;
                //Thread.Sleep(800);

                //if (ws.m_ExecutorOfCommandQueue.IsAlive) {
                //    // hardcore
                //    Thread.Sleep(5000);
                //    if (ws.m_ExecutorOfCommandQueue.IsAlive) {
                //        ws.m_ExecutorOfCommandQueue.Abort();
                //    }
                //}
                //break;

                //case Modes.Console:
                //ReadEvalPrintLoop.REPL();
                //break;

                //case Modes.SimpleConsole:
                //ReadEvalPrintLoop.REPL_Simple();
                //break;

                case Modes.Check:
                if(args.Length != 1) {
                    PrintUsage();
                    return int.MinValue;
                }
                InstallationChecker.CheckSetup();
                break;

                case Modes.OldFileUpgrade: {
                    string fileToOpen;
                    if(args.Length != 2) {
                        PrintUsage();
                        return int.MinValue;
                    }
                    fileToOpen = args[1];
                    OldFileToJupyter(fileToOpen);
                    break;
                }
                


                case Modes.JupyterBatch: {
                    string fileToOpen;
                    if(args.Length != 2) {
                        PrintUsage();
                        return int.MinValue;
                    }
                    fileToOpen = args[1];

                    RunJupyter(fileToOpen);
                    break;
                }

                case Modes.Batch:
                case Modes.TexBatch: {
                    string fileToOpen;
                    if(args.Length != 2) {
                        PrintUsage();
                        return int.MinValue;
                    }
                    fileToOpen = args[1];

                    Document doc;
                    if(fileToOpen.ToLowerInvariant().EndsWith(".tex")) {
                        LatexIO.SplitTexFile(fileToOpen, out _, out doc);
                    } else {
                        doc = Document.Deserialize(fileToOpen);
                    }
                    string OutDir = Path.GetDirectoryName(fileToOpen);
                    string DocNam = Path.GetFileNameWithoutExtension(fileToOpen) + ".texbatch";
                    InteractiveShell.CurrentDoc = doc;
                    InteractiveShell._CurrentDocFile = (new FileInfo(fileToOpen)).FullName;

                    // Which text boxes should be removed before 'restart' occurs
                    int f = 0;
                    if(mode == Modes.TexBatch) {

                        // bws was produced by Latex - some string replacements are necessary
                        for(int iEntry = 0; iEntry < doc.CommandAndResult.Count; iEntry++) {
                            var Entry = doc.CommandAndResult[iEntry];

                            // Check whether there are boxes before restart
                            if(Entry.Command.Equals("restart") || Entry.Command.Equals("restart;")) {
                                f = iEntry;
                            }

                            Entry.Command = LatexIO.Tex2Bws(Entry.Command);
                        }

                        BoSSSpadGnuplotExtensions.UseCairoLatex = true;
                    }

                    // All boxes before 'restart' should not be counted as error
                    int count = 0;
                    foreach(Document.Tuple dt in doc.CommandAndResult) {
                        Console.WriteLine(dt.Command);
                        bool success = dt.Evaluate();

                        if(!success && count >= f)
                            errCount++;

                        Console.WriteLine(Document.ResultStartMarker);
                        Console.WriteLine(dt.InterpreterTextOutput);
                        Console.WriteLine(Document.ResultEndMarker);

                        count++;
                    }

                    if(mode == Modes.TexBatch) {
                        LatexIO.Save_Texbatch(OutDir, DocNam, doc);
                    } else {
                        if(fileToOpen.EndsWith(".tex")) {

                        } else {
                            doc.Serialize(fileToOpen);
                        }
                    }
                    InteractiveShell.CurrentDoc = null;
                    break;
                }


                default:
                throw new NotImplementedException();
            }

            if (IinitializedMPI)
                BoSSS.Solution.Application.FinalizeMPI();

            return errCount;
        }

        private static void RunJupyter(string fileToOpen) {
            ProcessStartInfo psi = new ProcessStartInfo();
            psi.FileName = @"C:\Windows\System32\cmd.exe";

            psi.WorkingDirectory = Directory.GetCurrentDirectory();

            psi.RedirectStandardInput = true;
            //psi.RedirectStandardOutput = true;
            //psi.RedirectStandardError = true;

            var p = Process.Start(psi);

            //p.StandardInput.WriteLine("dir");
            p.StandardInput.WriteLine(@"C:\ProgramData\Anaconda3\Scripts\activate.bat");
            p.StandardInput.WriteLine("jupyter.exe nbconvert \"" + fileToOpen + "\" --to html --execute");
            p.StandardInput.WriteLine("exit");
            p.WaitForExit();

            Console.WriteLine("--------------------------------");
            Console.WriteLine("Done with fucking notebook");
            Console.WriteLine("Exshit code " + p.ExitCode);
        }

        static string GetStartupCode() {
            using(var stw = new StringWriter()) {
                stw.WriteLine("#r \"BoSSSpad.dll\"");
                stw.WriteLine("using System;");
                stw.WriteLine("using System.Collections.Generic;");
                stw.WriteLine("using System.Linq;");
                stw.WriteLine("using ilPSP;");
                stw.WriteLine("using ilPSP.Utils;");
                stw.WriteLine("using BoSSS.Platform;");
                stw.WriteLine("using BoSSS.Foundation;");
                stw.WriteLine("using BoSSS.Foundation.Grid;");
                stw.WriteLine("using BoSSS.Foundation.Grid.Classic;");
                stw.WriteLine("using BoSSS.Foundation.IO;");
                stw.WriteLine("using BoSSS.Solution;");
                stw.WriteLine("using BoSSS.Solution.Control;");
                stw.WriteLine("using BoSSS.Solution.GridImport;");
                stw.WriteLine("using BoSSS.Solution.Statistic;");
                stw.WriteLine("using BoSSS.Solution.Utils;");
                stw.WriteLine("using BoSSS.Solution.Gnuplot;");
                stw.WriteLine("using BoSSS.Application.BoSSSpad;");
                stw.WriteLine("using BoSSS.Application.XNSE_Solver;");
                stw.WriteLine("using static BoSSS.Application.BoSSSpad.BoSSSshell;");
                stw.WriteLine("Init();");

                return stw.ToString();
            }
        }


        private static string OldFileToJupyter(string fileToOpen) {
            Document doc;
            if(fileToOpen.ToLowerInvariant().EndsWith(".tex")) {
                LatexIO.SplitTexFile(fileToOpen, out _, out doc);
            } else {
                doc = Document.Deserialize(fileToOpen);
            }


            string DestFile = Path.GetFileNameWithoutExtension(fileToOpen) + ".ipynb";

            var cells = new List<NotebookCell>();
            foreach(var entry in doc.CommandAndResult) {
                string cmd = entry.Command;
                if(cmd.StartsWith("restart;"))
                    cmd = cmd.Replace("restart;", GetStartupCode());
                if(cmd.StartsWith("restart"))
                    cmd = cmd.Replace("restart", GetStartupCode());

                //cells.Add(new NotebookCell("C#", cmd));
                cells.Add(new NotebookCell("csharp", cmd));
            }


            var docNew = new NotebookDocument(cells.ToArray());

            var data = NotebookFileFormatHandler.Serialize(DestFile, docNew, System.Environment.NewLine);
            System.IO.File.WriteAllBytes(DestFile, data);
            return DestFile;
        }

        /// <summary>
        /// Write the usage statement to the console
        /// </summary>
        public static void PrintUsage()
        {
            //                 0        1         2         3         4         5         6         7         8
            //                 12345678901234567890123456789012345678901234567890123456789012345678901234567890
            Console.WriteLine();
            Console.WriteLine("BoSSSpad application. ");
            Console.WriteLine();
            Console.WriteLine("Usage is:");
            Console.WriteLine();
            //Console.WriteLine("Option 1: Worksheet mode (interactive/GUI):");
            //Console.WriteLine("-------------------------------------------");
            //Console.WriteLine("    BoSSSpad.exe                       Opens an empty worksheet");
            //Console.WriteLine(" or BoSSSpad.exe file.bws              Opens file.bws");
            //Console.WriteLine(" or BoSSSpad.exe --worksheet file.bws  Opens file.bws");
            //Console.WriteLine();
            //Console.WriteLine("Option 2: Console mode (interactive/text):");
            //Console.WriteLine("------------------------------------------");
            //Console.WriteLine("    BoSSSpad.exe --console             Starts a console session.");
            //Console.WriteLine("    BoSSSpad.exe --simpleconsole       Starts a simple (e.g. no autocompletion etc.) console session.");
            //Console.WriteLine();
            Console.WriteLine("Option 1: Batch mode (non-interactive/text):");
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --batch file.bws      Upgrade file.bws to Jupyter & batch exec");
            Console.WriteLine("    BoSSSpad.exe --texbatch file.tex   Upgrade file.tex to Jupyter & batch exec");
            Console.WriteLine("    BoSSSpad.exe --JupyterBatch file.ipynb     Batch execution of in file.ipynb");
            Console.WriteLine();
            Console.WriteLine("Option 2: Old file upgrade (non-interactive/text):");
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --OldFileUpgrade file.bws        Upgrade file.bws to Jupyter");
            Console.WriteLine();
            Console.WriteLine("Option 3: Installation check:");
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --check               Check the BoSSS installation.");
            Console.WriteLine();
            
        }

        /// <summary>
        /// Dummy function which ensures that certain referenced assemblies 
        /// </summary>
        private static void LinkEnforcer() {
            // If you remove these lines, this may break some worksheets and tutorials.
            Console.WriteLine(typeof(CNS.Program).FullName);
            Console.WriteLine(typeof(IBM_Solver.IBM_SolverMain).FullName);
            Console.WriteLine(typeof(XNSE_Solver.XNSE).FullName);
            Console.WriteLine(typeof(XNSERO_Solver.XNSERO).FullName);
        }
    }
}
