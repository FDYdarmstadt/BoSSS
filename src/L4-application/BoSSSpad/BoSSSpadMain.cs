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
using System.Linq;
using System.Reflection;
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
            JupyterBatch,


            /// <summary>
            /// Create a Jupyter notebook with BoSSS-init code
            /// </summary>
            Jupyterfile,

            /// <summary>
            /// application deployment and running on a cluster
            /// </summary>
            RunBatch

        }

        /*
        class KnownTypesBinder : System.Runtime.Serialization.SerializationBinder {

            Job m_owner;

            internal KnownTypesBinder(Job __owner) {
                m_owner = __owner;
            }

            Assembly a = typeof(BoSSS.Solution.Statistic.CellLocalization).Assembly;

            
            public override Type BindToType(string assemblyName, string typeName) {
                Console.WriteLine(a.FullName);
                var tts = a.GetExportedTypes();
                var tt = tts.First(t => t.FullName == typeName);
                return tt;
                //throw new NotImplementedException();

            }
        }
   */


        /// <summary>
        /// application entry point
        /// </summary>
        //[STAThread]
        public static int Main(string[] args) {

            /*
            string path = @"c:\Users\flori\AppData\Local\BoSSS-LocalJobs\Demo_BoundaryAndInitialData-ipPoisson2021Juni10_083737\control.obj";
            string text = File.ReadAllText(path);
            var obj = Solution.Control.AppControl.Deserialize(text, new KnownTypesBinder(null));
            Console.WriteLine("desez: " + obj.GetType());
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

            try
            {
                switch (mode)
                {

                    case Modes.Check:
                        if (args.Length != 1)
                        {
                            PrintUsage();
                            return int.MinValue;
                        }
                        InstallationChecker.CheckSetup();
                        break;

                    case Modes.OldFileUpgrade:
                        {
                            string fileToOpen;
                            if (args.Length != 2)
                            {
                                PrintUsage();
                                return int.MinValue;
                            }
                            fileToOpen = args[1];
                            OldFileToJupyter(fileToOpen);
                            break;
                        }



                    case Modes.JupyterBatch:
                        {
                            string fileToOpen;
                            if (args.Length != 2)
                            {
                                PrintUsage();
                                return int.MinValue;
                            }
                            fileToOpen = args[1];

                            errCount = RunJupyter(fileToOpen);
                            break;
                        }

                    case Modes.Batch:
                    case Modes.TexBatch:
                        {
                            string fileToOpen;
                            if (args.Length != 2)
                            {
                                PrintUsage();
                                return int.MinValue;
                            }
                            fileToOpen = args[1];
                            string ConvFile = OldFileToJupyter(fileToOpen);
                            errCount = RunJupyter(ConvFile);


                            break;
                        }

                    case Modes.Jupyterfile:
                        {
                            string fileToOpen;
                            if (args.Length != 2)
                            {
                                PrintUsage();
                                return int.MinValue;
                            }
                            fileToOpen = args[1];
                            Jupyterfile(fileToOpen);
                            break;
                        }

                    case Modes.RunBatch:
                        errCount = SubprogramRunbatch.RunBatch(args.Skip(1).ToArray());
                        break;

                    default:
                        throw new NotImplementedException();
                }
            } catch(Exception e) {
                Console.WriteLine(e.GetType().Name + ": " + e.Message);
                errCount = -666;
                //throw new AggregateException(e);
            }


            if (IinitializedMPI)
                BoSSS.Solution.Application.FinalizeMPI();

            return errCount;
        }


        private static Mutex JupyterMutex = new Mutex(false, "JupyterMutex");

        private static int RunJupyter(string fileToOpen) {
            ProcessStartInfo psi = new ProcessStartInfo();
            psi.FileName = @"C:\Windows\System32\cmd.exe";

            psi.WorkingDirectory = Directory.GetCurrentDirectory();

            psi.RedirectStandardInput = true;
            //psi.RedirectStandardOutput = true;
            //psi.RedirectStandardError = true;

            try {
                Console.WriteLine("Waiting for Jupyter mutex (can only use one Jupyter notebook at time) ...");
                JupyterMutex.WaitOne();
                Console.WriteLine("Mutex obtained!");

                var p = Process.Start(psi);

                //p.StandardInput.WriteLine("dir");
                p.StandardInput.WriteLine(@"C:\ProgramData\Anaconda3\Scripts\activate.bat");
                p.StandardInput.WriteLine("jupyter.exe nbconvert \"" + fileToOpen + "\" --to html --execute");
                p.StandardInput.WriteLine("exit");
                p.WaitForExit();

                Console.WriteLine("--------------------------------");
                Console.WriteLine("Done with notebook");
                Console.WriteLine("Exit code " + p.ExitCode);
                return p.ExitCode;
            } finally {
                JupyterMutex.ReleaseMutex();
            }
        }

        static string GetStartupCode() {
            using(var stw = new StringWriter()) {
                string path = typeof(BoSSSpadMain).Assembly.Location;
                stw.WriteLine("#r \"" + path + "\"");
                stw.WriteLine("using System;");
                stw.WriteLine("using System.Collections.Generic;");
                stw.WriteLine("using System.Linq;");
                stw.WriteLine("using ilPSP;");
                stw.WriteLine("using ilPSP.Utils;");
                stw.WriteLine("using BoSSS.Platform;");
                stw.WriteLine("using BoSSS.Platform.LinAlg;");
                stw.WriteLine("using BoSSS.Foundation;");
                stw.WriteLine("using BoSSS.Foundation.XDG;");
                stw.WriteLine("using BoSSS.Foundation.Grid;");
                stw.WriteLine("using BoSSS.Foundation.Grid.Classic;");
                stw.WriteLine("using BoSSS.Foundation.Grid.RefElements;");
                stw.WriteLine("using BoSSS.Foundation.IO;");
                stw.WriteLine("using BoSSS.Solution;");
                stw.WriteLine("using BoSSS.Solution.Control;");
                stw.WriteLine("using BoSSS.Solution.GridImport;");
                stw.WriteLine("using BoSSS.Solution.Statistic;");
                stw.WriteLine("using BoSSS.Solution.Utils;");
                stw.WriteLine("using BoSSS.Solution.AdvancedSolvers;");
                stw.WriteLine("using BoSSS.Solution.Gnuplot;");
                stw.WriteLine("using BoSSS.Application.BoSSSpad;");
                stw.WriteLine("using BoSSS.Application.XNSE_Solver;");
                stw.WriteLine("using BoSSS.Application.XNSFE_Solver;");
                stw.WriteLine("using static BoSSS.Application.BoSSSpad.BoSSSshell;");
                stw.WriteLine("Init();");

                return stw.ToString();
            }
        }


        private static void Jupyterfile(string fileToCreate) {

            var cells = new List<NotebookCell>();
            string cmd = GetStartupCode();
            cells.Add(new NotebookCell("csharp", cmd));

            var docNew = new NotebookDocument(cells.ToArray());

            var data = NotebookFileFormatHandler.Serialize(fileToCreate, docNew, System.Environment.NewLine);
            System.IO.File.WriteAllBytes(fileToCreate, data);
        }

        private static string OldFileToJupyter(string fileToOpen) {
            Document doc;
            if(fileToOpen.ToLowerInvariant().EndsWith(".tex")) {
                LatexIO.SplitTexFile(fileToOpen, out _, out doc);
            } else {
                doc = Document.Deserialize(fileToOpen);
            }


            string DestFile = Path.Combine(Path.GetDirectoryName(fileToOpen), Path.GetFileNameWithoutExtension(fileToOpen)) + ".ipynb";

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
            Console.WriteLine("Option 3: Jupyter notebook for BoSSS:");
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --Jupyterfile file.ipynb   Create file prepared for BoSSS");
            Console.WriteLine();
            Console.WriteLine("Option 4: Submit HPC job:");
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --RunBatch --help  deploy and run some solver on a ");
            Console.WriteLine("                                    cluster; Use --help for further info.");
            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine("Option 5: Installation check:");
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --check               Check the BoSSS installation.");
            
        }

        /// <summary>
        /// Dummy function which ensures that certain referenced assemblies 
        /// </summary>
        public static void LinkEnforcer() {
            // If you remove these lines, this may break some worksheets and tutorials.
            Console.WriteLine(typeof(CNS.Program).FullName);
            Console.WriteLine(typeof(IBM_Solver.IBM_SolverMain).FullName);
            Console.WriteLine(typeof(XNSE_Solver.XNSE).FullName);
            Console.WriteLine(typeof(XNSFE_Solver.XNSFE).FullName);
            Console.WriteLine(typeof(XNSERO_Solver.XNSERO).FullName);
            Console.WriteLine(typeof(ZwoLevelSetSolver.ZLS).FullName);
        }


        

    }
}
