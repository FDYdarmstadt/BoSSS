﻿/* =======================================================================
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

using BoSSS.Platform;
using Microsoft.DotNet.Interactive.Documents;
using Microsoft.DotNet.Interactive.Documents.Jupyter;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;
using ilPSP;
using System.IO.Pipes;
using System.Threading.Tasks;
using static ilPSP.Connectors.Matlab.BatchmodeConnector;
using System.Diagnostics.Metrics;
using BoSSS.Solution;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// entry point of the `BoSSSpad` application.
    /// </summary>
    public static class BoSSSpadMain {

        /// <summary>
        /// Modes of operation of BoSSSpad
        /// </summary>
        private enum Modes {
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
            int errCount = 0;
            

            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
           /*
            string path = @"\\fdygitrunner\ValidationTests\databases\bkup-2023Oct31_165410.LinslvPerfSer\sessions\dce38035-09a6-4f62-a82f-b7c089d27822\"; 
            var json = File.ReadAllText(Path.Combine(path, "profiling_bin.0.txt"));

            var p = OnlineProfiling.Deserialize(json);

            Console.Write(p.ToString());
            
            //*/

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
                if (args.Length <= 0) {
                    PrintUsage();
                    return int.MinValue;
                } else {
                    if (!args[0].StartsWith("--")) {
                        PrintUsage();
                        return int.MinValue;
                    }

                    parseModeSuccesfully = Enum<Modes>.TryParse(args[0].Substring(2), out mode);

                    if (!parseModeSuccesfully) {
                        PrintUsage();
                        return int.MinValue;
                    }
                }




            }

            // launch the app
            // ==============
            bool IinitializedMPI = BoSSS.Solution.Application.InitMPI(args);

            try {
                switch (mode) {

                    case Modes.Check:
                    if (args.Length != 1) {
                        PrintUsage();
                        return int.MinValue;
                    }
                    InstallationChecker.CheckSetup();
                    break;

                    case Modes.OldFileUpgrade: {
                        string fileToOpen;
                        if (args.Length != 2) {
                            PrintUsage();
                            return int.MinValue;
                        }
                        fileToOpen = args[1];
                        OldFileToJupyter(fileToOpen);
                        break;
                    }



                    case Modes.JupyterBatch: {
                        string fileToOpen;
                        if (args.Length != 2) { 
                            Console.Error.WriteLine($"Expecting exactly two arguments, but got {args.Length} (which are {args.ToConcatString("", ",", "")}");
                            PrintUsage();
                            return int.MinValue;
                        }
                        fileToOpen = args[1];

                        errCount = RunJupyter(fileToOpen);
                        break;
                    }

                    case Modes.Batch:
                    case Modes.TexBatch: {
                        string fileToOpen;
                        if (args.Length != 2) {
                            PrintUsage();
                            return int.MinValue;
                        }
                        fileToOpen = args[1];
                        string ConvFile = OldFileToJupyter(fileToOpen);
                        errCount = RunJupyter(ConvFile);


                        break;
                    }

                    case Modes.Jupyterfile: {
                        string fileToOpen;
                        if (args.Length != 2) {
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
            } catch (Exception e) {
                Console.Error.WriteLine(e.GetType().Name + ": " + e.Message);
                errCount = -666;
                //throw new AggregateException(e);
            }


            if (IinitializedMPI)
                BoSSS.Solution.Application.FinalizeMPI();

            return errCount;

            //*/
        }

        /// <summary>
        /// When worksheets are executed in batch mode, this 
        /// should prevent race-condition for a port during startup up dotnet-interactive.
        /// Jupyter worksheets can be batch-processed in parallel (using `jupyter nbconvert` or `papermill`), but
        /// during startup we might run into an exception 
        /// `System.IO.IOException: Failed to bind to address http://192.168.56.1:1004: address already in use.`
        /// 
        /// Therefore, this mutex is used, in combination with a named pipe named <see cref="BoSSSpadInitDone_PipeName"/> to prevent BoSSSpad from 
        /// starting another Jupyter instance before the port is acquired.
        /// </summary>
        static Mutex JupyterMutex;

        static BoSSSpadMain() {
            

            try {
                JupyterMutex = new Mutex(false, "JupyterMutex");
            } catch (Exception ee) {
                Console.Error.WriteLine("BoSSSpadMain, Exception in static ctor during obtaining Jupyter Mutex:" + ee);
                Console.Error.WriteLine("Terminating application.");
                System.Environment.Exit(-987);
           
            }

           
        }

        /// <summary>
        /// Name used for a synchronization pipe
        /// </summary>
        internal static string BoSSSpadInitDone_PipeName => "BoSSSpadInitDone" + System.Environment.UserName;


        private static int RunJupyter(string fileToOpen) {
            return RunPapermillAndNbconvert(fileToOpen);
            
            
            //int RetVal = RunJupyter(fileToOpen, false); // first try, don't allow errors
            /*if(RetVal == 0)
                return 0;

            // if we have an error, we have no output; then, the error might be hard to fix;
            // therefore, we have to re-run with `--allow-errors`;
            // Rem.: if --allow-errors is used, the return value of `jupyter nbconvert` is always 0; therefore we would not detect any error if we use it beforehand.
            //       as an alternative, one might convert to a notebook (`--to notebook --execute --allow-errors --output out.ipynb` 
            //       and then parse the `out.ipynb` for exceptions and errors.
            RunJupyter(fileToOpen, true);*/
            //return RetVal;
        }

        /*
        private static int RunJupyter(string fileToOpen, bool AllowErrors) {
            ProcessStartInfo psi = new ProcessStartInfo();
            string SuffixAllowErrors = AllowErrors ? "--allow-errors" : "";

            psi.WorkingDirectory = Directory.GetCurrentDirectory();

            psi.RedirectStandardInput = true;
            //psi.RedirectStandardOutput = true;
            //psi.RedirectStandardError = true;

            bool MutexReleased = false;
            try {
                Console.WriteLine("Waiting for Jupyter mutex (can only use one Jupyter notebook at time) ...");
                JupyterMutex.WaitOne();
                Console.WriteLine("Mutex obtained!");
                
                Process p;
                if(System.OperatingSystem.IsWindows()) {
                    psi.FileName = @"C:\Windows\System32\cmd.exe";

                    // wait here a bit to avoid a port conflict...
                    // when two notebooks are started simultaneously, we might run into the following:
                    //  ---> System.IO.IOException: Failed to bind to address http://192.168.56.1:1004: address already in use.
                    p = Process.Start(psi);
                    Random rnd = new Random();
                    Thread.Sleep(rnd.Next(1000, 5000) + Math.Abs(fileToOpen.GetHashCode()%2217));

                    //p.StandardInput.WriteLine("dir");
                    p.StandardInput.WriteLine(@"C:\ProgramData\Anaconda3\Scripts\activate.bat");
                    p.StandardInput.WriteLine("jupyter.exe nbconvert \"" + fileToOpen + "\" --to html --execute " + SuffixAllowErrors);
                    p.StandardInput.WriteLine("exit");

                    // wait here a bit more...
                    Thread.Sleep(rnd.Next(1000, 5000) + Math.Abs(fileToOpen.GetHashCode()%2217));

                    JupyterMutex.ReleaseMutex();
                    MutexReleased = true;
                    p.WaitForExit();

                } else {
                    psi.FileName = @"jupyter";
                    psi.Arguments = " nbconvert " + fileToOpen + " --to html --execute " + SuffixAllowErrors; 

                    p = Process.Start(psi);
                    p.WaitForExit();

                }
                Console.WriteLine("--------------------------------");
                Console.WriteLine("Done with notebook");
                Console.WriteLine("Exit code " + p.ExitCode);
                return p.ExitCode;
            } finally {
                if(!MutexReleased)
                    JupyterMutex.ReleaseMutex();
            }
        }
        */


        static bool UseAnacondaPython() {
            if(System.Environment.MachineName.Contains("hpccluster", StringComparison.InvariantCultureIgnoreCase))
                return false;
            if (System.Environment.MachineName.Contains("SHUBNIGGURATH", StringComparison.InvariantCultureIgnoreCase))
                return false;

            return System.OperatingSystem.IsWindows();
        }

        const bool UseMutexOnPapermill = true;
        const bool UseMutexOnNbconvert = true;

       
        private static int RunPapermillAndNbconvert(string fileToOpen) {
       
            string fileToOpen_out = Path.Combine(Path.GetDirectoryName(fileToOpen), Path.GetFileNameWithoutExtension(fileToOpen) + "_out.ipynb");

            string htmlResult = Path.Combine(Path.GetDirectoryName(fileToOpen), Path.GetFileNameWithoutExtension(fileToOpen) + ".html");
            string htmlResult_out = Path.Combine(Path.GetDirectoryName(fileToOpen), Path.GetFileNameWithoutExtension(fileToOpen) + "_out.html");


            bool MutexReleased = true;
            Random rnd = new Random();
            //Thread.Sleep(rnd.Next(1000, 5000) + Math.Abs(fileToOpen.GetHashCode() % 2217));
            try {



                void GetMutex(int RecDepth) {
                    if (MutexReleased) {
                        Console.WriteLine("Waiting for Jupyter mutex, " + DateTime.Now +" (can only start one Jupyter notebook at time) ...");
                        try {
                            JupyterMutex.WaitOne();
                        } catch (AbandonedMutexException ae) {
                            Console.WriteLine($"AbandonedMutexException caught during WaitOne() : ({ae}: {ae.Message})");
                            if (RecDepth > 10) {
                                Console.WriteLine("Stopping trying - unrecoverable");
                                Console.Error.WriteLine("Terminating application.");
                                System.Environment.Exit(-989);
                            } else {
                                Random rnd = new Random();
                                int msWait = rnd.Next(10000);
                                Console.WriteLine($"Retry No. {RecDepth} in {msWait} milliseconds.");
                                Thread.Sleep(msWait);


                                try {
                                    if (JupyterMutex != null) {
                                        JupyterMutex.ReleaseMutex();
                                        JupyterMutex.Dispose();
                                    }
                                    JupyterMutex = new Mutex(false, "JupyterMutex");
                                } catch (Exception eee) {
                                    Console.Error.WriteLine("BoSSSpadMain.RunPapermillAndNbconvert(...): Unrecoverable Exception during creation of JupyterMutex : " + eee);
                                    Console.Error.WriteLine("Terminating application.");
                                    System.Environment.Exit(-988);
                                }
                                GetMutex(RecDepth + 1);
                            }

                        } catch (Exception eee) {
                            Console.Error.WriteLine("BoSSSpadMain.RunPapermillAndNbconvert(...): Exception during WaitOne() :" + eee);
                            Console.Error.WriteLine("Terminating application.");
                            System.Environment.Exit(-988);
                        }
                        Console.WriteLine("Mutex obtained at " + DateTime.Now + ".");
                        MutexReleased = false;
                    }
                }

                void ReleaseMutex() {
                    //Thread.Sleep(rnd.Next(1000, 5000) + Math.Abs(fileToOpen.GetHashCode() % 2217));
                    Console.WriteLine("Releasing Jupyter mutex @ " + DateTime.Now + " ...");
                    if (!MutexReleased) {
                        try {
                            JupyterMutex.ReleaseMutex();
                        } catch (Exception eee) {
                            Console.Error.WriteLine("BoSSSpadMain.RunPapermillAndNbconvert(...): Exception (1) during ReleaseMutex():" + eee);
                            Console.Error.WriteLine("Terminating application.");
                            System.Environment.Exit(-989);
                        }

                    }
                    Console.WriteLine("Mutex released at " + DateTime.Now + ".");
                    MutexReleased = true;
                }


                void WaitForPipeConnection(string tempguid) {
                    try {


                        // Miss-used pipe for inter-process synchronization.
                        // An `EventWaitHandle` would be much nicer, but that works only on Windows-machines.
                        //static internal EventWaitHandle BoSSSpadInitDone = new EventWaitHandle(false, EventResetMode.ManualReset, "MyUniqueEventName");
                        using (NamedPipeClientStream BoSSSpadInitDone = new NamedPipeClientStream(".", tempguid, PipeDirection.InOut)) {
                            Console.WriteLine("Waiting for BoSSSpad to start up at " + DateTime.Now + ". Current directory is " + Directory.GetCurrentDirectory());
                            BoSSSpadInitDone.Connect(3 * 60 * 1000);
                            Console.WriteLine("BoSSSpad connected at " + DateTime.Now + "; now waiting for signal...");

                            using (var cts = new CancellationTokenSource()) {
                                Task t = new Task(delegate () {
                                    Console.WriteLine($" waiting for signal {DateTime.Now} ...");
                                    int str = BoSSSpadInitDone.ReadByte();
                                    while (str != 1) {
                                        Console.WriteLine($" received signal {DateTime.Now}, got {str}");
                                        str = BoSSSpadInitDone.ReadByte();
                                    }
                                    Console.WriteLine($" Finally received: {DateTime.Now}" + str);
                                }, cts.Token);

                                t.Start();

                                try {
                                    //{ 
                                    if (t.Wait(60 * 1000) == false) {
                                        cts.Cancel();
                                        Console.Error.WriteLine("Timeout waiting for Cancellation token pipe. " + DateTime.Now);
                                    }
                                } catch (Exception e) {
                                    Console.Error.WriteLine("Exception " + DateTime.Now + " while waiting for Cancellation token pipe: " + e);

                                }
                            }
                        }
                    } catch (Exception e) {
                        Console.Error.WriteLine("Exception " + DateTime.Now + " while waiting for pipe connection: " + e);
                        throw new AggregateException(e);
                    }
                }

                int papermill_exit, nbconvert_exit;
                if (UseAnacondaPython()) {
                    


                    int RunAnacondaShell(string command, bool useMutex, bool startupMutex) {
                        if (useMutex)
                            GetMutex(0);

                        ProcessStartInfo psi = new ProcessStartInfo();
                        psi.WorkingDirectory = Directory.GetCurrentDirectory();
                        psi.RedirectStandardInput = true;
                        psi.FileName = @"C:\Windows\System32\cmd.exe";

                        var tempguid = Guid.NewGuid().ToString();
                        Console.WriteLine("Temp GUID for named pipe = " + tempguid);
                        psi.EnvironmentVariables.Add(BoSSSpadInitDone_PipeName, tempguid);

                        // wait here a bit to avoid a port conflict...
                        // when two notebooks are started simultaneously, we might run into the following:
                        //  ---> System.IO.IOException: Failed to bind to address http://192.168.56.1:1004: address already in use.
                       
                        Process p = Process.Start(psi);
                        p.StandardInput.WriteLine(@"C:\ProgramData\Anaconda3\Scripts\activate.bat");
                        p.StandardInput.WriteLine(command);// "jupyter.exe nbconvert \"" + fileToOpen + "\" --to html ");
                        p.StandardInput.Flush();

                        if (startupMutex == true) {
                            // wait here a bit more...
                            //{ 
                            WaitForPipeConnection(tempguid);

                            // received a signal from worksheet that it is up and running, 
                            // so we can release the mutex **before** the external process exits
                            if (useMutex)
                                ReleaseMutex(); 
                        }

                        p.StandardInput.WriteLine("exit");
                        p.WaitForExit();

                        if(startupMutex == false) {
                            // the mutex should block the entire call to the external process.
                            if (useMutex)
                                ReleaseMutex();
                        }

                        return p.ExitCode;
                    }

                    papermill_exit = RunAnacondaShell($"papermill {fileToOpen} {fileToOpen_out}", UseMutexOnPapermill, true);
                    nbconvert_exit = RunAnacondaShell("jupyter.exe nbconvert \"" + fileToOpen_out + "\" --to html ", UseMutexOnNbconvert, false);
                    
                    //nbconvert_exit = RunAnacondaShell("jupyter.exe nbconvert \"" + fileToOpen_out + "\" --to html --execute");
                    //papermill_exit = nbconvert_exit;
                } else {


                    int RunExt(string executable, string arguments, bool useMutex, bool startupMutex) {
                        if(startupMutex)
                            GetMutex(0);

                        ProcessStartInfo psi = new ProcessStartInfo();
                        psi.WorkingDirectory = Directory.GetCurrentDirectory();
                        psi.RedirectStandardInput = true;

                        psi.FileName = executable;
                        psi.Arguments = arguments;

                        var tempguid = Guid.NewGuid().ToString();
                        Console.WriteLine("Temp GUID for named pipe = " + tempguid);
                        psi.EnvironmentVariables.Add(BoSSSpadInitDone_PipeName, tempguid);

                        Process p = Process.Start(psi);

                        // wait here a bit more...
                        if (startupMutex == true) {
                            // wait here a bit more...
                            //{ 
                            WaitForPipeConnection(tempguid);

                            // received a signal from worksheet that it is up and running, 
                            // so we can release the mutex **before** the external process exits
                            if (useMutex)
                                ReleaseMutex();
                        }

                        p.WaitForExit();

                        if (startupMutex == false) {
                            // the mutex should block the entire call to the external process.
                            if (useMutex)
                                ReleaseMutex();
                        }

                        return p.ExitCode;
                    }

                    papermill_exit = RunExt($"papermill", $"{fileToOpen} {fileToOpen_out}", UseMutexOnPapermill, true);
                    nbconvert_exit = RunExt("jupyter", "nbconvert \"" + fileToOpen_out + "\" --to html ", UseMutexOnNbconvert, false);
                    //nbconvert_exit = RunExt("jupyter.exe", "nbconvert \"" + fileToOpen_out + "\" --to html --execute");
                    //papermill_exit = nbconvert_exit;

                }
                Console.WriteLine("--------------------------------");
                Console.WriteLine("Done with notebook");
                Console.WriteLine("papermill code " + papermill_exit);
                Console.WriteLine("nbconvert code " + nbconvert_exit);

                if(File.Exists(htmlResult_out)) {
                    // try to restore the html result with the original name
                    File.Move(htmlResult_out, htmlResult, true);
                }


                return papermill_exit;
            } finally {
                if (!MutexReleased) {
                    try {
                        JupyterMutex.ReleaseMutex();
                    } catch (Exception eee) {
                        Console.Error.WriteLine("BoSSSpadMain.RunPapermillAndNbconvert(...): Exception (2) during ReleaseMutex():" + eee);
                        Console.Error.WriteLine("Terminating application.");
                        System.Environment.Exit(-990);
                    }
                }
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

            //var cells = new List<NotebookCell>();
            var cells = new List<InteractiveDocumentElement>();
            string cmd = GetStartupCode();
            cells.Add(new InteractiveDocumentElement("csharp", cmd));

            //var docNew = new NotebookDocument(cells.ToArray());
            var docNew = new InteractiveDocument(cells.ToArray());

            //var data = NotebookFileFormatHandler.Serialize(fileToCreate, docNew, System.Environment.NewLine);
            using(var stw = new StreamWriter(fileToCreate)) {
                Notebook.Write(docNew, System.Environment.NewLine, stw);
                //System.IO.File.WriteAllBytes(fileToCreate, data);
                stw.Flush();
            }
        }

        private static string OldFileToJupyter(string fileToOpen) {
            Document doc;
            if(fileToOpen.ToLowerInvariant().EndsWith(".tex")) {
                LatexIO.SplitTexFile(fileToOpen, out _, out doc);
            } else {
                doc = Document.Deserialize(fileToOpen);
            }


            string DestFile = Path.Combine(Path.GetDirectoryName(fileToOpen), Path.GetFileNameWithoutExtension(fileToOpen)) + ".ipynb";

            var cells = new List<InteractiveDocumentElement>();
            foreach(var entry in doc.CommandAndResult) {
                string cmd = entry.Command;
                if(cmd.StartsWith("restart;"))
                    cmd = cmd.Replace("restart;", GetStartupCode());
                if(cmd.StartsWith("restart"))
                    cmd = cmd.Replace("restart", GetStartupCode());

                //cells.Add(new NotebookCell("C#", cmd));
                cells.Add(new InteractiveDocumentElement("csharp", cmd));
            }


            var docNew = new InteractiveDocument(cells.ToArray());

            //var data = NotebookFileFormatHandler.Serialize(DestFile, docNew, System.Environment.NewLine);
            //System.IO.File.WriteAllBytes(DestFile, data);
            using(var stw = new StreamWriter(DestFile)) {
                Notebook.Write(docNew, System.Environment.NewLine, stw);
                //System.IO.File.WriteAllBytes(fileToCreate, data);
                stw.Flush();
            }
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
            Console.WriteLine(typeof(CNS.CNSProgram).FullName);
            Console.WriteLine(typeof(XNSE_Solver.XNSE).FullName);
            Console.WriteLine(typeof(XNSFE_Solver.XNSFE).FullName);
            Console.WriteLine(typeof(XNSERO_Solver.XNSERO).FullName);
            Console.WriteLine(typeof(ZwoLevelSetSolver.ZLS).FullName);
            Console.WriteLine(typeof(XNSEC.XNSEC).FullName);
        }


        

    }
}
