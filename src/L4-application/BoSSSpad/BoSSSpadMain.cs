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
using MPI.Wrappers;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
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

            /// <summary>
            /// Classic worksheet mode (with a GUI)
            /// </summary>
            Worksheet,

            /// <summary>
            /// Interactive console mode (without a GUI)
            /// </summary>
            Console,

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
            TexBatch
        }

        /// <summary>
        /// application entry point
        /// </summary>
        [STAThread]
        public static int Main(string[] args) { 
            
            int errCount = 0;

            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            // interpretation of command line options
            // ======================================
            string fileToOpen;
            Modes mode;
            {
                bool parseModeSuccesfully = true;
                if (args.Length == 0) {
                    // assuming the user wants to run the worksheet mode
                    mode = Modes.Worksheet;
                    fileToOpen = null;
                }
                else if (args.Length == 1) {
                    if (args[0].StartsWith("--")) {
                        parseModeSuccesfully = Enum<Modes>.TryParse(args[0].Substring(2), out mode);
                        fileToOpen = null;
                    }
                    else {
                        mode = Modes.Worksheet;
                        fileToOpen = args[0];
                    }
                }
                else if (args.Length == 2) {
                    parseModeSuccesfully = Enum<Modes>.TryParse(args[0].Substring(2), out mode);
                    fileToOpen = args[1];
                }
                else {
                    PrintUsage();
                    return int.MinValue;
                }

                if (!parseModeSuccesfully) {
                    PrintUsage();
                    return int.MinValue;
                }

                if (mode == Modes.Console && fileToOpen != null) {
                    PrintUsage();
                    return int.MinValue;
                }

                if ((mode == Modes.Batch || mode == Modes.TexBatch) && (fileToOpen == null)) {
                    PrintUsage();
                    return int.MinValue;
                }
            }

            // launch the app
            // ==============
            ilPSP.Environment.Bootstrap(
                new string[0],
                Utils.GetBoSSSInstallDir(),
                out bool mpiInitialized);

            switch (mode) {
                case Modes.Worksheet:
                var ws = new Worksheet(fileToOpen);
                ws.Shown += Worksheet.OnShown; // Workaround for wrong word-wrap on start-up of the application
                System.Windows.Forms.Application.Run(ws);

                ws.m_ExecutorOfCommandQueue_RegularTermination = false;
                Thread.Sleep(800);

                if (ws.m_ExecutorOfCommandQueue.IsAlive) {
                    // hardcore
                    Thread.Sleep(5000);
                    if (ws.m_ExecutorOfCommandQueue.IsAlive) {
                        ws.m_ExecutorOfCommandQueue.Abort();
                    }
                }
                break;

                case Modes.Console:
                ReadEvalPrintLoop.REPL();
                break;

                case Modes.Check:
                InstallationChecker.CheckSetup();
                break;

                case Modes.Batch:
                case Modes.TexBatch:
                Document doc;
                if (fileToOpen.ToLowerInvariant().EndsWith(".tex")) {
                    List<string> dummy;
                    LatexIO.SplitTexFile(fileToOpen, out dummy, out doc);
                } else {
                    doc = Document.Deserialize(fileToOpen);
                }
                string OutDir = Path.GetDirectoryName(fileToOpen);
                string DocNam = Path.GetFileNameWithoutExtension(fileToOpen) + ".texbatch";
                InteractiveShell.CurrentDoc = doc;
                InteractiveShell._CurrentDocFile = (new FileInfo(fileToOpen)).FullName;

                // Which text boxes should be removed before 'restart' occurs
                int f = 0;
                if (mode == Modes.TexBatch) {

                    // bws was produced by Latex - some string replacements are necessary
                    for (int iEntry = 0; iEntry < doc.CommandAndResult.Count; iEntry++) {
                        var Entry = doc.CommandAndResult[iEntry];

                        // Check whether there are boxes before restart
                        if (Entry.Command.Equals("restart") || Entry.Command.Equals("restart;")) {
                            f = iEntry;
                        }

                        Entry.Command = LatexIO.Tex2Bws(Entry.Command);
                    }

                    GnuplotExtensions.UseCairoLatex = true;
                }

                // All boxes before 'restart' should not be counted as error
                int count = 0;
                foreach (Document.Tuple dt in doc.CommandAndResult) {
                    Console.WriteLine(dt.Command);
                    bool success = dt.Evaluate();

                    if (!success && count >= f)
                        errCount++;

                    Console.WriteLine(Document.ResultStartMarker);
                    Console.WriteLine(dt.InterpreterTextOutput);
                    Console.WriteLine(Document.ResultEndMarker);

                    count++;
                }

                if (mode == Modes.TexBatch) {
                    LatexIO.Save_Texbatch(OutDir, DocNam, doc);
                } else {
                    if (fileToOpen.EndsWith(".tex")) {

                    } else {
                        doc.Serialize(fileToOpen);
                    }
                }
                InteractiveShell.CurrentDoc = null;
                break;

                default:
                throw new NotImplementedException();
            }

            if (mpiInitialized)
                csMPI.Raw.mpiFinalize();

            return errCount;
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
            Console.WriteLine("Option 1: Worksheet mode (interactive/GUI):");
            Console.WriteLine("-------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe                       Opens an empty worksheet");
            Console.WriteLine(" or BoSSSpad.exe file.bws              Opens file.bws");
            Console.WriteLine(" or BoSSSpad.exe --worksheet file.bws  Opens file.bws");
            Console.WriteLine();
            Console.WriteLine("Option 2: Console mode (interactive/text):");
            Console.WriteLine("------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --console             Starts a console session.");
            Console.WriteLine();
            Console.WriteLine("Option 3: Batch mode (non-interactive/text):");
            Console.WriteLine("--------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --batch file.bws      Batch execution of commands in file.bws");
            Console.WriteLine();
            Console.WriteLine("Option 4: Installation check:");
            Console.WriteLine("--------------------------------------------");
            Console.WriteLine("    BoSSSpad.exe --check               Check the BoSSS installation.");
            Console.WriteLine();
        }

        /// <summary>
        /// Dummy function which ensures that certain referenced assemblies 
        /// </summary>
        private static void LinkEnforcer()
        {
            // If you remove these lines, this may break some worksheets and tutorials.
            Console.WriteLine(typeof(CNS.Program).FullName);
            Console.WriteLine(typeof(IBM_Solver.IBM_SolverMain).FullName);
        }
    }
}
