using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;

namespace AllSpark {

    static class AllSpark {

        static List<int> CurrentlyRunning = new List<int>();
        static object synclock = new Object();

        class SomeTask {
            protected int i;  // task index
            protected string args = "";
            protected string exefile;
            public bool success = false;

            public void Run() {
                /*
                FileInfo exeFileInfo = new FileInfo(exefile);
                if(!exeFileInfo.Exists) {
                    Console.Error.WriteLine("unable to find file: '" + exeFileInfo.FullName + "'");
                    return;
                }
                psi.FileName = exeFileInfo.FullName;
                psi.WorkingDirectory = exeFileInfo.DirectoryName;
                */
                ProcessStartInfo psi = new ProcessStartInfo();
                psi.FileName = this.exefile;
                psi.Arguments = this.args;

                psi.UseShellExecute = false;
                psi.RedirectStandardOutput = true;
                psi.RedirectStandardError = true;

                Console.WriteLine("starting case {0} (exe '{1}' args '{2}')...", i + 1, psi.FileName, psi.Arguments);

                lock (synclock) {
                    CurrentlyRunning.Add(i);
                }

                using (var stdout = new StreamWriter("out_" + (i + 1) + ".txt")) {
                    using (var p = new Process()) {
                        p.StartInfo = psi;

                        DataReceivedEventHandler handler = delegate (object sender, DataReceivedEventArgs e) {
                            try {
                                stdout.WriteLine(e.Data);
                                stdout.Flush();
                            } catch (Exception ee) {
                                Console.WriteLine("FAILED " + (i + 1) + ": " + ee.Message + " (" + ee.GetType().FullName + ")");
                                throw;
                            }
                            return;
                        };
                        p.OutputDataReceived += handler;
                        p.ErrorDataReceived += handler;

                        try {
                            p.Start();
                            p.BeginErrorReadLine();
                            p.BeginOutputReadLine();
                            p.WaitForExit();
                            success = (p.ExitCode == 0);

                            stdout.Flush();
                            Console.WriteLine("finished " + (i + 1) + ", exit code " + p.ExitCode);
                        } catch (Exception e) {
                            Console.WriteLine("FAILED " + (i + 1) + ": " + e.Message + " (" + e.GetType().FullName + ")");
                            stdout.WriteLine("FAILED " + (i + 1) + ": " + e.Message + " (" + e.GetType().FullName + ")");
                            success = false;
                        }
                    }

                    stdout.Flush();
                }

                int[] running;
                lock (synclock) {
                    CurrentlyRunning.Remove(i);
                    running = CurrentlyRunning.ToArray();

                    Console.Write("Currently running:");
                    foreach (int j in running) {
                        Console.Write(" " + (j + 1));
                    }
                    Console.WriteLine(";");
                }

                Console.Out.Flush();
            }
        }

        class SweepTask : SomeTask {

            public SweepTask(int i, string _exefile, string cfile) {
                base.i = i;
                base.exefile = _exefile;
                base.args = " --control " + cfile + " --pstudy_case " + (i + 1).ToString();

                //psi.UseShellExecute = false;
                //psi.RedirectStandardOutput = true;
                //psi.RedirectStandardError = true;


            }
        }

        class GeneralTask : SomeTask {
            public GeneralTask(int i, string statement) {
                base.i = i;
                //Console.WriteLine(i + ": " + statement);

                var parts = Regex.Matches(statement, @"[\""].+?[\""]|[^ ]+")
                                .Cast<Match>()
                                .Select(m => m.Value)
                                .ToList();

                string exefile = parts.First();
                base.exefile = exefile.Replace("\"", "");

                if (parts.Count > 1) {
                    string args = "";
                    for (int k = 1; k < parts.Count; k++) {
                        args = args + parts[k];
                        if (k < parts.Count - 1)
                            args = args + " ";
                    }
                    base.args = args;
                }

            }
        }

        static int PrintUsage() {
            Console.WriteLine("Usage:");
            Console.WriteLine("AllSpark.exe sweep $num_threads $first_case $last_case $exe $control_file");
            //Console.WriteLine(" -or- ");
            //Console.WriteLine("AllSpark.exe sweep2 $num_threads $case_number $case_number ... $exe $control_file");
            Console.WriteLine(" -or- ");
            Console.WriteLine("AllSpark.exe general $num_threads $statement_1  $statement_1 ... $statement_N");
            Console.WriteLine();

            return -666;
        }

        static int Main(string[] args) {
            if (args.Length < 1) {
                return PrintUsage();
            }
            string mode = args[0];

            SomeTask[] tasks = null;
            int threads = 1;
            try {
                threads = int.Parse(args[1]); // arg 1: number of threads

                switch (mode) {
                    case "sweep": {
                            int IE = int.Parse(args[3]); // arg 3: last task (one-based)
                            int I0 = int.Parse(args[2]); // arg 2: first task

                            tasks = new SomeTask[IE - I0 + 1];
                            for (int i = I0; i <= IE; i++) {
                                tasks[i - I0] = new SweepTask(
                                    i - 1,
                                    args[4],  // arg 4: exe file
                                    args[5]   // arg 5: control file
                                    );
                            }

                            break;
                        }

                    case "general": {
                            tasks = new SomeTask[args.Length - 2];
                            for (int i = 2; i < args.Length; i++) {
                                //Console.WriteLine("#" + i + ":   " + args[i]);
                                tasks[i - 2] = new GeneralTask(i - 2, args[i].Replace("@@", "  "));
                            }

                            break;
                        }

                    default:
                        return PrintUsage();
                }
            } catch (Exception) {
                return PrintUsage();
            }

            var options = new ParallelOptions() { MaxDegreeOfParallelism = threads };
            Parallel.Invoke(options, tasks.Select(t => ((Action)t.Run)).ToArray());

            for (int i = 0; i < tasks.Length; i++) {
                if (!tasks[i].success) {
                    return i + 1;
                }
            }
            return 0;

            //tasks[0].RunExample();
        }
    }
}
