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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;


namespace MiniBatchProcessor {

    


    /// <summary>
    /// Server functionality.
    /// </summary>
    public class Server : ClientAndServer {

        /// <summary>
        /// %
        /// </summary>
        public Server(string BatchInstructionDir) : base(BatchInstructionDir) {

        }

        //static List<int> CurrentlyRunning = new List<int>();
        //static object synclock = new Object();

        FileStream ServerMutex;

        /// <summary>
        /// True, if the server is hosed by the current process and is running, false if not; 
        /// This state is determined based upon the existence, resp. a lock on the file
        /// `MiniBatchProcessor-ServerMutex.txt`.
        /// </summary>
        public static bool GetIsRunning(string __BatchInstructionDir) {
            __BatchInstructionDir ??= Configuration.GetDefaultBatchInstructionDir();

            {
                if(MiniBatchThreadIsMyChild) {
                    if(ServerInternal != null) {
                        return ServerInternal.IsAlive;
                    }

                    if(ServerExternal != null) {
                        return !ServerExternal.HasExited;
                    }
                }


                try {
                    

                    if(File.Exists(GetServerMutexPath(__BatchInstructionDir))) {
                        
                        // try to write exclusively; this should produce an exception in Linux and probably also Windows, if it is opened by any other process)
                        using(var ServerMutex2 = File.Open(GetServerMutexPath(__BatchInstructionDir), FileMode.Create, FileAccess.Write, FileShare.None)) {
                            ServerMutex2.WriteByte(1);
                            ServerMutex2.Close();
                        }
                                                
                        // try to delete; was sufficient to produce an exsption on Windows, (if opened by any other process), but unreliable on Linux
                        File.Delete(GetServerMutexPath(__BatchInstructionDir));

                        return false;
                    } else {
                        return false;
                    }

                } catch(IOException) {
                    return true;
                }
            }
        }

        /// <summary>
        /// If <see cref="StartIfNotRunning"/> has been used, the external process which is running the server.
        /// </summary>
        public static Process ServerExternal {
            get;
            private set;
        }

        /// <summary>
        /// If <see cref="StartIfNotRunning"/> has been used, the internal thread which is running the server.
        /// </summary>
        public static Thread ServerInternal {
            get;
            private set;
        }

        /// <summary>
        /// Starts the server in an external process, or in an internal background-thread if it is not already running.
        /// </summary>
        /// <param name="RunExternal">
        /// If true, an external process is used.
        /// </param>
        /// <returns>
        /// - true: the server was just started
        /// - false: the server is already running
        /// </returns>
        public static bool StartIfNotRunning(bool RunExternal = true, string BatchInstructionDir = null, bool Reset = false) {
            BatchInstructionDir ??= Configuration.GetDefaultBatchInstructionDir();
            
            if (ServerExternal != null && ServerExternal.HasExited) {
                ServerExternal.Dispose();
                ServerExternal = null;
            }

            if(ServerInternal != null && !ServerInternal.IsAlive) {
                ServerInternal = null;
            }

            if(ServerExternal != null || ServerInternal != null || GetIsRunning(BatchInstructionDir)) {
                Console.WriteLine("Mini batch processor is already running.");
                return false;
            }
            Random r = new Random();
            r.Next(100, 5000);
            if(GetIsRunning(BatchInstructionDir)) {
                Console.WriteLine("Mini batch processor is already running.");
                return false;
            }

            if(RunExternal) {
                Console.WriteLine("Starting mini batch processor in external process...");

                ProcessStartInfo psi = new ProcessStartInfo() { 
                    UseShellExecute = true // default value changed from .net framework to .net, set this true, so an external window is opened
                };
                psi.FileName = "dotnet";
                psi.Arguments = typeof(Server).Assembly.Location;

                Process p = Process.Start(psi);

                Console.WriteLine("Started mini batch processor on local machine, process id is " + p.Id + ".");

                ServerExternal = p;
                MiniBatchThreadIsMyChild = true;

            } else {
                Console.WriteLine("Starting mini batch processor in background thread...");

                Thread ServerThread = (new Thread(delegate () {
                    var s = new Server(null);
                    s._Main(new string[0], Reset);
                }));
                ServerThread.Start();

                ServerInternal = ServerThread;
                MiniBatchThreadIsMyChild = true;
            }

            int timeout = 0;
            while(!GetIsRunning(BatchInstructionDir)) {
                Console.WriteLine("waiting for startup...");
                Thread.Sleep(4875);
                timeout++;
                if(timeout > 20) {
                    Console.WriteLine("giving up.");
                    if(ServerExternal != null) {
                        if(!ServerExternal.HasExited) {
                            try {
                                ServerExternal.Kill();
                            } catch(Exception e) {
                                Console.Error.WriteLine($"Server startup fail: {e.GetType().Name}: {e.Message}");
                            }
                        }
                        ServerExternal = null;
                    }

                    if(ServerInternal != null) {
                        if(ServerInternal.IsAlive) {
                            try {
                                MiniBatchProcessor.Server.SendTerminationSignal();
                            } catch(Exception e) {
                                Console.Error.WriteLine($"Server startup fail: {e.GetType().Name}: {e.Message}");
                            }
                        }
                        ServerInternal = null;
                    }

                    return false;
                }
            }

           

            Console.WriteLine("started.");
            return true;
        }

        /// <summary>
        /// True, if the batch thread/process was started by this process;
        /// </summary>
        public static bool MiniBatchThreadIsMyChild {
            get;
            private set;
        }

        /// <summary>
        /// Sends a signal which should terminate a running server instance linked to the directory <paramref name="BatchInstructionDir"/>.
        /// </summary>
        public static void SendTerminationSignal(bool WaitForOtherJobstoFinish = true, int TimeOutInSeconds = 1800, string BatchInstructionDir = null) {
            DateTime st = DateTime.Now;

            if(BatchInstructionDir == null) {
                BatchInstructionDir = Configuration.GetDefaultBatchInstructionDir();
            }
            
            if (WaitForOtherJobstoFinish) {

                Client cln = new(BatchInstructionDir);


                Random rnd = new();

                JobData[] Q = cln.Queue.ToArray();
                JobData[] W = cln.Working.ToArray();

                while (Q.Length > 0 || W.Length > 0) {
                    Console.WriteLine($"Waiting for other jobs to finish; in queue: {Q.Length}, working: {W.Length}.");
                    foreach (var j in Q) {
                        Console.WriteLine(j);
                    }
                    foreach (var j in W) {
                        Console.WriteLine(j);
                    }
                    
                    Thread.Sleep(rnd.Next(60001));
                    if(TimeOutInSeconds > 0) {
                        var dur = DateTime.Now - st;
                        if(dur.TotalSeconds > TimeOutInSeconds) {
                            Console.WriteLine($" Waiting for {dur.TotalSeconds}; timeout of {TimeOutInSeconds} reached.");
                        }
                    }

                    Q = cln.Queue.ToArray();
                    W = cln.Working.ToArray();
                }
            
            }

            using (var s = File.Create(GetTerminationSignalPath(BatchInstructionDir))) {
                s.Flush();
                s.Close();
            }

            if(MiniBatchThreadIsMyChild) {
                while(Server.GetIsRunning(BatchInstructionDir)) {
                    
                    Thread.Sleep(1000 + ClientAndServer.IOwaitTime);
                }

                MiniBatchThreadIsMyChild = false;
            }

        }

        static string GetTerminationSignalPath(string BatchInstructionDir) {
            if (BatchInstructionDir == null) {
                BatchInstructionDir = Configuration.GetDefaultBatchInstructionDir();
            }
            return (Path.Combine(BatchInstructionDir, "MiniBatchProcessor-TerminationSignal.txt"));
        }


        static string GetServerMutexPath(string BatchInstructionDir) {
            if (BatchInstructionDir == null) {
                BatchInstructionDir = Configuration.GetDefaultBatchInstructionDir();
            }
            return (Path.Combine(BatchInstructionDir, "MiniBatchProcessor-ServerMutex.txt"));
        }

        string LogFilePath {
            get {
                return (Path.Combine(config.BatchInstructionDir, "MiniBatchProcessor-Log.txt"));
            }
        }

        bool Init(bool Reset) {

            
            string MutexFileName = Path.Combine(GetServerMutexPath(config.BatchInstructionDir));
            try {
                // Check that we are the only instance:
                if(File.Exists(MutexFileName)) {
                    // try to delete
                    File.Delete(MutexFileName);
                }
            } catch(IOException) {
                Console.WriteLine("Unable to obtain server mutex.");
                return false;
            }

            try { 
                if(File.Exists(GetTerminationSignalPath(config.BatchInstructionDir))) {
                    // try to delete
                    File.Delete(GetTerminationSignalPath(config.BatchInstructionDir));
                }
            } catch (IOException) {
                Console.WriteLine("Unable to delete termination signal");
                return false;
            }


            // create mutex file and share it with no one!
            try {
                ServerMutex = File.Open(MutexFileName, FileMode.Create, FileAccess.Write, FileShare.None);
                LogFile = new StreamWriter(File.Open(LogFilePath, FileMode.Append, FileAccess.Write, FileShare.ReadWrite));
                var ServerMutexS = new StreamWriter(ServerMutex);
                ServerMutexS.WriteLine("This file is used by the MiniBatchProcessor to ensure that only");
                ServerMutexS.WriteLine("one instance of the batch processor per computer/user is running.");
                ServerMutexS.Flush();
            } catch(IOException) {
                LogMessage("Unable to obtain server mutex (2).", error:true);
                return false;
            }

            // Rename present queue and working directory to start fresh
            if(Reset) {
                var date = DateTime.Now.ToString("ddMMMyyyy_HHmmss");
                foreach(string q in new[]{ 
                    Path.Combine(config.BatchInstructionDir, QUEUE_DIR),
                    Path.Combine(config.BatchInstructionDir, WORK_FINISHED_DIR) }) {

                    if(Directory.Exists(q)) {
                        /*
                        foreach(string f in Directory.GetFiles(q)) {
                            try {
                                File.Delete(f);
                            } catch(Exception) {
                                // ignore
                            }
                        }
                        */
                        Directory.Move(q, q + ".bkup-" + date);
                    }

                }
            }



            // see if there are any zombies left in the 'working' directory
            foreach (var J in Working) {
                //MoveWorkingToFinished(J);
                string exitTokenPath = Path.Combine(
                    config.BatchInstructionDir,
                    ClientAndServer.WORK_FINISHED_DIR,
                    J.ID.ToString() + "_exit.txt");
                if(!File.Exists(exitTokenPath)) {
                    File.WriteAllText(exitTokenPath, "-9876");
                }
            }

          

            return true;
        }
               

        void MoveQueueToWorking(JobData J) {
            string BaseDir = config.BatchInstructionDir;
            foreach(string nmn in new string[] { J.ID.ToString() }) {
                var Src = Path.Combine(BaseDir, ClientAndServer.QUEUE_DIR, nmn);
                var Dst = Path.Combine(BaseDir, ClientAndServer.WORK_FINISHED_DIR, nmn);


                if(File.Exists(Src)) {
                    int ReTryCount = 0;
                    while(true) {
                        try {
                            File.Copy(Src, Dst);
                            base.UpdateStatus(J.ID, JobStatus.Working);
                            return;
                        } catch(Exception e) {
                            if(ReTryCount < ClientAndServer.IO_OPS_MAX_RETRY_COUNT) {
                                ReTryCount++;
                                Thread.Sleep(ClientAndServer.IOwaitTime);
                            } else {
                                LogMessage(string.Format("{0} while trying to move file '{1}' to '{2}', message: {3}.", e.GetType().Name, Src, Dst, e.Message), error:true);
                                return;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Checks if job configuration is legal and can be started on this machine.
        /// </summary>
        /// <param name="J"></param>
        /// <param name="WriteInfo">
        /// If true, error messages are written to stderr.
        /// </param>
        bool CheckJob(JobData J, bool WriteInfo) {
            if (J.NoOfProcs > config.MaxProcessors) {
                if (WriteInfo) {
                    LogMessage(string.Format("Job #{0} is to big for this machine: configured to use {1} processors, max. allowed is {2}.", J.ID, J.NoOfProcs, config.MaxProcessors), error:true);
                }
                return false;
            }

            if (J.NoOfProcs <= 0) {
                if (WriteInfo) {
                    LogMessage(string.Format("Illegal configuration for job #{0}: cannot run with {1} processors.", J.ID, J.NoOfProcs), error:true);
                }
                return false;
            }

            return true;
        }

        private static readonly object padlock = new object();

        static string LastMessage = "";

        static StreamWriter LogFile;

        /// <summary>
        /// Writes a message to stdout, but only if it differs from the last message.
        /// </summary>
        internal static void LogMessage(string m, bool error = false) {
            lock(padlock) {
                if(!m.Equals(LastMessage)) {
                    string fm = DateTime.Now + ": " + m;
                    if(error)
                        fm = fm + "[ERROR] ";
                    if(ServerInternal == null) {
                        if(!error)
                            Console.WriteLine(fm);
                    }
                    if(error)
                        Console.Error.WriteLine(fm);
                    if (LogFile != null) {
                        LogFile.WriteLine(fm);
                        LogFile.Flush();
                    }
                    LastMessage = m;
                }
            }
        }


        /// <summary>
        /// Entry Point for server running as separate process
        /// </summary>
        public static int Main(string[] args) {

            string dir = null;
            if(args.Length > 0) {
                if (args[0].Equals("terminate")) {
                    Console.WriteLine("sending termination signal.");
                    SendTerminationSignal(WaitForOtherJobstoFinish:false);
                    return 0;
                }


                dir = args[0];
                if(!Directory.Exists(dir)) {
                    throw new IOException("specified batch directory '" + dir + "' does not exist.");
                }
            }

            if (GetIsRunning(args.Length > 0 ? args[0] : null)) {
                Console.WriteLine("MiniBatchProcessor server is already running -- only one instance is allowed, terminating this one.");
                return -777;
            }

            var s = new Server(dir);
            var r = s._Main(args, false);

            return r ? -777 : 0;
        }
       
        /// <summary>
        /// Main routine of the server; continuously checks the respective 
        /// directories (e.g. <see cref="ClientAndServer.QUEUE_DIR"/>) and runs jobs in external processes.
        /// </summary>
        bool _Main(string[] args, bool Reset) {
            
            if(!Init(Reset)) {
                Console.WriteLine("Terminating server init.");
                return false;
            }
            LogMessage("server started.");

            // infinity loop
            // =============

            //int AvailableProcs = Math.Max(1, Math.Min(Environment.ProcessorCount, config.MaxProcessors));
            int AvailableProcs = Math.Max(1, config.MaxProcessors);
            bool ExclusiveUse = false;

            var Running = new List<Tuple<Thread, ProcessThread>>();
            bool keepRunning = true;
            int PrevRunning = -1, PrevQueue = -1;
            while (keepRunning) {
                if(File.Exists(GetTerminationSignalPath(config.BatchInstructionDir))) {
                    // try to delete
                    File.Delete(GetTerminationSignalPath(config.BatchInstructionDir));
                    LogMessage("Received termination signal: server will be terminated, but jobs may continue.");
                    keepRunning = false;
                    ServerMutex.Close();
                    ServerMutex.Dispose();
                    ServerMutex = null;

                    LogFile.Close();
                    LogFile.Dispose();
                    LogFile = null;

                    File.Delete(GetServerMutexPath(config.BatchInstructionDir));
                    return true;
                }

                // see if any of the running jobs is finished
                for (int i = 0; i < Running.Count; i++) {
                    var TT = Running[i];
                    if (TT.Item1.IsAlive == false) {
                        Running.RemoveAt(i);
                        //MoveWorkingToFinished(TT.Item2.data);
                        AvailableProcs += TT.Item2.data.NoOfProcs;
                        ExclusiveUse = false;
                        i--;
                    }
                }

                var NextJobs = Queue.ToArray();

                if(NextJobs.Count() > 0 && NextJobs.Count() != PrevQueue) {
                    LogMessage("Number of jobs in queue: " + NextJobs.Count());
                }
                PrevQueue = NextJobs.Count();

                if(Running.Count > 0 && Running.Count != PrevRunning) {
                    LogMessage("Currently running " + Running.Count + " jobs.");
                }
                PrevRunning = Running.Count;

                // see if there are available processors
                if (AvailableProcs <= 0 || ExclusiveUse) {
                    Thread.Sleep(config.ServerPollingInSeconds*1000);
                    continue;
                }

                // sort out jobs which have problems
                NextJobs = NextJobs.Where(job => CheckJob(job, true) == true).ToArray();

                // sleep if there is nothing to do
                if (NextJobs.Count() <= 0) {
                    LogMessage(string.Format("No more jobs in queue. Running: {0}, Avail. procs.: {1}.", Running.Count, AvailableProcs));
                    Thread.Sleep(10000);
                    continue;
                } 

                // priorize (at the moment, only by ID)
                NextJobs = NextJobs.OrderBy(job => job.ID).ToArray();

                if (config.BackFilling) {
                    throw new NotImplementedException("todo");
                } else {
                    if ((NextJobs[0].NoOfProcs > AvailableProcs)
                        || (NextJobs[0].UseComputeNodesExclusive && AvailableProcs < config.MaxProcessors)
                        || (ExclusiveUse == true)) {
                        LogMessage($"Not enough available processors for job #{NextJobs[0].ID} - start is delayed; {NextJobs.Length} jobs are pending.");
                        Thread.Sleep(10000);
                        continue;
                    }
                }

                // Launch next job.
                {
                    ProcessThread task = new ProcessThread() {
                        data = NextJobs[0],
                        WorkDir = Path.Combine(config.BatchInstructionDir, ClientAndServer.WORK_FINISHED_DIR),
                        BatchInstructionDir = config.BatchInstructionDir,
                        TermAct = delegate() { base.UpdateStatus(NextJobs[0].ID, JobStatus.Finished);  }
                    };
                    Thread th = new Thread(task.Run);
                    th.Priority = ThreadPriority.Lowest;
                    Running.Add(new Tuple<Thread, ProcessThread>(th, task));
                    AvailableProcs -= task.data.NoOfProcs;
                    ExclusiveUse = task.data.UseComputeNodesExclusive;
                    MoveQueueToWorking(NextJobs[0]);

                    th.Start();
                }
            }

            return true;
        }

        /// <summary>
        /// Starts job in some external process and waits until it terminates.
        /// </summary>
        class ProcessThread {

            /// <summary>
            /// Job which is started.
            /// </summary>
            public JobData data;

            /// <summary>
            /// 
            /// </summary>
            public string WorkDir;

            /// <summary>
            /// 
            /// </summary>
            public string BatchInstructionDir;

            /// <summary>
            /// True if the process exited without errors.
            /// </summary>
            public bool success = false;

            /// <summary>
            /// Executed on termination
            /// </summary>
            public Action TermAct = null;

            /// <summary>
            /// Starts job in some external process and waits until it terminates.
            /// </summary>
            public void Run() {

                //FileInfo exeFileInfo = new FileInfo(exefile);
                //if(!exeFileInfo.Exists) {
                //    Console.Error.WriteLine("unable to find file: '" + exeFileInfo.FullName + "'");
                //    return;
                //}
                //psi.FileName = exeFileInfo.FullName;
                //psi.WorkingDirectory = exeFileInfo.DirectoryName;

                ProcessStartInfo psi = new ProcessStartInfo();
                switch (System.Environment.OSVersion.Platform) {
                    case System.PlatformID.Unix:
                        //psi.FileName = "mpirun";
                        //psi.Arguments = " -np " + data.NoOfProcs + "  " + data.exefile + " ";
                        if(data.NoOfProcs <= 1) {

                            // This is some "speciality" due to our CI workflow
                            // On our Linux server (Ubuntu 20.04.6 LTS, OpenMPI 4.0.3), 
                            // when we are using `mpirun` or `mpiexe` to execute validation tests
                            // (i.e., tests which execute a Jupyter worksheet)
                            // the worksheet crashes on MPI_Init;
                            // If executed without `mpirun` or `mpiexec`, the respective tests work.
                            //
                            // This is a somewhat known issue: https://github.com/open-mpi/ompi/issues/11063
                            // When we set the envrinoment variable "PMIX_MCA_gds" to "hash", this resolves the problem partly,
                            // (e.g., add `psi.EnvironmentVariables.Add("PMIX_MCA_gds", "hash");` here.)
                            // i.e. the worksheet test sucseeds, but the test runner seems to hang indefinitly during exit.
                            //
                            // However, since worksheets only use one processor anyway, 
                            // avoiding `mpirun` is a sufficient workaround for now (fk, 03aug23).
                            //
                            Server.LogMessage($"Skipping mpirun on {data.exefile} because MPI size is 1 anyway."); 
                            psi.FileName = data.exefile;

                        } else {
                            psi.FileName = "mpirun";
                            psi.Arguments = " -np " + data.NoOfProcs + "  " + data.exefile + " ";
                        }
                        break;
                    default:
                        psi.FileName = "mpiexec.exe";
                        psi.Arguments = " -n " + data.NoOfProcs + " " + data.exefile + " ";
                    break;
                }

                foreach (var a in data.Arguments)
                    psi.Arguments += a + " ";
                psi.WorkingDirectory = data.ExeDir;

                for (int i = 0; i < data.EnvVars.Length; i++) {
                    psi.EnvironmentVariables.Add(data.EnvVars[i].Item1, data.EnvVars[i].Item2);
                }

                psi.UseShellExecute = false;
                psi.RedirectStandardOutput = true;
                psi.RedirectStandardError = true;

                Server.LogMessage(string.Format("starting job #{2}, '{3}': {0} {1}", psi.FileName, psi.Arguments, data.ID, data.Name));

                // note: we create the stdout and stderr file on the output-directory to have a 
                //       constant path of these files for all times.
                string stdout_file = Path.Combine(WorkDir, "..", ClientAndServer.WORK_FINISHED_DIR, data.ID.ToString() + "_out.txt");
                string stderr_file = Path.Combine(WorkDir, "..", ClientAndServer.WORK_FINISHED_DIR, data.ID.ToString() + "_err.txt");

                try {
                    using(FileStream stdoutStream = new FileStream(stdout_file, FileMode.Create, FileAccess.Write, FileShare.Read),
                        stderrStream = new FileStream(stderr_file, FileMode.Create, FileAccess.Write, FileShare.Read)) {
                        using(StreamWriter stdout = new StreamWriter(stdoutStream),
                            stderr = new StreamWriter(stderrStream)) {
                            try {

                                var p = new Process();
                                p.StartInfo = psi;

                                p.OutputDataReceived += delegate (object sender, DataReceivedEventArgs e) {
                                    try {
                                        stdout.WriteLine(e.Data);
                                        stdout.Flush();
                                        return;
                                    } catch(Exception ex) {
                                        Server.LogMessage(string.Format("STDOUT file exception, unable to write " + e.Data + "; (job " + this.data.Name + ", " + ex.Message + ", " + ex.GetType().FullName + ")"));
                                        return;
                                    }
                                };
                                p.ErrorDataReceived += delegate (object sender, DataReceivedEventArgs e) {
                                    try {
                                        stdout.WriteLine(e.Data);
                                        stdout.Flush();
                                        return;
                                    } catch(Exception ex) {
                                        Server.LogMessage(string.Format("STDERR file exception, unable to write " + e.Data + "; (job " + this.data.Name + ", " + ex.Message + ", " + ex.GetType().FullName + ")"));
                                        return;
                                    }
                                };

                                p.Start();
                                p.BeginErrorReadLine();
                                p.BeginOutputReadLine();

                                while(!p.HasExited) {
                                    Thread.Sleep(5000);
                                }

                                success = (p.ExitCode == 0);
                                Server.LogMessage(string.Format("finished job #" + data.ID + " (" + data.Name + "), exit code " + p.ExitCode + "."));

                                using(var exit = new StreamWriter(Path.Combine(
                                    BatchInstructionDir,
                                    ClientAndServer.WORK_FINISHED_DIR,
                                    data.ID.ToString() + "_exit.txt"))) {
                                    exit.WriteLine(p.ExitCode);
                                    exit.Flush();
                                }

                                // wait a little bit longer before streams get closed - sometimes it seems
                                // stdout and stderr send data after process has exited.
                                Thread.Sleep(5000);

                            } catch(Exception e) {
                                Server.LogMessage(string.Format("FAILED " + psi.FileName + " " + psi.Arguments + ": " + e.Message + " (" + e.GetType().FullName + ")"));
                                stdout.WriteLine(DateTime.Now + ": FAILED " + psi.FileName + " " + psi.Arguments + ": " + e.Message + " (" + e.GetType().FullName + ")");
                                success = false;
                            }
                            stdout.Flush();
                            stderr.Flush();
                        } // writers
                    } // streams 

                    TermAct?.Invoke();

                } catch (Exception e) {

                    Console.WriteLine(e.GetType().Name + " in runner thread: " + e.Message);
                    Console.WriteLine(e.StackTrace);
                }

                Console.Out.Flush();
            } // Run() -- method
        }
    }
}
