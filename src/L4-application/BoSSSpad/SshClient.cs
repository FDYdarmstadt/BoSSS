using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Text.RegularExpressions;
using ilPSP.Tracing;
using System.Text;
using System.Threading;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.GridImport;
using log4net;
using ilPSP;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// internal wrapper around ssh/slurm commands
    /// </summary>
    abstract class SshClient : IDisposable {



        /// <summary>
        /// ctor
        /// </summary>
        public SshClient(string ServerName, string Username, PrivateKeyFile pkf, string sshclientExeToUse) {
            m_pkf = pkf;
            m_srvrname = ServerName;
            m_usrname = Username;
            m_sshclientExeToUse = sshclientExeToUse;
        }


        abstract public void Dispose();

        public SshClient(string ServerName, string Username, string pwd, string sshclientExeToUse) {
            throw new NotSupportedException("TO DO");
        }


        protected abstract void Connect();
       

        protected bool TestConnection() {
            string search = RunCommand("ls").stdout;
            bool connected = search.Contains(m_usrname);
            return connected;
        }

        readonly private PrivateKeyFile m_pkf;
        readonly private string m_usrname;
        readonly private string m_srvrname;
        readonly string m_sshclientExeToUse;


        public string SshclientExeToUse {
            get { 
                if (m_sshclientExeToUse.IsEmptyOrWhite()) {
                    return "ssh";
                } else {
                    return m_sshclientExeToUse;
                }
            }
        }

        public string KeyFilePath {
            get { return m_pkf.Path; }
        }

        public string UserName {
            get { return m_usrname; }
        }

        public string ServerName {
            get { return m_srvrname; }
        }

        public bool IsConnected {
            get {
                return TestConnection();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        internal string SubmitJob(string remotepath, out string _stdout, out string _stderr) {
            using (var tr = new FuncTrace()) {
                Connect();
                string sbatchCmd = "sbatch " + remotepath + "/batch.sh";
                var (resultString, err) = RunCommand(sbatchCmd);
                _stdout = resultString;
                _stderr = err;

                String SearchString = "Submitted batch job ";
                String jobId = Regex.Match(resultString, SearchString + "[0-9]*") // look for SearchString followed by a number (the Job ID)
                    .ToString() // convert to string
                    .Replace(SearchString, ""); // remove SearchString, leaving only the Job ID
                return jobId;
            }
        }

        abstract public (string stdout, string stderr) RunCommand(string command, bool verbose = false);


        

    }


    /// <summary>
    /// internal wrapper around ssh/slurm commands;
    /// slow version, which starts a separate ssh connection for every command.
    /// </summary>
    class SlowSshClient : SshClient {

        /// <summary>
        /// ctor
        /// </summary>
        public SlowSshClient(string ServerName, string Username, PrivateKeyFile pkf, string sshClientExeToUse) : base(ServerName, Username, pkf, sshClientExeToUse) {
         
            PlatformID CurrentSys = System.Environment.OSVersion.Platform;
            string shell = "";
            switch (CurrentSys) {
                case PlatformID.Unix:
                    shell = "bash";
                    break;
                case PlatformID.Win32NT:
                    shell = "cmd";
                    break;
                default:
                    throw new NotImplementedException("Unkonwn OS!");
            }

            Process cmd = new Process() {
                StartInfo = new ProcessStartInfo {
                    FileName = shell,
                    CreateNoWindow = true,
                    UseShellExecute = false,
                    RedirectStandardInput = true,
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                }
            };
            m_cmd = cmd;
        }


        public override void Dispose() {
            m_cmd.Dispose();
        }


        public SlowSshClient(string ServerName, string Username, string pwd, string sshClientExeToUse) : base(ServerName, Username, pwd, sshClientExeToUse) {
            throw new NotSupportedException("TO DO");
        }

        protected override void Connect() {
            if (!TestConnection()) {
                string std, err;
                ReadLines(out std, out err);
                Console.WriteLine(err);
                Console.WriteLine("Connection could not be established");
            }
        }

        private void ReadLines(out string std, out string err) {
            using (TextWriter stdoutW = new StringWriter(), stderrW = new StringWriter()) {
                while (m_cmd.StandardOutput.Peek() > -1) {
                    stdoutW.WriteLine(m_cmd.StandardOutput.ReadLine());
                }

                while (m_cmd.StandardError.Peek() > -1) {
                    stderrW.WriteLine(m_cmd.StandardError.ReadLine());
                }

                std = stdoutW.ToString();
                err = stderrW.ToString();
            }
        }

       

        private Process m_cmd;
               
        

        override public (string stdout, string stderr) RunCommand(string command, bool verbose = false) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = verbose;
                tr.Info("Starting test shell...");
                m_cmd.Start();
                tr.Info("started.");

                string sshCmd = "ssh " + this.UserName + "@" + this.ServerName + " -oStrictHostKeyChecking=no \"" + command + "\"";
                tr.Info("Command: " + sshCmd);
                m_cmd.StandardInput.WriteLine(sshCmd);
                tr.Info("command written; Waiting for completion...");
                //m_cmd.StandardInput.WriteLine(command);
                m_cmd.StandardInput.Flush();
                m_cmd.StandardInput.Close();
                m_cmd.WaitForExit();
                tr.Info("external shell terminated; exit code is " + m_cmd.ExitCode);
                string std, err;
                ReadLines(out std, out err);
                tr.Info("stdout received: " + std + "----- (end of stdout)");
                if (err.Length > 0)
                    tr.Error("stderr received: " + err + "----- (end of stderr)");


                return (std, err);
            }
        }



    }

    /// <summary>
    /// SSH client which tries to keep the SSH connection open,
    /// so only a single connection is required.
    /// Thus, it should be much faster than the <see cref="SlowSshClient"/>.
    /// </summary>
    class SingleSessionSshClient : SshClient {

        static Dictionary<string,int> instance_conter = new Dictionary<string,int>();

        public SingleSessionSshClient(string ServerName, string Username, PrivateKeyFile pkf, string sshClientExeToUse) : base(ServerName, Username, pkf, sshClientExeToUse) {
            ConstructorCommon(ServerName, Username);
        }

        private static void ConstructorCommon(string ServerName, string Username) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                string keyname = Username + "@" + ServerName;

                if (!instance_conter.ContainsKey(keyname)) {
                    instance_conter.Add(keyname, 0);
                }
                instance_conter[keyname]++;
                tr.Info($"Instance #{instance_conter[keyname]} of ssh client {keyname} instantiated.");

                if (instance_conter[keyname] > 10) {
                    tr.Warning($"Creating tons ssh-connections to same server ({keyname}) is bad.");
                }
            }
        }

        public SingleSessionSshClient(string ServerName, string Username, string pwd, string sshClientExeToUse) : base(ServerName, Username, pwd, sshClientExeToUse) {
            ConstructorCommon(ServerName, Username);
        }


        public override (string stdout, string stderr) RunCommand(string command, bool verbose = false) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = verbose;

                Connect();
                //Console.ForegroundColor = ConsoleColor.Yellow;
                //Console.BackgroundColor = ConsoleColor.Red;
                //Console.WriteLine(command);
                //Console.ResetColor();
                
                process.StandardInput.Write(command + "\n");
                var stdout = WaitForPrompt();

                return (stdout, "");
            }
        }

        Process process;
        AutoResetEvent outputWaitHandle;
        AutoResetEvent errorWaitHandle;
        Task processActiveTask;

        List<string> outCollector = new List<string>();

        void append(string r) {
            lock (outCollector) {
                outCollector.Add(r);
            }
        }

        void ClearOut() {
            lock (outCollector) {
                outCollector.Clear(); ;
            }
        }

        /// <summary>
        /// tests if the output of some command is finished, i.e. if we are waiting for a prompt
        /// </summary>
        bool IsAtPrompt() {
            lock (outCollector) {
                if (outCollector.Count <= 0)
                    return false;
                if (outCollector.Count > 1)
                    if (outCollector[outCollector.Count - 1].TrimEnd().EndsWith("]$"))
                        return true;
                if (outCollector.Count > 2)
                    if (outCollector[outCollector.Count - 2].TrimEnd().EndsWith("]$"))
                        return true;
                return false;
            }
        }


        void OnOutputDataReceived(object sender, DataReceivedEventArgs e) {
            if (e.Data == null) {
                outputWaitHandle.Set();
            } else {
                // in orderto remove escaped sequences regarding font color, etc., 
                // which can be present in the output, we use the following regex:
                // (from: https://social.msdn.microsoft.com/Forums/en-US/635f240a-e586-4a79-8350-b85111365366/how-can-i-strip-off-control-characters-from-a-string?forum=csharpgeneral)
                string line = Regex.Replace(e.Data, @"\e\[(\d+;)*(\d+)?[ABCDHJKfmsu]", "");
                m_Logger.Info("SSH stdout:" + line);
                append(line);
            }
        }
        void OnErrorDataReceived(object sender, DataReceivedEventArgs e) {
            if (e.Data == null) {
                errorWaitHandle.Set();
            } else {
                // in orderto remove escaped sequences regarding font color, etc., 
                // which can be present in the output, we use the following regex:
                // (from: https://social.msdn.microsoft.com/Forums/en-US/635f240a-e586-4a79-8350-b85111365366/how-can-i-strip-off-control-characters-from-a-string?forum=csharpgeneral)
                string line = Regex.Replace(e.Data, @"\e\[(\d+;)*(\d+)?[ABCDHJKfmsu]", "");
                m_Logger.Info("SSH stdout:" + line);
                append(line);
            }
        }


        string FormatAllReturn() {
            using (var stw = new StringWriter()) {
                lock (outCollector) {
                    foreach(var l in outCollector) {
                        stw.WriteLine(l);
                    }
                }
                return stw.ToString();
            }
        }

        const int timeout = 10 * 1000;

        public bool Disconnect(bool killImmediate) {
            void DisposeObjects() {
                if (process != null) {
                    process.Dispose();
                    process = null;
                }

                /*
                if (outputWaitHandle != null) {
                    outputWaitHandle.Dispose();
                    outputWaitHandle = null;
                }

                if (errorWaitHandle != null) {
                    errorWaitHandle.Dispose();
                    errorWaitHandle = null;
                }
                */

                if (processActiveTask != null) {
                    processActiveTask.Dispose();
                    processActiveTask = null;
                }
            }

            if (process != null && killImmediate == false) {
                //
                // gracefully closing the ssh connection
                //
                if (!process.HasExited) {

                    process.StandardInput.Write("exit\n");
                    process.StandardInput.Flush();

                    bool finished = processActiveTask.Wait(timeout);
                    bool outW = outputWaitHandle.WaitOne(timeout);
                    bool errW = errorWaitHandle.WaitOne(timeout);
                    bool all = finished && outW && errW;

                    //Console.Error.WriteLine("Exit: " + FormatAllReturn());
                    //Console.WriteLine($"{finished} {outW} {errW}");
                    if (all) {
                        DisposeObjects();
                        return true;
                    }

                }
            }

            if (process != null) {
                //
                // just kill everything
                //

                if (process.HasExited == false) {
                    process.Kill();

                }

                DisposeObjects();

            }

            return false;
        }

        ILog m_Logger = LogManager.GetLogger(typeof(SingleSessionSshClient));

        protected override void Connect() {
            using (var tr = new FuncTrace()) {
                if (process != null && process.HasExited == false)
                    return;

                Disconnect(true);


                outCollector = new List<string>();

                process = new Process();

                process.StartInfo.FileName = base.SshclientExeToUse;
                string keyfile;
                if(base.KeyFilePath != null) {
                    string delim;
                    if(System.OperatingSystem.IsWindows()) {
                        delim = "\"";
                    } else {
                        delim = "'";
                    }

                    keyfile = " -i " + delim + base.KeyFilePath + delim;
                } else {
                    keyfile = "";
                }

                process.StartInfo.Arguments = this.UserName + "@" + this.ServerName + keyfile + " -oStrictHostKeyChecking=no -t -t";

                //process.StartInfo.FileName = @"C:\Users\flori\source\repos\SshBashWrapper\StupidClint\bin\Debug\net6.0\StupidClint.exe";

                //process.StartInfo.FileName = "cmd";

                process.StartInfo.UseShellExecute = false;
                process.StartInfo.RedirectStandardOutput = true;
                process.StartInfo.RedirectStandardError = true;
                process.StartInfo.RedirectStandardInput = true;

                if (outputWaitHandle != null)
                    outputWaitHandle.Dispose();
                if (errorWaitHandle != null)
                    errorWaitHandle.Dispose();

                outputWaitHandle = new AutoResetEvent(false);
                errorWaitHandle = new AutoResetEvent(false);

                process.OutputDataReceived += OnOutputDataReceived;
                process.ErrorDataReceived += OnErrorDataReceived;


                process.Start();
                process.BeginOutputReadLine();
                process.BeginErrorReadLine();
                processActiveTask = process.WaitForExitAsync();

                process.StandardInput.Write("PS1=\"$PS1\\n\"" + "\n");
                process.StandardInput.Flush();

                string sshLogin;
                try {
                    sshLogin = WaitForPrompt();
                } catch (Exception) {
                    Disconnect(true);
                    throw;
                }

                tr.Info("SSH LOGIN: " + sshLogin);
            }
        }


        string WaitForPrompt() {
            string allOut;
            
            const int sleep = 100;
            int iMax = timeout * 10 / sleep;

            for (int i = 0; i < iMax; i++) {
                if (IsAtPrompt()) {
                    try {
                        allOut = FormatAllReturn();
                    } catch (Exception exc) {
                        allOut = $"During formatting output: {exc}";
                    }
                    ClearOut();
                    return allOut; ;
                }
                Thread.Sleep(sleep);
            }

            try {
                allOut = FormatAllReturn();
            } catch (Exception exc) {
                allOut = $"During formatting output: {exc}";
            }

            ClearOut();

            throw new IOException($"timeout in ssh connection; so far, stdout+stderr :::{allOut}<<<<<<<<<");
        }


        public override void Dispose() {
            Disconnect(true);
        }


        /*
        internal void Qnd() {

            const int timeout = 1000 * 3600;

            //var ret = MsHPC2012Client.ExecuteProcess("job", "", timeout);
            //Console.WriteLine("retO: " + ret.stdOut);
            //Console.WriteLine("retE: " + ret.stdErr);



            using (Process process = new Process()) {
                StringBuilder output = new StringBuilder();
                StringBuilder error = new StringBuilder();
                

                


                process.StartInfo.FileName = base.SshclientExeToUse;
                process.StartInfo.Arguments = m_usrname + "@" + m_srvrname + " -oStrictHostKeyChecking=no -t -t";

                //process.StartInfo.FileName = @"C:\Users\flori\source\repos\SshBashWrapper\StupidClint\bin\Debug\net6.0\StupidClint.exe";

                //process.StartInfo.FileName = "cmd";

                process.StartInfo.UseShellExecute = false;
                process.StartInfo.RedirectStandardOutput = true;
                process.StartInfo.RedirectStandardError = true;
                process.StartInfo.RedirectStandardInput = true;

                using (AutoResetEvent outputWaitHandle = new AutoResetEvent(false))
                using (AutoResetEvent errorWaitHandle = new AutoResetEvent(false)) {
                    process.OutputDataReceived += (sender, e) =>
                    {
                        if (e.Data == null) {
                            outputWaitHandle.Set();
                        } else {
                            string line = Regex.Replace(e.Data, @"\e\[(\d+;)*(\d+)?[ABCDHJKfmsu]", "");

                            output.AppendLine(line);
                            append(line);
                            //Console.WriteLine("rcv o: " + e.Data);
                            
                            for (int i = 0; i < e.Data.Length; i++) {
                                char ci = e.Data[i];
                                int cii = ci;
                                if (cii <= 31) {
                                    Console.BackgroundColor = ConsoleColor.Red;
                                    Console.ForegroundColor = ConsoleColor.Yellow;
                                    Console.Write($"&{cii}");

                                    Console.ResetColor();
                                } else if (char.IsWhiteSpace(ci)) {
                                    Console.BackgroundColor = ConsoleColor.Green;
                                    Console.ForegroundColor = ConsoleColor.Yellow;
                                    Console.Write(ci);
                                    Console.ResetColor();
                                } else {

                                    Console.Write(ci);

                                }
                            }
                            Console.WriteLine();
                            
                        }
                    };
                    process.ErrorDataReceived += (sender, e) =>
                    {
                        if (e.Data == null) {
                            errorWaitHandle.Set();
                        } else {
                            string line = Regex.Replace(e.Data, @"\e\[(\d+;)*(\d+)?[ABCDHJKfmsu]", "");
                            error.AppendLine(line);
                            append(line);
                            
                            //Console.Write("rcv E: ");
                            for (int i = 0; i < e.Data.Length; i++) {
                                char ci = e.Data[i];
                                int cii = ci;
                                if (cii <= 31) {
                                    Console.BackgroundColor = ConsoleColor.Red;
                                    Console.ForegroundColor = ConsoleColor.Yellow;
                                    Console.Write($"&{cii}");

                                    Console.ResetColor();
                                } else {
                                    Console.Write(ci);

                                }
                            }
                            Console.WriteLine();
                            
                        }
                    };


                    process.Start();
                    process.BeginOutputReadLine();
                    process.BeginErrorReadLine();
                    var task = process.WaitForExitAsync();

                    void w_cmd(string c) {
                        Console.WriteLine(" -----------------  " + c);
                        process.StandardInput.Write(c + "\n");
                        process.StandardInput.Flush();
                        Console.WriteLine(" *********************");
                        while (!IsAtPrompt()) {
                            Thread.Sleep(100);
                        }
                        writenFlush();
                    }


                    //w_cmd($"ssh {m_usrname}@{m_srvrname} -oStrictHostKeyChecking=no");
                    //w_cmd("csreport");
                    ////w_cmd("dir");
                    //w_cmd("exit");

                    process.StandardInput.Write("PS1=\"$PS1\\n\"" + "\n");
                    process.StandardInput.Flush();

                    Console.WriteLine("Waiting.....");
                    Thread.Sleep(20 * 1000);
                    writenFlush();
                    Console.WriteLine("Waiting.....");
                    Thread.Sleep(10 * 1000);

                    //w_cmd("PS1=\"$PS1\\n\"");


                    w_cmd("cd mono");
                    w_cmd("pwd");


                    //w_cmd("ls -l");
                    Thread.Sleep(10000);


                    Console.WriteLine("Schtopping...");

                    int i = 0;
                    while (!task.IsCompleted) {
                        Console.WriteLine("Waiting...");
                        Thread.Sleep(1000);
                        i++;
                        if (i == 1) {
                            process.StandardInput.Write("exit\n");
                            process.StandardInput.Flush();

                        }
                    }
                    Console.WriteLine("finished");


                    if (process.WaitForExit(timeout) //) { 
                        && outputWaitHandle.WaitOne(timeout) && errorWaitHandle.WaitOne(timeout)) {

                        Console.WriteLine("exit code = " + process.ExitCode);

                        string stdout = output.ToString();
                        string stderr = error.ToString();

                        Console.WriteLine("Stdout from ssh   :::::" + stdout + "<<<<<");
                        Console.WriteLine("Stderr from ssh   :::::" + stderr + "<<<<<");

                        // Process completed. Check process.ExitCode here.
                        return;
                    } else {
                        // Timed out.
                        throw new IOException("timeout waiting for ssh");
                    }
                }

            }
            



        }
        */
        /*
        private static void outputReader(Process process) {
            // Pass the standard output of the child to our standard output
            passThrough(process.StandardOutput.BaseStream, Console.OpenStandardOutput());
        }

        private static void passThrough(Stream instream, Stream outstream) {
            byte[] buffer = new byte[4096];
            while (true) {
                int len;
                while ((len = instream.Read(buffer, 0, buffer.Length)) > 0) {
                    outstream.Write(buffer, 0, len);
                    outstream.Flush();
                }
            }
        }
        */
    }


    class PrivateKeyFile {
        public PrivateKeyFile(string keypath) {
            m_path = keypath;
        }

        private string m_path;

        public string Path {
            get { return m_path; }
        }
    }

}
