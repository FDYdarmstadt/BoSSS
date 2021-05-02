using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Text.RegularExpressions;
using System.IO;

namespace BoSSS.Application.BoSSSpad {

    class SshClient {
        /// <summary>
        /// ctor
        /// </summary>
        public SshClient(string ServerName, string Username, PrivateKeyFile pkf) {
            m_pkf = pkf;
            m_srvrname = ServerName;
            m_usrname = Username;
            PlatformID CurrentSys = System.Environment.OSVersion.Platform;
            string shell = "";
            switch(CurrentSys) {
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

        public SshClient(string ServerName, string Username, string pwd) {
            throw new NotSupportedException("TO DO");
        }

        public void Connect() {
            if(!TestConnection()) {
                string std, err;
                ReadLines(out std, out err);
                Console.WriteLine(err);
                Console.WriteLine("Connection could not be established");
            }
        }

        private void ReadLines(out string std, out string err) {
            using(TextWriter stdoutW = new StringWriter(), stderrW = new StringWriter()) {
                while(m_cmd.StandardOutput.Peek() > -1) {
                    stdoutW.WriteLine(m_cmd.StandardOutput.ReadLine());
                }

                while(m_cmd.StandardError.Peek() > -1) {
                    stderrW.WriteLine(m_cmd.StandardError.ReadLine());
                }

                std = stdoutW.ToString();
                err = stderrW.ToString();
            }
        }

        private bool TestConnection() {
            string search = RunCommand("ls").stdout;
            bool connected = search.Contains(m_usrname);
            return connected;
        }

        private Process m_cmd;
        private PrivateKeyFile m_pkf;
        private string m_usrname;
        private string m_srvrname;

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
            get { return TestConnection(); }
        }


        public string SubmitJob(string remotepath) {

            Connect();
            string sbatchCmd = "sbatch " + remotepath + "/batch.sh";
            var (resultString, err) = RunCommand(sbatchCmd);

            String SearchString = "Submitted batch job ";
            String jobId = Regex.Match(resultString, SearchString + "[0-9]*") // look for SearchString followed by a number (the Job ID)
                .ToString() // convert to string
                .Replace(SearchString, ""); // remove SearchString, leaving only the Job ID
            Console.WriteLine(jobId);
            return jobId;
        }

        public (string stdout, string stderr) RunCommand(string command, bool verbose = false) {
            if(verbose)
                Console.WriteLine("Starting test shell...");
            m_cmd.Start();
            if(verbose)
                Console.WriteLine("started.");

            string sshCmd = "ssh " + m_usrname + "@" + m_srvrname + " \"" + command + "\"";
            if(verbose)
                Console.WriteLine("Command: " + sshCmd);
            m_cmd.StandardInput.WriteLine(sshCmd);
            if(verbose)
                Console.WriteLine("command written; Waiting for completion...");
            //m_cmd.StandardInput.WriteLine(command);
            m_cmd.StandardInput.Flush();
            m_cmd.StandardInput.Close();
            m_cmd.WaitForExit();
            if(verbose)
                Console.WriteLine("external shell terminated; exit code is " + m_cmd.ExitCode);
            string std, err;
            ReadLines(out std, out err);
            if(verbose) {
                Console.WriteLine("stdout received:");
                Console.WriteLine(std);
                Console.WriteLine("----- (end of stdout)");
                Console.WriteLine("stderr received:");
                Console.WriteLine(err);
                Console.WriteLine("----- (end of stderr)");
            }

            return (std, err);
        }



    }
    
    
    /*
    /// <summary>
    /// Wrapper around system SSH where the connection is opened only once.
    /// </summary>
    class SshClient_exp {

        /// <summary>
        /// ctor
        /// </summary>
        public SshClient_exp(string ServerName, string Username, PrivateKeyFile pkf) {
            m_pkf = pkf;
            m_srvrname = ServerName;
            m_usrname = Username;
            PlatformID CurrentSys = System.Environment.OSVersion.Platform;
            string shell = "";
            switch(CurrentSys) {
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
                    FileName = @"C:\Program Files\Git\usr\bin\ssh.exe",
                    Arguments = " " + m_usrname + "@" + m_srvrname + " ",
                    CreateNoWindow = true,
                    UseShellExecute = false,
                    RedirectStandardInput = true,
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                }
            };
            m_cmd = cmd;

            m_cmd.Start();
            //m_cmd.StandardInput.WriteLine("ssh " + m_usrname + "@" + m_srvrname + " ");
            

            m_cmd.StandardInput.Flush();
            m_cmd.StandardInput.WriteLine("exit");
            m_cmd.StandardInput.Flush();
            m_cmd.StandardInput.Close();
            m_cmd.WaitForExit();
            
            Console.WriteLine("Opening SSH connection: ");
            ReadLines(out var stdout, out var stderr);
            Console.WriteLine(stdout);
            Console.WriteLine(stderr);
        }
        private void ReadLines(out string std, out string err) {
            using(TextWriter stdoutW = new StringWriter(), stderrW = new StringWriter()) {
                while(m_cmd.StandardOutput.Peek() > -1) {
                    stdoutW.WriteLine(m_cmd.StandardOutput.ReadLine());
                }

                while(m_cmd.StandardError.Peek() > -1) {
                    stderrW.WriteLine(m_cmd.StandardError.ReadLine());
                }

                std = stdoutW.ToString();
                err = stderrW.ToString();
            }
        }

        private bool TestConnection() {
            string search = RunCommand("ls").stdout;
            bool connected = search.Contains(m_usrname);
            return connected;
        }

        private Process m_cmd;
        private PrivateKeyFile m_pkf;
        private string m_usrname;
        private string m_srvrname;

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
            get { return TestConnection(); }
        }


        public string SubmitJob(string remotepath) {

            //Connect();

            string sbatchCmd = "sbatch " + remotepath + "/batch.sh";
            var (resultString, err) = RunCommand(sbatchCmd);
            
            String SearchString = "Submitted batch job ";
            String jobId = Regex.Match(resultString, SearchString + "[0-9]*") // look for SearchString followed by a number (the Job ID)
                .ToString() // convert to string
                .Replace(SearchString, ""); // remove SearchString, leaving only the Job ID
            Console.WriteLine(jobId);
            return jobId;
        }

        public (string stdout, string stderr) RunCommand(string command) {
            //m_cmd.Start();
            //m_cmd.StandardInput.WriteLine("ssh " + m_usrname + "@" + m_srvrname + " \"" + command + "\"");
            m_cmd.StandardInput.WriteLine(command);
            
            m_cmd.StandardInput.Flush();
            //m_cmd.StandardInput.Close();
            //m_cmd.WaitForExit();
            string std, err;
            ReadLines(out std, out err);
            return (std, err);
        }

        public void Final() {
            m_cmd.StandardInput.Close();
            m_cmd.WaitForExit();
            string std, err;
            ReadLines(out std, out err);
            Console.WriteLine(std);
            Console.WriteLine(err);

        }

    }

    //*/
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
