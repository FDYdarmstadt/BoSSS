using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Text.RegularExpressions;

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
            var stdout = new List<string>();
            while(m_cmd.StandardOutput.Peek() > -1) {
                stdout.Add(m_cmd.StandardOutput.ReadLine());
            }

            var stderr = new List<string>();
            while(m_cmd.StandardError.Peek() > -1) {
                stderr.Add(m_cmd.StandardError.ReadLine());
            }
            std = String.Join("\n", stdout);
            err = String.Join("\n", stderr);
            //Console.WriteLine(std);
            //Console.WriteLine(err);
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

        public (string stdout, string stderr) RunCommand(string command) {
            m_cmd.Start();
            m_cmd.StandardInput.WriteLine("ssh " + m_usrname + "@" + m_srvrname + " \"" + command + "\"");
            //m_cmd.StandardInput.WriteLine(command);
            m_cmd.StandardInput.Flush();
            m_cmd.StandardInput.Close();
            m_cmd.WaitForExit();
            string std, err;
            ReadLines(out std, out err);
            return (std, err);
        }



    }
    public class PrivateKeyFile {
        public PrivateKeyFile(string keypath) {
            m_path = keypath;
        }

        private string m_path;

        public string Path {
            get { return m_path; }
        }
    }
}
