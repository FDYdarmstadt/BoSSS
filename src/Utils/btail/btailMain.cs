using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace btail {
    
    /// <summary>
    /// This is a loose clone of the 'tail' shell command, showing the rolling output of a text file;
    /// usually some stdout file of some job.
    /// </summary>
    public class TailMain {
        /// <summary>
        /// Tests if a string is empty or contains only whitespaces.
        /// </summary>
        public static bool IsEmptyOrWhite(string  s) {

            if (s == null)
                return true;

            int L = s.Length;
            for (int l = 0; l < L; l++) {
                if (!char.IsWhiteSpace(s[l]))
                    return false;
            }

            return true;
        }

        enum State {
            WaitingForFile = 1,

            Tailing = 2,

            Terminating = 3
        }

        /// <summary>
        /// Performs a encoding of the file names into base64, to avoid escape-char bullshit.
        /// </summary>
        public static void SetArgs(ProcessStartInfo psi, string StdOutFileName, string StdErrFileName) {
            // Environment variables would be nicer, but we want to use shell execution...
            //psi.EnvironmentVariables.Add("TAIL_STDOUT", StdOutFileName);
            //psi.EnvironmentVariables.Add("TAIL_STDERR", StdOutFileName);

            var StdOutFileBytes = Encoding.UTF8.GetBytes(StdOutFileName);
            var StdOutFileCode = Convert.ToBase64String(StdOutFileBytes);

            var StdErrFileBytes = Encoding.UTF8.GetBytes(StdErrFileName);
            var StderrFileCode = Convert.ToBase64String(StdErrFileBytes);

            psi.Arguments = StdOutFileCode + "  " + StderrFileCode;

        }


        static void Main(string[] args) {
            string StdOutFileCode = args[0];
            var StdOutBytes = Convert.FromBase64String(StdOutFileCode);
            string StdOutFileName = Encoding.UTF8.GetString(StdOutBytes);

            string StdErrFileCode = args[1];
            var StdErrBytes = Convert.FromBase64String(StdErrFileCode);
            string StdErrFileName = Encoding.UTF8.GetString(StdErrBytes);

            
            //string StdOutFileName = System.Environment.GetEnvironmentVariable("TAIL_STDOUT");
            //string StdErrFileName = System.Environment.GetEnvironmentVariable("TAIL_STDERR");

            if (IsEmptyOrWhite(StdOutFileName))
                throw new ApplicationException("Missing valid value for environment variable TAIL_STDOUT.");
            if (IsEmptyOrWhite(StdErrFileName))
                throw new ApplicationException("Missing valid value for environment variable TAIL_STDERR.");

            StreamReader StdOut = null;
            StreamReader StdErr = null;


            int i = 0;
            State s = State.WaitingForFile;
            while(s != State.Terminating) {
                i++;

                if (File.Exists(StdErrFileName) && StdErr == null) {
                    StdErr = new StreamReader(
                        new FileStream(StdErrFileName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                }


                if(s == State.WaitingForFile) {
                    

                    if(File.Exists(StdOutFileName)) {
                        StdOut = new StreamReader(
                            new FileStream(StdOutFileName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                        s = State.Tailing;
                    } else {
                        // keep waiting

                        if (i <= 1)
                            Console.WriteLine("Waiting for Job to launch...");
                    }


                } else if(s == State.Tailing) {

                    Debug.Assert(StdOut != null);
                    for(string L = StdOut.ReadLine(); L != null; L = StdOut.ReadLine()) {
                        Console.WriteLine(L);
                    }

                    if (StdErr != null) {
                        bool ch = false;
                        for (string L = StdErr.ReadLine(); L != null; L = StdErr.ReadLine()) {
                            ch = true;
                            Console.ForegroundColor = ConsoleColor.Yellow;
                            Console.BackgroundColor = ConsoleColor.Red;
                            Console.WriteLine(L);
                        }
                        if(ch)
                            Console.ResetColor();
                    }

                }


                Thread.Sleep(10000);
            }


            // termination
            // ===========

            {
                for (int j = 0; j < 10; j++)
                    Console.WriteLine();

                Console.WriteLine("Press any key to close this window.");
                Console.ReadKey(true);
                Console.WriteLine("Bye!");
                Thread.Sleep(5000);
            }
        }
    }
}
