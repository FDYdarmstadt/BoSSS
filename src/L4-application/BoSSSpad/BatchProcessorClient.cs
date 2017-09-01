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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {
    /// <summary>
    /// Abstraction for all kind of batch systems (aka. job managers).
    /// </summary>
    abstract public class BatchProcessorClient {

        /// <summary>
        /// Base directory where the executables should be deployed.
        /// </summary>
        public string DeploymentBaseDirectory {
            get;
            protected set;
        }

        /// <summary>
        /// True, if it is required to deploy the runtime on the destination system,
        /// in order to execute a job.
        /// </summary>
        public bool DeployRuntime {
            get;
            set;
        }

        /// <summary>
        /// Returns the directory where the assemblies for <paramref name="myJob"/> 
        /// are deployed if <paramref name="myJob"/> is assigned to this batch processor.
        /// </summary>
        virtual public string GetNewDeploymentDir(Job myJob) {
            string ShortName = Path.GetFileNameWithoutExtension(myJob.EntryAssembly.Location);
            string DeployDir;
            int Counter = 0;
            do {
                string Suffix = Counter > 0 ? "-" + Counter : "";
                string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HH-mm-ss");
                DeployDir = Path.Combine(DeploymentBaseDirectory, ShortName + DateNtime + Suffix);
                Counter++;
            } while (Directory.Exists(DeployDir) == true);

            return DeployDir;
        }

        /// <summary>
        /// Submits the job to the batch system.
        /// </summary>
        /// <param name="myJob">Job to submit.</param>
        /// <returns>
        /// An optional identifier token (<see cref="Job.BatchProcessorIdentifierToken"/>).
        /// </returns>
        abstract public object Submit(Job myJob);

        /// <summary>
        /// Try to get some information about a job from the job manager.
        /// </summary>
        /// <param name="myJob"></param>
        /// <param name="SubmitCount">
        /// Number of jobs in the batch processor which match <paramref name="myJob"/>.
        /// </param>
        /// <param name="isRunning">
        /// True, if <paramref name="myJob"/> is currently running.
        /// </param>
        /// <param name="isFailed">
        /// True, if <paramref name="myJob"/> has terminated with an error state.
        /// </param>
        /// <param name="wasSuccessful">
        /// True, if <paramref name="myJob"/> has terminated without an error state.
        /// </param>
        /// <param name="DeployDir">
        /// If job is already known, the directory of the executable.
        /// </param>
        public abstract void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir);

        /// <summary>
        /// Path to standard output file, if present - otherwise null.
        /// </summary>
        public abstract string GetStdoutFile(Job myJob);

        /// <summary>
        /// Path to standard error file, if present - otherwise null.
        /// </summary>
        public abstract string GetStderrFile(Job myJob);

        /// <summary>
        /// Copies the executable files to the <see cref="DeploymentBaseDirectory"/>, 
        /// but does not submit the job.
        /// </summary>
        virtual public void DeployExecuteables(Job myJob, IEnumerable<Tuple<byte[],string>> AdditionalFiles) {

            Console.WriteLine("Deploying executables and additional files ...");

            // Collect files
            List<string> files = new List<string>();
            {
                //string SystemPath = Path.GetDirectoryName(typeof(object).Assembly.Location);
                string MainAssemblyDir = Path.GetDirectoryName(myJob.EntryAssembly.Location);
                foreach (var a in myJob.AllDependentAssemblies) {
                    if(Path.GetDirectoryName(a.Location).Equals(MainAssemblyDir)
                        || Path.GetFileName(a.Location).StartsWith("BoSSS")
                        || Path.GetFileName(a.Location).StartsWith("ilPSP")
                        || !a.GlobalAssemblyCache) { 
                        files.Add(a.Location);
                    }
                }
            }

            // create deployment directory.
            string DeployDir = myJob.DeploymentDirectory;
            if (!Directory.Exists(DeployDir))
                Directory.CreateDirectory(DeployDir);

            Console.WriteLine("Deployment directory: " + DeployDir);
            string OriginDir = null;
                        
            // copy files
            foreach(var fOrg in files) {
                string fNmn = Path.GetFileName(fOrg);
                if (fNmn.Equals("mscorlib.dll"))
                    throw new ApplicationException("internal error - something went wrong during filtering.");

                string fTarget = Path.Combine(DeployDir, fNmn);
                if(File.Exists(fTarget)) {
                    throw new IOException("File '" + fTarget + "' already exists - wont over write.");
                }

                File.Copy(fOrg, fTarget);
                if(OriginDir == null || !OriginDir.Equals(Path.GetDirectoryName(fOrg))) {
                    if (OriginDir != null)
                        Console.WriteLine();
                    OriginDir = Path.GetDirectoryName(fOrg);
                    Console.WriteLine("Source directory: " + OriginDir);
                    Console.Write("   copied: ");
                }

                Console.Write(Path.GetFileName(fOrg) + " ");
            }

            // additional files
            if(AdditionalFiles != null) {
                foreach(var t in AdditionalFiles) {
                    string fTarget = Path.Combine(DeployDir, t.Item2);
                    if (File.Exists(fTarget)) {
                        throw new IOException("File '" + fTarget + "' already exists - wont over write.");
                    }
                    File.WriteAllBytes(fTarget, t.Item1);
                    Console.WriteLine("   writing file: " + t.Item2);
                }
            }
                        
            //
            if (DeployRuntime) {
                string BosssInstall =  BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                string BosssBinNative = Path.Combine(BosssInstall, "bin", Path.Combine("native", "win"));
                CopyDirectoryRec(BosssBinNative, DeployDir, "amd64");
                Console.WriteLine("   copied 'amd64' runtime.");
            }
            
            // finally
            Console.WriteLine("deployment finished.");
        }

        /// <summary>
        /// (tries to) do a recursive copy of a directory
        /// </summary>
        /// <param name="srcDir"></param>
        /// <param name="dstDir"></param>
        /// <param name="filter">search pattern/filter</param>
        static void CopyDirectoryRec(string srcDir, string dstDir, string filter) {
            string[] srcFiles = Directory.GetFiles(srcDir);
            
            foreach (string srcFile in srcFiles) {
                TryCopy(srcFile, Path.Combine(dstDir, Path.GetFileName(srcFile)));
            }

            string[] subDirs;
            if (filter == null)
                subDirs = Directory.GetDirectories(srcDir);
            else
                subDirs = Directory.GetDirectories(srcDir, filter);
            foreach (string srcAbsDir in subDirs) {
                string srcRelDir = Path.GetFileName(srcAbsDir);
                string dstSubDir = Path.Combine(dstDir, srcRelDir);
                if (!Directory.Exists(dstSubDir))
                    Directory.CreateDirectory(dstSubDir);
                CopyDirectoryRec(srcAbsDir, dstSubDir, null);
            }
        }

        /// <summary>
        /// Utility function which tries to copy a file from
        /// <paramref name="sourceFileName"/> to
        /// <paramref name="destFileName"/> overwriting existing files if
        /// required. Issues a warning (but proceeds as normal) if the copy
        /// process fails.
        /// </summary>
        /// <param name="sourceFileName">
        /// The path to the file to be copied
        /// </param>
        /// <param name="destFileName">The path to the destination</param>
        static void TryCopy(string sourceFileName, string destFileName) {
            try {
                File.Copy(sourceFileName, destFileName, true);
            } catch (Exception e) {
                Console.WriteLine("WARNING: Unable to copy to: '"
                    + destFileName + "': " + e.GetType().Name + " says:'" + e.Message + "'");
            }
        }
    }
}

