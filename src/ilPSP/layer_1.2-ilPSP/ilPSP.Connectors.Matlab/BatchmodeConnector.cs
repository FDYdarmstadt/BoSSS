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

using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;

namespace ilPSP.Connectors.Matlab {

    /// <summary>
    /// provides a bride to the MATLAB system.
    /// </summary>
    /// <remarks>
    /// There is a .NET binding provided by MATLAB. However, this primitive one
    /// has no external dependencies and work with octave too, even under Cygwin.
    /// note that only MPI process rank 0 actually creates a MATLAB instance.
    /// </remarks>
    public class BatchmodeConnector : IDisposable {

        /// <summary>
        /// to distinct between MATLAB and octave.
        /// </summary>
        public enum Flavor {

            /// <summary>
            /// the clone.
            /// </summary>
            Octave_cygwin,

            /// <summary>
            /// the clone.
            /// </summary>
            Octave,

            /// <summary>
            /// the original.
            /// </summary>
            Matlab
        }

        /// <summary>
        /// the rank of the current process in the MPI-WORLD communicator
        /// </summary>
        int Rank {
            get {
                var comm = csMPI.Raw._COMM.WORLD;
                int rnk;
                csMPI.Raw.Comm_Rank(comm, out rnk);
                return rnk;
            }
        }


        private string get_program_path(string pname) {
            char[] chArr = new char[] { Path.PathSeparator };
            string path = System.Environment.GetEnvironmentVariable("PATH");
            if (path == null) {
                return null;
            } else {
                IEnumerable<string> pathDirs = path.Split(chArr, StringSplitOptions.RemoveEmptyEntries);

                foreach (string s4 in pathDirs) {
                    try {
                        string matlab_exe = Path.Combine(s4, pname);
                        if (File.Exists(matlab_exe))
                            return matlab_exe;
                    } catch (Exception) {
                    }
                }
            }
            return null;
        }

        static Flavor s_Flav;

        /// <summary>
        ///  octave or MATLAB
        /// </summary>
        static public Flavor Flav {
            get {
                return s_Flav;
            }
            set {
                s_Flav = value;
            }
        }

        /// <summary>
        /// Where to find the executable on the current system.
        /// If NULL, the standard installation path is assumed.
        /// In the case of Cygwin/octave, the path to the Cygwin bash.exe;
        /// </summary>
        static public string MatlabExecuteable {
            get;
            set;
        }

        /// <summary>
        /// Static ctor.
        /// </summary>
        static BatchmodeConnector() {
            Flav = Flavor.Matlab;
            //Flav = Flavor.Octave;
            MatlabExecuteable = null;  //"D:\\cygwin64\\bin\\bash.exe";

            try {
                TempDirMutex = new Mutex(false, "BoSSSbatchmodeconnector_IOmutex-" + System.Environment.UserName); // temp dir path is separate for each user => exclusivity within a user is sufficient
            } catch(Exception ee) {
                Console.Error.WriteLine("BatchmodeConnector, static ctor:" + ee);
                TempDirMutex = null; 
            }

            random = new Random();
        }

        /// <summary>
        /// Inter-process synchronization of file IO
        /// </summary>
        static Mutex TempDirMutex;

        /// <summary>
        /// a random timer, initialized at app startup, to obtain random wait times on IO problems
        /// </summary>
        static Random random;

        /// <summary>
        /// creates a new instance of the MATLAB connector.
        /// </summary>
        /// <param name="WorkingPath">
        /// working directory of the MATLAB instance;
        /// if NULL, a temporary directory is created.
        /// </param>
        public BatchmodeConnector(string WorkingPath = null) {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            this.m_Flav = Flav;

            // create/check working path
            // =========================

            if (Rank == 0) {
                if (WorkingPath == null) {
                    try {
                        if(TempDirMutex != null) // in the batch environment, Mutex sometimes throws UnauthorizedAccessException;
                            //                      then, we have no mutex and must rely on something else
                            //                      regarding IO, nothing seems reliable.
                            TempDirMutex.WaitOne();

                        bool Exists = false;
                        do {
                            var tempPath = Path.GetTempPath();
                            var tempDir = Guid.NewGuid().ToString(); // GUIDs should be unlikely to create collisions
                            WorkingDirectory = new DirectoryInfo(Path.Combine(tempPath, tempDir));
                            try {
                                Exists = WorkingDirectory.Exists;
                                if (!Exists) {
                                    WorkingDirectory.Create();
                                    DelWorkingDir = true;
                                }
                            } catch (Exception) {
                                // some unfortunate IO collision
                                Exists = false;
                                Thread.Sleep(random.Next(1000, 10000)); // wait somewhat between 1 and 10 seconds
                            }
                        } while(Exists == true);
                    
                    
                    
                    } finally {
                        if(TempDirMutex != null)
                            TempDirMutex.ReleaseMutex();
                    }
                } else {
                    WorkingDirectory = new DirectoryInfo(WorkingPath);
                    if (!WorkingDirectory.Exists)
                        throw new ArgumentException("Given working directory is inexistent.");
                }

            }

            var wd_FullName = MPIExtensions.MPIBroadcast(this.WorkingDirectory?.FullName, 0, csMPI.Raw._COMM.WORLD);
            this.WorkingDirectory = new DirectoryInfo(wd_FullName);

            // more checks
            // ===========
            if (MatlabExecuteable != null) {
                if (!File.Exists(MatlabExecuteable))
                    throw new ArgumentException("Unable to find file '" + MatlabExecuteable + "' on this system.");
            }


            // setup MATLAB process
            // ====================

            if (Rank == 0) {


                psi = new ProcessStartInfo();
                psi.WorkingDirectory = WorkingDirectory.FullName;
                psi.UseShellExecute = false;
                //psi.RedirectStandardOutput = true;
                //psi.RedirectStandardError = true;
                //psi.RedirectStandardInput = true;


                PlatformID CurrentSys = System.Environment.OSVersion.Platform;
                switch (CurrentSys) {
                    case PlatformID.Win32NT:
                    case PlatformID.Win32S:
                    case PlatformID.Win32Windows: {
                            if (m_Flav == Flavor.Matlab) {
                                if (MatlabExecuteable == null) {
                                    MatlabExecuteable = get_program_path("matlab.exe");
                                    if (MatlabExecuteable == null)
                                        throw new ApplicationException("Unable to find 'matlab.exe' in your PATH environment; please provide path to 'matlab.exe'.");
                                }

                                psi.FileName = MatlabExecuteable;
                                psi.Arguments = "-nosplash -nodesktop -minimize -wait -r " + CMDFILE + " -logfile " + LOGFILE;
                            } else if (m_Flav == Flavor.Octave) {
                                if (MatlabExecuteable == null) {
                                    MatlabExecuteable = get_program_path("octave-cli.exe");
                                    if (MatlabExecuteable == null)
                                        throw new ApplicationException("Unable to find 'octave-cli.exe' in your PATH environment; please provide path to 'octave-cli.exe'.");
                                }

                                psi.FileName = MatlabExecuteable;
                                psi.Arguments = " --no-gui " + CMDFILE + ".m > " + LOGFILE;
                                //psi.UseShellExecute = false;
                            } else if (m_Flav == Flavor.Octave_cygwin) {
                                this.Cygwin = true;

                                if (MatlabExecuteable == null) {
                                    if (File.Exists("c:\\cygwin64\\bin\\octave")) {
                                        psi.FileName = "c:\\cygwin64\\bin\\bash.exe";
                                    } else if (File.Exists("c:\\cygwin\\bin\\octave")) {
                                        psi.FileName = "c:\\cygwin\\bin\\bash.exe";
                                    } else {
                                        throw new NotSupportedException("Cygwin/Octave are expected to be in the default path: C:\\cygwin or c:\\cygwin64");
                                    }
                                } else {
                                    //throw new NotSupportedException("Cygwin/Octave are expected to be in the default path: C:\\cygwin or c:\\cygwin64");
                                    if (!MatlabExecuteable.EndsWith("bash.exe")) {
                                        throw new NotSupportedException("For Cygwin/Octave, the 'MatlabExecuteable' is expected to point to 'bash.exe'.");
                                    }

                                    psi.FileName = MatlabExecuteable;
                                }
                                psi.Arguments = "--login -c \"cd " + TranslatePath(WorkingDirectory.FullName) + " "
                                    + "&& octave --no-gui " + CMDFILE + ".m" + " > " + LOGFILE + "  \"";
                                //+ "pwd && ls - l && pwd";



                            } else {
                                throw new NotImplementedException();
                            }

                            break;
                        }



                    case PlatformID.Unix:
                        // if (this.m_Flav != Flavor.Octave){
                        //     throw new NotImplementedException("Use Octave instead");
                        // }
                        if (m_Flav == Flavor.Octave) {
                            if (MatlabExecuteable == null) {
                                MatlabExecuteable = get_program_path("octave-cli");
                                if (MatlabExecuteable == null)
                                    throw new ApplicationException("Unable to find 'octave-cli' in your PATH environment");
                            }

                            psi.FileName = MatlabExecuteable;
                            psi.Arguments = " --no-gui " + CMDFILE + ".m > " + LOGFILE;
                            //psi.UseShellExecute = false;

                        } else if (m_Flav == Flavor.Matlab) {
                            if (MatlabExecuteable == null) {
                                MatlabExecuteable = get_program_path("matlab");
                                if (MatlabExecuteable == null)
                                    throw new ApplicationException("Unable to find 'matlab' in your PATH environment");
                            }

                            psi.FileName = MatlabExecuteable;
                            psi.Arguments = "-nosplash -nodesktop -minimize -wait -r " + CMDFILE + " -logfile " + LOGFILE;
                        } else { 
                            throw new NotImplementedException();
                        }
                        break;
                    case PlatformID.MacOSX: {
                        throw new NotImplementedException("will implement on request");
                    }

                    default:
                        throw new NotSupportedException("unable to use MATLAB on " + CurrentSys.ToString());
                }


                CreatedFiles.Add(Path.Combine(WorkingDirectory.FullName, LOGFILE));


                var ScriptsToWrite = new List<Tuple<string, string>>();
                ScriptsToWrite.Add(new Tuple<string, string>("ReadMsr.m", Resource1.ReadMsr));
                ScriptsToWrite.Add(new Tuple<string, string>("SaveVoronoi.m", Resource1.SaveVoronoi));

                foreach (var t in ScriptsToWrite) {
                    string name = t.Item1;
                    string script = t.Item2;
                    var rmPath = Path.Combine(WorkingDirectory.FullName, name);
                    CreatedFiles.Add(rmPath);
                    File.WriteAllText(rmPath, script);
                }

            }

            // create command file
            // ===================
            if (Rank == 0) {
                var p = Path.Combine(WorkingDirectory.FullName, CMDFILE + ".m");
                CommandFile = new StreamWriter(p, false);
                CreatedFiles.Add(p);

            }
        }

        bool Cygwin;

        TextWriter CommandFile;

        /// <summary>
        /// MATLAB or octave?
        /// </summary>
        public Flavor m_Flav {
            get;
            private set;
        }

        /// <summary>
        /// the working directory of the MATLAB process
        /// </summary>
        public DirectoryInfo WorkingDirectory {
            get;
            private set;
        }

        /// <summary>
        /// List of all Files that should be deleted during disposing
        /// </summary>
        private List<string> CreatedFiles = new List<string>();

        /// <summary>
        /// delete working directory after finishing!
        /// </summary>
        private bool DelWorkingDir;

        private string TranslatePath(string p) {
            if (Cygwin) {
                if (p[1] != ':')
                    throw new ArgumentException();

                p = "/cygdrive/" + p.Substring(0, 1).ToLower() + "/" + p.Substring(3);
                p = p.Replace('\\', '/');
                return p;

            } else {
                return p;
            }
        }

        /// <summary>
        /// Transfers a sparse matrix to MATLAB.
        /// </summary>
        /// <param name="M">
        /// Matrix to transfer;
        /// If the matrix lives only on a sub-communicator of <see cref="IMPI_CommConstants.WORLD"/>, 
        /// this should be null on those processors which do not belong to the world-communicator.
        /// </param>
        /// <param name="MatlabName">the name which <paramref name="M"/> should have in the MATLAB session</param>
        public void PutSparseMatrix(IMutableMatrixEx M, string MatlabName) {
            var comm = csMPI.Raw._COMM.WORLD;
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No Data can be put after Execute() has been called..");

            string filepath;
            if (Rank == 0)
                filepath = Path.Combine(WorkingDirectory.FullName, MatlabName);
            else
                filepath = null;
            filepath = filepath.MPIBroadcast(0);
            if(M != null) // if M is on a sub-communicator of WORLD, this should be null
                M.SaveToTextFileSparse(filepath);
            
            if (Rank == 0)
                CreatedFiles.Add(filepath);

            Cmd(MatlabName + " = ReadMsr('" + TranslatePath(filepath) + "');");
        }

        /// <summary>
        /// Transfers a full matrix to Matlab
        /// </summary>
        /// <param name="M">matrix to transfer</param>
        /// <param name="MatlabName">
        /// The name which <paramref name="M"/> should have in the Matlab
        /// session
        /// </param>
        public void PutMatrix(IMatrix M, string MatlabName) {
            var comm = csMPI.Raw._COMM.WORLD;
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed) {
                throw new InvalidOperationException("No Data can be put after Execute() has been called..");
            }

            string filepath;
            if (Rank == 0)
                filepath = Path.Combine(WorkingDirectory.FullName, MatlabName);
            else
                filepath = null;
            M.SaveToTextFile(filepath, FileMode.Create);
            if (Rank == 0)
                CreatedFiles.Add(filepath);

            Cmd(MatlabName + " = dlmread('" + TranslatePath(filepath) + "');");
        }

        

        /// <summary>
        /// transfers a vector <paramref name="vec"/> to MATLAB;
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="vec">vector to transfer</param>
        /// <param name="MatlabName">the name which <paramref name="vec"/> should have in the MATLAB session</param>
        public void PutVector<T>(T vec, string MatlabName) where T : IEnumerable<double> {
            var comm = csMPI.Raw._COMM.WORLD;
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No commands can be added after Execute() has been called.");

            string filepath;
            if (Rank == 0)
                filepath = Path.Combine(WorkingDirectory.FullName, MatlabName);
            else
                filepath = null;
            vec.SaveToTextFile(filepath, comm);
            if (Rank == 0)
                CreatedFiles.Add(filepath);

            Cmd(MatlabName + " = dlmread('" + TranslatePath(filepath) + "');");
        }

        /// <summary>
        /// Transfers a vector without MPI-gathering <paramref name="vec"/> on rank 0, i.e. multiple vectors, one per MPI process.
        /// <paramref name="MatlabName"/>is extended with a suffix `_rank`, where `rank` is the respective MPI rank.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="vec"></param>
        /// <param name="MatlabName"></param>
        public void PutVectorPerMPIrank<T>(T vec, string MatlabName) where T : IEnumerable<double> {
            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No commands can be added after Execute() has been called.");

            string Mname = String.Concat(MatlabName + "_" + rank);
            string workdir;
            if (Rank == 0)
                workdir = WorkingDirectory.FullName;
            else
                workdir = null;
            workdir = workdir.MPIBroadcast(0);
            string filepath = Path.Combine(workdir, Mname);

            if (vec != null)
                VectorIO.SaveToTextFile(vec, filepath, csMPI.Raw._COMM.SELF);

            string[] Mnames = Mname.MPIGatherO(0);
            if (Rank == 0) {
                foreach (string Mn in Mnames) {
                    CreatedFiles.Add(Path.Combine(workdir, Mn));
                }
            }
            Mnames = Mnames.MPIBroadcast(0);
            foreach (string Mn in Mnames)
                Cmd(Mn + " = dlmread('" + TranslatePath(Path.Combine(workdir, Mn)) + "');");
        }


        /// <summary>
        /// Transfers a matrix without MPI-gathering <paramref name="M"/> on rank 0, i.e. multiple matrices, one per MPI process.
        /// <paramref name="MatlabName"/>is extended with a suffix `_rank`, where `rank` is the respective MPI rank.
        /// </summary>
        /// <param name="M"></param>
        /// <param name="MatlabName"></param>
        public void PutMatrixPerMPIrank(IMatrix M, string MatlabName) {
            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No commands can be added after Execute() has been called.");

            string Mname = String.Concat(MatlabName + "_" + rank);
            string workdir;
            if (Rank == 0)
                workdir = WorkingDirectory.FullName;
            else
                workdir = null;
            workdir = workdir.MPIBroadcast(0);
            string filepath = Path.Combine(workdir, Mname);

            if (M != null)
                M.SaveToTextFile(filepath, FileMode.Create);

            string[] Mnames = Mname.MPIGatherO(0);
            if (Rank == 0) {
                foreach (string Mn in Mnames) {
                    CreatedFiles.Add(Path.Combine(workdir, Mn));
                }
            }
            Mnames = Mnames.MPIBroadcast(0);
            foreach (string Mn in Mnames)
                Cmd(Mn + " = dlmread('" + TranslatePath(Path.Combine(workdir, Mn)) + "');");
        }


        /// <summary>
        /// Transfers a sparse matrix without MPI-gathering <paramref name="M"/> on rank 0, i.e. multiple matrices, one per MPI process.
        /// <paramref name="MatlabName"/>is extended with a suffix `_rank`, where `rank` is the respective MPI rank.
        /// </summary>
        /// <param name="M"></param>
        /// <param name="MatlabName"></param>
        public void PutSparseMatrixPerMPIrank(IMutableMatrixEx M, string MatlabName) {
            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No Data can be put after Execute() has been called..");
            if (comm == M.MPI_Comm || M == null)
                throw new NotSupportedException("has to be on sub comm of world and not empty");

            string Mname = String.Concat(MatlabName + "_" + rank);
            string workdir;
            if (Rank == 0)
                workdir = WorkingDirectory.FullName;
            else
                workdir = null;
            workdir = workdir.MPIBroadcast(0);
            string filepath = Path.Combine(workdir, Mname);

            if (M != null)
                M.SaveToTextFileSparse(filepath);

            string[] Mnames = Mname.MPIGatherO(0);
            if (Rank == 0) {
                foreach (string Mn in Mnames) {
                    CreatedFiles.Add(Path.Combine(workdir, Mn));
                }
            }
            Mnames = Mnames.MPIBroadcast(0);
            foreach (string Mn in Mnames)
                Cmd(Mn + " = ReadMsr('" + TranslatePath(Path.Combine(workdir, Mn)) + "');");
        }

        /// <summary>
        /// Imports a matrix from MATLAB.
        /// </summary>
        /// <param name="M">
        /// Output memory, can be null; in this case, a <see cref="MultidimensionalArray"/> will 
        /// be allocated which can be obtained from <see cref="OutputObjects"/> under the name <paramref name="MatlabName"/>.
        /// </param>
        /// <param name="MatlabName">Name of the matrix in MATLAB.</param>
        /// <remarks>
        /// Must be called before <see cref="Execute(bool)"/>, the size of the matrix must be known.
        /// 1. (Optional) Allocate output memory <paramref name="M"/>. E.g. if the size of the matrix is not known
        ///    in advance, use null. In this case, a <see cref="MultidimensionalArray"/> will be allocated.
        /// 2. Write all your MATLAB-commands using <see cref="Cmd(string)"/>, resp. <see cref="Cmd(string, object[])"/>.
        /// 3. Call <see cref="GetMatrix(IMatrix, string)"/>.
        /// 4. Call <see cref="Execute(bool)"/>.
        ///    After the execution of MATLAB statements in <see cref="Execute(bool)"/> is finished, 
        ///    <see cref="Execute(bool)"/> tries to load the matrix from a text file and stores it in <paramref name="M"/>.
        /// 5. Optionally, the result can be obtained from <see cref="OutputObjects"/>.
        /// </remarks>
        public void GetMatrix(IMatrix M, string MatlabName) {
            var comm = csMPI.Raw._COMM.WORLD;
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No commands can be added after Execute() has been called.");

            string filepath;
            if (Rank == 0)
                filepath = Path.Combine(WorkingDirectory.FullName, MatlabName);
            else
                filepath = null;

            if (Rank == 0)
                CreatedFiles.Add(filepath);
            if (M != null)
                m_OutputObjects.Add(MatlabName, M);
            else
                m_OutputObjects.Add(MatlabName, typeof(MultidimensionalArray));

            Cmd(string.Format("dlmwrite('{0}',{1},'precision',16);", TranslatePath(filepath), MatlabName));
        }

        /// <summary>
        /// 
        /// </summary>
        public void Cmd() {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Quick hack to return the 2nd output argument from output from `"[V, C] = voronoin(...);"`.
        /// </summary>
        public void GetStaggeredIntArray(int[][] A, string MatlabName) {
            var comm = csMPI.Raw._COMM.WORLD;
            ilPSP.MPICollectiveWatchDog.Watch(comm);
            if (Executed == true)
                throw new InvalidOperationException("No commands can be added after Execute() has been called.");
            if (A == null)
                throw new ArgumentNullException();

            string filepath;
            if (Rank == 0)
                filepath = Path.Combine(WorkingDirectory.FullName, MatlabName);
            else
                filepath = null;

            if (Rank == 0)
                CreatedFiles.Add(filepath);
            m_OutputObjects.Add(MatlabName, A);


            Cmd(string.Format("SaveVoronoi({1},'{0}');", TranslatePath(filepath), MatlabName));
        }




        Dictionary<string, object> m_OutputObjects = new Dictionary<string, object>();

        /// <summary>
        /// Dictionary of objects (matrices, vectors, etc.) which should be returned from the MATLAB script to 
        /// this application.
        /// - key: name of the object in MATLAB
        /// - value: matrix, vector, etc. where the return value is stored after <see cref="Execute(bool)"/> has been called.
        /// </summary>
        public IDictionary<string, object> OutputObjects {
            get {
                return new System.Collections.ObjectModel.ReadOnlyDictionary<string, object>(m_OutputObjects);
            }
        }


        /// <summary>
        /// imports a vector form MATLAB; 
        /// </summary>
        /// <param name="MatlabName">name of the vector, </param>
        public void GetVector(string MatlabName) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// executes some MATLAB command
        /// </summary>
        /// <param name="MatlabCommand">the MATLAB command</param>
        public void Cmd(string MatlabCommand) {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            if (Executed == true)
                throw new InvalidOperationException("No commands can be added after Execute has been called.");

            if (Rank == 0) {
                CommandFile.WriteLine(MatlabCommand);
            }
        }

        /// <summary>
        /// executes some MATLAB command
        /// </summary>
        /// <param name="MatlabCommand">the MATLAB command</param>
        /// <param name="formatParams"></param>
        public void Cmd(string MatlabCommand, params object[] formatParams) {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            if (Executed == true) {
                throw new InvalidOperationException("No commands can be added after Execute has been called.");
            }

            if (Rank == 0) {
                CommandFile.WriteLine(MatlabCommand, formatParams);
            }
        }


        const string LOGFILE = "matlab_connector_logfile.log";
        const string CMDFILE = "matlab_connector_commands";
        bool Executed = false;
        ProcessStartInfo psi;


        bool SuccessfulExe = false;

        /// <summary>
        /// executes all commands
        /// </summary>
        /// <returns>
        /// a reader to the standard output of the MATLAB process.
        /// </returns>
        public void Execute(bool PrintOutput = true) {
            SuccessfulExe = false;
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            if (Executed == true)
                throw new InvalidOperationException("Execute can be called only once.");

            int Rank = this.Rank;

            // run MATLAB
            // ==========
            Cmd("exit"); // be sure to exit!

            if (Rank == 0) {

                CommandFile.Flush();
                CommandFile.Close();

                //psi.RedirectStandardOutput = true;
                //psi.RedirectStandardError = true;
                //psi.RedirectStandardInput = true;

                var proc = Process.Start(psi);
                proc.WaitForExit();
            }


            // return 
            // ======


            if (Rank == 0) {
                var p = Path.Combine(WorkingDirectory.FullName, LOGFILE);

                if (PrintOutput && File.Exists(p)) {
                    //var stdout = new StreamReader(p);
                    //string line = stdout.ReadLine();
                    //while (line != null) {
                    //    Console.WriteLine(line);
                    //    line = stdout.ReadLine();
                    //}
                    //stdout.Dispose();
                    var text = File.ReadAllText(p);
                    Console.WriteLine(text);
                }

                CreatedFiles.Add(p);
            }

            foreach (string key in m_OutputObjects.Keys.ToArray()) {
                object outputObj = m_OutputObjects[key];

                string filepath;
                if (Rank == 0)
                    filepath = Path.Combine(WorkingDirectory.FullName, key);
                else
                    filepath = null;

                if (outputObj is IMatrix) {
                    // ++++++++++++++++++++++++++
                    // load pre-allocated matrix
                    // ++++++++++++++++++++++++++
                    IMatrix outputMtx = (IMatrix)outputObj;
                    if (Rank == 0)
                        outputMtx.LoadFromTextFile(filepath);


                    var _outputMtx = MPIExtensions.MPIBroadcast(outputMtx, 0, csMPI.Raw._COMM.WORLD);

                    if (Rank != 0) {
                        outputMtx.Clear();
                        outputMtx.Acc(1.0, _outputMtx);
                    }

                    if (!object.ReferenceEquals(outputObj, outputMtx)) {

                    }

                } else if (outputObj is Type && ((Type)outputObj) == typeof(MultidimensionalArray)) {
                    // ++++++++++++++++++++++++++
                    // load matrix which is NOT pre-allocated, unknown size
                    // ++++++++++++++++++++++++++

                    IMatrix outputMtx = null;
                    if (Rank == 0)
                        outputMtx = IMatrixExtensions.LoadFromTextFile(filepath);

                    var _outputMtx = MPIExtensions.MPIBroadcast(outputMtx, 0, csMPI.Raw._COMM.WORLD);
                    m_OutputObjects[key] = _outputMtx;

                } else if (outputObj is int[][]) {
                    int[][] outputStAry = (int[][])outputObj;

                    if (Rank == 0)
                        LoadStaggeredArray(outputStAry, filepath);

                    var _outputStAry = MPIExtensions.MPIBroadcast(outputStAry, 0, csMPI.Raw._COMM.WORLD);

                    if (Rank != 0) {
                        for (int i = 0; i < Math.Max(_outputStAry.Length, outputStAry.Length); i++) { // we use max to ensure an index-out-of-range if something is fishy
                            outputStAry[i] = _outputStAry[i];
                        }
                    }
                } else {
                    throw new NotImplementedException("output object type not implemented.");
                }
            }

            SuccessfulExe = true;
        }


        static void LoadStaggeredArray(int[][] StAry, string filepath) {
            using (var str = new StreamReader(filepath)) {
                int iLine = -1;

                for (string line = str.ReadLine(); line != null; line = str.ReadLine()) {
                    iLine++;

                    string[] NumbersStr = line.Split(default(char[]), StringSplitOptions.RemoveEmptyEntries);
                    int[] NumbersInt = new int[NumbersStr.Length];
                    for (int j = 0; j < NumbersStr.Length; j++) {
                        NumbersInt[j] = int.Parse(NumbersStr[j]);
                    }

                    StAry[iLine] = NumbersInt;
                }
            }
        }

        /// <summary>
        /// Kills the MATLAB process
        /// </summary>
        public void Dispose() {
            //using(new FuncTrace()) {
            { 
                if(SuccessfulExe) {
                    try {
                        if (TempDirMutex != null)
                            TempDirMutex.WaitOne();
                        foreach(var f in CreatedFiles) {
                            try {
                                File.Delete(f);
                            } catch(Exception e) {
                                Console.Error.WriteLine($"BatchModeConnector: {e.GetType()} during deletion of file {f}: {e.Message}");
                            }
                        }
                        if(DelWorkingDir) {
                            try {
                                Directory.Delete(WorkingDirectory.FullName, true);
                            } catch(Exception e) {
                                Console.Error.WriteLine($"BatchModeConnector: {e.GetType()} during deletion of directory {WorkingDirectory.FullName}: {e.Message}");
                            }
                        }
                    } finally {
                        if(TempDirMutex != null)
                            TempDirMutex.ReleaseMutex();
                    }
                } else {
                    // keeping files for diagnostic purposes
                }
            }
        }
    }
}
