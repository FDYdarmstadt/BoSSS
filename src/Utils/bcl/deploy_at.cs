using Renci.SshNet;
using Renci.SshNet.Common;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;

namespace bcl.deploy_at {

    /// <summary>
    /// Utility program to copy a program to a remote location while taking
    /// along all dependencies. Required assemblies will be taken from
    /// <see cref="MyEnvironment.BOSSS_BIN"/> and
    /// <see cref="MyEnvironment.BOSSS_BIN_NATIVE"/>.
    /// </summary>
    class Program : IProgram {

        /// <summary>
        /// The (BoSSS-)executable to be deployed.
        /// </summary>
        FileInfo assembly;

        /// <summary>
        /// The directory to which <see cref="assembly"/> should be deployed.
        /// </summary>
        DestinationDriver targetDirectory;

        /// <summary>
        /// 
        /// </summary>
        bool SurpressNative = false;

        /// <summary>
        /// Tries to find assembly <paramref name="dependencyName"/> in the
        /// global bin directory.
        /// </summary>
        /// <param name="dependencyName">
        /// The BoSSS assembly to be copied
        /// </param>
        /// <param name="searchPath">
        /// Directory in which we search for assemblies
        /// </param>
        private FileInfo GetBoSSSAssemblyLocation(AssemblyName dependencyName, DirectoryInfo searchPath) {
            string fileName = dependencyName.Name + ".dll";
            FileInfo[] files = searchPath.GetFiles(fileName);
            if (files.Length == 0) {
                fileName = dependencyName.Name + ".exe";
                files = searchPath.GetFiles(fileName);
            }

            if (files.Length > 1) {
                throw new EnviromentException("Found more than one file with name \"" + fileName + "\". Hell has probably just frozen over.");
            } else if (files.Length == 0) {
                // do recursion
                DirectoryInfo[] subDirectories = searchPath.GetDirectories();
                foreach (DirectoryInfo subDirectory in subDirectories) {
                    FileInfo file = GetBoSSSAssemblyLocation(dependencyName, subDirectory);
                    if (file != null) {
                        return file;
                    }
                }
                return null;
            } else {
                return files[0];
            }
        }
        

        
        /// <summary>
        /// Recursively walks through the dependencies of
        /// <paramref name="assemblyFile"/> and collects them in
        /// <paramref name="results"/>
        /// </summary>
        /// <param name="assemblyFile">The root of the recursion</param>
        /// <param name="results">A container for the determined files</param>
        private void CollectAllDependencies(FileInfo assemblyFile, List<string> results) {
            CollectAllDependenciesRecursive(assemblyFile, results, assemblyFile.Directory);
            //Assembly assembly = Assembly.LoadFile(assemblyFile.FullName);
            //List<Assembly> assis = new List<Assembly>();
            //GetAllAssembliesRecursive(assembly, assis, assemblyFile.Directory);
            //results.AddRange(assis.Select(a => a.CodeBase));
        }

        
        /*
        /// <summary>
        /// Recursive collection of all dependencies of some assembly.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="searchPath"></param>
        /// <param name="assiList">
        /// Output, list where all dependent assemblies are collected.
        /// </param>
        void GetAllAssembliesRecursive(Assembly a, List<Assembly> assiList, DirectoryInfo searchPath) {
            if (a.GlobalAssemblyCache)
                return;

            if (assiList.Contains(a))
                return;

            assiList.Add(a);

            foreach (var _b in a.GetReferencedAssemblies()) {
                Assembly b;
                if (_b.FullName.Contains("Hpc"))
                    Debugger.Break();

                try {
                    FileInfo fi_b = GetBoSSSAssemblyLocation(_b, searchPath);
                    b = Assembly.LoadFrom(fi_b.FullName);
                } catch (Exception e_outer) {
                    try {
                        b = Assembly.Load(_b);
                    } catch (Exception) {
                        throw e_outer;
                    }
                }


                if (b.GlobalAssemblyCache)
                    continue;
                GetAllAssembliesRecursive(b, assiList, searchPath);
            }
        }
        */
        
        /// <summary>
        /// Recursively walks through the dependencies of
        /// <paramref name="assemblyFile"/> and collects them in
        /// <paramref name="results"/>
        /// </summary>
        /// <param name="assemblyFile">The root of the recursion</param>
        /// <param name="results">A container for the determined files;
        /// relative paths to <paramref name="searchPath"/>
        /// </param>
        /// <param name="searchPath">
        /// Directory in which we search for assemblies
        /// </param>
        private void CollectAllDependenciesRecursive(FileInfo assemblyFile, List<string> results, DirectoryInfo searchPath) {
            string relativeAssemblyPath;
            string FullAssemblyPath = assemblyFile.FullName;
            FullAssemblyPath = Path.GetFullPath(FullAssemblyPath);

            string StartPath = searchPath.FullName;

            if (!FullAssemblyPath.StartsWith(StartPath)) {
                throw new ApplicationException("internal error");
            }

            relativeAssemblyPath = FullAssemblyPath.Substring(
                StartPath.Length + 1, // Remove directory seperator
                FullAssemblyPath.Length - StartPath.Length - 1);

            if (results.Contains(relativeAssemblyPath)) {
                // Already added, end of recursion
                return;
            } else {
                results.Add(relativeAssemblyPath);
            }

            Assembly assembly = Assembly.LoadFile(assemblyFile.FullName);

            // Collect all directly referenced assemblies
            foreach (AssemblyName dependencyName in assembly.GetReferencedAssemblies()) {
                FileInfo requiredFile = GetBoSSSAssemblyLocation(dependencyName, searchPath);
                if (requiredFile != null) {
                    CollectAllDependenciesRecursive(requiredFile, results, searchPath);
                }
            }
        }
        //*/
        #region IProgram Member

        /// <summary>
        /// Reads <see cref="assembly"/> and copies it with all its
        /// dependencies to <see cref="targetDirectory"/>.
        /// </summary>
        public void Execute() {
            
            // Load executable and find all dependencies
            List<string> requiredFiles = new List<string>();
            CollectAllDependencies(assembly, requiredFiles);

            // Copy the files (managed code)
            bool requiresNative = false;
            foreach (string requiredFile in requiredFiles) {
                string pathToOrigin = Path.Combine(assembly.DirectoryName, requiredFile);
                FileInfo origin = new FileInfo(pathToOrigin);

                if (!origin.Exists) {
                    throw new EnviromentException("Could not find required file " + origin.FullName);
                }

                if (requiredFile.ToLowerInvariant().Contains("bosss.platform.dll")) {
                    requiresNative = true;
                }

                try {
                    targetDirectory.PutFile(origin, requiredFile);
                } catch (Exception e) {
                    throw new EnviromentException(
                        "Could not write file \"" + requiredFile + "\" (Message: " + e.Message + ")");
                }
            }

            // Copy unmanaged code, if necessary
            if (requiresNative
                && (targetDirectory is LocalFSDriver)
                && Path.DirectorySeparatorChar == '\\'
                && !SurpressNative) {
                nativelib_inst.Program.CopyDirectoryRec(bcl.myEnv.BOSSS_BIN_NATIVE, ((LocalFSDriver)targetDirectory).target, "amd64");
            }

            targetDirectory.Dispose();
        }

        /// <summary>
        /// Tries to extract the path to the executable and the target
        /// directory from the first two entries of <paramref name="args"/>.
        /// </summary>
        /// <param name="args">
        /// An array containing one string representing the path to the
        /// executable and the deployment location.
        /// </param>
        public void DecodeArgs(string[] args) {
            if (args.Length > 3) {
                throw new UserInputException("Too many arguments");
            } else if (args.Length < 2) {
                throw new UserInputException("Missing option");
            } else {
                assembly = new FileInfo(args[0]);
                if (!assembly.Exists) {
                    throw new UserInputException(
                        "Could not find file " + assembly.FullName);
                } else if (!assembly.Name.EndsWith(".exe") && !assembly.Name.EndsWith(".dll")) {
                    throw new UserInputException(
                        "Given file is not an executable (.exe) or library (.dll) file");
                }

                targetDirectory = DestinationDriver.CreateDriver(args[1]);
            }
            if (args.Length == 3) {
                SurpressNative = !(args[2] == "0" || args[2].Equals("false", StringComparison.InvariantCultureIgnoreCase));
            }
        }

        /// <summary>
        /// See <see cref="IProgram.PrintUsage"/>.
        /// </summary>
        public void PrintUsage() {
            Console.WriteLine("Usage: bcl deploy-at $binaryFile $target {surpressNative}");
            Console.WriteLine("  Arguments:");
            Console.WriteLine("      binaryFile       An executable, e.g. NSE2b.exe,");
            Console.WriteLine("                       which should be deployed to e.g. a ");
            Console.WriteLine("                       compute cluster.");
            Console.WriteLine("      target           destination directory; either local");
            Console.WriteLine("                       or 'sftp://...'.");
            Console.WriteLine("      surpressNative   0 (default) or 1; if 1, the native");
            Console.WriteLine("                       libs, i.e. the 'amd64' and 'x86' ");
            Console.WriteLine("                       directories are not created (applies");
            Console.WriteLine("                       to MS windows only.");
            Console.WriteLine();
        }

        #endregion

        /// <summary>
        /// methods to copy files to the destination (either in file system or SFTP)
        /// </summary>
        abstract class DestinationDriver : IDisposable {

            public static DestinationDriver CreateDriver(string target) {
                if (target.ToLowerInvariant().StartsWith("sftp://")) {
                    return new SFTPDriver(target);
                } else {
                    return new LocalFSDriver(target);
                }
            }

            public abstract void PutFile(FileInfo Origin, params string[] RelativeDestinationPath);

            abstract protected string CombinePath(string[] RelPath);

            public abstract void Dispose();
        }

        /// <summary>
        /// used when copying to local file system
        /// </summary>
        class LocalFSDriver : DestinationDriver {

            public DirectoryInfo target;

            public LocalFSDriver(string path) {
                if (path.ToLowerInvariant().StartsWith("file://"))
                    path = path.Substring(7);

                target = new DirectoryInfo(path);
                if (!target.Exists)
                    target.Create();
            }

            override protected string CombinePath(string[] RelPath) {
                string ret = target.FullName;
                for (int i = 0; i < RelPath.Length; i++)
                    ret = Path.Combine(ret, RelPath[i]);
                return ret;
            }

            public override void PutFile(FileInfo Origin, params string[] RelativeDestinationPath) {
                Origin.CopyTo(CombinePath(RelativeDestinationPath), true);
            }

            public override void Dispose() {
            }
        }

        /// <summary>
        /// used when copying to a SFTP server
        /// </summary>
        class SFTPDriver : DestinationDriver {

            public SFTPDriver(string sftpPath) {
                if (!sftpPath.ToLowerInvariant().StartsWith("sftp://")) {
                    throw new ArgumentException("no sftp path");
                }

                string restPath = sftpPath.Substring(7); // remove sftp://
                string[] user_N_address_N_Path = restPath.Split('@', ':');

                if (user_N_address_N_Path.Length != 3) {
                    throw new ArgumentException("unable to parse sftp address: required format: 'sftp://user@adress:path-on-server';");
                }

                string server_password;
                {
                    Console.Write("Enter password for '" + user_N_address_N_Path[0] + "@" + user_N_address_N_Path[1] + "': ");
                    List<char> Password = new List<char>();
                    ConsoleKeyInfo key;
                    do {
                        key = Console.ReadKey(true);
                        Password.Add(key.KeyChar);
                    } while (key.Key != ConsoleKey.Enter);

                    Password.RemoveAt(Password.Count - 1);

                    server_password = new string(Password.ToArray());
                    Console.WriteLine();
                }

                m_Sftp = new SftpClient(user_N_address_N_Path[1], user_N_address_N_Path[0], server_password);
                m_Sftp.Connect();

                {
                    m_WorkingPath = user_N_address_N_Path[2].Split(new char[] { '/' }, StringSplitOptions.RemoveEmptyEntries);
                    for (int i = 0; i < m_WorkingPath.Length; i++) {
                        if (m_WorkingPath[i].Equals("~")) {
                            throw new ArgumentException(
                                "Your remote path should not start with '~'; just use a relative path"
                                    + "(which will automatically be relative to your home directory)",
                                nameof(sftpPath));
                        }
                    }

                    if (user_N_address_N_Path[2].StartsWith("/")) {
                        // absolute path
                        m_WorkingPath[0] = "/" + m_WorkingPath[0];
                    }
                }

                CreateDirectoryIfNotExistent();
            }

            /// <summary>
            /// the path where the application should be deployed (target directory),
            /// each directory is a separate entry of the array, i.e. the complete path reads as<br/>
            /// <see cref="m_WorkingPath"/>[0]/m_WorkingPath[1] ...
            /// </summary>
            string[] m_WorkingPath;

            /// <summary>
            /// combines the <see cref="m_WorkingPath"/>, up to directory level <paramref name="DirectoryLevel"/>,
            /// into one string
            /// </summary>
            /// <param name="DirectoryLevel"></param>
            /// <returns></returns>
            string PartOfWorkingPath(int DirectoryLevel) {
                StringWriter stw = new StringWriter();
                for (int i = 0; i < DirectoryLevel; i++) {
                    stw.Write(m_WorkingPath[i]);
                    if ((i + 1) < DirectoryLevel)
                        stw.Write("/");
                }
                return stw.ToString();
            }

            /// <summary>
            /// combines the <see cref="m_WorkingPath"/> into one string
            /// </summary>
            /// <returns></returns>
            string WorkingPath() {
                return PartOfWorkingPath(m_WorkingPath.Length);
            }

            /// <summary>
            /// appends the <paramref name="RelPath"/> to the current <see cref="m_WorkingPath"/>
            /// </summary>
            protected override string CombinePath(string[] RelPath) {
                StringWriter stw = new StringWriter();
                stw.Write(WorkingPath());
                if (RelPath.Length > 0)
                    stw.Write("/");
                for (int i = 0; i < RelPath.Length; i++) {
                    stw.Write(RelPath[i]);
                    if ((i + 1) < RelPath.Length)
                        stw.Write("/");
                }
                return stw.ToString();
            }

            /// <summary>
            /// checks if the <see cref="m_WorkingPath"/>-directory exists and creates it,
            /// if neccessary
            /// </summary>
            void CreateDirectoryIfNotExistent() {
                for (int i = 1; i <= m_WorkingPath.Length; i++) {
                    try {
                        string path = PartOfWorkingPath(i);
                        try {
                            var dummy = m_Sftp.ListDirectory(path);
                        } catch (SshException) {
                            m_Sftp.CreateDirectory(PartOfWorkingPath(i));
                            var dummy = m_Sftp.ListDirectory(path);
                        }
                    } catch (SshException sftpEx) {
                        throw new ApplicationException("SFTP error: " + sftpEx.Message, sftpEx);
                    }
                }
            }

            /// <summary>
            /// copies the file <paramref name="Origin"/>;
            /// </summary>
            /// <param name="Origin"></param>
            /// <param name="RelativeDestinationPath"></param>
            public override void PutFile(FileInfo Origin, params string[] RelativeDestinationPath) {
                using (var file = File.OpenRead(Origin.FullName)) {
                    m_Sftp.UploadFile(file, CombinePath(RelativeDestinationPath));
                }
            }

            SftpClient m_Sftp;

            /// <summary>
            /// closes the SFTP connection
            /// </summary>
            public override void Dispose() {
                if (m_Sftp != null) {
                    m_Sftp.Dispose();
                    m_Sftp = null;
                }
            }
        }
    }
}
