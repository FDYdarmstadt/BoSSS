using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace bcl {

    /// <summary>
    /// Valid modes for storing the binaries
    /// </summary>
    public enum BinMode {

        /// <summary>
        /// Binaries are stored in the "Release" subfolder
        /// </summary>
        Release,

        /// <summary>
        /// Binaries are stored in the "Debug" subfolder
        /// </summary>
        Debug
    }

    /// <summary>
    /// Utility class to encapsulate information about the installation
    /// environment.
    /// </summary>
    class MyEnvironment {

        /// <summary>
        /// See <see cref="FindSystemPath"/>, <see cref="SearchGit"/> and
        /// <see cref="InitAndVerify"/>.
        /// </summary>
        public MyEnvironment() {
            FindSystemPath();
            SearchGit();
            InitAndVerify();
        }

        /// <summary>
        /// reads the PATH variable of the system
        /// </summary>
        private void FindSystemPath() {
            // get Path
            {
                char[] pathSep = new char[1];
                switch (Path.DirectorySeparatorChar) {
                    case '/': // Unix
                        pathSep[0] = ':';
                        break;
                    case '\\': // Win
                        pathSep[0] = ';';
                        break;
                    default:
                        throw new NotSupportedException("Unknown system; directory separator seems to be '" + Path.DirectorySeparatorChar + "' ?");
                }

                // try to be as robust as possible...
                IDictionary allEnvVars = Environment.GetEnvironmentVariables();
                foreach (string EnvVar in allEnvVars.Keys) {
                    if (EnvVar.Equals("path", StringComparison.InvariantCultureIgnoreCase)) {
                        string full = (string)allEnvVars[EnvVar];
                        string[] dirs = full.Split(pathSep, StringSplitOptions.RemoveEmptyEntries);

                        foreach (string d in dirs) {
                            try {
                                DirectoryInfo dirInfo = new DirectoryInfo(d);
                                if (dirInfo.Exists)
                                    _Path.Add(dirInfo);
                            } catch (Exception) {
                            }
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Looks for a git executable and stores its location
        /// </summary>
        private void SearchGit() {
            // search for git tool
            foreach (DirectoryInfo dir in _Path) {

                FileInfo[] files = dir.GetFiles();

                foreach (FileInfo file in files) {
                    if (file.Name == "git" || file.Name == "git.exe") {
                        // bingo

                        GitExe = file;
                        return;
                    }
                }
            }

            Console.WriteLine("Err: Could not locate git executable");
            //throw new EnviromentException("Could not locate git executable");
        }

        /// <summary>
        /// Recursively walks through all subdirectories of
        /// <see cref="BOSSS_ROOT"/> and adds every git repository to a list.
        /// </summary>
        /// <remarks>
        /// Private git repositories (i.e. repositories without an origin) will
        /// be omitted.
        /// </remarks>
        /// <returns>A list of git repositories</returns>
        public List<DirectoryInfo> GetGitRepositories() {
            List<DirectoryInfo> resultList = new List<DirectoryInfo>();
            GetGitRepositories(BOSSS_ROOT, resultList);

            // Remove private repositories
            List<DirectoryInfo> newList = new List<DirectoryInfo>();
            foreach (DirectoryInfo repository in resultList) {
                if (!GitRepositoryIsPrivate(repository)) {
                    newList.Add(repository);
                }
            }

            return newList;
        }

        /// <summary>
        /// Recursively walks through all subdirectories of
        /// <paramref name="currentDirectory"/> and adds every git repository
        /// to a list.
        /// </summary>
        /// <param name="currentDirectory">The directory to start from</param>
        /// <param name="resultList">The list storing the result</param>
        private void GetGitRepositories(
            DirectoryInfo currentDirectory, List<DirectoryInfo> resultList) {

            foreach (DirectoryInfo directory in currentDirectory.GetDirectories()) {
                if (directory.Name.Equals(".git", StringComparison.InvariantCultureIgnoreCase)) {
                    // Add $currentDirectory and break (since subdirectories of
                    // a repository should not contain more repositories)
                    resultList.Add(currentDirectory);
                    break;
                } else if (directory.Name.EndsWith(".git", StringComparison.InvariantCultureIgnoreCase)) {
                    // This is a bare clone of a repository. Add this
                    // _sub_directory but do not break since there still might
                    // be other relevant subdirectories.
                    resultList.Add(currentDirectory);
                } else {
                    // Do recursion
                    GetGitRepositories(directory, resultList);
                }
            }
        }

        /// <summary>
        /// Determines whether <paramref name="repository"/> contains a private
        /// git repository (i.e. a repository without an origin)
        /// </summary>
        /// <param name="repository">
        /// The directory containing the repository
        /// </param>
        /// <returns>
        /// True if the repository is private, false otherwise
        /// </returns>
        private bool GitRepositoryIsPrivate(DirectoryInfo repository) {
            ProcessStartInfo psi = new ProcessStartInfo();
            psi.FileName = GitExe.FullName;
            psi.Arguments = " config --get remote.origin.url";
            psi.UseShellExecute = false;
            psi.RedirectStandardOutput = true;
            psi.WorkingDirectory = repository.FullName;

            Process git = Process.Start(psi);
            git.WaitForExit();
            return (git.StandardOutput.ReadToEnd().Trim().Length == 0);
        }

        /// <summary>
        /// a heuristic, which verifies that some directory <paramref name="_d"/>
        /// contains BoSSS;
        /// </summary>
        /// <param name="_d">The directory in question</param>
        /// <returns>
        /// True if it _could_ be a BoSSS source code directory, false otherwise
        /// </returns>
        private bool VerifyBoSSSDir(DirectoryInfo _d) {

            var src = _d.GetDirectories("src").FirstOrDefault();
            if (src == null)
                return false;

            var pub = src.GetDirectories("public").FirstOrDefault();
            if (pub == null)
                return false;
            var DBE = src.GetDirectories("DBEv2").FirstOrDefault();
            if (DBE == null)
                return false;
            var ilPSP = src.GetDirectories("ilPSP").FirstOrDefault();
            if (ilPSP == null)
                return false;
            var libs = _d.GetDirectories("libs").FirstOrDefault();
            if (libs == null)
                return false;

            return true;
        }

        /// <summary>
        /// all directories in Path;
        /// contains olny directories that exist.
        /// </summary>
        public List<DirectoryInfo> _Path = new List<DirectoryInfo>();

        /// <summary>
        /// Git Executable
        /// </summary>
        public FileInfo GitExe;

        /// <summary>
        /// User home directory
        /// </summary>
        public DirectoryInfo HOME;

        /// <summary>
        /// path to the BoSSS root folder
        /// </summary>
        public DirectoryInfo BOSSS_ROOT;

        /// <summary>
        /// path to the directory which contains the native libraries,
        /// i.e. the 'amd64' and 'x86' subdirectories
        /// </summary>
        public DirectoryInfo BOSSS_BIN_NATIVE;

        /// <summary>
        /// path to the global binary folder
        /// </summary>
        public DirectoryInfo BOSSS_BIN;

        /// <summary>
        /// path to the documentation file folder ("doc");
        /// </summary>
        public DirectoryInfo BOSSS_DOC;

        /// <summary>
        /// path to the user configuration file folder ("$HOME/.BoSSS/etc");
        /// </summary>
        public DirectoryInfo BOSSS_ETC;

        /// <summary>
        /// path to the user folder ("$HOME/.BoSSS");
        /// </summary>
        public DirectoryInfo BOSSS_USER;

        /// <summary>
        /// path to the external libraries folder
        /// </summary>
        public DirectoryInfo BOSSS_LIBS;

        /// <summary>
        /// bosss L0-L4 repository: %BOSSS_ROOT%/src/public
        /// </summary>
        public DirectoryInfo BOSSS_SRC_PUBLIC;

        /// <summary>
        /// %BOSSS_ROOT/src
        /// </summary>
        public DirectoryInfo BOSSS_SRC;

        /// <summary>
        /// database explorer source code repository: %BOSSS_ROOT/src/DBE
        /// </summary>
        public DirectoryInfo BOSSS_SRC_DBE;

        /// <summary>
        /// Path to the database explorer application
        /// </summary>
        public DirectoryInfo BOSSS_BIN_DBE;

        /// <summary>
        /// Path to BoSSS installation directory.
        /// </summary>
        public DirectoryInfo BOSSS_INSTALL;

        /// <summary>
        /// prints some paths
        /// </summary>
        public void PrintInfo() {
            Console.WriteLine("\n\nFolders:");
            Console.WriteLine("Installation (binaries): " + BOSSS_INSTALL.FullName);
            Console.WriteLine("Source code:             " + ((BOSSS_ROOT != null) ? BOSSS_ROOT.FullName : "null"));
            Console.WriteLine("User configuration:      " + BOSSS_USER.FullName);
            Console.WriteLine("Native libraries:        " + BOSSS_BIN_NATIVE.FullName);
        }

        /// <summary>
        /// tries to initialize <see cref="BOSSS_ROOT"/>
        /// </summary>
        /// <param name="d"></param>
        void FindBoSSSRootRec(DirectoryInfo d) {
            if (VerifyBoSSSDir(d)) {
                BOSSS_ROOT = d;
                return;
            }

            DirectoryInfo _d = d.Parent;
            if (_d == null) {
                //Console.WriteLine("ERROR: unable to find BoSSS root directory - unvalid installation?");
                //throw new EnviromentException("ERROR: unable to find BoSSS root directory - unvalid installation?");
                return;
            }
            FindBoSSSRootRec(_d);
        }

        /// <summary>
        /// checks the BoSSS installation and imports environment variables;
        /// </summary>
        void InitAndVerify() {
            FindBoSSSRootRec(new DirectoryInfo(Directory.GetCurrentDirectory()));

            // Define ETC, BIN folder
            // ======================
            {
                if (BOSSS_ROOT != null)
                    BOSSS_DOC = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName, "doc"));

                if (BOSSS_ROOT != null)
                    BOSSS_SRC = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName, "src"));

                if (BOSSS_ROOT != null)
                    BOSSS_SRC_DBE = new DirectoryInfo(Path.Combine(BOSSS_SRC.FullName, "DBE"));

                if (BOSSS_ROOT != null)
                    BOSSS_SRC_PUBLIC = new DirectoryInfo(Path.Combine(BOSSS_SRC.FullName, "public"));

                BOSSS_INSTALL = GetBoSSSInstallDir();

                if (Path.DirectorySeparatorChar == '\\') {
                    // native libs currently only on windows
                    BOSSS_BIN_NATIVE = new DirectoryInfo(Path.Combine(BOSSS_INSTALL.FullName, Path.Combine("bin", Path.Combine("native", "win"))));
                }

                BOSSS_BIN_DBE = new DirectoryInfo(Path.Combine(BOSSS_INSTALL.FullName, Path.Combine("bin", "DBE")));

                BOSSS_BIN = new DirectoryInfo(Path.Combine(BOSSS_INSTALL.FullName, "bin"));

                if (BOSSS_ROOT != null)
                    BOSSS_LIBS = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName, "libs"));
            }

            // user config folder
            // ==================
            {
                if (Path.DirectorySeparatorChar == '\\') {
                    // Win -> search for USERPROFILE
                    string si1 = System.Environment.GetEnvironmentVariable("USERPROFILE");
                    if (si1 == null)
                        throw new ApplicationException("missing USERPROFILE environment variable");
                    HOME = new DirectoryInfo(si1);
                } else {
                    // Unix -> search for HOME
                    string si1 = System.Environment.GetEnvironmentVariable("HOME");
                    if (si1 == null)
                        throw new ApplicationException("missing HOME environment variable");
                    HOME = new DirectoryInfo(si1);
                }

                BOSSS_USER = new DirectoryInfo(Path.Combine(HOME.FullName, ".BoSSS"));
                BOSSS_ETC = new DirectoryInfo(Path.Combine(BOSSS_USER.FullName, "etc"));

                if (!BOSSS_USER.Exists) {
                    if (BOSSS_INSTALL.Exists) {
                        // create only if there is a full bosss installation on this machine
                        BOSSS_USER.Create();
                        BOSSS_ETC.Create();

                    }
                }
            }
        }

        /// <summary>
        /// searches for the User- or Machine-environment variable 'BOSSS_INSTALL'
        /// and verifies the existence of this directory.
        /// </summary>
        /// <returns></returns>
        static DirectoryInfo GetBoSSSInstallDir() {
            string si1 = System.Environment.GetEnvironmentVariable("BOSSS_INSTALL", EnvironmentVariableTarget.User);
            string si2 = System.Environment.GetEnvironmentVariable("BOSSS_INSTALL", EnvironmentVariableTarget.Machine);

            string si = null;
            if (si1 != null && si1.Length > 0) {
                si = si1;
            } else if (si2 != null && si2.Length > 0) {
                si = si2;
            }


            if (si != null && si.Length > 0) {
                bool success = false;
                try {
                    success = !Directory.Exists(si);
                } catch (Exception) {
                }

                if (success) {
                    Console.Error.WriteLine("WARNING: Environment variable 'BOSSS_INSTALL' is defined as '" + si + "', but directory seems to be non-existent.");
                }
            } else {
                Console.Error.WriteLine("WARNING: Unable to find environment variable 'BOSSS_INSTALL'; Hopefully the directory of the executable contains 'amd64' and 'x86', otherwise we're screwed.");
            }
            try {
                return new DirectoryInfo(si);
            } catch (Exception) {
                return null;
            }

        }
    }
}

