using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Collections;

namespace bcl {
    
    /// <summary>
    /// 
    /// </summary>
    public enum BinMode {
        Release,

        Debug
    }


    /// <summary>
    /// 
    /// </summary>
    class MyEnviroment {
        public MyEnviroment() {

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
                        throw new NotSupportedException("Unknown system; directory seperator seems to be '" + Path.DirectorySeparatorChar + "' ?");
                }

                // try to be as robust as possible...
                IDictionary allEnvVars = Environment.GetEnvironmentVariables();
                foreach (string EnvVar in allEnvVars.Keys) {
                    if (EnvVar.Equals("path", StringComparison.InvariantCultureIgnoreCase)) {
                        string full = (string)allEnvVars[EnvVar];
                        string[] dirs = full.Split(pathSep, StringSplitOptions.RemoveEmptyEntries);

                        foreach (string d in dirs) {
                            DirectoryInfo dirInfo = new DirectoryInfo(d);
                            if (dirInfo.Exists)
                                _Path.Add(dirInfo);
                        }
                    }
                }
            }
        }
               

        /// <summary>
        /// 
        /// </summary>
        void SearchGit() {
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
        }

        /// <summary>
        /// a heuristic, which verifies that some directory <paramref name="_d"/>
        /// contains BoSSS;
        /// </summary>
        /// <param name="_d"></param>
        /// <returns></returns>
        bool VerifyBoSSSDir(DirectoryInfo _d) {
            DirectoryInfo[] subdirs = _d.GetDirectories();
            DirectoryInfo bin = null, doc = null, etc = null;

            foreach (DirectoryInfo d in subdirs) {
                if (d.Name == "bin") bin = d;
                if (d.Name == "doc") doc = d;
                if (d.Name == "etc") etc = d;
            }

            if (bin == null || doc == null || etc == null)
                return false;

            DirectoryInfo bin_Release = null, bin_native = null;
            foreach (DirectoryInfo d in bin.GetDirectories()) {
                if (d.Name == "Release") bin_Release = d;
                if (d.Name == "native") bin_native = d;
            }
            if (bin_Release == null || bin_native == null)
                return false;

            //DirectoryInfo doc_notes = null;
            //foreach (DirectoryInfo d in doc.GetDirectories()) {
            //    if (d.Name == "notes.git") doc_notes = d;
            //}
            //if (doc_notes == null)
            //    return false;

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
        /// path to the BoSSS root folder
        /// </summary>
        public DirectoryInfo BOSSS_ROOT;

        /// <summary>
        /// path to the directory which contains the native libraries,
        /// e.g. %BOSSS_ROOT%\bin\native\win32
        /// </summary>
        public DirectoryInfo BoSSS_BIN_NATIVE;
        
        /// <summary>
        /// path to the global binary folder
        /// </summary>
        public DirectoryInfo BOSSS_BIN;

        /// <summary>
        /// path to the global documentation file folder ("doc");
        /// </summary>
        public DirectoryInfo BOSSS_DOC;

        /// <summary>
        /// path to the global confoguration file folder ("etc");
        /// </summary>
        public DirectoryInfo BOSSS_ETC;

        /// <summary>
        /// bosss L0-L4 repository: %BOSSS_ROOT%/src/public
        /// </summary>
        public DirectoryInfo BoSSS_SRC_PUBLIC;

        /// <summary>
        /// %BOSSS_ROOT/src
        /// </summary>
        public DirectoryInfo BoSSS_SRC;

        /// <summary>
        /// database explorer repository: %BOSSS_ROOT/src/DBE
        /// </summary>
        public DirectoryInfo BoSSS_SRC_DBE;

        /// <summary>
        /// tries to iniialize <see cref="BOSSS_ROOT"/>
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
                throw new EnviromentException("ERROR: unable to find BoSSS root directory - unvalid installation?");
            }
            FindBoSSSRootRec(_d);
        }


        /// <summary>
        /// checks the BoSSS installation and imports enviroment variables;
        /// </summary>
        void InitAndVerify() {

            //// find the directory which contains the executable
            //// ================================================
            {
                System.Reflection.Assembly a = System.Reflection.Assembly.GetEntryAssembly();
                //BOSSS_BIN = System.IO.Path.GetDirectoryName(a.Location);

                DirectoryInfo di = new DirectoryInfo(Path.GetDirectoryName( a.Location));
                
                FindBoSSSRootRec(di);
                
                if (di.Name.Equals("Debug"))
                    M = BinMode.Debug;
                else if (di.Name.Equals("Release"))
                    M = BinMode.Release;
                else
                    throw new EnviromentException("unable to determine build mode - Debug or Release?");
            }
            //// find the BoSSS root directory
            //// =============================
            //{
            //    //System.Diagnostics.Debugger.Break();

            //    string root = BOSSS_BIN
            //        + Path.DirectorySeparatorChar
            //        + ".."
            //        + Path.DirectorySeparatorChar
            //        + "..";

            //    string msg;
            //    DirectoryInfo root_dir = new DirectoryInfo(root);
            //    if (!VerivyRootDir(root_dir, out msg)) {
            //        // Maybe we are running from build directory

            //        string msg2;
            //        // build dir is $BOSSS_ROOT/src/Utils/bcl/bin/{Debug|Release}
            //        root_dir = root_dir.Parent.Parent.Parent;
            //        if (!VerivyRootDir(root_dir, out msg2)) {
            //            throw new ApplicationException("incomplete BoSSS installation: " + msg + " -- OR -- " + msg2);
            //        }
            //    }

            //    BOSSS_ROOT = root_dir.FullName;

            //    Console.WriteLine("ROOT DIR = " + BOSSS_ROOT);
            //}

            // Define ETC, BIN folder
            // ======================
            {
                BOSSS_DOC = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName, "doc"));

                BOSSS_ETC = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName, "etc"));

                BoSSS_SRC = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName, "src"));

                BoSSS_SRC_DBE = new DirectoryInfo(Path.Combine(BoSSS_SRC.FullName, "DBE"));

                BoSSS_SRC_PUBLIC = new DirectoryInfo(Path.Combine(BoSSS_SRC.FullName, "public"));

                BoSSS_BIN_NATIVE = new DirectoryInfo(Path.Combine(BOSSS_ROOT.FullName,
                                                     Path.Combine("bin",
                                                     Path.Combine("native", "win32"))));
                switch (M) {
                    case BinMode.Debug:
                        BOSSS_BIN = new DirectoryInfo(Path.Combine(Path.Combine(BOSSS_ROOT.FullName, "bin"), "Debug"));
                        break;
                    case BinMode.Release:
                        BOSSS_BIN = new DirectoryInfo(Path.Combine(Path.Combine(BOSSS_ROOT.FullName, "bin"), "Release"));
                        break;
                    default:
                        throw new NotImplementedException("internal error");
                }
            }

        }

        /// <summary>
        /// Debug or Release configuration?
        /// </summary>
        public BinMode M;

    }
}

