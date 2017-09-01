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

using log4net;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;

namespace MPI.Wrappers.Utils {

    /// <summary>
    /// Tries to do the impossible: platform independent loading of dynamic libraries.<br/> 
    /// </summary>
    public abstract class DynLibLoader : IDisposable {

        /// <summary>
        /// the logger
        /// </summary>
        static ILog logger = LogManager.GetLogger(typeof(DynLibLoader));

        /// <summary>
        /// converts the function name <paramref name="s"/> into the desired 
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public delegate string GetNameMangling(string s);

        DynamicLibraries.DynLibHandle m_LibHandle;

        /// <summary>
        /// 'dll'/'so'--handle of the main library underlying this wrapper
        /// </summary>
        public DynamicLibraries.DynLibHandle LibHandle {
            get {
                return m_LibHandle;
            }
        }

        private DynamicLibraries.DynLibHandle[] m_PreqLibHandle;

        /// <summary>
        /// 'dll'/'so'--handles of prerequisite libraries onderlyuing this wrapper
        /// </summary>
        public IEnumerable<DynamicLibraries.DynLibHandle> PreqLibHandle {
            get {
                return m_PreqLibHandle.ToList().AsReadOnly();
            }
        }

        static bool IsDelegate(Type t) {
            if (t == null)
                return false;

            if (t.Equals(typeof(Delegate)))
                return true;
            if (t.Equals(typeof(MulticastDelegate)))
                return true;
            if (t.Equals(typeof(object)))
                return false;

            return IsDelegate(t.BaseType);
        }

        // <summary>
        // all directories in which libraries are searched;
        // </summary>
        bool ContainesWildCards(string name) {
            if (name.Contains("["))
                return true;
            if (name.Contains("]"))
                return true;
            if (name.Contains("("))
                return true;
            if (name.Contains(")"))
                return true;
            if (name.Contains("{"))
                return true;
            if (name.Contains("}"))
                return true;
            if (name.Contains("|"))
                return true;
            if (name.Contains("?"))
                return true;
            if (name.Contains("*"))
                return true;
            if (name.Contains("$"))
                return true;
            return false;
        }

        /// <summary>
        /// diagnostic output, for debugging;
        /// </summary>
        protected TextWriter m_debugoutput;


		[DllImport ("libc")]
		static extern int uname (IntPtr buf);

        static bool IsRunningOnMac ()
        {
            IntPtr buf = IntPtr.Zero;
            try {
                buf = Marshal.AllocHGlobal (8192);
                // This is a hacktastic way of getting sysname from uname ()
                if (uname (buf) == 0) {
                    string os = Marshal.PtrToStringAnsi (buf);
                    if (os == "Darwin")
                        return true;
                }
            } catch {
            } finally {
                if (buf != IntPtr.Zero)
                    Marshal.FreeHGlobal (buf);
            }
            return false;
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_LibNames">
        /// a list of library names; e.g. different names of the library on various systems, 
        /// or alternative implementations (e.g. various BLAS - implementations) on the same system
        /// </param>
        /// <param name="OsFilter">
        /// Operating system filter
        /// </param>
        /// <param name="NameMangling">
        /// name mangling, suggestions are: <see cref="LeadingUnderscore_SmallLetters"/>, <see cref="SmallLetters_TrailingUnderscore"/>
        /// or <see cref="CAPITAL_LETTERS"/>
        /// </param>
        /// <param name="PointerSizeFilter">
        /// library is excluded from search if the size of pointers (usually 4 or 8 bytes) in not equal; 
        /// This is for libs which require a different name mangling in 64 and 32 bit mode;
        /// Ignored, when negative;
        /// </param>
        /// <param name="PrequesiteLibraries">
        /// 1st index: correlates with library names
        /// 2nd index: chain of required prequesites 
        /// 3rd index: for each required prequesite, a list of alternatives.
        /// </param>
        protected DynLibLoader(string[] _LibNames, string[][][] PrequesiteLibraries, GetNameMangling[] NameMangling, PlatformID[] OsFilter, int[] PointerSizeFilter) {
            m_debugoutput = new StringWriter();
            if (_LibNames.Length != OsFilter.Length || OsFilter.Length != NameMangling.Length || PointerSizeFilter.Length != NameMangling.Length || _LibNames.Length != PrequesiteLibraries.Length)
                throw new ApplicationException("all arrays must have the same length.");
            PlatformID CurrentSys = System.Environment.OSVersion.Platform;

            // MacOS tweak
            if (CurrentSys == PlatformID.Unix && IsRunningOnMac ())
                CurrentSys = PlatformID.MacOSX;


            List<DirectoryInfo> LibrarySearchPath = GetLibrarySearchDirs(CurrentSys);

            Info("Pointer size: " + IntPtr.Size);

            for (int k = 0; k < _LibNames.Length; k++) {

                // apply OS and byteness-filters
                // =============================
                {
                    Info("Trying with '" + _LibNames[k] + "' ... ");
                    if (CurrentSys != OsFilter[k]) {
                        Info("WRONG OS");
                        continue; // don't search for the library 'LibNames[i]' on this system
                    }
                    if (PointerSizeFilter[k] >= 0 && PointerSizeFilter[k] != IntPtr.Size) {
                        Info("MISMATCH IN POINTER BYTENESS");
                        continue; // don't search for the library 'LibNames[i]' on this system
                    }
                    Info("FILERS PASSED");
                }

                // load prequesite library
                // =======================
                bool PreqSuccess = true;
                if (PrequesiteLibraries[k] != null) {
                    Info("prerequisite libraries required; trying to load these...");
                    string[][] PreqChain = PrequesiteLibraries[k];
                    this.m_PreqLibHandle = new DynamicLibraries.DynLibHandle[PreqChain.Length];

                    int cnt = 0;
                    foreach (string[] PreqAlt in PreqChain) {
                        foreach (var LibName in GetLibraryFiles(PreqAlt, LibrarySearchPath)) {
                            Info("trying to load '" + LibName + "' (prerequisite) ... ");
                            string errstr;
                            m_PreqLibHandle[cnt] = DynamicLibraries.LoadDynLib(LibName, out errstr, true);
                            if (m_PreqLibHandle[cnt].val == IntPtr.Zero) {
                                Info("FAILED: >" + errstr + "<");
                            } else {
                                Info("SUCCESS. going on to next prerequisite.");
                                break;
                            }
                        }
                        if (m_PreqLibHandle[cnt].val == IntPtr.Zero) {
                            Info("FAILED on loading any of the prerequisite libraries " + CatStrings(PreqAlt) + ".");
                            for (int i = 0; i < cnt; i++) {
                                DynamicLibraries.UnloadDynLib(m_PreqLibHandle[i]);
                            }
                            PreqSuccess = false;
                            break;
                        }
                        cnt++;
                    }
                } else {
                    this.m_PreqLibHandle = new DynamicLibraries.DynLibHandle[0];
                }

                if (!PreqSuccess) {
                    Info("UNABLE TO LOAD MAIN LIBRARY. Failed on some of the prerequesites.");
                    continue;
                }


                // load the main library
                // =====================
                bool mainLibsuccess = false;
                foreach (var LibName in GetLibraryFiles(_LibNames[k], LibrarySearchPath)) {

                    Info("trying to load '" + LibName + "' ... ");
                    string errstr;
                    m_LibHandle = DynamicLibraries.LoadDynLib(LibName, out errstr, false);
                    if (m_LibHandle.val == IntPtr.Zero) {
                        Info("FAILED: >" + errstr + "<");
                        continue;
                    } else {
                        Info("SUCCESS.");
                    }

                    bool success = LoadAllDelegates(NameMangling[k], LibName);

                    if (success) {
                        // successfully loaded all library functions
                        //Console.WriteLine("success in : " + LibNames[i]);
                        //Console.WriteLine("success.");
                        CurrentLibraryName = LibName;
                        mainLibsuccess = true;
                        break;
                    } else {
                        UnloadAllDelegates();
                    }
                }

                if (mainLibsuccess)
                    return;
            }

            // error
            {
                StringWriter stw = new StringWriter();
                int cnt = 1;
                foreach (string s in _LibNames) {
                    stw.Write(s);
                    if (cnt < _LibNames.Length)
                        stw.Write(", ");
                    cnt++;
                }
                throw new ApplicationException("unable to find/load dynamic library " + stw.ToString() + " on current system."
                    + "Extended info:" + System.Environment.NewLine
                    + "------------------------------------------------------------------------" + System.Environment.NewLine
                    + m_debugoutput.ToString()
                    + "------------------------------------------------------------------------" + System.Environment.NewLine);
            }
        }

        static private string CatStrings(IEnumerable<string> s) {
            StringWriter ret = new StringWriter();
            var t = s.ToArray();
            for (int i = 0; i < t.Length; i++) {
                ret.Write(t[i]);
                if (i < t.Length - 1)
                    ret.Write(", ");
            }
            return ret.ToString();
        }

        private List<string> GetLibraryFiles(IEnumerable<string> LibNames, List<DirectoryInfo> LibrarySearchPath) {
            List<string> ret = new List<string>();
            foreach (var ln in LibNames) {
                ret.AddRange(GetLibraryFiles(ln, LibrarySearchPath));
            }
            return ret;
        }

        private List<string> GetLibraryFiles(string _LibName, List<DirectoryInfo> LibrarySearchPath) {
			List<string> LibNames = new List<string>();

            if(Path.IsPathRooted(_LibName)) {
                Info ("Got absolute library path '" + _LibName + "'. ");
                if(File.Exists (_LibName)) {
                    Info ("Absolute path file '" + _LibName + "' exists. ");
                    LibNames.Add (_LibName);
                } else {
                    Info ("Absolute path file '" + _LibName + "' is missing. ");
                }
				return LibNames;
            }



            if (!ContainesWildCards(_LibName)) {
                LibNames.Add(_LibName); // NO wildcard -- add without a path, so that the operation system has a chance to pick
                                        // the system-default version.
            }

            // now, we search for the lib in different paths:
            foreach (DirectoryInfo di in LibrarySearchPath) {
                Info("searching for pattern '" + _LibName + "' in '" + di.FullName + "' ... ");
                if (!di.Exists) {
                    Info("Skipping: " + di.FullName + " does not exist.");
                    continue;
                }

                FileInfo[] files;
                try {
                    files = di.GetFiles(_LibName);
                } catch (DirectoryNotFoundException) {
                    // may happen if _LibName is a relative path, e.g. "openmpi/libmpi_mpifh.so"
                    Info("Skipping: " + di.FullName + " does not exist.");
                    continue;
                }

                foreach (var f in files) {
                    Info("Found '" + f.FullName + "'");
                    LibNames.Add(f.Name);
                }
            }

            return LibNames;
        }

        private bool LoadAllDelegates(GetNameMangling NameMangling, string LibName) {
            Type myType = this.GetType();
            FieldInfo[] fields = myType.GetFields(BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public);
            m_LoadedSymbols.Clear();
            Delegate2FunctionPointer.Clear();

            // loop over all delegates ....
            bool success = true;
            foreach (FieldInfo fld in fields) {
                if (IsDelegate(fld.FieldType)) {
                    m_LoadedSymbols.Add(fld);

                    // get function name in DLL
                    string UnmanagedName = NameMangling(fld.Name);
                    Info("Trying to load '" + fld.Name + "' as '" + UnmanagedName + "' ...");

                    // get function pointer
                    string errstr2;
                    IntPtr FuncPtr = DynamicLibraries.LoadSymbol(m_LibHandle, UnmanagedName, out errstr2);
                    if (FuncPtr == IntPtr.Zero) {
                        //throw new ApplicationException("Library '" + LibNames[i] + "' not working - missing function '" + UnmanagedName + "';");
                        // try next library ..
                        DynamicLibraries.UnloadDynLib(m_LibHandle);
                        m_LibHandle.val = IntPtr.Zero;
                        success = false;
                        Info("FAILED: >" + errstr2 + "<");
                        break;
                    } else {
                        Info("SUCCESS.");
                    }


                    // create delegate
                    Delegate del = Marshal.GetDelegateForFunctionPointer(FuncPtr, fld.FieldType);
                    if (del == null) {
                        Warn("unable to get delegate from function pointer: library '" + LibName + "', symbol '" + UnmanagedName + "';");
                        throw new ApplicationException("unable to get delegate from function pointer: library '" + LibName + "', symbol '" + UnmanagedName + "';");
                    }
                    fld.SetValue(this, del);
                    Delegate2FunctionPointer.Add(del, FuncPtr);
                }
            }
            return success;
        }

        /// <summary>
        /// all directories in "PATH" or "LD_LIBRARY_PATH"
        /// </summary>
        private List<DirectoryInfo> GetLibrarySearchDirs(PlatformID CurrentSys) {
            List<DirectoryInfo> LibrarySearchPath = new List<DirectoryInfo>();
            {
                switch (CurrentSys) {
                    case PlatformID.MacOSX:
                    case PlatformID.Unix: {

                            AddDirsFromSystemEnvironment(LibrarySearchPath, "LD_LIBRARY_PATH", new char[] { ':' });

                            foreach (var dir in new string[] { "/usr/lib", "/lib" }) {
                                try {
                                    var di = new DirectoryInfo(dir);
                                    if (di.Exists) {
                                        LibrarySearchPath.Add(di);
                                        m_debugoutput.WriteLine("Adding search directory: " + dir);
                                    }
                                } catch (Exception) {

                                }
                            }
                            break;
                        }
                    case PlatformID.Win32S:
                    case PlatformID.Win32Windows:
                    case PlatformID.Win32NT: {
                            string varName = "PATH";
                            AddDirsFromSystemEnvironment(LibrarySearchPath, varName, new char[] { ';' });
                            break;
                        }

                    default:
                        break;
                }
            }
            return LibrarySearchPath;
        }

        private void AddDirsFromSystemEnvironment(List<DirectoryInfo> LibrarySearchPath, string varName, char[] splitChars) {
            string ld_library_path = System.Environment.GetEnvironmentVariable(varName);
            Info("PATH = " + ld_library_path);
            if (ld_library_path == null) {
                ld_library_path = "";
            }
            string[] dirs = ld_library_path.Split(splitChars, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 0; i < dirs.Length; i++) {

                try {
                    var di = new DirectoryInfo(dirs[i]);
                    if (di.Exists) {
                        Info("Adding search directory: " + dirs[i]);
                        LibrarySearchPath.Add(di);
                    } else {
                        Warn("Nonexistent entry in " + varName + ": directory: " + dirs[i]);
                    }
                } catch (Exception) {
                    Warn("Strange entry in " + varName + ": directory: " + dirs[i]);
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public Dictionary<Delegate, IntPtr> Delegate2FunctionPointer = new Dictionary<Delegate, IntPtr>();


        /// <summary>
        /// all members of this object that were initialized by the constructor
        /// </summary>
        protected List<FieldInfo> m_LoadedSymbols = new List<FieldInfo>();

        /// <summary>
        /// gets the actual name of the loaded library
        /// </summary>
        public string CurrentLibraryName {
            get;
            private set;
        }

        /// <summary>
        /// sets all fields in <see cref="m_LoadedSymbols"/> to null;
        /// </summary>
        protected void UnloadAllDelegates() {
            foreach (FieldInfo f in m_LoadedSymbols) {
                f.SetValue(this, null);
            }
        }

        /// <summary>
        /// Logs a message if compiled in debug mode.
        /// </summary>
        /// <param name="obj">
        /// The message to be logged
        /// </param>
        private void Info(object obj) {
#if DEBUG
            logger.Info(obj);
#endif
            m_debugoutput.WriteLine(obj.ToString());
        }

        /// <summary>
        /// Issues a warning if compiled in debug mode.
        /// </summary>
        /// <param name="obj">
        /// The message to be logged
        /// </param>
        private void Warn(object obj) {
#if DEBUG
            logger.Warn(obj);
#endif
            m_debugoutput.WriteLine(obj.ToString());
        }


        /// <summary>
        /// DLL function name mangling: no change
        /// </summary>
        public static string Identity(string Name) {
            return Name;
        }

        /// <summary>
        /// DLL function name mangling: capital letters
        /// </summary>
        public static string CAPITAL_LETTERS(string Name) {
            return Name.ToUpperInvariant();
        }

        /// <summary>
        /// DLL function name mangling: lower case names, with leading underscore
        /// </summary>
        public static string SmallLetters_TrailingUnderscore(string Name) {
            return Name.ToLowerInvariant() + "_";
        }

        /// <summary>
        /// DLL function name mangling: lower case names, with trailing underscore
        /// </summary>
        public static string LeadingUnderscore_SmallLetters(string Name) {
            return "_" + Name.ToLowerInvariant();
        }

        ///// <summary>
        ///// destructor - calls dispose
        ///// </summary>
        //~DynLibLoader() {
        //    this.Dispose();
        //}

        #region IDisposable Members

        /// <summary>
        /// unloads the library
        /// </summary>
        virtual public void Dispose() {
            if (m_LibHandle.val != IntPtr.Zero) {
                DynamicLibraries.UnloadDynLib(m_LibHandle);
                m_LibHandle.val = IntPtr.Zero;
                UnloadAllDelegates();
            }
        }

        #endregion
    }
}
