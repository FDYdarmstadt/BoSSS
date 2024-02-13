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
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;
using ilPSP;
using ilPSP.Tracing;
using log4net;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Collection of random utility functions.
    /// </summary>
    public static class Utils {

        /// <summary>
        /// BoSSS install directory, e.g. 'C:\Program Files\FDY\BoSSS'.
        /// Within this directory, one would find e.g. the subdirectory 'bin\native\win\amd64'
        /// where native libraries are included.
        /// </summary>
        public const string BOSSS_INSTALL = "BOSSS_INSTALL";

        /// <summary>
        /// Name of environment variable which 
        /// overrides for location of native libraries if defined.
        /// </summary>
        public const string BOSSS_NATIVE_OVERRIDE = "BOSSS_NATIVE_OVERRIDE";

        /// <summary>
        /// Returns the installation path of the BoSSS native libraries.
        /// The search priority should be:
        /// 1. if defined, the <see cref="BOSSS_NATIVE_OVERRIDE"/> environment variable
        /// 2. if existent, a local 'amd64' subdirectory
        /// 3. if <see cref="BOSSS_INSTALL"/>, the respective subdirectory, e.g. `C:\Program Files\FDY\BoSSS\bin\native\win\amd64`
        /// </summary>
        public static string GetNativeLibraryDir(ILog logger = null) {
            {
                string overr = System.Environment.GetEnvironmentVariable(BOSSS_NATIVE_OVERRIDE);
                if (overr != null) {
                    Console.WriteLine($"Note: {BOSSS_NATIVE_OVERRIDE} variable is defined: libraries will be searched at: '{overr}'.");
                    return overr;
                }
            }

            string ilPSP_Directory = GetBoSSSInstallDir(logger);

            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                // ++++++++++
                // MS windows 
                // ++++++++++

                // search for "amd64"  subdirectories locally (in application directory)
                // =====================================================================
                DirectoryInfo nativeDir = null; // directory where we search for native lib's
                {
                    DirectoryInfo di = null;

                    System.Reflection.Assembly a = System.Reflection.Assembly.GetEntryAssembly();
                    if (a == null) {
                        // Entry assembly might be null if called from
                        // unmanaged code; fall back to executing assembly in
                        // this case
                        a = System.Reflection.Assembly.GetExecutingAssembly();
                    }

                    di = new DirectoryInfo(Path.GetDirectoryName(a.Location));

                    if (IntPtr.Size == 8) {
                        // seem to run on 64-Bit windows
                        //Console.WriteLine("running 64 bit");
                        nativeDir = di.GetDirectories("amd64").FirstOrDefault();
                    } else if (IntPtr.Size == 4) {
                        // seem to run on 32-Bit windows
                        nativeDir = di.GetDirectories("x86").FirstOrDefault();

                    } else {
                        throw new ApplicationException("something very strange: IntPtr.Size == " + IntPtr.Size + ".");
                    }

                    if (nativeDir != null) {
                        return nativeDir.FullName;
                    }
                }


                // if not found, search for a BoSSS-Installation globally
                // ======================================================

                if (nativeDir == null || !nativeDir.Exists) {

                    DirectoryInfo _di = null;

                    if (ilPSP_Directory != null && ilPSP_Directory.Length > 0) {
                        // search in the path of the optional ilPSP directory
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++

                        DirectoryInfo d = new DirectoryInfo(ilPSP_Directory);
                        _di = new DirectoryInfo(Path.Combine(d.FullName, "bin", "native", "win"));
                        if (!_di.Exists) {
                            Console.Error.WriteLine("WARNING: illegal BoSSS installation; missing directory '" + _di.FullName + "';");
                        }
                    } else {
                        Console.Error.WriteLine("WARNING: Native libraries: local search failed, missing BOSSS_INSTALL - depending on system settings.");
                        _di = new DirectoryInfo(".");
                    }

                    if (IntPtr.Size == 8) {
                        // seem to run on 64-Bit windows
                        //Console.WriteLine("running 64 bit");
                        nativeDir = _di.GetDirectories("amd64").FirstOrDefault();
                    } else if (IntPtr.Size == 4) {
                        // seem to run on 32-Bit windows
                        Console.Error.WriteLine("Warning: seem to run on 32 bit - unless you build native libraries yourself, this is not supported; expecting libraries in directory 'x86'.");
                        nativeDir = _di.GetDirectories("x86").FirstOrDefault();
                    } else {
                        throw new ApplicationException("something very strange: IntPtr.Size == " + IntPtr.Size + ".");
                    }

                    if (nativeDir != null)
                        return nativeDir.FullName;
                }

            } else if (System.Environment.OSVersion.Platform == PlatformID.Unix || System.Environment.OSVersion.Platform == PlatformID.MacOSX) {
                // ++++
                // Unix
                // ++++

                // search for "amd64"  subdirectories locally (in application directory)
                // =====================================================================
                DirectoryInfo nativeDir = null; // directory where we search for native lib's
                {
                    DirectoryInfo di = null;

                    System.Reflection.Assembly a = System.Reflection.Assembly.GetEntryAssembly();
                    if (a == null) {
                        // Entry assembly might be null if called from
                        // unmanaged code; fall back to executing assembly in
                        // this case
                        a = System.Reflection.Assembly.GetExecutingAssembly();
                    }

                    di = new DirectoryInfo(Path.GetDirectoryName(a.Location));

                    if (IntPtr.Size == 8) {
                        nativeDir = di.GetDirectories("amd64-openmpi").FirstOrDefault();
                    } else if (IntPtr.Size == 4) {
                        // seem to run on 32-Bit windows
                        Console.Error.WriteLine("Warning: seem to run on 32 bit - unless you build native libraries yourself, this is not supported.");

                    } else {
                        throw new ApplicationException("something very strange: IntPtr.Size == " + IntPtr.Size + ".");
                    }
                }

                if (nativeDir != null)
                    return nativeDir.FullName;

                // if not found, search for a BoSSS-Installation globally
                // ======================================================

                if (nativeDir == null || !nativeDir.Exists) {

                    DirectoryInfo _di = null;

                    if (ilPSP_Directory != null && ilPSP_Directory.Length > 0) {
                        // search in the path of the optional ilPSP directory
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++

                        DirectoryInfo d = new DirectoryInfo(ilPSP_Directory);
                        _di = new DirectoryInfo(Path.Combine(d.FullName, "bin", "native", "linux"));
                        if (!_di.Exists) {
                            Console.Error.WriteLine("WARNING: illegal BoSSS installation; missing directory '" + _di.FullName + "';");
                        }
                    } else {
                        Console.Error.WriteLine("WARNING: Native libraries: local search failed, missing BOSSS_INSTALL - depending on system settings.");
                        _di = new DirectoryInfo(".");
                    }

                    if (IntPtr.Size == 8) {
                        // seem to run on 64-Bit windows
                        //Console.WriteLine("running 64 bit");
                        nativeDir = _di.GetDirectories("amd64-openmpi").FirstOrDefault();
                    } else if (IntPtr.Size == 4) {
                        // seem to run on 32-Bit linux
                        nativeDir = _di.GetDirectories("x86").FirstOrDefault();
                    } else {
                        throw new ApplicationException("something very strange: IntPtr.Size == " + IntPtr.Size + ".");
                    }

                }

                if (nativeDir != null)
                    return nativeDir.FullName;
            }

            Console.Error.WriteLine("Unable to find directory for native libraries; BoSSS will most likely crash.");
            return "";
        }

        /// <summary>
        /// searches for the User- or Machine-environment variable 'BOSSS_INSTALL'
        /// and verifies the existence of this directory.
        /// </summary>
        /// <returns></returns>
        public static string GetBoSSSInstallDir(ILog logger = null) {
            logger = logger ?? NullLog.Instance;

            void Print(string s) {
                logger.Info(s);
            }

            string si1 = System.Environment.GetEnvironmentVariable(
                BOSSS_INSTALL, EnvironmentVariableTarget.User);
            if (si1 != null && si1.Length > 0)
                Print($"found USER variable '{BOSSS_INSTALL}': '" + si1 + "'.");
            else
                Print($"unable to find a USER variable '{BOSSS_INSTALL}'.");

            string si2 = System.Environment.GetEnvironmentVariable(BOSSS_INSTALL, EnvironmentVariableTarget.Machine);
            if (si2 != null && si2.Length > 0)
                Print($"found MACHINE variable '{BOSSS_INSTALL}': '" + si2 + "'.");
            else
                Print($"unable to find a MACHINE variable '{BOSSS_INSTALL}'.");

            string si3 = System.Environment.GetEnvironmentVariable (BOSSS_INSTALL);
            if (si3 != null && si3.Length > 0)
                Print ($"found variable '{BOSSS_INSTALL}': '" + si3 + "'.");
            else
                Print($"unable to find a variable '{BOSSS_INSTALL}'.");


            string si = null;
            if (si1 != null && si1.Length > 0) {
                si = si1;
                Print($"Picking USER setting.");
            } else if (si2 != null && si2.Length > 0) {
                si = si2;
                Print($"Picking MACHINE setting.");
            } else if (si3 != null && si3.Length > 0) {
                si = si3;
                Print($"Picking general setting.");
            }

            if (si != null && si.Length > 0) {
                bool success = false;
                try {
                    success = !Directory.Exists(si);
                } catch (Exception) {
                }

                if (success) {
                    //logger.Error(Err);
                    //throw new ApplicationException(Err);
                    Console.WriteLine($"   !!!  WARNING  !!! Environment variable '{BOSSS_INSTALL}' is defined as '" + si + "', but directory seems to be non-existent. Could lead to Problems!");
                }
            } else {
                logger.Info($"Unable to find environment variable '{BOSSS_INSTALL}'; Hopefully the directory of the executable contains 'amd64' and 'x86', otherwise we're screwed.");
            }
            return si;
        }

        /// <summary>
        /// Finds the settings directory (%USERPROFILE%/.BoSSS);
        /// If not existent (1st startup), the directory and a dummy config is created.
        /// </summary>
        public static string GetBoSSSUserSettingsPath() {
            using (var tr = new FuncTrace()) {
                string UserProfile = System.Environment.GetEnvironmentVariable("USERPROFILE")
                    ?? System.Environment.GetEnvironmentVariable("HOME");

                if (UserProfile == null) {
                    return "";
                }

                string settingsPath = Path.Combine(UserProfile, ".BoSSS");
                if (!Directory.Exists(settingsPath)) {
                    try {
                        Console.WriteLine("Creating: " + settingsPath);
                        Directory.CreateDirectory(settingsPath);
                    } catch (IOException e) {
                        tr.Info(String.Format(
                            "Creating user settings path failed with message '{0}';"
                                + " proceeding without user settings",
                            e.Message));
                        return "";
                    }
                }

                string etcpath = Path.Combine(settingsPath, "etc");
                if (!Directory.Exists(etcpath)) {
                    try {
                        Console.WriteLine("Creating: " + etcpath);
                        Directory.CreateDirectory(etcpath);
                    } catch (IOException e) {
                        tr.Info(String.Format(
                            "Creating user settings 'etc' sub-path failed with message '{0}';"
                                + " proceeding without user settings",
                            e.Message));
                        return "";
                    }
                }

                string batchpath = Path.Combine(settingsPath, "batch");
                if (!Directory.Exists(batchpath)) {
                    try {
                        Console.WriteLine("Creating: " + batchpath);
                        Directory.CreateDirectory(batchpath);
                    } catch (IOException e) {
                        tr.Info(String.Format(
                            "Creating user settings 'batch' sub-path failed with message '{0}';"
                                + " proceeding without user settings",
                            e.Message));
                        return "";
                    }
                }


                string dbeConfigFilePath = Path.Combine(etcpath, "DBE.xml");
                if (!File.Exists(dbeConfigFilePath)) {
                    try {
                        XmlDocument doc = new XmlDocument();
                        doc.LoadXml(Properties.Resources.DBE_empty);
                        doc.Save(dbeConfigFilePath);
                    } catch (IOException e) {
                        tr.Info(String.Format(
                            "Creating default DBE.xml failed with message '{0}';"
                                + " proceeding without DBE.xml",
                            e.Message));
                        return "";
                    }
                }

                return settingsPath;
            }
        }

        /// <summary>
        /// If not existent (1st time plot), the directory is created
        /// </summary>
        public static string GetExportOutputPath() {
            string TempDir;
            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                TempDir = System.Environment.GetEnvironmentVariable("LOCALAPPDATA");
            } else {
                TempDir = System.Environment.GetEnvironmentVariable("HOME");
            }

            if(TempDir.IsEmptyOrWhite()) {
                TempDir = Directory.GetCurrentDirectory();
            }

            DirectoryInfo BoSSSTempDir = new DirectoryInfo(Path.Combine(TempDir, "BoSSS"));
            if (!BoSSSTempDir.Exists)
                BoSSSTempDir.Create();
            DirectoryInfo plot = new DirectoryInfo(Path.Combine(BoSSSTempDir.FullName, "plots"));
            if (!plot.Exists)
                plot.Create();
            return plot.FullName;
        }

        /// <summary>
        /// Retrieves the directory where the exports for the selected
        /// <paramref name="session"/> are stored.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static string GetExportDirectory(ISessionInfo session) {
            string dirname =
                (session.ProjectName.IsEmptyOrWhite() ? "NO-PROJ" : session.ProjectName)
                + "__" +
                (session.Name.IsEmptyOrWhite() ? "NO-NAME" : session.Name)
                + "__" +
                session.ID.ToString();

            foreach (char c in System.IO.Path.GetInvalidFileNameChars()) {
                if (dirname.Contains(c))
                    dirname = dirname.Replace(c, '_');
            }

            return Path.Combine(
                Utils.GetExportOutputPath(),
                "sessions",
                   dirname);
        }

        /// <summary>
        /// Retrieves the directory where the plots for the selected
        /// <paramref name="grid"/> are stored.
        /// </summary>
        /// <param name="grid">
        /// The selected grid.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static string GetExportDirectory(IGridInfo grid) {
            string path = Path.Combine(
                Utils.GetExportOutputPath(),
                StandardFsDriver.GridsDir, grid.ID.ToString());
            return (path);
        }

        /// <summary>
        /// Formats an XmlDocment in a pretty manner
        /// </summary>
        /// <param name="xmlAsString">The xml to be formatted</param>
        /// <returns>
        /// A pretty version of <paramref name="xmlAsString"/>
        /// </returns>
        public static string PrettyFormatXMLAsString(string xmlAsString) {
            StringBuilder stringBuilder = new StringBuilder();
            StringWriter stringWriter = new StringWriter(stringBuilder);
            XmlTextWriter xmlWriter = new XmlTextWriter(stringWriter);
            xmlWriter.Formatting = Formatting.Indented;

            try {
                XmlDocument doc = new XmlDocument();
                doc.LoadXml(xmlAsString.Trim(new char[] { (char)0, ' ', '\r', '\n', '\t' }));
                doc.Save(xmlWriter);

                return stringBuilder.ToString();
            } catch {
                return null;
            } finally {
                xmlWriter.Close();
                stringWriter.Dispose();
            }
        }

        /// <summary>
        /// Formats a XmlDocment in a pretty manner and writes the result into
        /// a document.
        /// </summary>
        /// <param name="xmlAsString">The xml to be formatted</param>
        /// <returns>A document consisting of the pretty xml</returns>
        public static XmlDocument PrettyFormatXMLAsXmlDocument(string xmlAsString) {
            try {
                XmlDocument doc = new XmlDocument();
                doc.LoadXml(xmlAsString.Trim(new char[] { (char)0, ' ', '\r', '\n', '\t' }));
                return doc;
            } catch {
                return null;
            }
        }

        /// <summary>
        /// Concatenate two strings with a combination string in between
        /// </summary>
        /// <param name="s1">First string</param>
        /// <param name="s2">Second string</param>
        /// <param name="maxLength1">
        /// Maximum length of the first string in the output
        /// </param>
        /// <param name="maxLength2">
        /// Maximum length of the second string in the output
        /// </param>
        /// <param name="combineString">
        /// String appearing in the middle between (the cutted)
        /// <paramref name="s1"/> and <paramref name="s2"/>.
        /// </param>
        /// <returns>
        /// <paramref name="s1"/> + <paramref name="combineString"/>
        /// + <paramref name="s2"/> with the additional constraint that
        /// <paramref name="s1"/> and <paramref name="s2"/> will be chopped
        /// according to <paramref name="maxLength1"/> and
        /// <paramref name="maxLength2"/> respectively.
        /// </returns>
        /// <remarks>
        /// <paramref name="combineString"/> will be omitted if either of the
        /// strings to be combined is null.
        /// </remarks>
        public static string CreateCombinedString(string s1, string s2, int maxLength1, int maxLength2, string combineString) {
            s1 = ChopString(s1, maxLength1);
            s2 = ChopString(s2, maxLength2);

            if (s1 == null) {
                return s2;
            } else if (s2 == null) {
                return s1;
            } else {
                return s1 + combineString + s2;
            }
        }

        /// <summary>
        /// Returns true if the specified path contains only valid characters.
        /// </summary>
        /// <param name="path">Any path. May be null or empty.</param>
        static public bool IsValidPath(string path) {
            Regex r = new Regex(@"^(([a-zA-Z]\:)|(\\))(\\{1}|((\\{1})[^\\]([^/:*?<>""|]*))+)$");
            return r.IsMatch(path);
        }

        /// <summary>
        /// Chops a string after <paramref name="maxLength"/> characters.
        /// </summary>
        /// <param name="s">The string to be chopped</param>
        /// <param name="maxLength">The cut-off length</param>
        /// <returns>
        /// The chopped string; null, if <paramref name="s"/> was null
        /// </returns>
        private static string ChopString(string s, int maxLength) {
            if (s == null) {
                return null;
            } else {
                return s.Substring(0, Math.Min(s.Length, maxLength));
            }
        }

        /// <summary>
        /// Gathers paths to all files in the specified subdirectory
        /// whose file names match any Guid out of a given set of Guids.
        /// </summary>
        /// <param name="uids">A collection of file Guids</param>
        /// <param name="dirPath">The path to the directory containing the files</param>
        /// <param name="extension">Optional specification of the file name
        /// extension (without the preceding dot)</param> 
        /// <returns>A collection of paths to the files associated with the IDs</returns>
        public static IEnumerable<string> GetPathsFromGuids(IEnumerable<Guid> uids, string dirPath, string extension = "*") {
            return uids.SelectMany(uid => GetPathsFromGuid(uid, dirPath, extension));
        }

        /// <summary>
        /// Retrieves the paths to the files in the specified subdirectory
        /// whose file names match a given Guid.
        /// </summary>
        /// <param name="uid">A file Guid</param>
        /// <param name="dirPath">The path to the directory containing the file</param>
        /// <param name="extension">Optional specification of the file name 
        /// extension (without the preceding dot).</param>
        /// <returns>A path to the file associated with the Guid</returns>
        public static IEnumerable<string> GetPathsFromGuid(Guid uid, string dirPath,
            string extension = "*") {
            IEnumerable<string> foundFilePaths = Directory.GetFiles(dirPath, uid.ToString() + "*")
                .Where(f => extension.Equals("*") || f.EndsWith("." + extension));
            // small dirty hack here ^ to make the wild card work
            int count = foundFilePaths.Count();

            if (count == 0) {
                throw new Exception("No file with matching Guid and extension found.");
            } else {
                return foundFilePaths;
            }
        }

      
        /// <summary>
        /// Retrieves the write time of a physical file associated with an
        /// ISessionInfo object.
        /// </summary>
        /// <param name="session">The session in question.</param>
        /// <returns>The last time the file has been written to disk.</returns>
        public static DateTime GetSessionFileWriteTime(ISessionInfo session) {
            
            try {
                
                // way to slow: searching for the dir every time when we check if the session info is up-to-date results in a unresponsive database
                //string sessFolderPath = Directory.GetDirectories(
                //    Path.Combine(session.Database.Controller.DBDriver.FsDriver.BasePath, StandardFsDriver.SessionsDir),
                //    session.ID.ToString()).FirstOrDefault();
                //if (sessFolderPath == null) {

                string sessFolderPath = Path.Combine(session.Database.Controller.DBDriver.FsDriver.BasePath, StandardFsDriver.SessionsDir, session.ID.ToString());
                if(!Directory.Exists(sessFolderPath)) {
                    // Session has probably been deleted; Setting it to max
                    // value makes sure cached data disappears
                    return DateTime.MaxValue;
                }

                string entityFilePath = Path.Combine(sessFolderPath, "Session.info");
                var t = File.GetLastWriteTime(entityFilePath);
                return t;
            } catch (Exception) {
                Console.WriteLine("Error at loading of session:"+ session.ID);
                throw;
            }

            
        }

        /// <summary>
        /// Retrieves the write time of a physical file associated with an
        /// IGridInfo object.
        /// </summary>
        /// <param name="grid">The grid in question.</param>
        /// <returns>The last time the file has been written to disk.</returns>
        public static DateTime GetGridFileWriteTime(IGridInfo grid) {
            string gridFolderPath = Path.Combine(
                grid.Database.Path, StandardFsDriver.GridsDir);
            string gridFileName = grid.ID.ToString() + ".grid";
            string gridFilePath = Path.Combine(gridFolderPath, gridFileName);
            return File.GetLastWriteTime(gridFilePath);
        }

        /// <summary>
        /// Retrieves the write time of a physical file associated with an
        /// ITimestepInfo object.
        /// </summary>
        /// <param name="timestep">The timestep in question.</param>
        /// <returns>The last time the file has been written to disk.</returns>
        public static DateTime GetTimestepFileWriteTime(ITimestepInfo timestep) {
            string timestepFolderPath = Path.Combine(timestep.Database.Path,
                    StandardFsDriver.TimestepDir);
            string timestepFileName = timestep.ID.ToString() + ".ts";
            string timestepFilePath = Path.Combine(timestepFolderPath, timestepFileName);
            return File.GetLastWriteTime(timestepFilePath);
        }
    }
}
