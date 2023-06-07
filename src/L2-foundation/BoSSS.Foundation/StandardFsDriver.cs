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
using System.IO;
using System.IO.Compression;
using System.Linq;
using ilPSP.Tracing;
using MPI.Wrappers;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// IO driver for standard file systems;
    /// </summary>
    public class StandardFsDriver : IFileSystemDriver, IDisposable {

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="databasePath">
        /// - if pointing to a directory, a normal BoSSS database
        /// - if pointing to a file, a ZIP archive of a BoSSS database
        /// </param>
        public StandardFsDriver(string databasePath) {
            if(!File.Exists(databasePath) && !Directory.Exists(databasePath))
                throw new FileNotFoundException($"Database directory/file {databasePath} does not exist.");

            int Rank, Size;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out Rank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out Size);

            // A hack to support zip-files
            // ---------------------------
            if(File.Exists(databasePath)) {
                // Database is provided a ZIP-file, we need to uncompress it
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                if(Rank == 0) {
                    // we use the application directory (and not some temp-directory)
                    // because the application directory is most likely to be accessible from all
                    // MPI processes.
                    System.Reflection.Assembly a = System.Reflection.Assembly.GetEntryAssembly();

                    if(a == null) {
                        // Entry assembly might be null if called from
                        // unmanaged code; fall back to executing assembly in
                        // this case
                        a = System.Reflection.Assembly.GetExecutingAssembly();
                    }

                    //string appdir = Path.GetDirectoryName(a.Location);
                    string appdir = Directory.GetCurrentDirectory();

                    UnZippedDirectory = Path.Combine(appdir, "bosss_db_tmp");
                    if(Directory.Exists(UnZippedDirectory)) {
                        // this may be some artifact from a previous run
                        Directory.Delete(UnZippedDirectory, true);
                    }
                    if(File.Exists(UnZippedDirectory)) {
                        throw new IOException("Can't create temporary database directory because there is already a file with the same name ('" + UnZippedDirectory + "').");
                    }

                    /*
                    ZipStorer zip = ZipStorer.Open(databasePath, FileAccess.Read);
                    // Read the central directory collection
                    List<ZipStorer.ZipFileEntry> dir = zip.ReadCentralDir();

                    // Look for the desired file
                    foreach(ZipStorer.ZipFileEntry entry in dir) {
                        var out_file = Path.Combine(UnZippedDirectory, entry.FilenameInZip);
                        zip.ExtractFile(entry, out_file);
                    }
                    zip.Close();
                    */

                    System.IO.Compression.ZipFile.ExtractToDirectory(databasePath, UnZippedDirectory);


                    DirectoryInfo di = new DirectoryInfo(UnZippedDirectory);
                    try {
                        VerifyDirectoryStructure(di.FullName);
                    } catch(ArgumentException) {
                        di = di.GetDirectories()[0]; // brutal hack.
                    }

                    databasePath = MPIExtensions.MPIBroadcast(di.FullName, 0, csMPI.Raw._COMM.WORLD);
                } else {
                    databasePath = MPIExtensions.MPIBroadcast(default(string), 0, csMPI.Raw._COMM.WORLD);
                }
            }

            // check base path
            VerifyDirectoryStructure(databasePath);
            BasePath = databasePath;

        }

        /// <summary>
        /// If database was extracted from zip file, this is the path, and on exit/dispose we have to delete it.
        /// </summary>
        string UnZippedDirectory = null;

        /// <summary>
        /// checks whether a path <paramref name="path"/> is the root
        /// of a valid BoSSS database (i.e. it contains all
        /// necessary subfolders);
        /// If it isn't, a <see cref="ArgumentException"/> is thrown;
        /// </summary>
        /// <param name="path"></param>
        public static void VerifyDirectoryStructure(string path) {
            if(!Directory.Exists(path)) {
                throw new ArgumentException(
                    "Database Error: base directory (" + path + ") does not exist.");
            }

            if(!Directory.Exists(Path.Combine(path, TimestepDir))) {
                throw new ArgumentException(
                    "Database Error: field state directory (" + Path.Combine(path, TimestepDir) + ") does not exist.");
            }

            if(!Directory.Exists(Path.Combine(path, DistVectorDataDir))) {
                throw new ArgumentException(
                    "Database Error: distributed vector data directory (" + Path.Combine(path, DistVectorDataDir) + ") does not exist.");
            }

            if(!Directory.Exists(Path.Combine(path, SessionsDir))) {
                throw new ArgumentException(
                    "Database Error: session directory (" + Path.Combine(path, SessionsDir) + ") does not exist.");
            }

            if(!Directory.Exists(Path.Combine(path, GridsDir))) {
                throw new ArgumentException(
                    "Database Error: grid directory (" + Path.Combine(path, GridsDir) + ") does not exist.");
            }
        }

        /// <summary>
        /// base directory of all sessions
        /// </summary>
        public static string SessionsDir {
            get {
                return "sessions";
            }
        }

        /// <summary>
        /// base directory for grids
        /// </summary>
        public static string GridsDir {
            get {
                return "grids";
            }
        }

        /// <summary>
        /// relative path to data chunks of distributed vectors
        /// </summary>
        public static string DistVectorDataDir {
            get {
                return "data";
            }
        }

        /// <summary>
        /// relative path to timestep state objects
        /// </summary>
        public static string TimestepDir {
            get {
                return "timesteps";
            }
        }

        /// <summary>
        /// Base path to a database
        /// </summary>
        public string BasePath {
            get;
            private set;
        }


        /// <summary>
        /// This is intended for files which store log-files: 
        /// Opens some file in the database **cooperatively**, i.e. if it is opened for writing (<paramref name="create"/> == true)
        /// this method will allow reading by any other process.
        /// Simultaneous writing, however, should yield to an exception.
        /// </summary>
        /// <param name="create">
        /// - true: creates a new file for writing;
        /// - false: open a file for reading;
        /// </param>
        /// <param name="RelPath">
        /// relative path within the base paths.
        /// </param>
        /// <returns>
        /// a file stream
        /// </returns>
        /// <param name="ForceOverride">
        /// when opening a stream for writing (<paramref name="create"/>=true), this argument
        /// toggles how existing files should be treated:
        /// - false: an exception is thrown if the file already exists;
        /// - true: an existing file would be overwritten
        /// </param>
        /// <returns></returns>
        private Stream OpenFileCooperatively(bool create, string RelPath, bool ForceOverride) {
            using(var tr = new FuncTrace()) {
                tr.Info("opening file '" + RelPath + "', create='" + create + "'");
                if(create) {
                    // ++++++++++++++++++++++++++++++++++
                    // create new/overwrite existing file
                    // ++++++++++++++++++++++++++++++++++
                    string fullpath = Path.Combine(BasePath, RelPath);

                    FileStream fs = null;

                    fs = new FileStream(fullpath,
                        ForceOverride ? FileMode.Create : FileMode.CreateNew,
                        FileAccess.ReadWrite, FileShare.Read);

                    return fs;

                } else {
                    // ++++++++++++++++++++++++++++
                    // try to open file for reading
                    // ++++++++++++++++++++++++++++

                    FileStream fs = new FileStream(Path.Combine(BasePath, RelPath),
                                FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                    return fs;

                }
            }

        }


        /// <summary>
        /// This is intended for files which store objects: 
        /// Opens some file in the database **exclusively**, i.e. if it is opened for writing (<paramref name="create"/> == true)
        /// this method will block on any other process 
        /// for write **and read** operations until the file 
        /// is released by the first process.
        /// Simultaneous reading should be allowed, however.
        /// </summary>
        /// <param name="create">
        /// - true: creates a new file for writing;
        /// - false: open a file for reading;
        /// </param>
        /// <param name="RelPath">
        /// relative path within the base paths.
        /// </param>
        /// <returns>
        /// a file stream
        /// </returns>
        /// <param name="ForceOverride">
        /// when opening a stream for writing (<paramref name="create"/>=true), this argument
        /// toggles how existing files should be treated:
        /// - false: an exception is thrown if the file already exists;
        /// - true: an existing file would be overwritten
        /// </param>
        /// <returns></returns>
        private Stream OpenFileExclusiveBlocking(bool create, string RelPath, bool ForceOverride) {
            using(var tr = new FuncTrace()) {
                tr.Info("opening file '" + RelPath + "', create='" + create + "'");

                string fullpath = Path.Combine(BasePath, RelPath);
               
                if(create) {
                    // ++++++++++++++++++++++++++++++++++
                    // create new/overwrite existing file
                    // ++++++++++++++++++++++++++++++++++

                    if(ForceOverride == false && File.Exists(fullpath))
                        throw new IOException($"Unable to open file for writing: not allowed to overwrite, but file '{fullpath}' already exists;");

                    FileStream fs = null;
                    int i = 0;
                    while(fs == null) {
                        try {
                            fs = new FileStream(fullpath,
                                ForceOverride ? FileMode.Create : FileMode.CreateNew,
                                FileAccess.ReadWrite, FileShare.None);
                        } catch(IOException ioe) {
                            fs = null;
                            tr.Info($"Write access to file {fullpath} is delayed by {ioe.GetType().Name}: {ioe.Message}.");
                            i++;
                            if(i > 1000) {
                                tr.Logger.Error($"Write access to file {fullpath} TIMEOUT: {ioe.GetType().Name}: {ioe.Message}.");
                                throw new IOException("File write open failed more than 1000 times. ", ioe);
                            }
                            System.Threading.Thread.Sleep(System.DateTime.Now.Millisecond + 501);
                        }
                    }
                    //if(i > 0)
                    //    Console.WriteLine();
                    return fs;

                } else {
                    // ++++++++++++++++++++++++++++
                    // try to open file for reading
                    // ++++++++++++++++++++++++++++

                    FileNotFoundException exc = null;

                    FileStream fs = null;
                    int i = 0;
                    while(fs == null) {
                        try {
                            fs = new FileStream(fullpath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                        } catch(FileNotFoundException fnf) {
                            exc = fnf;
                            throw exc;
                        } catch(IOException ioe) {
                            fs = null;
                            tr.Info($"Read access to file {fullpath} is delayed by {ioe.GetType().Name}: {ioe.Message}.");
                            i++;
                            if(i > 1000) {
                                tr.Logger.Error($"Read access to file {fullpath} TIMEOUT: {ioe.GetType().Name}: {ioe.Message}.");
                                throw new IOException("File read open failed more than 1000 times. ", ioe);
                            }
                            System.Threading.Thread.Sleep(System.DateTime.Now.Millisecond * 3 + 27);
                        }
                    }
                    //if(i > 0)
                    //    Console.WriteLine();
                    return fs;

                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="create">
        /// true: creates a new file for writing, an exception is thrown if the file already exists;
        /// false: open a file for reading;
        /// </param>
        /// <param name="vecGuid"></param>
        /// <param name="part"></param>
        /// <returns></returns>
        public Stream GetDistVectorDataStream(bool create, Guid vecGuid, int part) {
            string filename = Path.Combine(
                DistVectorDataDir, vecGuid.ToString() + "." + (part + 1) + ".data");
            return OpenFileExclusiveBlocking(create, filename, true);
        }

        /// <summary>
        /// See <see cref="IFileSystemDriver"/>
        /// </summary>
        /// <param name="create"></param>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public Stream GetTimestepStream(bool create, Guid sessionGuid) {
            string filename = Path.Combine(TimestepDir, sessionGuid.ToString() + ".ts");
            return OpenFileExclusiveBlocking(create, filename, false);
        }

        /// <summary>
        /// closes logfiles, ...
        /// </summary>
        public void Dispose() {
            if(IsDisposed)
                return;

            if(UnZippedDirectory != null) {
                if(Directory.Exists(UnZippedDirectory)) {
                    Directory.Delete(UnZippedDirectory, true);
                }
            }
            IsDisposed = true;
        }

        /// <summary>
        /// true if <see cref="Dispose"/> has been called
        /// </summary>
        public bool IsDisposed {
            get;
            private set;
        }

        /// <summary>
        /// creates a new session directory.
        /// </summary>
        public void CreateSessionDirectory(Guid sessionGuid) {
            if(sessionGuid.Equals(Guid.Empty))
                throw new ArgumentException();

            // create session directory
            string SessionDir = Path.Combine(BasePath, SessionsDir, sessionGuid.ToString());
            Directory.CreateDirectory(SessionDir);
        }

        /// <summary>
        /// <see cref="IFileSystemDriver.GetGridStream"/>
        /// </summary>
        public Stream GetGridStream(bool create, Guid GridGuid) {
            string filename = Path.Combine(GridsDir, GridGuid.ToString() + ".grid");
            return OpenFileVersionized(create, filename);
        }

        private Stream OpenFileVersionized(bool create, string filename) {
            if(create == true) {
                return OpenFileExclusiveBlocking(create, filename, true);
            } else {
                return OpenFileExclusiveBlocking(create, filename, false);
            }
        }

        /// <summary>
        /// See <see cref="IFileSystemDriver"/>
        /// </summary>
        /// <param name="create">
        /// - true for a write stream (also overwriting)
        /// - false for a read stream;
        /// </param>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public Stream GetSessionInfoStream(bool create, Guid sessionGuid) {
            string dirname = Path.Combine(SessionsDir, sessionGuid.ToString());
            string fullDirName = Path.Combine(BasePath, dirname);
            if(!Directory.Exists(fullDirName))
                Directory.CreateDirectory(fullDirName);

            string filename = Path.Combine(dirname, "Session.info");

            return OpenFileExclusiveBlocking(create, filename, true);
        }

        /// <summary>
        /// returns a stream for logging information
        /// </summary>
        /// <param name="logName">
        /// name for the logfile;
        /// this name identifies the log within the session;
        /// </param>
        /// <param name="sessionGuid">
        /// the session in which the log should be created.
        /// </param>
        /// <returns></returns>
        public Stream GetNewLogStream(string logName, Guid sessionGuid) {
            string filename = Path.Combine(
                SessionsDir, sessionGuid.ToString(), logName + ".txt");

            Stream fs = OpenFileCooperatively(true, filename, true);

            return fs;
        }

        

        /// <summary>
        /// returns a text writer for logging information
        /// </summary>
        /// <param name="logName">
        /// name for the logfile;
        /// this name identifies the log within the session;
        /// </param>
        /// <param name="sessionGuid">
        /// the session in which the log should be created.
        /// </param>
        /// <returns></returns>
        public TextWriter GetNewLog(string logName, Guid sessionGuid) {
            return new StreamWriter(GetNewLogStream(logName, sessionGuid));
        }

        /// <summary>
        /// Returns the path to the time-step log.
        /// </summary>
        /// <param name="sessionGuid">The ID of the session.</param>
        /// <returns>The path to the session's time-step log.</returns>
        public string GetTimestepLogPath(Guid sessionGuid) {
            return Path.Combine(
                BasePath, SessionsDir, sessionGuid.ToString(), "TimestepLog.txt");
        }

        /// <summary>
        /// Returns a text reader for reading the timestep log.
        /// </summary>
        /// <param name="sessionGuid">The ID of the session.</param>
        /// <returns>A text reader for the session's timestep log.</returns>
        public Stream GetTimestepLogStream(Guid sessionGuid) {
            return new FileStream(
                GetTimestepLogPath(sessionGuid),
                FileMode.Open,
                FileAccess.Read,
                FileShare.ReadWrite);
        }

        /// <summary>
        /// Tests whether a certain file exists in the database or not.
        /// </summary>
        /// <param name="RelPath"></param>
        /// <returns>
        /// true: if file exists;
        /// false: if file not exists;
        /// </returns>
        public bool FileExists(string RelPath) {
            return File.Exists(Path.Combine(BasePath, RelPath));
        }

        /// <summary>
        /// Gathers the <see cref="Guid"/>s from the names of the files and 
        /// subdirectories in the specified <paramref name="subdir"/>.
        /// </summary>
        /// <param name="subdir">
        /// The directory to search, relative from the root path of the database.
        /// </param>
        /// <param name="searchpattern">
        /// A filter to search for certain files or directories.
        /// </param>
        /// <param name="searchDirs">
        /// false, to parse file names into GUIDs;
        /// true, to parse subdirectory names.
        /// </param>
        /// <returns>
        /// A list of unique identifiers from the file or directory names.
        /// </returns>
        private IEnumerable<Guid> ParseDirectory(string subdir, string searchpattern, bool searchDirs) {
            HashSet<Guid> ret = new HashSet<Guid>();
            DirectoryInfo ddir = new DirectoryInfo(Path.Combine(BasePath, subdir));


            for (int i = 0; i < 2; i++) {

                FileSystemInfo[] gridFiles;
                if (searchDirs == false) {
                    if (i == 0)
                        gridFiles = ddir.GetFiles(searchpattern);
                    else
                        gridFiles = ddir.GetFiles(searchpattern + ".V2");
                } else {
                     var _subdirs = ddir.GetDirectories(searchpattern);


                    gridFiles = _subdirs.Where(dir => dir.GetFiles().Length > 0).ToArray(); // ignore empty session directories
                    // (sometime, when a session is deleted, the session directory cannot be deleted and an empty directory remains;
                    // this causes lots of error messages in BoSSSpad; therefore, we ignore these directories.)

                    i = 3;
                }

                foreach (var grdFile in gridFiles) {
                    string GuidStr = grdFile.Name.Substring(0, 36);
                    Guid grdGuid = new Guid(GuidStr);

                    if (!ret.Contains(grdGuid))
                        // we need to check, because there can be "old" and ".V2" - files, and
                        // because there can be the same file in multiple pa
                        ret.Add(grdGuid);
                }
            }

            return ret;
        }

        ///// <summary>
        ///// all Guid's in the 'sessions' - subdirectory
        ///// </summary>
        //public ICollection<Guid> GetAllSessions() {

        //    List<Guid> ret = new List<Guid>();

        //    DirectoryInfo di = new DirectoryInfo(Path.Combine(BasePath, SessionsDir));
        //    var allSessions = di.GetDirectories("*", SearchOption.TopDirectoryOnly);

        //    foreach (var d in allSessions) {
        //        try {
        //            Guid g = new Guid(d.Name);
        //            ret.Add(g);
        //        } catch (Exception) {
        //        }
        //    }

        //    return ret;
        //}

        /// <summary>
        /// all Guid's in the 'grids' - subdirectory
        /// </summary>
        public IEnumerable<Guid> GetAllGridGUIDs() {
            return ParseDirectory(GridsDir, "*.grid", false);
        }

        /// <summary>
        /// all guids in the 'data\headers' - subdirectoy
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Guid> GetAllDataVectorGUIDs() {
            return ParseDirectory(DistVectorDataDir, "*.1.data", false);
        }

        /// <summary>
        /// A list of the GUIDS of all sessions
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Guid> GetAllSessionGUIDs() {
            return ParseDirectory(SessionsDir, "*", true);
        }

    }
}
