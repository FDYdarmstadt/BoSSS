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
using System.Linq;
using System.Xml;
using ilPSP;
using BoSSS.Platform;
using Renci.SshNet;
using Renci.SshNet.Common;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Manages operations on higher-level objects of a database,
    /// e.g. session and grid info objects, as opposed to the file-based
    /// operations in <see cref="BoSSS.Foundation.IO.IFileSystemDriver "/>
    /// </summary>
    public class DatabaseController : IDatabaseController {

        private static IDictionary<string, IFileSystemDriver> fsDrivers =
            new Dictionary<string, IFileSystemDriver>();

        /// <summary>
        /// Allows registering a callback function that can ask for passwords
        /// if the connection to the database is password-proected
        /// </summary>
        public static Func<string, string> PasswordCallback;

        /// <summary>
        /// Creates the appropriate file system driver for accessing the given
        /// <paramref name="path"/>
        /// </summary>
        /// <param name="path">
        /// A fully qualified path.
        /// </param>
        /// <returns>
        /// If <paramref name="path"/> starts with 'sftp', an
        /// <see cref="SFTPFileSystemDriver"/> will be used. If the path is
        /// null, <see cref="NullFileSystemDriver"/> will be used. Otherwise,
        /// a local path will be assumed, which thus makes use of
        /// <see cref="StandardFsDriver"/>
        /// </returns>
        public static IFileSystemDriver GetFileSystemDriver(string path) {
            if (fsDrivers.ContainsKey(path)) {
                if (fsDrivers[path].IsDisposed) {
                    fsDrivers.Remove(path);
                } else {
                    return fsDrivers[path];
                }
            }

            IFileSystemDriver fsDriver;
            string sftpString = "sftp://";
            if (path.StartsWith(sftpString)) {
                try {
                    string s = path.Substring(sftpString.Length);
                    string user = s.Substring(0, s.IndexOf('@'));

                    s = s.Substring(s.IndexOf('@') + 1);
                    int hostPartEnd = s.IndexOf(':');

                    string host;
                    string relPath;
                    if (hostPartEnd >= 0) {
                        host = s.Substring(0, hostPartEnd);
                        relPath = s.Substring(s.IndexOf(':') + 1);
                    } else {
                        host = s;
                        relPath = "/";
                    }

                    KeyboardInteractiveConnectionInfo ci =
                        new KeyboardInteractiveConnectionInfo(host, user);
                    ci.AuthenticationPrompt +=
                            delegate(object sender, AuthenticationPromptEventArgs arg) {
                                var prompt = arg.Prompts.SingleOrDefault(p => p.Request == "Password: ");
                                if (prompt != null) {
                                    prompt.Response = PasswordCallback(user + "@" + host);
                                }
                            };
                    ci.Timeout = TimeSpan.FromSeconds(15);

                    fsDriver = new SFTPFileSystemDriver(ci, relPath);
                } catch (Exception) {
                    throw new ArgumentException(String.Format(
                        "Given path '{0}' has an invalid format", path));
                }
            } else if (File.Exists(path) || Directory.Exists(path)) {
                fsDriver = new StandardFsDriver(path);
            } else if (path == null) {
                fsDriver = NullFileSystemDriver.Instance;
            } else {
                throw new IOException(String.Format(
                    "Could not find or access a database located at '{0}'", path));
            }

            fsDrivers.Add(path, fsDriver);
            return fsDriver;
        }

        /// <summary>
        /// Loads the list of databases from the user's DBE.xml.
        /// </summary>
        /// <returns>A list of databases as listed in DBE.xml.</returns>
        public static IList<IDatabaseInfo> LoadDatabaseInfosFromXML() {
            string path = Path.Combine(Utils.GetBoSSSUserSettingsPath(), "etc", "DBE.xml");
            if (!File.Exists(path)) {
                return new List<IDatabaseInfo>();
            }

            Stream fs = new FileStream(
                path,
                FileMode.Open,
                FileAccess.Read);

            XmlDocument xmlDoc = new XmlDocument();
            xmlDoc.Load(fs);
            XmlNodeList xmlDatabases = xmlDoc.SelectNodes("/DBEControl/Databases/Database");

            IList<IDatabaseInfo> databases = new List<IDatabaseInfo>();


            foreach (XmlNode xmlDatabase in xmlDatabases) {
                XmlElement xmlPath = xmlDatabase.SelectSingleNode("path") as XmlElement;
                string dbpath = xmlPath.GetAttribute("value");
                databases.Add(new DatabaseInfo(dbpath));
            }

            return databases;
        }

        /// <summary>
        /// Creates a new instance of <see cref="DatabaseController"/>
        /// </summary>
        public DatabaseController(DatabaseInfo database) {
            // Creation of StandardFsDriver automatically verifies the
            // database path, so it does not need to be checked in here.
            DBDriver = new DatabaseDriver(GetFileSystemDriver(database.Path));
            Database = database;
        }

        /// <summary>
        /// The driver's database
        /// </summary>
        public IDatabaseInfo Database {
            get;
            private set;
        }

        /// <summary>
        /// IO DatabaseDriver for lower-level operations
        /// </summary>
        public IDatabaseDriver DBDriver {
            get;
            private set;
        }

        /// <summary>
        /// The number of sessions within this database. Does not require the
        /// actual sessions to be loaded.
        /// </summary>
        public int SessionCount {
            get {
                return DBDriver.FsDriver.GetAllSessionGUIDs().Count();
            }
        }

        /// <summary>
        /// The number of grids within this database. Does not require the
        /// actual grids to be loaded.
        /// </summary>
        public int GridCount {
            get {
                return DBDriver.FsDriver.GetAllGridGUIDs().Count();
            }
        }

        /// <summary>
        /// The sessions associated with this driver's database
        /// </summary>
        public IEnumerable<ISessionInfo> Sessions {
            get {
                HashSet<Guid> currentGuids = new HashSet<Guid>(
                    DBDriver.FsDriver.GetAllSessionGUIDs());

                // Remove obsolete entries
                Guid[] cachedGuids = m_Sessions.Keys.ToArray();
                foreach (Guid cachedGuid in cachedGuids) {
                    if (!currentGuids.Contains(cachedGuid)) {
                        m_Sessions.Remove(cachedGuid);
                    }
                }


                // Add new entries
                foreach (Guid guid in currentGuids) {
                    if (!m_Sessions.ContainsKey(guid)) {
                        m_Sessions[guid] = new SessionProxy(guid, this.Database);
                    }
                }

                return m_Sessions.Values;
            }
        }

        /// <summary>
        /// Cache for the database's session.
        /// </summary>
        private IDictionary<Guid, ISessionInfo> m_Sessions = new Dictionary<Guid, ISessionInfo>();

        /// <summary>
        /// Retrieves all sessions using a particular <paramref name="grid"/>.
        /// </summary>
        /// <param name="grid">The grid in question.</param>
        /// <returns>A collection of sessions using the grid.</returns>
        public IEnumerable<ISessionInfo> GetSessionInfos(IGridInfo grid) {
            IList<ISessionInfo> gridSessions = new List<ISessionInfo>();

            // Sweep all sessions to find the sessions using the grid
            foreach (ISessionInfo session in Sessions) {
                foreach (IGridInfo sessionGrid in GetGridInfos(session)) {
                    if (grid.ID.Equals(sessionGrid.ID)) {
                        gridSessions.Add(session);
                    }
                }
            }

            return gridSessions;
        }

        /// <summary>
        /// Obtains the session information associated with the given
        /// <paramref name="sessionID"/> through deserialization.
        /// </summary>
        /// <param name="sessionID">
        /// The id of the requested session
        /// </param>
        /// <returns>A lightweight representation of the session.</returns>
        public ISessionInfo GetSessionInfo(Guid sessionID) {
            if (!m_Sessions.ContainsKey(sessionID)) {
                m_Sessions.Add(sessionID, new SessionProxy(sessionID, Database));
            }
            return m_Sessions[sessionID];
        }

        /// <summary>
        /// Retrieves all session info objects of this database through
        /// deserialization. Strictly intended for populating
        /// <see cref="m_Sessions"/>.
        /// </summary>
        /// <returns>All session infos stored in the database.</returns>
        private IEnumerable<ISessionInfo> LoadSessionInfos() {
            return DBDriver.FsDriver.GetAllSessionGUIDs()
                .Select(sessID => GetSessionInfo(sessID));
        }

        /// <summary>
        /// Saves a session info object to a file on the disk.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        public void SaveSessionInfo(ISessionInfo session) {
            DBDriver.SaveSessionInfo(session);
        }

        private FileManager GetFileManager() {
            return new FileManager();
        }

        /// <summary>
        /// Deletes a session from the database and its files from the filesystem.
        /// </summary>
        /// <param name="session">The session to be deleted.</param>
        public void DeleteSession(ISessionInfo session) {
            FileManager fileManager = GetFileManager();

            // Lists of file and directory paths marked for deletion.
            // The delete operation will be executed at the very end, in
            // order to avoid reading from files that have already been removed.
            List<string> filesToDelete = new List<string>();
            IList<string> dirsToDelete = new List<string>();

            // Sessions subdirectory
            string sessionsSubDir = Path.Combine(session.Database.Path,
                StandardFsDriver.SessionsDir, session.ID.ToString());
            dirsToDelete.Add(sessionsSubDir);

            // try/catch block to delete a session with missing time-step log
            // (that happens when a session crashes before writing it's first timestep)
            try {
                // Timesteps
                foreach (string timestepFile in GetTimestepFiles(session)) {
                    filesToDelete.Add(timestepFile);
                }

                // Data subdirectory: distance data vectors
                foreach (ITimestepInfo tmstp in GetTimestepInfos(session)) {
                    List<string> storageVectorFiles = Utils.GetPathsFromGuid(
                        tmstp.StorageID,
                        Path.Combine(session.Database.Path, StandardFsDriver.DistVectorDataDir))
                        .ToList();
                    filesToDelete.AddRange(storageVectorFiles);
                }
#if DEBUG
                } catch (FileNotFoundException fnf) {
                    if (fnf.Message.Contains("TimestepLog.txt")) {
                        Console.WriteLine("No timestep log file found. Ignoring timestep files.");
                    }
                }
#else
            } catch (FileNotFoundException) {
                //Swallow
            }
#endif

            // Delete all the files marked for deletion
            foreach (string file in filesToDelete) {
                try {
                    fileManager.Delete(file);
                } catch (FileNotFoundException) {
                }
            }
            foreach (string dir in dirsToDelete) {
                try {
                    fileManager.DeleteDirectory(dir);
                } catch (FileNotFoundException) {
                }
            }
        }

        /// <summary>
        /// Copies a session to another database.
        /// </summary>
        /// <param name="session">The session to be copied.</param>
        /// <param name="dest">The destination database.</param>
        /// <returns>An info object for the session in the destination database.</returns>
        public ISessionInfo CopySession(ISessionInfo session, IDatabaseInfo dest) {
            return CopyOrMoveSession(session, dest, true);
        }

        /// <summary>
        /// Moves a session to another database.
        /// </summary>
        /// <param name="session">The session to be moved.</param>
        /// <param name="dest">The destination database.</param>
        /// <returns>An info object for the session in the destination database.</returns>
        public ISessionInfo MoveSession(ISessionInfo session, IDatabaseInfo dest) {
            return CopyOrMoveSession(session, dest, false);
        }

        /// <summary>
        /// Copies or moves a session to another database.
        /// </summary>
        /// <param name="session">Session to be moved.</param>
        /// <param name="dest">Destination database.</param>
        /// <param name="copy">true, if the session is to be copied; 
        /// false if the session should be moved.</param>
        /// <returns>An ISessionInfo object of the copied/moved session.</returns>
        private ISessionInfo CopyOrMoveSession(ISessionInfo session,
            IDatabaseInfo dest, bool copy) {
            // Lists of file paths marked for copy/move: First entry is source,
            // second entry is destination.
            // The operation will be executed at the very end, in order to avoid
            // reading from files that have already been moved
            Dictionary<string, string> filesToCopyOrMove =
                new Dictionary<string, string>();

            // Directories to be created (both copy/move) and deleted
            // (after move operation only)
            IList<string> dirsToCreate = new List<string>();
            IList<string> dirsToDelete = new List<string>();

            // sessions subfolder
            string sessionsFolderSource = Path.Combine(session.Database.Path,
                StandardFsDriver.SessionsDir, session.ID.ToString());
            string sessionsFolderDest = Path.Combine(dest.Path,
                StandardFsDriver.SessionsDir, session.ID.ToString());

            dirsToCreate.Add(sessionsFolderDest);

            // if in move mode, mark source directory for deletion
            // after all the files have been moved
            if (!copy) {
                dirsToDelete.Add(sessionsFolderSource);
            }

            foreach (string sessionFilePathSource in Directory.GetFiles(sessionsFolderSource)) {
                string fileName = Path.GetFileName(sessionFilePathSource);
                filesToCopyOrMove.Add(sessionFilePathSource,
                    Path.Combine(sessionsFolderDest, fileName));
            }

            // grids subdirectory
            foreach (IGridInfo grid in GetGridInfos(session)) {
                // Copy the grid in any case. Even if the session is moved there
                // may be other sessions in the source database using the grid.
                CopyGrid(grid, dest);
            }

            // timesteps subdirectory
            IEnumerable<Guid> timestepUids = DBDriver.GetTimestepGuids(session.ID);

            foreach (Guid timestepID in timestepUids) {
                string tmstpPathSource = Utils.GetPathsFromGuid(
                    timestepID,
                    Path.Combine(session.Database.Path, StandardFsDriver.TimestepDir),
                    "ts").Single(); // will throw exception if there is more than one .ti for this timestep
                string tmstpPathDest = Path.Combine(
                    dest.Path,
                    StandardFsDriver.TimestepDir,
                    Path.GetFileName(tmstpPathSource));

                filesToCopyOrMove.Add(tmstpPathSource, tmstpPathDest);
            }

            // data subdirectory: distance data vectors
            foreach (ITimestepInfo tmstp in GetTimestepInfos(session)) {
                IList<string> storageVectorSrcPaths = Utils.GetPathsFromGuid(
                    tmstp.StorageID,
                    Path.Combine(session.Database.Path, StandardFsDriver.DistVectorDataDir))
                    .ToList();

                IList<string> storageVectorDestPaths = storageVectorSrcPaths.
                    Select(srcPath =>
                        Path.Combine(dest.Path, StandardFsDriver.DistVectorDataDir,
                        Path.GetFileName(srcPath)))
                    .ToList();

                for (int i = 0; i < storageVectorSrcPaths.Count; i++) {
                    filesToCopyOrMove.Add(storageVectorSrcPaths[i],
                        storageVectorDestPaths[i]);
                }
            }

            // Execute the copy or move operations
            FileManager fileManager = GetFileManager();

            foreach (string dirToCreate in dirsToCreate) {
                fileManager.CreateDirectory(dirToCreate);
            }

            foreach (KeyValuePair<string, string> fileToCopyOrMove
                in filesToCopyOrMove) {

                if (copy) {
                    fileManager.Copy(fileToCopyOrMove.Key,
                        fileToCopyOrMove.Value, false);
                } else {
                    fileManager.Move(fileToCopyOrMove.Key,
                        fileToCopyOrMove.Value);
                }
            }

            // because there is no fileManager.Move(Directory) available, we
            // need to delete the original directories by hand after copying them
            // this especially applies to the subfolder of db_root/sessions
            foreach (string dirToDelete in dirsToDelete) {
                fileManager.DeleteDirectory(dirToDelete);
            }

            return session.CopyFor(dest);
        }

        /// <summary>
        /// Deletes a time-step from the database and its files from the file
        /// system.
        /// </summary>
        /// <param name="timestep">The time-step to be deleted.</param>
        public void DeleteTimestep(ITimestepInfo timestep) {
            FileManager fileManager = GetFileManager();
            string path = Path.Combine(
                Path.Combine(timestep.Database.Path, "timesteps"),
                timestep.ID + ".ts");

            List<string> filesToDelete = new List<string>();
            filesToDelete.Add(path);

            // Data subdirectory: data vectors
            string dataPath = Path.Combine(
                timestep.Database.Path, StandardFsDriver.DistVectorDataDir);
            filesToDelete.AddRange(Utils.GetPathsFromGuid(timestep.StorageID, dataPath));

            // Delete all the files marked for deletion
            DBDriver.RemoveTimestepGuid(timestep.Session.ID, timestep.ID);
            foreach (string file in filesToDelete) {
                try {
                    fileManager.Delete(file);
                } catch (FileNotFoundException) {
                }
            }
        }

        private IDictionary<Guid, IDictionary<Guid, ITimestepInfo>> m_Timesteps =
            new Dictionary<Guid, IDictionary<Guid, ITimestepInfo>>();

        /// <summary>
        /// Retrieves information about all time-steps in a session.
        /// </summary>
        /// <param name="session">The session in question.</param>
        /// <returns>All information about the session's time-steps.</returns>
        public IList<ITimestepInfo> GetTimestepInfos(ISessionInfo session) {
            if (!m_Timesteps.ContainsKey(session.ID)) {
                m_Timesteps.Add(session.ID, new Dictionary<Guid, ITimestepInfo>());
            }

            HashSet<Guid> currentGuids = new HashSet<Guid>(
                DBDriver.GetTimestepGuids(session.ID));

            // Remove obsolete entries
            Guid[] cachedGuids = m_Timesteps[session.ID].Keys.ToArray();
            foreach (Guid cachedGuid in cachedGuids) {
                if (!currentGuids.Contains(cachedGuid)) {
                    m_Timesteps[session.ID].Remove(cachedGuid);
                }
            }

            // Add new entries
            foreach (Guid timestepGuid in currentGuids) {
                if (!m_Timesteps[session.ID].ContainsKey(timestepGuid)) {
                    m_Timesteps[session.ID].Add(
                        timestepGuid, new TimestepProxy(timestepGuid, session));
                }
            }

            return m_Timesteps[session.ID].Values.ToList();
        }

        private Dictionary<Guid, IGridInfo> m_Grids = new Dictionary<Guid, IGridInfo>();

        /// <summary>
        /// A list of all grids within this database
        /// </summary>
        public IEnumerable<IGridInfo> Grids {
            get {
                HashSet<Guid> currentGuids = new HashSet<Guid>(
                    DBDriver.FsDriver.GetAllGridGUIDs());

                // Remove obsolete entries
                Guid[] cachedGuids = m_Grids.Keys.ToArray();
                foreach (Guid cachedGuid in cachedGuids) {
                    if (!currentGuids.Contains(cachedGuid)) {
                        m_Grids.Remove(cachedGuid);
                    }
                }

                // Add new entries
                foreach (Guid guid in currentGuids) {
                    if (!m_Grids.ContainsKey(guid)) {
                        m_Grids[guid] = GetGridInfo(guid);
                    }
                }

                return m_Grids.Values;
            }
        }

        /// <summary>
        /// Retrieves a grid information object by reading the grid file(s)
        /// stored in the specified database
        /// </summary>
        /// <param name="gridID">The Guid of the grid to be read.</param>
        /// <returns>
        /// The grid info object for the selected grid.
        /// </returns>
        public IGridInfo GetGridInfo(Guid gridID) {
            if (!m_Grids.ContainsKey(gridID)) {
                m_Grids.Add(gridID, new GridProxy(gridID, Database));
            }

            return m_Grids[gridID];
        }

        private Dictionary<Guid, GridInitializationContext> m_gridInitializationContexts
            = new Dictionary<Guid, GridInitializationContext>();


        // Saving the timestep-initialization contexts should not be necessay anymore,
        // since the DG basis is now a leightweigt object -- all expensive caching is done in 'ChefBasis',
        // see GridData.BasisData.
        // fk, 05sep16
        //private Dictionary<Guid, TimestepInitializationContext> m_timestepInitializationContexts
        //    = new Dictionary<Guid, TimestepInitializationContext>();

        /// <summary>
        /// Adds an initialization context for a particular grid object, i.e.
        /// such the grid data object and grid-global objects (like a basis)
        /// don't have to be created over and over again.
        /// </summary>
        /// <param name="gridData"></param>
        public void AddGridInitializationContext(Grid.Classic.GridData gridData) {
            if (m_gridInitializationContexts.ContainsKey(gridData.Grid.ID)) {
                throw new ArgumentException(
                    "An initialization context for the given grid already exists");
            }

            m_gridInitializationContexts[gridData.Grid.ID] = new GridInitializationContext(gridData);
        }

        /// <summary>
        /// Retrieves the initialization for a particular time-step.
        /// </summary>
        /// <param name="ts"></param>
        /// <returns></returns>
        public IInitializationContext GetInitializationContext(ITimestepInfo ts) {
            // Saving the timestep-initialization contexts should not be necessay anymore,
            // since the DG basis is now a leightweigt object -- all expensive caching is done in 'ChefBasis',
            // see GridData.BasisData.
            // fk, 05sep16
            
            //if (!m_timestepInitializationContexts.ContainsKey(ts.ID)) {
            if (!m_gridInitializationContexts.ContainsKey(ts.GridID)) {
                Grid.Classic.GridCommons grid = DBDriver.LoadGrid(ts.GridID, Database);
                Grid.Classic.GridData gridData = new Grid.Classic.GridData(grid);
                m_gridInitializationContexts[ts.GridID] =
                    new GridInitializationContext(gridData);
            }

            //m_timestepInitializationContexts[ts.ID] =
            //    new TimestepInitializationContext(m_gridInitializationContexts[ts.GridID]);
            //}
            //return m_timestepInitializationContexts[ts.ID];

            return new TimestepInitializationContext(m_gridInitializationContexts[ts.GridID]);
        }

        /// <summary>
        /// Saves a grid info object to the persistent memory.
        /// </summary>
        /// <param name="grid">The grid to be saved.</param>
        public void SaveGridInfo(IGridInfo grid) {
            if (grid is Grid.Classic.GridCommons) {
                DBDriver.SaveGrid((Grid.Classic.GridCommons)grid);
            } else if (grid is GridProxy) {
                DBDriver.SaveGrid(grid.Cast<GridProxy>().RealGrid.Cast<Grid.Classic.GridCommons>());
            } else {
                throw new NotSupportedException("As of now, only GridCommons objects can be saved.");
            }
        }

        /// <summary>
        /// Retrieves a grid information object by reading the grid file(s)
        /// associated with the specified <paramref name="session"/>.
        /// </summary>
        /// <param name="session">The session in question.</param>
        /// The collection of grids associated to <paramref name="session"/>.
        public IEnumerable<IGridInfo> GetGridInfos(ISessionInfo session) {
            List<IGridInfo> grids = new List<IGridInfo>();

            foreach (ITimestepInfo ti in session.Timesteps) {
                if (!grids.Contains(ti.Grid)) {
                    grids.Add(ti.Grid);
                }
            }

            return grids;
        }

        /// <summary>
        /// Deletes a grid from the database,
        /// </summary>
        /// <param name="grid">The grid to be deleted.</param>
        /// <param name="safelyDelete">
        /// true if should only be deleted if it is not in use by any session;
        /// false if it should be deleted regardless of its usage status.
        /// </param>
        /// <returns>true if the deletion was successful, false otherwise</returns>
        public bool DeleteGrid(IGridInfo grid, bool safelyDelete = true) {
            // Sweep all sessions to find out if the grid is used by any of them
            if (safelyDelete && GetSessionInfos(grid).Count() > 0) {
                // grid in use => no deletion
                return false;
            }

            // At this point, the grid is not in use (or safelyDelete is false)
            // and will be deleted.
            // collect all the files associated with the grid
            IEnumerable<string> filesToDelete = GetGridFiles(grid);

            // Delete the files
            FileManager fileManager = GetFileManager();
            foreach (string file in filesToDelete) {
                fileManager.Delete(file);
            }

            m_Grids.Remove(grid.ID);
            return true;
        }

        /// <summary>
        /// Copies the grid to a new database without overwriting,
        /// i.e. the copy operation is only performed 
        /// if the grid does not yet exist in the target database
        /// </summary>
        /// <param name="grid">The grid to copy</param>
        /// <param name="dest">The target database</param>
        /// <returns>
        /// An IGridInfo object referring to the grid in the target database
        /// </returns>
        public IGridInfo CopyGrid(IGridInfo gridInfo, IDatabaseInfo dest) {
            Grid.Classic.GridCommons grid;
            if (gridInfo is GridProxy) {
                grid = gridInfo.As<GridProxy>().RealGrid;
            } else if (gridInfo is Grid.Classic.GridCommons) {
                grid = (Grid.Classic.GridCommons)gridInfo;
            } else {
                throw new NotSupportedException();
            }

            FileManager fileManager = GetFileManager();

            IList<string> gridDirSrcFullPaths = Utils.GetPathsFromGuid(
                grid.ID,
                Path.Combine(grid.Database.Path, StandardFsDriver.GridsDir),
                "grid").ToList();
            IList<string> gridDirDestFullPaths = gridDirSrcFullPaths
                .Select(srcPath =>
                Path.Combine(dest.Path, StandardFsDriver.GridsDir, Path.GetFileName(srcPath)))
                .ToList();

            for (int i = 0; i < gridDirSrcFullPaths.Count; i++) {
                // copy only if grid file in grids subdirectory is
                // nonexistent in destination database
                if (!File.Exists(gridDirDestFullPaths[i])) {
                    List<Guid> dataGuids = new List<Guid>();

                    // load GridCommons object to retrieve the storage guid
                    // for the grid data
                    Guid gridStorageID = GetGridStorageID(grid.ID);
                    dataGuids.Add(gridStorageID);

                    // Don't forget optional custom partitionings!
                    foreach (var s in grid.m_PredefinedGridPartitioning) {
                        Guid partitioningGuid = s.Value.Guid;
                        dataGuids.Add(partitioningGuid);
                    }

                    // gather paths
                    string gridDataSrcBasePath =
                        Path.Combine(grid.Database.Path, StandardFsDriver.DistVectorDataDir);
                    string gridDataDestBasePath =
                        Path.Combine(dest.Path, StandardFsDriver.DistVectorDataDir);
                    IList<string> gridDataSrcFullPaths =
                        Utils.GetPathsFromGuids(dataGuids, gridDataSrcBasePath).ToList();
                    
                    // copy from grids subdirectory
                    fileManager.Copy(gridDirSrcFullPaths[i], gridDirDestFullPaths[i], false);

                    // copy data files
                    foreach (string dataFile in gridDataSrcFullPaths) {
                        fileManager.Copy(dataFile,
                            Path.Combine(gridDataDestBasePath,
                            Path.GetFileName(dataFile)), false);
                    }
                }
            }

            // return IGridInfo object of new grid in destination database
            return grid.CopyFor(dest);
        }

        /// <summary>
        /// Clears the entire database, leaving only the basic folder
        /// structure behind.
        /// </summary>
        public void ClearDatabase() {
            FileManager fileManager = GetFileManager();

            // cycle through all subfolders
            // i.e. "data", "grids", "sessions" and "timesteps"
            foreach (string subDirectoryPath
                in Directory.GetDirectories(Database.Path)) {

                //special case: sessions subdirectory
                if (Path.GetFileName(subDirectoryPath)
                    .Equals(StandardFsDriver.SessionsDir)) {
                    foreach (string sessionSubDirectory
                        in Directory.GetDirectories(subDirectoryPath)) {
                        fileManager.DeleteDirectory(sessionSubDirectory);
                    }
                } else {
                    foreach (string subDirectoryFile
                        in Directory.GetFiles(subDirectoryPath)) {
                        fileManager.Delete(subDirectoryFile);
                    }
                }
            }

            // Verify that the folder structure is still correct
            if (!IsDBFolderStructureValid(Database)) {
                throw new Exception("Database folder structure not valid anymore.");
            }
        }

        /// <summary>
        /// Disposes of all unused objects in the database.
        /// </summary>
        public void CleanDatabase() {
            // All files that are used by any session loaded from the database
            HashSet<string> requiredFiles = new HashSet<string>();

            foreach (ISessionInfo session in Sessions) {
                // Session timesteps
                IEnumerable<Guid> timestepUids = DBDriver.GetTimestepGuids(session.ID);
                try {
                    timestepUids = DBDriver.GetTimestepGuids(session.ID);
                } catch (IOException) {
                    // Timestep log not readable, session is probably useless anway
                    continue;
                }

                foreach (string timestepFile in GetTimestepFiles(session)) {
                    requiredFiles.Add(timestepFile);
                }

                // Data vector
                foreach (ITimestepInfo tmstp in GetTimestepInfos(session)) {
                    IEnumerable<string> storageVectorFiles = Utils.GetPathsFromGuid(
                        tmstp.StorageID,
                        Path.Combine(session.Database.Path, StandardFsDriver.DistVectorDataDir),
                        "data");
                    requiredFiles.AddRange(storageVectorFiles);
                }
            }

            // Don't delete any grids during clean-up: this means all files
            // of all grids have to be in allFilesInDatabase
            IEnumerable<Guid> gridIDs = DBDriver.FsDriver.GetAllGridGUIDs();

            foreach (Guid gridID in gridIDs) {
                requiredFiles.AddRange(GetGridFiles(gridID, Database));
            }

            // Gather all files and directories to be deleted.
            // [*] excludes the sessions subdirectory because it contains
            // folders, not files that will be removed.
            // [**] selects all files that are not in allFilesInDatabase
            IEnumerable<string> filesToDelete = Directory
                .GetFiles(Database.Path, "*.*", SearchOption.AllDirectories)
                .Where(f => !Path.GetDirectoryName(f).Contains(StandardFsDriver.SessionsDir)) //[*]
                .Where(f => !requiredFiles.Contains(f)); //[**]

            IEnumerable<string> dirsToDelete = Directory.
                GetDirectories(Path.Combine(Database.Path, StandardFsDriver.SessionsDir), "*",
                SearchOption.TopDirectoryOnly).
                Where(d => !Sessions.Any(d2 => d2.ID.ToString() == Path.GetFileName(d)));

            // Execute delete operations
            FileManager fileManager = GetFileManager();

            foreach (string fileToDelete in filesToDelete) {
                fileManager.Delete(fileToDelete);
            }

            foreach (string dirToDelete in dirsToDelete) {
                fileManager.DeleteDirectory(dirToDelete);
            }
        }

        /// <summary>
        /// Checks whether the basic folder structure of a database is existing
        /// </summary>
        private bool IsDBFolderStructureValid(IDatabaseInfo database) {
            try {
                StandardFsDriver.VerifyDirectoryStructure(database.Path);
            } catch (ArgumentException) {
                return false;
            }

            return true;
        }

        /// <summary>
        /// Returns the storage ID of a grid by deserializing the
        /// <see cref="BoSSS.Foundation.Grid.GridCommons"/> object.
        /// </summary>
        /// <param name="gridId">Id of the grid</param>
        /// <returns>The storage ID of the grid</returns>
        private Guid GetGridStorageID(Guid gridId) {
            Grid.Classic.GridCommons gridComm = DBDriver.LoadGrid(gridId, Database);
            return gridComm.StorageGuid;
        }

        /// <summary>
        /// Retrieves all files associated with a grid.
        /// </summary>
        /// <param name="grid">The info object of the grid</param>
        /// <returns>Paths to all files associated with <paramref name="grid"/></returns>
        public IEnumerable<string> GetGridFiles(IGridInfo grid) {
            return GetGridFiles(grid.ID, grid.Database);
        }

        /// <summary>
        /// Retrieves all files associated with a grid.
        /// </summary>
        /// <param name="gridID">The grid ID</param>
        /// <param name="database">The database where the grid is located</param>
        /// <returns>Paths to all files associated with the grid</returns>
        private IEnumerable<string> GetGridFiles(Guid gridID, IDatabaseInfo database) {
            // Add main file
            string gridMainFile = Utils.GetPathsFromGuid(gridID,
                Path.Combine(database.Path, StandardFsDriver.GridsDir), "grid")
                .Single();
            IList<string> gridFiles = new List<string> { gridMainFile };

            // Add data files
            Guid gridStorageID = GetGridStorageID(gridID);
            foreach (string gridDataFile in Directory.GetFiles(
                Path.Combine(database.Path, StandardFsDriver.DistVectorDataDir),
                gridStorageID.ToString() + ".*", SearchOption.AllDirectories)) {
                gridFiles.Add(gridDataFile);
            }

            return gridFiles.Distinct(); // remove duplicates
        }

        /// <summary>
        /// Returns a collection of paths to all files associated with a session's timesteps.
        /// </summary>
        private IEnumerable<string> GetTimestepFiles(ISessionInfo session) {
            IEnumerable<Guid> timestepUids = DBDriver.GetTimestepGuids(session.ID);
            return Utils.GetPathsFromGuids(
                timestepUids, Path.Combine(session.Database.Path, "timesteps"))
                .Where(pth => pth.EndsWith(".ts"));
        }
    }
}
