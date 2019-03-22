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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Linq;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using log4net.Appender;
using log4net.Config;
using log4net.Layout;
using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using log4net;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// The (formerly called) IO-master provides storage access on an
    /// object-level (in comparison, <see cref="IFileSystemDriver"/>-
    /// implementations provide IO on a stream-level), mainly by using
    /// serialization.
    /// </summary>
    public partial class DatabaseDriver : MPIProcess, IDatabaseDriver {

        IVectorDataSerializer standardVectorSerializer;
        GridDatabaseDriver gridDatabaseDriver;
        SessionDatabaseDriver sessionsDatabaseDriver;
        TimeStepDatabaseDriver timestepDatabaseDriver;
        IFileSystemDriver fsDriver; 

        /// <summary>
        /// </summary>
        /// <param name="fsDriver">
        /// the IO Driver; can be null;
        /// </param>
        public DatabaseDriver(IFileSystemDriver fsDriver)
        {
            this.fsDriver = fsDriver;
            ISerializer standardSerializer = new SerializerVersion0();
            ISerializer objectTypeSerializer = new SerializerVersion1();

            standardVectorSerializer = new VectorDataSerializer(fsDriver, standardSerializer);
            IVectorDataSerializer objectTypeVectorSerializer = new VectorDataSerializer(fsDriver, objectTypeSerializer);
            IVectorDataSerializer versionedVectorSerializer = new VersionManager(objectTypeVectorSerializer, standardVectorSerializer);

            gridDatabaseDriver = new GridDatabaseDriver(versionedVectorSerializer, fsDriver);
            sessionsDatabaseDriver = new SessionDatabaseDriver(standardSerializer, fsDriver);
            timestepDatabaseDriver = new TimeStepDatabaseDriver(standardVectorSerializer, fsDriver);
        }

        /// <summary>
        /// need to close on dispose
        /// </summary>
        TextWriterAppender logger_output = null;

        public void Dispose()
        {
            if (this.FsDriver is IDisposable)
            {
                ((IDisposable)this.FsDriver).Dispose();

            }
            sessionsDatabaseDriver.Dispose();

            if (logger_output != null)
                logger_output.Close();
        }

        /// <summary>
        /// Returns a write-stream for some new log file.
        /// </summary>
        public Stream GetNewLogStream(SessionInfo si, string name) {
            var id = si.ID;

            int Rank = MyRank;

            Stream file = null;
            if (fsDriver != null && id != Guid.Empty)
                file = this.FsDriver.GetNewLogStream(name + "." + Rank, id);

            if (file == null) {
                // create trace file in local directory

                string tracefilename = name + "." + Rank + ".txt";
                file = new FileStream( tracefilename, FileMode.Create, FileAccess.Write, FileShare.Read);
            }

            return file;
        }

        /// <summary>
        /// some fucking fake.
        /// </summary>
        static bool configAllreadyDone = false;

        /// <summary>
        /// the file system driver
        /// </summary>
        public IFileSystemDriver FsDriver => fsDriver;

        /// <summary>
        /// Tracing setup.
        /// </summary>
        public void InitTraceFile(SessionInfo si) {

            if (logger_output != null)
                throw new ApplicationException("Already called."); // is seems this object is designed so that it stores at max one session per lifetime

            var id = si.ID;

            // create tracer
            // =============

            int Rank = MyRank;

            Stream tracerfile = null;
            if (fsDriver != null && id != Guid.Empty)
                tracerfile = this.FsDriver.GetNewLogStream("trace." + Rank, id);

            TextWriter tracertxt = null;
            if (tracerfile == null && Tracer.NamespacesToLog.Length > 0 && !configAllreadyDone) {
                // create trace file in local directory

                string tracefilename = "trace." + Rank + ".txt";
                tracerfile = new FileStream(tracefilename, FileMode.Create, FileAccess.Write, FileShare.Read);
            }

            if (tracerfile != null) {
                var zipper = new System.IO.Compression.GZipStream(tracerfile, System.IO.Compression.CompressionMode.Compress);
                tracertxt = new StreamWriter(zipper);
            }

            if (tracertxt != null) {
                TextWriterAppender fa = new TextWriterAppender();
                fa.ImmediateFlush = true;
                fa.Writer = tracertxt;
                fa.Layout = new PatternLayout("%date %-5level %logger: %message%newline");
                fa.ActivateOptions();
                BasicConfigurator.Configure(fa);
                logger_output = fa;
                configAllreadyDone = true;
            }

        }

        /// <summary>
        /// Tracing setup.
        /// </summary>
        void CloseTraceFile() {
            if (logger_output != null)
                logger_output.Close();
        }

        /// <summary>
        /// Creates a new session;
        /// </summary>
        public SessionInfo CreateNewSession(IDatabaseInfo database)
        {
            return sessionsDatabaseDriver.CreateNewSession(database);
        }

        /// <summary>
        /// Saves a vector to the database, allocating a GUID automatically.
        /// </summary>
        /// <param name="vector"></param>
        /// <returns>
        /// The GUID that was allocated to identify the vector within the storage system.
        /// </returns>
        public Guid SaveVector<T>(IList<T> vector)
        {
            return standardVectorSerializer.SaveVector(vector);
        }

        /// <summary>
        /// saves a vector to the database, under a specified Guid
        /// </summary>
        /// <param name="vector">
        /// the part of the vector which is stored on the local process.
        /// </param>
        /// <param name="id">
        /// the Guid under which the vector should be stored; must be the same
        /// on all MPI processes.
        /// </param>
        public void SaveVector<T>(IList<T> vector, Guid id)
        {
            standardVectorSerializer.SaveVector(vector,id);
        }

        /// <summary>
        /// Loads a vector from the database
        /// </summary>
        /// <param name="id"></param>
        /// <param name="part">
        /// Optional partition of the vector among MPI processors: if null, a
        /// partition is defined by the loader logic.
        /// </param>
        public IList<T> LoadVector<T>(Guid id, ref Partitioning part)
        {
            return standardVectorSerializer.LoadVector<T>(id, ref part);
        }

        /// <summary>
        /// tests whether a grid with GUID <paramref name="g"/> exists in database, or not;
        /// </summary>
        public bool GridExists(Guid g)
        {
            return gridDatabaseDriver.GridExists(g);
        }

        /// <summary>
        /// Loads the grid info object for the given
        /// <paramref name="gridId"/> from the given
        /// <param name="database"></param>
        /// <returns></returns>
        public IGridInfo LoadGridInfo(Guid gridId, IDatabaseInfo database)
        {
            return gridDatabaseDriver.LoadGridInfo(gridId, database);
        }

        /// <summary>
        /// loads the grid identified by <paramref name="uid"/> from the
        /// given <paramref name="database"/>
        /// </summary>
        /// <param name="gridId">The unique identifier of the grid.</param>
        /// <param name="database">
        /// The database that is associated with the grid.
        /// </param>
        /// <returns>
        /// The loaded grid
        /// </returns>
        public IGrid LoadGrid(Guid gridId, IDatabaseInfo database)
        {
            return gridDatabaseDriver.LoadGrid(gridId, database);
        }

        /// <summary>
        /// Loads the actual grid data for the given <paramref name="grid"/>.
        /// That is, loads the actual cell data.
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public IGrid LoadGridData(IGrid grid)
        {
            return gridDatabaseDriver.LoadGridData(grid);
        }

        /// <summary>
        /// Loads the given <paramref name="sessionId"/> from the given
        /// <paramref name="database"/>.
        /// </summary>
        /// <param name="sessionId"></param>
        /// <param name="database"></param>
        /// <returns></returns>
        public SessionInfo LoadSession(Guid sessionId, IDatabaseInfo database)
        {
            return sessionsDatabaseDriver.LoadSession(sessionId, database);
        }

        /// <summary>
        /// Retrieves the directory where the files for the selected
        /// <paramref name="session"/> are stored.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static string GetSessionDirectory(ISessionInfo session) {
            return SessionDatabaseDriver.GetSessionDirectory(session);
        }

        /// <summary>
        /// Saves a session info object to a file on the disk.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        public void SaveSessionInfo(ISessionInfo session)
        {
            sessionsDatabaseDriver.SaveSessionInfo(session);
        }

        /// <summary>
        /// Searches for an equivalent grid in the database and, if none is found
        /// saves a grid object to the database.
        /// </summary>
        /// <param name="grid">
        /// On entry, the grid which should be saved to the database.
        /// On exit, either unchanged, or the equivalent grid.
        /// </param>
        /// <param name="EquivalentGridFound">
        /// Inidicates that an equivalent grid was found.
        /// </param>
        /// <param name="database"></param>
        public Guid SaveGridIfUnique(ref IGrid grid, out bool EquivalentGridFound, IDatabaseInfo database)
        {
            return gridDatabaseDriver.SaveGridIfUnique(ref grid, out EquivalentGridFound, database);
        }

        /// <summary>
        /// Saves the given grid object to the database;
        /// </summary>
        /// <returns>
        /// the Guid of the <see cref="IGrid"/>-object that was saved
        /// (equal to the <see cref="IDatabaseEntityInfo{T}.ID"/>-property).
        /// </returns>
        /// <param name="grid">
        /// The grid to save.
        /// </param>
        /// <param name="database">
        /// chaos
        /// </param>
        public Guid SaveGrid(IGrid grid, IDatabaseInfo database)
        {
            return gridDatabaseDriver.SaveGrid(grid, database);
        }

        /// <summary>
        /// loads a single <see cref="TimestepInfo"/>-object from the database.
        /// </summary>
        public TimestepInfo LoadTimestepInfo(Guid timestepGuid, ISessionInfo session, IDatabaseInfo database) {
            return timestepDatabaseDriver.LoadTimestepInfo<TimestepInfo>(timestepGuid, session, database);
        }

        /// <summary>
        /// loads a single <see cref="TimestepInfo"/>-object from the database.
        /// </summary>
        public T LoadTimestepInfo<T>(Guid timestepGuid, ISessionInfo session, IDatabaseInfo database)
            where T : TimestepInfo
        {
            return timestepDatabaseDriver.LoadTimestepInfo<T>(timestepGuid, session, database);
        }

        /// <summary>
        /// Gathers all time-step IDs of a session.
        /// </summary>
        /// <param name="sessionGuid">ID of the session.</param>
        /// <returns>A collection of th session's timestep IDs.</returns>
        public IEnumerable<Guid> GetTimestepGuids(Guid sessionGuid)
        {
            return timestepDatabaseDriver.GetTimestepGuids(sessionGuid);
        }

        /// <summary>
        /// Removes the given <paramref name="timestepGuid"/> from the
        /// time-step log for the given <paramref name="sessionGuid"/>
        /// </summary>
        /// <param name="sessionGuid"></param>
        /// <param name="timestepGuid"></param>
        public void RemoveTimestepGuid(Guid sessionGuid, Guid timestepGuid)
        {
            timestepDatabaseDriver.RemoveTimestepGuid(sessionGuid, timestepGuid);
        }

        /// <summary>
        /// Loads a time-step from the database into previously allocated
        /// DG-fields (<paramref name="PreAllocatedFields"/>).
        /// </summary>
        public void LoadFieldData(ITimestepInfo info, IGridData grdDat, IEnumerable<DGField> PreAllocatedFields)
        {
            timestepDatabaseDriver.LoadFieldData(info, grdDat, PreAllocatedFields);
        }

        /// <summary>
        /// Loads a time-step from the database.
        /// </summary>
        /// <remarks>
        /// By using this method, it is ensured that the loaded/returned fields
        /// have the same DG polynomial degree as in the database.
        /// </remarks>
        public IEnumerable<DGField> LoadFields(ITimestepInfo info, IGridData grdDat, IEnumerable<string> NameFilter = null)
        {
            return timestepDatabaseDriver.LoadFields(info, grdDat, NameFilter);
        }

        /// <summary>
        /// Saves a time-step to the database's persistent memory.
        /// </summary>
        /// <param name="_tsi">Contains Id etc.</param>
        public void SaveTimestep(TimestepInfo _tsi)
        {
            timestepDatabaseDriver.SaveTimestep( _tsi);
        }
    }
}
