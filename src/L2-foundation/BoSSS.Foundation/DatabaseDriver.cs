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
    public partial class DatabaseDriver : IDatabaseDriver {

        /// <summary>
        /// Indicates whether the content of the database should be serialized
        /// in a human-readable (debugging) format, or in a significantly
        /// smaller binary format
        /// </summary>
        private static readonly bool DebugSerialization = false;

        /// <summary>
        /// the file system driver
        /// </summary>
        public IFileSystemDriver FsDriver {
            get {
                return m_fsDriver;
            }
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
            string path = Path.Combine(
                session.Database.Path,
                StandardFsDriver.SessionsDir,
                session.ID.ToString());
            return path;
        }


        /// <summary>
        /// MPI rank of actual process within the MPI world communicator
        /// </summary>
        public int MyRank {
            get {
                int rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                return rank;
            }
        }

        /// <summary>
        /// Number of MPI processes within the MPI world communicator
        /// </summary>
        public int Size {
            get {
                int size;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                return size;
            }
        }

        /// <summary>
        /// need to close on dispose
        /// </summary>
        TextWriterAppender logger_output = null;

        /// <summary>
        /// Creates a new session;
        /// </summary>
        public SessionInfo CreateNewSession(IDatabaseInfo database) {

            Guid id = Guid.NewGuid();
            if (m_fsDriver == null)
                id = Guid.Empty;
            id = id.MPIBroadcast(0);

            // init driver
            // ===========
            if (m_fsDriver != null) {

                if (MyRank == 0)
                    m_fsDriver.CreateSessionDirectory(id);

                // ensure that the session directory is available, before ANY mpi process continues.
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                System.Threading.Thread.Sleep(1000);
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }

            SessionInfo si = new SessionInfo(id, database);

            // create copy of stdout and stderr
            // ================================
            if (this.FsDriver != null) {
                m_stdout = this.FsDriver.GetNewLog("stdout." + MyRank, id);
                m_stderr = this.FsDriver.GetNewLog("stderr." + MyRank, id);

                ilPSP.Environment.StdOut.WriterS.Add(m_stdout);
                ilPSP.Environment.StdErr.WriterS.Add(m_stderr);
            }
            return si;
        }

        TextWriter m_stdout = null;
        TextWriter m_stderr = null;

        /// <summary>
        /// Returns a write-stream for some new log file.
        /// </summary>
        public Stream GetNewLogStream(SessionInfo si, string name) {
            var id = si.ID;

            int Rank = MyRank;

            Stream file = null;
            if (m_fsDriver != null && id != Guid.Empty)
                file = this.FsDriver.GetNewLogStream(name + "." + Rank, id);

            if (file == null) {
                // create trace file in local directory

                string tracefilename = name + "." + Rank + ".txt";
                file = new FileStream(tracefilename, FileMode.Create, FileAccess.Write, FileShare.Read);
            }

            return file;
        }

        /// <summary>
        /// some fucking fake.
        /// </summary>
        static bool configAllreadyDone = false;

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
            if (m_fsDriver != null && id != Guid.Empty)
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
        /// 
        /// </summary>
        /// <param name="driver">
        /// the IO Driver; can be null;
        /// </param>
        public DatabaseDriver(IFileSystemDriver driver) {
            m_fsDriver = driver;
        }

        /// <summary>
        /// serialization formatter used for all bigger (data) objects
        /// </summary>
        JsonSerializer m_Formatter = new JsonSerializer() {
            NullValueHandling = NullValueHandling.Ignore,
            TypeNameHandling = TypeNameHandling.Auto,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
            Binder = new MySerializationBinder()
        };

        class MySerializationBinder : Newtonsoft.Json.Serialization.DefaultSerializationBinder {

            public override Type BindToType(string assemblyName, string typeName) {

                if (assemblyName.Equals("BoSSS.Foundation") && typeName.Equals("BoSSS.Foundation.Grid.Cell[]")) {
                    typeName = "BoSSS.Foundation.Grid.Classic.Cell[]";
                }

                if (assemblyName.Equals("BoSSS.Foundation") && typeName.Equals("BoSSS.Foundation.Grid.BCElement[]")) {
                    typeName = "BoSSS.Foundation.Grid.Classic.BCElement[]";
                }


                Type T = base.BindToType(assemblyName, typeName);
                return T;
            }

        }

        /// <summary>
        /// 
        /// </summary>
        IFileSystemDriver m_fsDriver;

        /// <summary>
        /// header file for distributed stored vectors
        /// </summary>
        [Serializable]
        public class DistributedVectorHeader {

            /// <summary>
            /// GUID to identify the vector in the storage system.
            /// </summary>
            public Guid m_Guid;

            /// <summary>
            /// partition of the vector;
            /// the i-th entry is the first index of the vector that is stored in the i-th part.
            /// the length of this array is equal to the number of parts plus 1, the last entry
            /// is the total length of the vector.
            /// </summary>
            public long[] Partitioning;

            /// <summary>
            /// A hack, which indicates that we are using the most modern version.
            /// </summary>
            public bool UseGzip;
        }

        /// <summary>
        /// Saves a vector to the database, allocating a GUID automatically.
        /// </summary>
        /// <param name="vector"></param>
        /// <returns>
        /// The GUID that was allocated to identify the vector within the storage system.
        /// </returns>
        public Guid SaveVector<T>(IList<T> vector) {
            // allocate GUID
            // =============
            Guid id;
            id = Guid.NewGuid();
            id = id.MPIBroadcast(0);

            SaveVector(vector, id);

            return id;
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
        public void SaveVector<T>(IList<T> vector, Guid id) {
            using(new FuncTrace()) {
                if(m_fsDriver == null)
                    throw new NotSupportedException("Can't save data when not connected to database.");

#if DEBUG
                // never trust the user
                Guid id0 = id.MPIBroadcast(0);
                if(id0 != id) {
                    throw new ArgumentException(
                        "Guid differs at least on one MPI process.",
                        nameof(id));
                }
#endif
                Exception e = null;
                try {

                    // save parts
                    using(Stream s = m_fsDriver.GetDistVectorDataStream(true, id, MyRank)) {
                        using(var s2 = GetJsonWriter(new GZipStream(s, CompressionMode.Compress))) {
                            // Use a tuple since the Json format expects one object per
                            // file (with some tricks, we could avoid this, but it's
                            // not really worth the effort)
                            var tuple = new Tuple<DistributedVectorHeader, IList<T>>(null, vector);

                            // save header (only proc 0)
                            // =========================
                            var _part = (new Partitioning(vector.Count));
                            if(MyRank == 0) {
                                DistributedVectorHeader header = new DistributedVectorHeader();
                                header.m_Guid = id;
                                header.UseGzip = !DebugSerialization;
                                long[] part = new long[_part.MpiSize + 1];
                                for(int i = 0; i < _part.MpiSize; i++) {
                                    part[i + 1] = part[i] + _part.GetLocalLength(i);
                                }
                                header.Partitioning = part;
                                tuple = new Tuple<DistributedVectorHeader, IList<T>>(header, vector);
                            }

                            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                            m_Formatter.Serialize(s2, tuple);

                            s2.Close();
                            s.Close();
                        }
                    }
                } catch(Exception ee) {
                    Console.Error.WriteLine(ee.GetType().Name + " on rank " + this.MyRank + " saving vector " + id + ": " + ee.Message);
                    Console.Error.WriteLine(ee.StackTrace);
                    e = ee;
                }

                e.ExceptionBcast();
            }
        }

        /// <summary>
        /// Loads a vector from the database
        /// </summary>
        /// <param name="id"></param>
        /// <param name="part">
        /// Optional partition of the vector among MPI processors: if null, a
        /// partition is defined by the loader logic.
        /// </param>
        public IList<T> LoadVector<T>(Guid id, ref Partitioning part) {
            using (new FuncTrace()) {

                // load header
                // -----------
                DistributedVectorHeader header = null;

                // load on p0
                Exception e = null;
                try {
                    if (MyRank == 0) {
                        // load header
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, 0))
                        using (var reader = GetJsonReader(new GZipStream(s, CompressionMode.Decompress))) {
                            header = m_Formatter.Deserialize<Tuple<DistributedVectorHeader, IList<T>>>(reader).Item1;
                            reader.Close();
                            s.Close();
                        }

                    }
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                // broadcast
                header = header.MPIBroadcast(0);

                // define partition, if necessary
                if (part == null) {
                    int L = (int)header.Partitioning.Last();
                    int i0 = (L * this.MyRank) / this.Size;
                    int ie = (L * (this.MyRank + 1)) / this.Size;
                    part = new Partitioning(ie - i0);
                }


                // check size
                if (part.TotalLength != header.Partitioning.Last())
                    throw new ArgumentException("IO error: length of vector to load differs from total length of list.");


                // find first part
                // ===============

                long myI0;
                int myLen;
                myI0 = part.i0;
                myLen = part.LocalLength;

                int p = 0;
                while (header.Partitioning[p] <= myI0) {
                    p++;
                }
                p--; // now p is the index of the first vector part

                // load parts
                // ==========
                List<T> ret = new List<T>(part.LocalLength);
                try {
                    int d = 0;
                    do {
                        IList<T> vecP;
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, p))
                        using (var reader = GetJsonReader(new GZipStream(s, CompressionMode.Decompress))) {
                            vecP = m_Formatter.Deserialize<Tuple<DistributedVectorHeader, IList<T>>>(reader).Item2;
                            reader.Close();
                            s.Close();
                        }

                        int srdInd = (int)((myI0 + d) - header.Partitioning[p]);
                        int CopyLen = Math.Min(vecP.Count - srdInd, myLen - d);

                        for (int i = 0; i < CopyLen; i++) {
                            ret.Add(vecP[srdInd + i]);
                            d++;
                        }

                        p++; // next part;

                    } while (d < myLen);
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                // return
                // ======
                return ret;
            }
        }

        /// <summary>
        /// Loads a vector from the database
        /// </summary>
        /// <param name="id"></param>
        /// <param name="part">
        /// Optional partition of the vector among MPI processors: if null, a
        /// partition is defined by the loader logic.
        /// </param>
        public IList<T> LoadVectorSpeziale<T>(Guid id, ref Partitioning part, out Partitioning verySpecial) {
            using (new FuncTrace()) {

                // load header
                // -----------
                DistributedVectorHeader header = null;

                // load on p0
                Exception e = null;
                try {
                    if (MyRank == 0) {
                        // load header
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, 0))
                        using (var reader = GetJsonReader(new GZipStream(s, CompressionMode.Decompress))) {
                            header = m_Formatter.Deserialize<Tuple<DistributedVectorHeader, IList<T>>>(reader).Item1;
                            reader.Close();
                            s.Close();
                        }

                    }
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                // broadcast
                header = header.MPIBroadcast(0);

                // define partition, if necessary
                if (part == null) {
                    int L = (int)header.Partitioning.Last();
                    int i0 = (L * this.MyRank) / this.Size;
                    int ie = (L * (this.MyRank + 1)) / this.Size;
                    part = new Partitioning(ie - i0);
                }

                verySpecial = new Partitioning((int)(header.Partitioning[this.MyRank + 1] - header.Partitioning[this.MyRank]));

                // check size
                if (part.TotalLength != header.Partitioning.Last())
                    throw new ArgumentException("IO error: length of vector to load differs from total length of list.");


                // find first part
                // ===============

                long myI0;
                int myLen;
                myI0 = 0;
                myLen = (int)header.Partitioning.Last();

                int p = 0;
                while (header.Partitioning[p] <= myI0) {
                    p++;
                }
                p--; // now p is the index of the first vector part

                // load parts
                // ==========
                List<T> ret = new List<T>(part.LocalLength);
                try {
                    int d = 0;
                    do {
                        IList<T> vecP;
                        using (var s = m_fsDriver.GetDistVectorDataStream(false, id, p))
                        using (var reader = GetJsonReader(new GZipStream(s, CompressionMode.Decompress))) {
                            vecP = m_Formatter.Deserialize<Tuple<DistributedVectorHeader, IList<T>>>(reader).Item2;
                            reader.Close();
                            s.Close();
                        }

                        int srdInd = (int)((myI0 + d) - header.Partitioning[p]);
                        int CopyLen = Math.Min(vecP.Count - srdInd, myLen - d);

                        for (int i = 0; i < CopyLen; i++) {
                            ret.Add(vecP[srdInd + i]);
                            d++;
                        }

                        p++; // next part;

                    } while (d < myLen);
                } catch (Exception ee) {
                    e = ee;
                }
                e.ExceptionBcast();

                // return
                // ======
                return ret;
            }
        }

        /// <summary>
        /// Saves a session info object to a file on the disk.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        public void SaveSessionInfo(ISessionInfo session) {
            using (Stream s = FsDriver.GetSessionInfoStream(true, session.ID))
            using (var writer = GetJsonWriter(s)) {
                m_Formatter.Serialize(writer, session);
                writer.Close();
                s.Close();
            }
        }

        /// <summary>
        /// finalize
        /// </summary>
        public void Dispose() {
            if (this.FsDriver is IDisposable) {
                ((IDisposable)this.FsDriver).Dispose();

                if (m_stdout != null) {
                    Debug.Assert(ilPSP.Environment.StdOut.WriterS.Contains(m_stdout));
                    Debug.Assert(ilPSP.Environment.StdErr.WriterS.Contains(m_stderr));

                    Console.Out.Flush();
                    Console.Error.Flush();

                    ilPSP.Environment.StdOut.WriterS.Remove(m_stdout);
                    ilPSP.Environment.StdErr.WriterS.Remove(m_stderr);

                    m_stderr.Close();
                    m_stdout.Close();
                    m_stderr.Dispose();
                    m_stdout.Dispose();
                }

                if (logger_output != null)
                    logger_output.Close();
            }
        }

        /// <summary>
        /// Saves a time-step to the database's persistent memory.
        /// </summary>
        /// <param name="_tsi">
        /// (input/output)
        /// A time-step-info object that has not yet been saved in the database;
        /// object state (e.g. <see cref="ITimestepInfo.StorageID"/>) will be modified on output;
        /// </param>
        public void SaveTimestep(TimestepInfo _tsi) {
            using (var tr = new FuncTrace()) {

                if (!(_tsi.ID.Equals(Guid.Empty) && _tsi.StorageID.Equals(Guid.Empty)))
                    throw new ArgumentException("Timestep is already saved in database");
                var fields = _tsi.Fields.ToArray();
                var GridDat = fields[0].GridDat;

                {
                    List<DGField> FieldsFlatten = new List<DGField>();
                    TimestepInfo.FlattenHierarchy(FieldsFlatten, fields);
                    foreach (var f in FieldsFlatten) {
                        if (!object.ReferenceEquals(f.GridDat, GridDat))
                            throw new ArgumentException("mismatch in GridData object.");

                        if (!fields.Contains(f, (a, b) => object.ReferenceEquals(a, b))) {
                            // here, we ensure that the 'fields' -- list is complete, i.e.
                            // that the flatten hierarchy contains no field which is not already a memeber of 'fields'.
                            // The purpose is e.g. to prevent saving an XDG field without the required level-set field.
                            throw new ArgumentException(
                                "Unable to save timestep: field '" + f.Identification
                                    + "', which is required by at least one of the"
                                    + " given fields, must also be contained in the"
                                    + " given list of fields.",
                                "_tsi");
                        }
                    }
                }

                // build vector
                // ============
                int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                var vec = new CellFieldDataSet[J];
                var _fields = fields.ToArray();
                int NF = _fields.Length;
                var Permutation = GridDat.CurrentGlobalIdPermutation.Values;
                for (int j = 0; j < J; j++) { // loop over cells
                    vec[j] = new CellFieldDataSet();
                    vec[j].GlobalID = Permutation[j];
                    //vec[j].DGCoordinateData = new CellFieldDataSet.CellFieldData[NF];
                    for (int idxF = 0; idxF < NF; idxF++) { // loop over fields
                        var field = _fields[idxF];
                        int N = field.Basis.GetLength(j);
                        double[] Coords = new double[N];
                        for (int n = 0; n < N; n++) {
                            Coords[n] = field.Coordinates[j, n];
                        }
                        //vec[j].DGCoordinateData[idxF] = new CellFieldDataSet.CellFieldData() {
                        //    Data = Coords
                        //};
                        vec[j].AppendDGCoordinates(Coords);
                        Debug.Assert(ArrayTools.ListEquals(Coords, vec[j].GetDGCoordinates(idxF)));
                    }
                }

                // Save dg coordinates
                // ===================
                Guid VectorGuid = SaveVector(vec);
                _tsi.StorageID = VectorGuid;
                

                // Save state object 
                // =================
                _tsi.ID = Guid.NewGuid().MPIBroadcast(0);
                Exception e = null;
                if (MyRank == 0) {
                    try {
                        //tsi = new TimestepInfo(physTime, currentSession, TimestepNo, fields, VectorGuid);
                        using(var s = FsDriver.GetTimestepStream(true, _tsi.ID))
                        using(var writer = GetJsonWriter(s)) {
                            m_Formatter.Serialize(writer, _tsi);
                            writer.Close();
                            s.Close();
                        }
                    } catch (Exception ee) {
                        e = ee;
                        Console.Error.WriteLine(ee.GetType().Name + " on rank " + MyRank + " saving time-step " + _tsi.TimeStepNumber + ": " + ee.Message);
                        Console.Error.WriteLine(ee.StackTrace);
                    }
                }
                e.ExceptionBcast();

                // log session
                // ===========
                SessionInfo currentSession = (SessionInfo)( _tsi.Session); // hack
                if (MyRank == 0) {

                    try {
                        currentSession.LogTimeStep(_tsi.ID);
                    } catch (Exception ee) {
                        e = ee;
                        Console.Error.WriteLine(ee.GetType().Name + " on rank " + MyRank + " saving time-step " + _tsi.TimeStepNumber + ": " + ee.Message);
                        Console.Error.WriteLine(ee.StackTrace);
                    }
                }
                e.ExceptionBcast();

                _tsi.Database = currentSession.Database;
            }
        }

        /// <summary>
        /// Loads the given <paramref name="sessionId"/> from the given
        /// <paramref name="database"/>.
        /// </summary>
        /// <param name="sessionId"></param>
        /// <param name="database"></param>
        /// <returns></returns>
        public SessionInfo LoadSession(Guid sessionId, IDatabaseInfo database) {
            using (var tr = new FuncTrace()) {
                tr.Info("Loading session " + sessionId);

                using (Stream s = FsDriver.GetSessionInfoStream(false, sessionId))
                using (var reader = GetJsonReader(s)) {
                    SessionInfo loadedSession = m_Formatter.Deserialize<SessionInfo>(reader);
                    loadedSession.Database = database;
                    loadedSession.WriteTime = Utils.GetSessionFileWriteTime(loadedSession);

                    reader.Close();
                    s.Close();

                    return loadedSession;
                }
            }
        }

        /// <summary>
        /// loads a single <see cref="TimestepInfo"/>-object from the database.
        /// </summary>
        public TimestepInfo LoadTimestepInfo(Guid timestepGuid, ISessionInfo session, IDatabaseInfo database) {
            using (var tr = new FuncTrace()) {
                tr.Info("Loading time-step " + timestepGuid);

                TimestepInfo tsi = null;
                if (MyRank == 0) {
                    using (Stream s = FsDriver.GetTimestepStream(false, timestepGuid))
                    using (var reader = GetJsonReader(s)) {
                        tsi = m_Formatter.Deserialize<TimestepInfo>(reader);
                        tsi.Session = session;
                        reader.Close();
                        s.Close();
                    }
                    tsi.ID = timestepGuid;
                }
                tsi = tsi.MPIBroadcast(0);
                tsi.Database = database;
                tsi.WriteTime = Utils.GetTimestepFileWriteTime(tsi);
                return tsi;
            }
        }

        /// <summary>
        /// Gathers all time-step IDs of a session.
        /// </summary>
        /// <param name="sessionGuid">ID of the session.</param>
        /// <returns>A collection of th session's timestep IDs.</returns>
        public IEnumerable<Guid> GetTimestepGuids(Guid sessionGuid) {
            IList<Guid> timestepUids = new List<Guid>();

            try {
                using (StreamReader timestepLogReader =
                new StreamReader(FsDriver.GetTimestepLogStream(sessionGuid))) {

                    while (!timestepLogReader.EndOfStream) {
                        timestepUids.Add(Guid.Parse(timestepLogReader.ReadLine()));
                    }
                }
            } catch (FileNotFoundException) {
                return new Guid[0];
            }

            return timestepUids;
        }

        /// <summary>
        /// Removes the given <paramref name="timestepGuid"/> from the
        /// time-step log for the given <paramref name="sessionGuid"/>
        /// </summary>
        /// <param name="sessionGuid"></param>
        /// <param name="timestepGuid"></param>
        public void RemoveTimestepGuid(Guid sessionGuid, Guid timestepGuid) {
            string logPath = FsDriver.GetTimestepLogPath(sessionGuid);
            string[] lines = File.ReadAllLines(logPath);

            bool match = false;
            List<string> reducedLines = new List<string>();
            foreach (string line in lines) {
                Guid guid = Guid.Parse(line);
                if (guid.Equals(timestepGuid)) {
                    match = true;
                } else {
                    reducedLines.Add(line);
                }
            }

            if (!match) {
                throw new IOException(String.Format(
                    "Time-step guid {0} was not present in the time-step log for session {1}",
                    timestepGuid,
                    sessionGuid));
            }

            File.WriteAllLines(logPath, reducedLines);
        }

        /// <summary>
        /// Loads a time-step from the database into previously allocated
        /// DG-fields (<paramref name="PreAllocatedFields"/>).
        /// </summary>
        public void LoadFieldData(ITimestepInfo info, IGridData grdDat, IEnumerable<DGField> PreAllocatedFields) {
            using (var tr = new FuncTrace()) {
                DGField[] Fields = PreAllocatedFields.ToArray(); // enforce 'evaluation' of the enum (in the case it is some delayed linq-expr).
                List<DGField> FieldsFlatten = new List<DGField>();
                TimestepInfo.FlattenHierarchy(FieldsFlatten, Fields);
                foreach (var f in FieldsFlatten) {
                    if (!Fields.Contains(f, (a, b) => object.ReferenceEquals(a, b)))
                        throw new ArgumentException("Unable to load timestep: field '" + f.Identification + "', which is required by at least one of the given fields, must also be contained in the given list of fields.", "PreAllocatedFields");
                }

                // Load data vector
                // ================
                var partition = grdDat.CellPartitioning;
                var DataVec = this.LoadVector<CellFieldDataSet>(info.StorageID, ref partition);

                // Permute data vector
                // ===================


                var SortedDataVec = new CellFieldDataSet[DataVec.Count];

                {
                    // tau   is the GlobalID-permutation that we have for the loaded vector
                    // sigma is the current GlobalID-permutation of the grid
                    var sigma = grdDat.CurrentGlobalIdPermutation;
                    var tau = new Permutation(DataVec.Select(cd => cd.GlobalID).ToArray(), csMPI.Raw._COMM.WORLD);

                    // compute resorting permutation
                    Permutation invSigma = sigma.Invert();
                    Permutation Resorting = invSigma * tau;
                    tau = null;
                    invSigma = null;

                    // put dg coordinates into right order
                    Resorting.ApplyToVector(DataVec, SortedDataVec);
                }


                // Load the fields
                // ===============
                HashSet<object> loadedObjects = new HashSet<object>(ReferenceComparer.Instance);

                foreach (var Field in Fields) {
                    Field.LoadData(info, SortedDataVec, loadedObjects);
                }
            }
        }

        /// <summary>
        /// Loads a time-step from the database.
        /// </summary>
        /// <remarks>
        /// By using this method, it is ensured that the loaded/returned fields
        /// have the same DG polynomial degree as in the database.
        /// </remarks>
        public IEnumerable<DGField> LoadFields(ITimestepInfo info, IGridData grdDat, IEnumerable<string> NameFilter = null) {
            using (var tr = new FuncTrace()) {
                // check
                // =====
                if (!info.Grid.ID.Equals(grdDat.GridID))
                    throw new ArgumentException("Mismatch in Grid.");

                // Instantiate
                // ==========
                IEnumerable<DGField.FieldInitializer> F2LoadInfo;
                if (NameFilter != null) {
                    F2LoadInfo = info.FieldInitializers.Where(fi => NameFilter.Contains(fi.Identification, (a, b) => a.Equals(b)));
                } else {
                    F2LoadInfo = info.FieldInitializers;
                }

                IInitializationContext ic = info.Database.Controller.GetInitializationContext(info);
                var fields = F2LoadInfo.Select(fi => fi.Initialize(ic)).ToArray();
                List<DGField> fieldsFlattened = new List<DGField>();
                TimestepInfo.FlattenHierarchy(fieldsFlattened, fields);

                this.LoadFieldData(info, grdDat, fieldsFlattened);

                return fieldsFlattened;
            }
        }

        private static JsonReader GetJsonReader(Stream s) {
            if (DebugSerialization) {
                return new JsonTextReader(new StreamReader(s));
            } else {
                return new BsonReader(s);
            }
        }

        private static JsonWriter GetJsonWriter(Stream s) {
            if (DebugSerialization) {
                return new JsonTextWriter(new StreamWriter(s));
            } else {
                return new BsonWriter(s);
            }
        }
    }
}
