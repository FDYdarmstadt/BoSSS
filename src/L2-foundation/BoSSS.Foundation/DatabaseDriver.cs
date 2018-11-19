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

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// The (formerly called) IO-master provides storage access on an
    /// object-level (in comparison, <see cref="IFileSystemDriver"/>-
    /// implementations provide IO on a stream-level), mainly by using
    /// serialization.
    /// </summary>
    public class DatabaseDriver : IDatabaseDriver {

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
        /// tests whether a grid with GUID <paramref name="g"/> exists in database, or not;
        /// </summary>
        public bool GridExists(Guid g) {
            try {
                Stream strm = m_fsDriver.GetGridStream(false, g);
                strm.Close();
            } catch (Exception) {
                return false;
            }
            return true;
        }

        /// <summary>
        /// Loads the grid info object for the given
        /// <paramref name="gridGuid"/> from the given
        /// <paramref name="database"/>
        /// </summary>
        /// <param name="gridGuid"></param>
        /// <param name="database"></param>
        /// <returns></returns>
        public IGridInfo LoadGridInfo(Guid gridGuid, IDatabaseInfo database) {
            using (var tr = new FuncTrace()) {
                tr.Info("Loading grid " + gridGuid);

                Grid.Classic.GridCommons grid = null;
                if (MyRank == 0) {
                    using (Stream s = FsDriver.GetGridStream(false, gridGuid))
                    using (var reader = GetJsonReader(s)) {
                        grid = m_Formatter.Deserialize<Grid.Classic.GridCommons>(reader);
                    }
                }

                grid = grid.MPIBroadcast(0);
                grid.Database = database;
                grid.WriteTime = Utils.GetGridFileWriteTime(grid);
                return grid;
            }
        }

        /// <summary>
        /// Loads the actual grid data for the given <paramref name="grid"/>.
        /// That is, loads the actual cell data.
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public IGrid LoadGridData(GridCommons grid) {
            // load grid data
            Partitioning p = null;
            grid.Cells = LoadVector<Cell>(grid.StorageGuid, ref p).ToArray();

            p = null;
            if (!grid.BcCellsStorageGuid.Equals(Guid.Empty))
                grid.BcCells = LoadVector<BCElement>(grid.BcCellsStorageGuid, ref p).ToArray();

            grid.InitNumberOfCells();

            // return
            // ------
            return grid;
        }

        /// <summary>
        /// loads the grid identified by <paramref name="uid"/> from the
        /// given <paramref name="database"/>
        /// </summary>
        /// <param name="uid">The unique identifier of the grid.</param>
        /// <param name="database">
        /// The database that is associated with the grid.
        /// </param>
        /// <returns>
        /// The loaded grid
        /// </returns>
        public IGrid LoadGrid(Guid uid, IDatabaseInfo database) {
            return LoadGridData((Grid.Classic.GridCommons)LoadGridInfo(uid, database));
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

        static bool GridCommons_CustomEquality(Grid.Classic.GridCommons A, Grid.Classic.GridCommons B) {
            if (object.ReferenceEquals(A, B))
                return true;
            if ((A == null) != (B == null))
                return false;

            ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

            //Debugger.Launch();
            //int[] locNoOfCells = new int[] {
            //            A.Cells.Length, B.Cells.Length,
            //            (A.BcCells != null ? A.BcCells.Length : 0), (B.BcCells != null ? B.BcCells.Length : 0) };
            //int[] glbNoOfCells = new int[4];
            //unsafe
            //{
            //    fixed (int* pLocNoOfCells = locNoOfCells, pGlbNoOfCells = glbNoOfCells)
            //    {
            //        csMPI.Raw.Allreduce((IntPtr)(pLocNoOfCells), (IntPtr)(pGlbNoOfCells), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            //    }
            //}

            int glbNoOfCells_A = A.NumberOfCells;
            int glbNoOfCells_B = B.NumberOfCells;
            int glbNoOfBcCells_A = A.NumberOfBcCells;
            int glbNoOfBcCells_B = B.NumberOfBcCells;

            if (glbNoOfCells_A != glbNoOfCells_B)
                return false;

            if (glbNoOfBcCells_A != glbNoOfBcCells_B)
                return false;

            if (!ArrayTools.ListEquals(A.RefElements, B.RefElements, (a, b) => object.ReferenceEquals(a, b)))
                return false;
            if (!ArrayTools.ListEquals(A.EdgeRefElements, B.EdgeRefElements, (a, b) => object.ReferenceEquals(a, b)))
                return false;
            if (!ArrayTools.ListEquals(A.EdgeTagNames, B.EdgeTagNames, (a, b) => (a.Key == b.Key && a.Value.Equals(b.Value))))
                return false;
            if (!ArrayTools.ListEquals(A.PeriodicTrafo, B.PeriodicTrafo, (a, b) => a.ApproximateEquals(b)))
                return false;

            return true;
        }

        bool GridCommons_CellEquality(Grid.Classic.GridCommons A, Grid.Classic.GridCommons B) {
            if (object.ReferenceEquals(A, B))
                return true;
            if ((A == null) != (B == null))
                return false;


            if (A.Cells == null)
                throw new ArgumentException();
            int A_NumberOfBcCells = A.NumberOfBcCells;

            int match = 1;

            {

                // load cells of grid B, if required
                // ---------------------------------

                Cell[] B_Cells;
                if (B.Cells == null) {
                    Partitioning p = A.CellPartitioning;
                    B_Cells = this.LoadVector<Cell>(B.StorageGuid, ref p).ToArray();
                } else {
                    B_Cells = B.Cells;
                }

                if (A.Cells.Length != B_Cells.Length)
                    throw new ApplicationException();

                // put the cells of B into the same order as those of A
                // ----------------------------------------------------

                {
                    // tau   is the GlobalID-permutation that we have for the loaded vector
                    // sigma is the current GlobalID-permutation of the grid
                    var sigma = new Permutation(A.Cells.Select(cell => cell.GlobalID).ToArray(), csMPI.Raw._COMM.WORLD);
                    var tau = new Permutation(B_Cells.Select(cell => cell.GlobalID).ToArray(), csMPI.Raw._COMM.WORLD);

                    if (sigma.TotalLength != tau.TotalLength)
                        // should have been checked already
                        throw new ArgumentException();

                    // compute resorting permutation
                    Permutation invSigma = sigma.Invert();
                    Permutation Resorting = invSigma * tau;
                    tau = null;      // Werfen wir sie dem GC zum Fraﬂe vor!
                    invSigma = null;

                    // put dg coordinates into right order
                    Resorting.ApplyToVector(B_Cells.CloneAs(), B_Cells);
                }

                // compare cells
                // -------------

                for (int j = 0; j < A.Cells.Length; j++) {
                    Cell Ca = A.Cells[j];
                    Cell Cb = B_Cells[j];

                    Debug.Assert(Ca.GlobalID == Cb.GlobalID);

                    if (!ArrayTools.ListEquals(Ca.NodeIndices, Cb.NodeIndices, (ia, ib) => ia == ib)) {
                        match = 0;
                        break;
                    }

                    if (Ca.Type != Cb.Type) {
                        match = 0;
                        break;
                    }

                    if (Ca.CellFaceTags != null || Cb.CellFaceTags != null) {

                        CellFaceTag[] CFTA = Ca.CellFaceTags != null ? Ca.CellFaceTags : new CellFaceTag[0];
                        CellFaceTag[] CFTB = Cb.CellFaceTags != null ? Cb.CellFaceTags : new CellFaceTag[0];

                        if (CFTA.Length != CFTB.Length) {
                            match = 0;
                            break;
                        }

                        bool setMatch = true;
                        for (int i1 = 0; i1 < CFTA.Length; i1++) {
                            bool b = false;
                            for (int j1 = 0; j1 < CFTB.Length; j1++) {
                                if (CFTA[i1].Equals(CFTB[j1])) {
                                    b = true;
                                    break;
                                }
                            }

                            if (b == false) {
                                setMatch = false;
                                break;
                            }
                        }

                        if (!setMatch) {
                            match = 0;
                            break;
                        }
                    }


                    double h = Math.Min(Ca.TransformationParams.MindistBetweenRows(), Cb.TransformationParams.MindistBetweenRows());
                    double L2Dist = Ca.TransformationParams.L2Dist(Cb.TransformationParams);
                    if (L2Dist > h * 1.0e-9) {
                        match = 0;
                        break;
                    }

                }
            }


            if (A_NumberOfBcCells > 0) {
                BCElement[] B_BcCells;
                if (B.BcCells == null && !B.BcCellsStorageGuid.Equals(Guid.Empty)) {
                    Partitioning p = A.BcCellPartitioning;
                    B_BcCells = this.LoadVector<BCElement>(B.StorageGuid, ref p).ToArray();
                } else {
                    B_BcCells = B.BcCells;
                }

                if (A.BcCells.Length != B_BcCells.Length)
                    throw new ApplicationException("Internal error.");


                // put the cells of B into the same order as those of A
                // ----------------------------------------------------

                {
                    long Offset = A.NumberOfCells_l;

                    // tau   is the GlobalID-permutation that we have for the loaded vector
                    // sigma is the current GlobalID-permutation of the grid
                    var sigma = new Permutation(A.BcCells.Select(cell => cell.GlobalID - Offset).ToArray(), csMPI.Raw._COMM.WORLD);
                    var tau = new Permutation(B_BcCells.Select(cell => cell.GlobalID - Offset).ToArray(), csMPI.Raw._COMM.WORLD);

                    if (sigma.TotalLength != tau.TotalLength)
                        // should have been checked already
                        throw new ArgumentException();

                    // compute resorting permutation
                    Permutation invSigma = sigma.Invert();
                    Permutation Resorting = invSigma * tau;
                    tau = null;      // Werfen wir sie dem GC zum Fraﬂe vor!
                    invSigma = null;

                    // put dg coordinates into right order
                    Resorting.ApplyToVector(B_BcCells.CloneAs(), B_BcCells);
                }


                // compare cells
                // -------------

                for (int j = 0; j < A.BcCells.Length; j++) {
                    BCElement Ca = A.BcCells[j];
                    BCElement Cb = B_BcCells[j];

                    Debug.Assert(Ca.GlobalID == Cb.GlobalID);

                    if (!ArrayTools.ListEquals(Ca.NodeIndices, Cb.NodeIndices, (ia, ib) => ia == ib)) {
                        match = 0;
                        break;
                    }

                    if (Ca.Type != Cb.Type) {
                        match = 0;
                        break;
                    }

                    if (Ca.Conformal != Cb.Conformal) {
                        match = 0;
                        break;
                    }

                    if (Ca.EdgeTag != Cb.EdgeTag) {
                        match = 0;
                        break;
                    }


                    if (Ca.NeighCell_GlobalIDs != null || Cb.NeighCell_GlobalIDs != null) {

                        long[] NgA = Ca.NeighCell_GlobalIDs != null ? Ca.NeighCell_GlobalIDs : new long[0];
                        long[] NgB = Cb.NeighCell_GlobalIDs != null ? Cb.NeighCell_GlobalIDs : new long[0];

                        if (NgA.Length != NgB.Length) {
                            match = 0;
                            break;
                        }

                        bool setMatch = true;
                        for (int i1 = 0; i1 < NgA.Length; i1++) {
                            bool b = false;
                            for (int j1 = 0; j1 < NgB.Length; j1++) {
                                if (NgA[i1] == NgB[j1]) {
                                    b = true;
                                    break;
                                }
                            }

                            if (b == false) {
                                setMatch = false;
                                break;
                            }
                        }

                        if (!setMatch) {
                            match = 0;
                            break;
                        }
                    }


                    double h = Math.Min(Ca.TransformationParams.MindistBetweenRows(), Cb.TransformationParams.MindistBetweenRows());
                    double L2Dist = Ca.TransformationParams.L2Dist(Cb.TransformationParams);
                    if (L2Dist > h * 1.0e-9) {
                        match = 0;
                        break;
                    }

                }
            }


            match = match.MPIMin();
            return (match > 0);
        }


        /// <summary>
        /// Searches for an equivalent grid in the database and, if none is found
        /// saves a grid object to the database.
        /// </summary>
        /// <param name="grd">
        /// On entry, the grid which should be saved to the database.
        /// On exit, either unchanged, or the equivalent grid.
        /// </param>
        /// <param name="EquivalentGridFound">
        /// Inidicates that an equivalent grid was found.
        /// </param>
        /// <param name="database"></param>
        public Guid SaveGridIfUnique(ref Grid.Classic.GridCommons grd, out bool EquivalentGridFound, IDatabaseInfo database) {
            using (new FuncTrace()) {

                var Grids = database.Grids;
                foreach (var GrdInf in Grids) {
                    Grid.Classic.GridCommons GrdInDb = (Grid.Classic.GridCommons)this.LoadGridInfo(GrdInf.ID, database);

                    if (GridCommons_CustomEquality(grd, GrdInDb) == false)
                        continue;

                    if (GridCommons_CellEquality(grd, GrdInDb)) {
                        grd = (Grid.Classic.GridCommons)LoadGridData(GrdInDb);
                        EquivalentGridFound = true;
                        //Console.WriteLine("Found equivalent grid: " + grd.GridGuid);
                        return grd.ID;
                    }
                }

                //Console.WriteLine("Grid is unique - saving...");

                EquivalentGridFound = false;
                var g = SaveGrid(grd, database);
                //Console.WriteLine("done: " + g.ToString());
                return g;
            }
        }



        /// <summary>
        /// Saves the given grid object to the database;
        /// </summary>
        /// <returns>
        /// the Guid of the <see cref="GridCommons"/>-object that was saved
        /// (equal to the <see cref="GridCommons.GridGuid"/>-property).
        /// </returns>
        /// <param name="grd">
        /// The grid to save.
        /// </param>
        /// <param name="database">
        /// chaos
        /// </param>
        public Guid SaveGrid(Grid.Classic.GridCommons grd, IDatabaseInfo database) {
            using (new FuncTrace()) {
                if (grd.GridGuid.Equals(Guid.Empty)) {
                    throw new ApplicationException("cannot save grid with empty Guid (Grid Guid is " + Guid.Empty.ToString() + ");");
                }
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                // save required grid data
                // =======================
                if (grd.StorageGuid != null && grd.StorageGuid != Guid.Empty) {
                    SaveVector(grd.Cells, grd.StorageGuid);
                } else {
                    grd.StorageGuid = SaveVector(grd.Cells);
                }

                grd.InitNumberOfCells();

                // save opt. data
                // ==============
                foreach (var s in grd.m_PredefinedGridPartitioning) {
                    int[] cellToRankMap = s.Value.CellToRankMap;
                    if (cellToRankMap == null) {
                        // Partitioning has not been loaded; do it now
                        Partitioning currentPartitioning = grd.CellPartitioning;
                        cellToRankMap = LoadVector<int>(s.Value.Guid, ref currentPartitioning).ToArray();
                    }

                    SaveVector(cellToRankMap, s.Value.Guid);
                }

                if (grd.BcCells != null && grd.BcCells.Length > 0)
                    grd.BcCellsStorageGuid = SaveVector(grd.BcCells);

                // save header data
                // ================
                if (MyRank == 0) {
                    using (Stream s = m_fsDriver.GetGridStream(true, grd.GridGuid))
                    using (var writer = GetJsonWriter(s)) {
                        m_Formatter.Serialize(writer, grd);
                        writer.Close();
                        s.Close();
                    }
                }

                // return
                // ======

                grd.Database = database;

                return grd.GridGuid;
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
        /// <param name="physTime">Physical time of the time-step.</param>
        /// <param name="TimestepNo">Time-step number.</param>
        /// <param name="currentSession">
        /// The session associated with the time-step.
        /// </param>
        /// <param name="fields">The fields of the time-step.</param>
        /// <param name="GridDat">grid data object (required if <paramref name="fields"/> is empty)</param>
        /// <returns>
        /// An object containing information about the time-step.
        /// </returns>
        public TimestepInfo SaveTimestep(double physTime, TimestepNumber TimestepNo, SessionInfo currentSession, IGridData GridDat, IEnumerable<DGField> fields) {
            using (var tr = new FuncTrace()) {

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
                                "fields");
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

                // Save state object 
                // =================
                TimestepInfo tsi = null;
                Exception e = null;
                if (MyRank == 0) {
                    try {
                        tsi = new TimestepInfo(physTime, currentSession, TimestepNo, fields, VectorGuid);
                        using(var s = FsDriver.GetTimestepStream(true, tsi.ID))
                        using(var writer = GetJsonWriter(s)) {
                            m_Formatter.Serialize(writer, tsi);
                            writer.Close();
                            s.Close();
                        }
                    } catch (Exception ee) {
                        e = ee;
                        Console.Error.WriteLine(ee.GetType().Name + " on rank " + MyRank + " saving time-step " + TimestepNo + ": " + ee.Message);
                        Console.Error.WriteLine(ee.StackTrace);
                    }
                }
                e.ExceptionBcast();
                tsi = tsi.MPIBroadcast(0);

                // return
                // ======
                if (MyRank == 0) {
                    try {
                        currentSession.LogTimeStep(tsi.ID);
                    } catch (Exception ee) {
                        e = ee;
                        Console.Error.WriteLine(ee.GetType().Name + " on rank " + MyRank + " saving time-step " + TimestepNo + ": " + ee.Message);
                        Console.Error.WriteLine(ee.StackTrace);
                    }
                }
                e.ExceptionBcast();

                tsi.Database = currentSession.Database;
                return tsi;
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
