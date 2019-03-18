using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using ilPSP.Utils;

namespace BoSSS.Foundation.IO
{
    class TimeStepDatabaseDriver : MPIProcess
    {

        readonly IVectorDataSerializer Driver;
        IFileSystemDriver fsDriver;
        public TimeStepDatabaseDriver(IVectorDataSerializer driver, IFileSystemDriver FsDriver)
        {
            fsDriver = FsDriver;
            Driver = driver;
        }

        /// <summary>
        /// loads a single <see cref="TimestepInfo"/>-object from the database.
        /// </summary>
        public T LoadTimestepInfo<T>(Guid timestepGuid, ISessionInfo session, IDatabaseInfo database)
            where T : TimestepInfo //
        {
            using (var tr = new FuncTrace())
            {
                tr.Info("Loading time-step " + timestepGuid);

                T tsi = null;
                if (MyRank == 0)
                {
                    using (Stream s = fsDriver.GetTimestepStream(false, timestepGuid))
                    {
                        tsi = (T)Driver.Deserialize(s, typeof(T));
                        tsi.Session = session;
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
        /// Saves a time-step to the database's persistent memory.
        /// </summary>
        /// <param name="_tsi">Contains Id etc.</param>
        public void SaveTimestep(TimestepInfo _tsi)
        {
            using (var tr = new FuncTrace())
            {

                if (!(_tsi.ID.Equals(Guid.Empty) && _tsi.StorageID.Equals(Guid.Empty)))
                    throw new ArgumentException("Timestep is already saved in database");
                var fields = _tsi.Fields.ToArray();
                var GridDat = fields[0].GridDat;

                {
                    List<DGField> FieldsFlatten = new List<DGField>();
                    TimestepInfo.FlattenHierarchy(FieldsFlatten, fields);
                    foreach (var f in FieldsFlatten)
                    {
                        if (!object.ReferenceEquals(f.GridDat, GridDat))
                            throw new ArgumentException("mismatch in GridData object.");

                        if (!fields.Contains(f, (a, b) => object.ReferenceEquals(a, b)))
                        {
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
                for (int j = 0; j < J; j++)
                { // loop over cells
                    vec[j] = new CellFieldDataSet();
                    vec[j].GlobalID = Permutation[j];
                    //vec[j].DGCoordinateData = new CellFieldDataSet.CellFieldData[NF];
                    for (int idxF = 0; idxF < NF; idxF++)
                    { // loop over fields
                        var field = _fields[idxF];
                        int N = field.Basis.GetLength(j);
                        double[] Coords = new double[N];
                        for (int n = 0; n < N; n++)
                        {
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
                Guid VectorGuid = Driver.SaveVector(vec);
                _tsi.StorageID = VectorGuid;


                // Save state object 
                // =================
                _tsi.ID = Guid.NewGuid().MPIBroadcast(0);
                Exception e = null;
                if (MyRank == 0)
                {
                    try
                    {
                        //tsi = new TimestepInfo(physTime, currentSession, TimestepNo, fields, VectorGuid);
                        using (var s = fsDriver.GetTimestepStream(true, _tsi.ID))
                        {
                            Driver.Serialize(s, _tsi, typeof(TimestepInfo));
                            s.Close();
                        }
                    }
                    catch (Exception ee)
                    {
                        e = ee;
                        Console.Error.WriteLine(ee.GetType().Name + " on rank " + MyRank + " saving time-step " + _tsi.TimeStepNumber + ": " + ee.Message);
                        Console.Error.WriteLine(ee.StackTrace);
                    }
                }
                e.ExceptionBcast();

                // log session
                // ===========
                SessionInfo currentSession = (SessionInfo)(_tsi.Session); // hack
                if (MyRank == 0)
                {

                    try
                    {
                        currentSession.LogTimeStep(_tsi.ID);
                    }
                    catch (Exception ee)
                    {
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
        /// loads a single <see cref="TimestepInfo"/>-object from the database.
        /// </summary>
        public TimestepInfo LoadTimestepInfo(Guid timestepGuid, ISessionInfo session, IDatabaseInfo database)
        {
            using (var tr = new FuncTrace())
            {
                tr.Info("Loading time-step " + timestepGuid);

                TimestepInfo tsi = null;
                if (MyRank == 0)
                {
                    using (Stream s = fsDriver.GetTimestepStream(false, timestepGuid))
                    {
                        tsi = (TimestepInfo)Driver.Deserialize(s, typeof(TimestepInfo));
                        tsi.Session = session;
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
        public IEnumerable<Guid> GetTimestepGuids(Guid sessionGuid)
        {
            IList<Guid> timestepUids = new List<Guid>();

            try
            {
                using (StreamReader timestepLogReader =
                new StreamReader(fsDriver.GetTimestepLogStream(sessionGuid)))
                {

                    while (!timestepLogReader.EndOfStream)
                    {
                        timestepUids.Add(Guid.Parse(timestepLogReader.ReadLine()));
                    }
                }
            }
            catch (FileNotFoundException)
            {
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
        public void RemoveTimestepGuid(Guid sessionGuid, Guid timestepGuid)
        {
            string logPath = fsDriver.GetTimestepLogPath(sessionGuid);
            string[] lines = File.ReadAllLines(logPath);

            bool match = false;
            List<string> reducedLines = new List<string>();
            foreach (string line in lines)
            {
                Guid guid = Guid.Parse(line);
                if (guid.Equals(timestepGuid))
                {
                    match = true;
                }
                else
                {
                    reducedLines.Add(line);
                }
            }

            if (!match)
            {
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
        public void LoadFieldData(ITimestepInfo info, IGridData grdDat, IEnumerable<DGField> PreAllocatedFields)
        {
            using (var tr = new FuncTrace())
            {
                DGField[] Fields = PreAllocatedFields.ToArray(); // enforce 'evaluation' of the enum (in the case it is some delayed linq-expr).
                List<DGField> FieldsFlatten = new List<DGField>();
                TimestepInfo.FlattenHierarchy(FieldsFlatten, Fields);
                foreach (var f in FieldsFlatten)
                {
                    if (!Fields.Contains(f, (a, b) => object.ReferenceEquals(a, b)))
                        throw new ArgumentException("Unable to load timestep: field '" + f.Identification + "', which is required by at least one of the given fields, must also be contained in the given list of fields.", "PreAllocatedFields");
                }

                // Load data vector
                // ================
                var partition = grdDat.CellPartitioning;
                var DataVec = this.Driver.LoadVector<CellFieldDataSet>(info.StorageID, ref partition);

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

                foreach (var Field in Fields)
                {
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
        public IEnumerable<DGField> LoadFields(ITimestepInfo info, IGridData grdDat, IEnumerable<string> NameFilter = null)
        {
            using (var tr = new FuncTrace())
            {
                // check
                // =====
                if (!info.Grid.ID.Equals(grdDat.GridID))
                    throw new ArgumentException("Mismatch in Grid.");

                // Instantiate
                // ==========
                IEnumerable<DGField.FieldInitializer> F2LoadInfo;
                if (NameFilter != null)
                {
                    F2LoadInfo = info.FieldInitializers.Where(fi => NameFilter.Contains(fi.Identification, (a, b) => a.Equals(b)));
                }
                else
                {
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
    }
}
