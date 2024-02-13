﻿using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;
using System.Diagnostics;

namespace BoSSS.Foundation.IO {
    class GridDatabaseDriver : MPIProcess {
        IFileSystemDriver fsDriver;
        readonly IVectorDataSerializer Driver;
        readonly ReflectionVectorDataSerializer dynamicDriver;

        public GridDatabaseDriver(IVectorDataSerializer driver, IFileSystemDriver FsDriver) {
            fsDriver = FsDriver;
            Driver = driver;
            dynamicDriver = new ReflectionVectorDataSerializer(driver);
        }

        Stream GetGridStream(bool create, Guid id) {
            return fsDriver.GetGridStream(create, id);
        }

        /// <summary>
        /// tests whether a grid with GUID <paramref name="g"/> exists in database, or not;
        /// </summary>
        public bool GridExists(Guid g) {
            try {
                Stream strm = fsDriver.GetGridStream(false, g);
                strm.Close();
            } catch (Exception) {
                return false;
            }
            return true;
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
        /// Indicates that an equivalent grid was found.
        /// </param>
        /// <param name="database"></param>
        public Guid SaveGridIfUnique(ref IGrid grid, out bool EquivalentGridFound, IDatabaseInfo database) {
            using (new FuncTrace("SaveGridIfUnique")) {
                var Grids = database.Grids;
                foreach (var GrdInf in Grids) {
                    IGrid gridInDatabase = (IGrid)this.LoadGridInfo(GrdInf.ID, database);

                    IEqualityComparer<IGrid> cellComparer = grid.GridSerializationHandler.CellComparer;
                    IEqualityComparer<IGrid> referenceComparer = grid.GridSerializationHandler.BasePropertiesComparer;

                    if (referenceComparer.Equals(grid, gridInDatabase)) {
                        gridInDatabase = LoadGridData(gridInDatabase);
                        if (cellComparer.Equals(grid, gridInDatabase)) {
                            grid = gridInDatabase;
                            EquivalentGridFound = true;
                            return grid.ID;
                        }
                    }
                }
                EquivalentGridFound = false;
                var g = SaveGrid(grid, database);
                return g;
            }
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
        public Guid SaveGrid(IGrid grid, IDatabaseInfo database) {
            using (new FuncTrace()) {
                CheckIfSuitableForSaving(grid);
                SaveGridData(grid.GridSerializationHandler, database);
                if (MyRank == 0) {
                    SaveGridInfo(grid);
                }
                Guid id = grid.ID;
                return id;
            }
        }

        void CheckIfSuitableForSaving(IGrid grid) {
            if (grid.ID.Equals(Guid.Empty)) {
                throw new ApplicationException("cannot save grid with empty Guid (Grid Guid is " + Guid.Empty.ToString() + ");");
            }
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
        }

        void SaveGridData(IGridSerializationHandler grdHandler, IDatabaseInfo database) {
            SaveVectorData(grdHandler);
            grdHandler.Database = database;
        }

        void SaveVectorData(IGridSerializationHandler grid) {
            using (var tr = new FuncTrace()) {
                object[][] vectorData = grid.GetVectorData();
                Guid[] vectorGuids = grid.GetVectorGuids();
                Type[] vectorTypes = grid.GetVectorTypes();

                int numberOfVectors = vectorData.Length;
                for (int i = 0; i < numberOfVectors; ++i) {
                    Guid guid = vectorGuids[i];
                    Type vectorType = vectorTypes[i];
                    object[] vector = vectorData[i];
                    if (vector != null) {
                        if (guid != Guid.Empty) {
                            dynamicDriver.SaveVector(vectorData[i], vectorType, guid);
                        } else {
                            vectorGuids[i] = dynamicDriver.SaveVector(vectorData[i], vectorType);
                        }

                    }
                }
                grid.SetVectorGuids(vectorGuids);
            }
        }

        void SaveGridInfo(IGrid grid) {
            using (new FuncTrace()) {
                using (Stream stream = GetGridStream(true, grid.ID)) {
                    Driver.Serialize(stream, grid, typeof(IGrid));
                }
            }
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
            using (new FuncTrace()) {
                IGridInfo gridInfo = LoadGridInfo(uid, database);
                IGrid grid = LoadGridData((IGrid)gridInfo);
                return grid;
            }
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

                IGrid grid = null;
                grid = DeserializeGrid(gridGuid);
                /*
                if (MyRank == 0)
                {
                    grid = DeserializeGrid(gridGuid);
                }
                
                grid = grid.MPIBroadcast(0);
                */
                grid.GridSerializationHandler.Database = database;
                grid.WriteTime = Utils.GetGridFileWriteTime(grid);
                return grid;
            }
        }

        public IGrid DeserializeGrid(Guid gridGuid) {
            using (var tr = new FuncTrace()) {
                // dbg_launch();
                using (Stream s = GetGridStream(false, gridGuid)) {
                    IGrid grid = (IGrid)Driver.Deserialize(s, typeof(IGrid));
                    return grid;
                }
            }
        }

        /// <summary>
        /// Loads the actual grid data for the given <paramref name="grid"/>.
        /// That is, loads the actual cell data.
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public IGrid LoadGridData(IGrid grid) {
            LoadGridData(grid.GridSerializationHandler);
            return grid;
        }

        public void LoadGridData(IGridSerializationHandler gridSerializationHandler) {
            using (new FuncTrace()) {
                Type[] vectorTypes = gridSerializationHandler.GetVectorTypes();
                Guid[] guids = gridSerializationHandler.GetVectorGuids();
                Partitioning p = null;

                int numberOfVectors = vectorTypes.Length;
                object[][] vectors = new object[numberOfVectors][];
                for (int i = 0; i < numberOfVectors; ++i) {
                    Guid guid = guids[i];
                    if (guid != Guid.Empty) {
                        Type vectorType = vectorTypes[i];
                        vectors[i] = (object[])dynamicDriver.LoadVector(guid, vectorType, ref p);
                    }
                }
                gridSerializationHandler.SetVectorData(vectors);
            }
        }
    }
}
