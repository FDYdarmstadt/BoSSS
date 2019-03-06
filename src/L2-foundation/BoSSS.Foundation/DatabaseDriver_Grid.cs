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

namespace BoSSS.Foundation.IO
{
    public partial class DatabaseDriver
    {
        /// <summary>
        /// tests whether a grid with GUID <paramref name="g"/> exists in database, or not;
        /// </summary>
        public bool GridExists(Guid g)
        {
            try
            {
                Stream strm = m_fsDriver.GetGridStream(false, g);
                strm.Close();
            }
            catch (Exception)
            {
                return false;
            }
            return true;
        }

        /// <summary>
        /// Searches for an equivalent grid in the database and, if none is found
        /// saves a grid object to the database.
        /// </summary>
        /// <param name="_grd">
        /// On entry, the grid which should be saved to the database.
        /// On exit, either unchanged, or the equivalent grid.
        /// </param>
        /// <param name="EquivalentGridFound">
        /// Inidicates that an equivalent grid was found.
        /// </param>
        /// <param name="database"></param>
        public Guid SaveGridIfUnique(ref IGrid _grd, out bool EquivalentGridFound, IDatabaseInfo database)
        {
            using (new FuncTrace())
            {
                GridCommons grd = (GridCommons)_grd;

                var Grids = database.Grids;
                foreach (var GrdInf in Grids)
                {
                    GridCommons GrdInDatabase = (Grid.Classic.GridCommons)this.LoadGridInfo(GrdInf.ID, database);

                    if (GrdInDatabase.HasEqualReferences(grd) == true)
                    {
                        GrdInDatabase = (GridCommons)LoadGridData(GrdInDatabase);
                        if (grd.HasEqualCells(GrdInDatabase))
                        {
                            _grd = GrdInDatabase;
                            EquivalentGridFound = true;
                            return grd.ID;
                        }
                    }
                }
                EquivalentGridFound = false;
                var g = SaveGrid(grd, database);
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
        public Guid SaveGrid(IGrid grid, IDatabaseInfo database)
        {
            using (new FuncTrace())
            {
                Guid id = SaveGridData(grid, database);
                SaveGridInfo(grid);
                return id;
            }
        }

        Guid SaveGridData(IGrid _grd, IDatabaseInfo database)
        {
            GridCommons grd = _grd as GridCommons;
            if (grd == null)
            {
                throw new InvalidCastException("Grid not supported by save method");
            }

            if (_grd.ID.Equals(Guid.Empty))
            {
                throw new ApplicationException("cannot save grid with empty Guid (Grid Guid is " + Guid.Empty.ToString() + ");");
            }
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            // save required grid data
            // =======================
            if (grd.StorageGuid != null && grd.StorageGuid != Guid.Empty)
            {
                SaveVector(grd.Cells, grd.StorageGuid);
            }
            else
            {
                grd.StorageGuid = SaveVector(grd.Cells);
            }

            grd.InitNumberOfCells();

            // save opt. data
            // ==============

            foreach (var s in grd.m_PredefinedGridPartitioning)
            {
                int[] cellToRankMap = s.Value.CellToRankMap;
                if (cellToRankMap == null)
                {
                    // Partitioning has not been loaded; do it now
                    throw new Exception("Should have been already loaded ");
                }

                SaveVector(cellToRankMap, s.Value.Guid);
            }


            if (grd.BcCells != null && grd.BcCells.Length > 0)
                grd.BcCellsStorageGuid = SaveVector(grd.BcCells);

            // save header data (save GridInfo)
            // ==========================


            // return
            // ======

            grd.Database = database;

            return grd.ID;
        }

        void SaveGridInfo(IGrid grid)
        {
            if (MyRank == 0)
            {
                using (Stream s = m_fsDriver.GetGridStream(true, grid.ID))
                using (var writer = GetJsonWriter(s))
                {
                    m_Formatter.Serialize(writer, grid);
                    writer.Close();
                    s.Close();
                }
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
        public IGridInfo LoadGridInfo(Guid gridGuid, IDatabaseInfo database)
        {
            using (var tr = new FuncTrace())
            {
                tr.Info("Loading grid " + gridGuid);

                Grid.Classic.GridCommons grid = null;
                if (MyRank == 0)
                {
                    using (Stream s = FsDriver.GetGridStream(false, gridGuid))
                    using (var reader = GetJsonReader(s))
                    {
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
        public IGrid LoadGridData(GridCommons grid)
        {
            // load grid data
            Partitioning p = null;
            grid.Cells = LoadVector<Cell>(grid.StorageGuid, ref p).ToArray();

            p = null;
            if (!grid.BcCellsStorageGuid.Equals(Guid.Empty))
                grid.BcCells = LoadVector<BCElement>(grid.BcCellsStorageGuid, ref p).ToArray();

            for (int i = 0; i < grid.m_PredefinedGridPartitioning.Count; ++i)
            {
                var s = grid.m_PredefinedGridPartitioning.ElementAt(i);
                if (s.Value.CellToRankMap == null)
                {
                    // Partitioning has not been loaded; do it now
                    Partitioning currentPartitioning = grid.CellPartitioning;
                    int[] cellToRankMap = LoadVector<int>(s.Value.Guid, ref currentPartitioning).ToArray();
                    grid.m_PredefinedGridPartitioning[s.Key] = new GridCommons.GridPartitioningVector { Guid = s.Value.Guid, CellToRankMap = cellToRankMap };
                }
            }

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
        public IGrid LoadGrid(Guid uid, IDatabaseInfo database)
        {
            IGridInfo grid = LoadGridInfo(uid, database);
            return LoadGridData((Grid.Classic.GridCommons)grid);
        }
    }
}
