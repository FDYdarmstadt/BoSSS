using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using NUnit.Framework;
using BoSSS.Platform;

namespace BoSSS.Application.DatabaseTests
{
    class GridDBDriverTests : TestDatabase
    {
        (IGrid, IGridInfo) LoadGrid()
        {
            var gridInfo = databaseWithFiles.Controller.Grids.First();
            Assert.IsNotNull(gridInfo);
            IGrid grid;
            if (gridInfo is GridProxy)
            {
                grid = gridInfo.As<GridProxy>().RealGrid;
            }
            else
            {
                grid = (IGrid)gridInfo;              
            }

            return (grid, gridInfo);
        }

        [Test]
        public void SaveGridIfUnique()
        {
            (IGrid grid , IGridInfo gridInfo) = LoadGrid();
            bool isNotUnique;
            databaseWithFiles.Controller.DBDriver.SaveGridIfUnique(ref grid, out isNotUnique, databaseWithFiles);
            Assert.IsTrue(isNotUnique == true, "Same grid was not recognized.");

            var grid2 = databaseWithFiles.Controller.CopyGrid(gridInfo, emptyDatabase);
            databaseWithFiles.Controller.DBDriver.SaveGridIfUnique(ref grid, out isNotUnique, emptyDatabase);
            Assert.IsTrue(isNotUnique == true, "Copied grid was not recognized.");
        }

        [Test]
        public void LoadAndSaveGrid()
        {
            (IGrid grid, IGridInfo gridInfo) = LoadGrid();

            Guid id = emptyDatabase.Controller.DBDriver.SaveGrid(grid, emptyDatabase);
            IGrid newGrid = emptyDatabase.Controller.DBDriver.LoadGrid(id, emptyDatabase);
            Assert.IsTrue(newGrid.Equals(grid));
        }
    }
}
