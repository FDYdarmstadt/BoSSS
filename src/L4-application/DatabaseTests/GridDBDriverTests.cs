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
        [Test]
        public void TestSaveGridIfUnique()
        {
            (IGrid grid , IGridInfo gridInfo) = LoadGrid();
            bool isNotUnique;
            databaseWithFiles.Controller.DBDriver.SaveGridIfUnique(ref grid, out isNotUnique, databaseWithFiles);
            Assert.IsTrue(isNotUnique == true, "Same grid was not recognized.");

            var grid2 = databaseWithFiles.Controller.CopyGrid(gridInfo, emptyDatabase);
            databaseWithFiles.Controller.DBDriver.SaveGridIfUnique(ref grid, out isNotUnique, emptyDatabase);
            Assert.IsTrue(isNotUnique == true, "Copied grid was not recognized.");
        }

        public (IGrid, IGridInfo) LoadGrid()
        {
            var gridInfo = databaseWithFiles.Controller.Grids.First();
            Assert.IsNotNull(gridInfo);
            IGrid grid;
            if (gridInfo is GridProxy)
            {
                grid = gridInfo.As<GridProxy>().RealGrid;
            }
            else if (gridInfo is Foundation.Grid.Classic.GridCommons)
            {
                grid = (Foundation.Grid.Classic.GridCommons)gridInfo;
            }
            else
            {
                throw new NotSupportedException();
            }
            return (grid, gridInfo);
        }

        [Test]
        void LoadGridInfo()
        {

        }

        [Test]
        void LoadGridData()
        {

        }

        [Test]
        public void SaveGrid()
        {

        }

        [Test]
        void SaveGridInfo()
        {

        }

        [Test]
        public void SaveGridData()
        {

        }
    }
}
