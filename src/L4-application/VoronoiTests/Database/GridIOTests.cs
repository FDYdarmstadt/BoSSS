using System;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using NUnit.Framework;

namespace VoronoiTests.Database
{
    class GridIOTests : DatabaseTest
    {       
        static GridMethods gridMethods;

        static GridMethods GridMethods {
            get {
                return gridMethods ?? (gridMethods = new GridMethods(Database)); 
            }
        }

        static VoronoiGrid CreateVoronoiGrid()
        {
            Vector[] DomainBndyPolygon = new[] {
                        new Vector(-1,1),
                        new Vector(1,1),
                        new Vector(1,-1),
                        new Vector(-1,-1)
                    };
            VoronoiGrid grid = VoronoiGrid2D.FromPolygonalDomain(DomainBndyPolygon, 0, 400);
            return grid;
        }

        [Test]
        public static void SaveAndLoad()
        {
            IGrid originalGrid = CreateVoronoiGrid();
            Guid gridId = GridMethods.SaveGrid(originalGrid);
            IGrid gridFromDatabase = GridMethods.LoadGrid(gridId);

            bool gridsAreEqual = GridMethods.AreEqual(originalGrid, gridFromDatabase);
            Assert.IsTrue(gridsAreEqual);
        }
        
        [Test]
        public static void LoadGridProxy()
        {
            IGrid grid = CreateVoronoiGrid();
            GridProxy proxyGrid = GridMethods.GetGridProxy(grid); 
            IGrid gridFromProxy = proxyGrid.RealGrid;
        }
        
        [Test]
        public static void ProxyToString()
        {
            IGrid grid = CreateVoronoiGrid();
            GridProxy proxyGrid = GridMethods.GetGridProxy(grid);
            string proxyName = proxyGrid.ToString();
        }

        [Test]
        public static void TestLoadGrids()
        {
            IGrid grid = CreateVoronoiGrid();
            Guid gridId = GridMethods.SaveGrid(grid);
            IGrid databaseGrid = GridMethods.LoadGrid(gridId);
            IDatabaseInfo database = databaseGrid.Database;
            IEnumerable<IGridInfo> grids = database.Grids;
        }
    }
}
