using System;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.SipPoisson;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using NUnit.Framework;

namespace VoronoiTests.Database
{
    class GridIOTests : DatabaseTest
    {       
        static GridDatabaseMethods gridMethods;

        static GridDatabaseMethods GridMethods {
            get {
                if(gridMethods == null)
                {
                    gridMethods = new GridDatabaseMethods(Database);
                }
                return gridMethods;
            }
        }

        static IGrid CreateVoronoiGrid()
        {
            AppControl VoronoiTest = SipHardcodedControl.TestVoronoi_LDomain(200);
            IApplication application = new SipPoissonMain();
            RunApplication(application, VoronoiTest);
            return application.Grid;
        }

        [Test]
        public static void SaveAndLoad()
        {
            IGrid grid = CreateVoronoiGrid();
            Guid gridId = GridMethods.SaveGrid(grid);
            IGrid databaseGrid = GridMethods.LoadGrid(gridId);
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
