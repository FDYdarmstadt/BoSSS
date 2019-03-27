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

namespace VoronoiTests
{
    class DatabaseTests : TestProgram
    {
        static void Main(string[] args)
        {
            SetUp();
            LoadGridProxy();
        }
        
        static IDatabaseInfo GetDatabase()
        {
            string dataBasePath = @"C:\Users\beck\Documents\VoronoiVersuche";
            IDatabaseInfo database = new DatabaseInfo(dataBasePath);
            return database;
        }

        public static IGrid SaveAndLoad()
        {
            IDatabaseInfo database = GetDatabase();
            IGrid grid = GetVoronoiGrid();
            Guid gridId = SaveGrid(grid, database);
            IGrid databaseGrid = LoadGrid(gridId, database);
            return databaseGrid;
        }

        static IGrid InitializeGridWithDatabase()
        {
            return SaveAndLoad();
        }

        static IGrid GetVoronoiGrid()
        {
            AppControl VoronoiTest = SipHardcodedControl.TestVoronoi_LDomain(200);
            IApplication application = new SipPoissonMain();
            RunApplication(application, VoronoiTest);
            return application.Grid;
        }

        static void RunApplication(IApplication app, AppControl ctrl)
        {
            app.Init(ctrl);
            app.RunSolverMode();
        }

        static Guid SaveGrid(IGrid grid, IDatabaseInfo database)
        {
            Guid id = database.Controller.DBDriver.SaveGrid( grid, database);
            return id;
        }

        static IGrid LoadGrid(Guid gridId, IDatabaseInfo database)
        {
            IGrid grid = database.Controller.DBDriver.LoadGrid(gridId, database);
            return grid;
        }

        static IGridInfo LoadGridInfo(Guid gridId, IDatabaseInfo database)
        {
            IGridInfo gridInfo = database.Controller.DBDriver.LoadGridInfo(gridId, database);
            return gridInfo;
        }

        static GridProxy InitializeGridProxy()
        {
            IGrid grid = GetVoronoiGrid();
            IDatabaseInfo database = GetDatabase();
            Guid gridId = SaveGrid(grid, database);
            GridProxy proxyGrid = new GridProxy(grid.ID, grid.Database);
            return proxyGrid;
        }

        [Test]
        public static void LoadGridProxy()
        {
            GridProxy proxyGrid = InitializeGridProxy(); 
            IGrid gridFromProxy = proxyGrid.RealGrid;
        }

        [Test]
        public static void ProxyToString()
        {
            GridProxy proxyGrid = InitializeGridProxy();
            string proxyName = proxyGrid.ToString();
        }

        [Test]
        public static void TestLoadGrids()
        {
            IGrid grid = InitializeGridWithDatabase();
            IDatabaseInfo database = grid.Database;
            IEnumerable<IGridInfo> grids = database.Grids;
        }

    }
}
