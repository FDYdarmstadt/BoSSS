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
            GridProxy();
        }
        
        public static IGrid SaveAndLoad()
        {
            string dataBasePath = @"C:\Users\beck\Documents\VoronoiVersuche";
            IDatabaseInfo database = new DatabaseInfo(dataBasePath);
            IGrid grid = GetVoronoiGrid();
            Guid gridId = SaveGrid(grid, database);
            IGrid databaseGrid = LoadGrid(gridId, database);
            return databaseGrid;
        }

        static IGrid GetVoronoiGrid()
        {
            AppControl VoronoiTest = SipHardcodedControl.TestVoronoi_LDomain(200);
            IApplication application = new SipPoissonMain();
            Run(application, VoronoiTest);
            return application.Grid;
        }

        static void Run(IApplication app, AppControl ctrl)
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

        public static void GridProxy()
        {
            IGrid grid = SaveAndLoad();
            GridProxy proxyGrid = new GridProxy(grid.ID, grid.Database);
            IGrid gridFromProxy = proxyGrid.RealGrid;
        }
    }

    /// <summary>
    /// Abstract base class for tests.
    /// </summary>
    [TestFixture]
    public abstract class TestProgram
    {
        /// <summary>
        /// Performs bootstrapping.
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp()
        {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                Application.GetBoSSSInstallDir(),
                out dummy);
            Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CurrentCulture;
        }
    }
}
