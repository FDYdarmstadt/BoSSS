using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.SipPoisson;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;

namespace VoronoiTests.Database
{
    class GridMethods
    {
        readonly IDatabaseInfo database;

        public GridMethods(IDatabaseInfo database)
        {
            this.database = database;
        }      

        public Guid SaveGrid(IGrid grid)
        {
            Guid id = database.Controller.DBDriver.SaveGrid(grid, database);
            return id;
        }

        public IGrid LoadGrid(Guid gridId)
        {
            IGrid grid = database.Controller.DBDriver.LoadGrid(gridId, database);
            return grid;
        }

        IGridInfo LoadGridInfo(Guid gridId)
        {
            IGridInfo gridInfo = database.Controller.DBDriver.LoadGridInfo(gridId, database);
            return gridInfo;
        }

        public GridProxy GetGridProxy(IGrid grid)
        {
            Guid gridId = SaveGrid(grid);
            GridProxy proxyGrid = new GridProxy(grid.ID, grid.Database);
            return proxyGrid;
        }

        public bool AreEqual(IGrid gridA, IGrid gridB)
        {
            IEqualityComparer<IGrid> cellComparer = gridA.GridSerializationHandler.CellComparer;
            IEqualityComparer<IGrid> referenceComparer = gridA.GridSerializationHandler.ReferenceComparer;
            bool areEqual = cellComparer.Equals(gridA, gridB);
            areEqual = areEqual && referenceComparer.Equals(gridA, gridB);
            return areEqual;
        }
    }
}
