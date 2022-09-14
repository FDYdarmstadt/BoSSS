using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.DatabaseTests {
    
    /// <summary>
    /// 
    /// </summary>
    class MiscTests : DatabaseTest {


        [Test]
        public void GridEquivalenceTest() {
            var dbinfo = emptyDatabase;
            var driver = dbinfo.Controller.DBDriver;

            /*
            var grid1 = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, +1, 10), GenericBlas.Linspace(-1, +1, 10));
            var grid1_equiv = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, +1, 10), GenericBlas.Linspace(-1, +1, 10));

            var grid2 = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-0.1, +0.1, 10), GenericBlas.Linspace(-0.1, +0.1, 10));
            */

            IGrid CreateGrid(int Res, double scale) {
                double[] xNodes = GenericBlas.Linspace(0, 3 * scale, Res + 1);
                double[] yNodes = xNodes;
                double[] zNodes = GenericBlas.Linspace(-3 * scale, 3 * scale, Res * 2 + 1);
                var g = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                return g;
            }

            var grid1 = CreateGrid(6, 1);
            var grid1_equiv = CreateGrid(6, 1);
            var grid2 = CreateGrid(6, 0.1);


            IGrid _grid1 = grid1;
            driver.SaveGridIfUnique(ref _grid1, out bool eq1, dbinfo);
            Assert.IsFalse(eq1, "No grid should be in the database, but strangely some equivalent was found");
            Assert.IsTrue(object.ReferenceEquals(_grid1, grid1), "No grid should be in the database, but strangely some equivalent was found");


            IGrid _grid1_equiv = grid1_equiv;
            driver.SaveGridIfUnique(ref _grid1_equiv, out bool eq2, dbinfo);
            Assert.IsTrue(eq2, "Equivalent grid was not identified as equivalent.");
            Assert.IsFalse(object.ReferenceEquals(_grid1_equiv, grid1_equiv), "Equivalent grid was not identified as equivalent.");
            Assert.IsTrue(_grid1.ID == _grid1_equiv.ID, "Equivalent grid was not correctly identified.");


            IGrid _grid2 = grid2;
            driver.SaveGridIfUnique(ref _grid2, out bool eq3, dbinfo);
            Assert.IsFalse(eq3, "NON-Equivalent grid was identified as equivalent.");
            Assert.IsTrue(object.ReferenceEquals(_grid2, grid2), "NON-Equivalent grid was identified as equivalent.");
            


        }

    }
}
