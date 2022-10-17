using BoSSS.Platform.Utils.Geom;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.GridGen {
    
    /// <summary>
    /// some examples for grid creation
    /// </summary>
    static public class HardcodedControl {


        public static GridGenControl Cart3D(int xRes = 8, int yRes = 8, int zRes = 8) {
            // --control cs:BoSSS.Application.GridGen.HardcodedControl.Cart3D()
            var C = new GridGenControl();

            C.savetodb = true;
            C.DbPath = @"C:\Users\flori\default_bosss_db";

            C.GridName = "Hello from GridGen";

            C.GridBlocks = new GridGenControl.MeshBlock[] {
                new GridGenControl.Cartesian3D() {
                    xNodes = GenericBlas.Linspace(-1, 1, xRes+1),
                    yNodes = GenericBlas.Linspace(-1, 1, yRes+1),
                    zNodes = GenericBlas.Linspace(-1, 1, zRes+1)
                }
            };

            C.BoundaryRegions.Add((new BoundingBox(new double[] { -2, -2, -1 }, new double[] { +2, +2, -1 }), "xyPlane"));
            C.BoundaryRegions.Add((null, "allOthers"));



            return C;
        }


    }
}
