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


        public static GridGenControl Cart3D() {
            // --control cs:BoSSS.Application.GridGen.HardcodedControl.Cart3D()
            var C = new GridGenControl();

            C.savetodb = true;
            C.DbPath = @"C:\Users\flori\default_bosss_db";

            C.GridName = "Hello from GridGen";

            C.GridBlocks = new GridGenControl.MeshBlock[] {
                new GridGenControl.Cartesian3D() {
                    xNodes = GenericBlas.Linspace(-1, 1, 4),
                    yNodes = GenericBlas.Linspace(-1, 1, 4),
                    zNodes = GenericBlas.Linspace(-1, 1, 4)
                }
            };

            C.BoundaryRegions.Add((new BoundingBox(new double[] { -2, -2, -1 }, new double[] { +2, +2, -1 }), "xyPlane"));
            C.BoundaryRegions.Add((null, "allOthers"));
            


            return C;
        }


    }
}
