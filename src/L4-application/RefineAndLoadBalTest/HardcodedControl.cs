
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;

namespace BoSSS.Application.RefineAndLoadBal {
    public static class HardcodedControl {

         static RefineAndLoadBalControl OnlyBalancing() {
            RefineAndLoadBalControl C = new RefineAndLoadBalControl();

            //C.InitialValues.Add("LevelSet", )


            C.GridFunc = delegate () {
                double[] nodes = GenericBlas.Linspace(-5, 5, 21);
                var grd = Grid2D.Cartesian2DGrid(nodes, nodes);
                return grd;
            };

            double s = 1.0;
            C.LevelSet = (double[] X, double t) => -(X[0] - s*t).Pow2() - X[1].Pow2() + (2.4).Pow2();
           

            C.GridPartType = GridPartType.none;
            C.DynamicLoadBalancing_ImbalanceThreshold = 0.0;
            C.DynamicLoadBalancing_Period = 3;

            C.NoOfTimesteps = 20;
            C.dtFixed = 0.1;

            return C;
        }

        static RefineAndLoadBalControl Refine() {
            RefineAndLoadBalControl C = new RefineAndLoadBalControl();

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-5, 5, 4);
                double[] Ynodes = GenericBlas.Linspace(-5, 5, 4);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                return grd;
            };

            C.LevelSet = (double[] X, double t) => X[0];
           

            C.NoOfTimesteps = 1;
            C.dtFixed = 0.0001;


            return C;
        }

    }
}
