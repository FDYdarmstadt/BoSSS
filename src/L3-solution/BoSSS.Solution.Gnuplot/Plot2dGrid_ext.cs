/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Platform;
//using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Gnuplot {

    /// <summary>
    /// Grid plotting
    /// </summary>
    public static class Plot2dGrid_ext {


        public static void PlotCoordinateLabels(this MultidimensionalArray X, string filename) {
            using (var gp = new Gnuplot()) {
                int L = X.NoOfRows;
                if (X.NoOfCols != 2)
                    throw new NotSupportedException("works only for 2D grid");

                gp.Cmd("set terminal png");
                gp.Cmd("set output \"{0}\"", filename);

                gp.SetXRange(-4, 4);
                gp.SetYRange(-4, 4);

                for (int l = 0; l < L; l++) {
                    gp.Cmd("set label \"{0}\" at {1},{2}", l, X[l, 0], X[l, 1]);
                }

                gp.PlotSlope(0, 0);
                gp.Execute();
            }
        }



        /// <summary>
        /// plot a 2D grid
        /// </summary>
        public static void Plot2DGrid(this GridCommons grd) {
            using (var gp = new Gnuplot()) {
                if (grd.SpatialDimension != 2)
                    throw new NotSupportedException("works only for 2D grid");

                int J = grd.Cells.Length;
                for (int j = 0; j < J; j++) {
                    var Cell_j = grd.Cells[j];
                    //var Kref = grd.RefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));

                    //var Vtx = Kref.Vertices;
                    //var _Vtx = Kref.GetSubDivTree(3).GlobalVertice;
                    //var Vtx = MultidimensionalArray.Create(_Vtx.GetLength(0), _Vtx.GetLength(1));
                    //Vtx.SetA2d(_Vtx);
                    //Vtx.Scale(1);

                    var Vtx_glob = Cell_j.TransformationParams;


                    double[] xNodes = Vtx_glob.ExtractSubArrayShallow(-1, 0).To1DArray();
                    double[] yNodes = Vtx_glob.ExtractSubArrayShallow(-1, 1).To1DArray();
                    double xC = xNodes.Sum() / xNodes.Length;
                    double yC = yNodes.Sum() / yNodes.Length;
                    for (int k = 0; k < xNodes.Length; k++) {
                        double dx = xNodes[k] - xC;
                        dx *= 0.95;
                        xNodes[k] = xC + dx;

                        double dy = yNodes[k] - yC;
                        dy *= 0.95;
                        yNodes[k] = yC + dy;
                    }

                    gp.PlotXY(xNodes, yNodes, title: j.ToString(), format: (new PlotFormat(Style: Styles.LinesPoints, lineColor: ((LineColors)j))));



                }

                gp.Execute();


                Console.ReadKey();
            }

        }

    }
}
