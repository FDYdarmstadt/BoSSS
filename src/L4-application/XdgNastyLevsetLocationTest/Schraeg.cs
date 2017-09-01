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
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XdgNastyLevsetLocationTest {

    /// <summary>
    /// a level-set, from the top-left corner to the lower-right corner of the grid
    /// </summary>
    class Schraeg : ITest {

        public Schraeg(double slope = 0.0, double intercept = 0.0)
            : this(new double[] { slope }, new double[] { intercept })
        { }

        public Schraeg(IEnumerable<double> slopes, IEnumerable<double> intercepts) {

            Testcases = new List<Tuple<double, double>>(slopes.Count()*intercepts.Count());

            for (int i = 0; i < slopes.Count(); i++) {
                for (int k = 0; k < intercepts.Count(); k++) {
                    Testcases.Add(new Tuple<double, double>(slopes.ElementAt(i), intercepts.ElementAt(k)));
                }
            }

        }

        List<Tuple<double, double>> Testcases;
        int TestCaseCount = -1;

        double slope = double.NaN;
        double intercept = double.NaN;

        public bool NextTestCase() {
            if (TestCaseCount >= Testcases.Count - 1) {
                return false;
            } else {
                TestCaseCount++;
                slope = Testcases[TestCaseCount].Item1;
                intercept = Testcases[TestCaseCount].Item2;
                //Console.WriteLine("slope = {0}, intercept = {1}", this.slope, this.intercept);
                return true;
            }
        }

        public void ResetTest() {
            TestCaseCount = -1;
            slope = double.NaN;
            intercept = double.NaN;
        }

        double XMIN = -3;
        double XMAX = 3;
        double YMIN = -2;
        double YMAX = 2;

        int NN = 3;

        double dx {
            get {
                return (XMAX - XMIN)/((double)(NN - 1));
            }
        }
        double dy {
            get {
                return (YMAX - YMIN)/((double)(NN - 1));
            }
        }

        public GridCommons GetGrid() {
            double[] Xnodes = GenericBlas.Linspace(XMIN, XMAX, NN);
            double[] Ynodes = GenericBlas.Linspace(YMIN, YMAX, NN);
            return Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
        }

        public bool VolumeTestSupported {
            get { return true; }
        }

        public bool EdgeTestSupported {
            get { return true; }
        }

        public bool LevelsetTestSupported {
            get { return true; }
        }

        public double GetLevelSet(double x, double y) {
            return (x*(1.0 + this.slope)/XMAX + (y + this.intercept)/YMAX);
        }

        double SpeciesSign(string sp) {
            switch (sp) {
                case "A": return -1;
                case "B": return +1;
                default: throw new ArgumentException();
            }
        }

        public double CellVolume(int jCell, string species, GridData g) {
            Debug.Assert(g.SpatialDimension == 2);
            MultidimensionalArray CellCenter = MultidimensionalArray.Create(1,2);
            g.TransformLocal2Global(g.Cells.GetRefElement(jCell).Center, CellCenter, jCell);

            double phi = this.GetLevelSet(CellCenter[0,0], CellCenter[0,1]);
            double sign = phi*SpeciesSign(species);

            if(Math.Abs(phi) <= 1.0e-5)
                return 0.5*dx*dy;

            return (dx*dy)*(sign.Heaviside());
        }

        public double EdgeArea(int iEdge, string species, GridData g) {
            Debug.Assert(g.SpatialDimension == 2);
            MultidimensionalArray EdgeVerticesGlobal = MultidimensionalArray.Create(2, 2);
            MultidimensionalArray EdgeCenterGlobal = MultidimensionalArray.Create(1, 2);

            NodeSet EdgeCenter = g.Edges.GetRefElement(iEdge).Center;
            NodeSet EdgeVertices = g.Edges.GetRefElement(iEdge).Vertices;
            g.TransformLocal2Global(EdgeCenter.GetVolumeNodeSet(g, g.Edges.Edge2CellTrafoIndex[iEdge, 0]), EdgeCenterGlobal, g.Edges.CellIndices[iEdge, 0]);
            g.TransformLocal2Global(EdgeVertices.GetVolumeNodeSet(g, g.Edges.Edge2CellTrafoIndex[iEdge, 0]), EdgeVerticesGlobal, g.Edges.CellIndices[iEdge, 0]);


            double phi = this.GetLevelSet(EdgeCenterGlobal[0, 0], EdgeCenterGlobal[0, 1]);
            double sign = phi*SpeciesSign(species);

            if (sign <= 0)
                return 0;

            double vecX = EdgeVerticesGlobal[1, 0] - EdgeVerticesGlobal[0, 0];
            double vecY = EdgeVerticesGlobal[1, 1] - EdgeVerticesGlobal[0, 1];

            if (Math.Abs(vecX) > Math.Abs(vecY))
                // horizontal edge
                return dx;
            else
                // vertical edge
                return dy;
        }

        public double LevelsetArea(int jCell, GridData g) {
            Debug.Assert(g.SpatialDimension == 2);
            MultidimensionalArray CellCenter = MultidimensionalArray.Create(1, 2);
            g.TransformLocal2Global(g.Cells.GetRefElement(jCell).Center, CellCenter, jCell);

            double phi = this.GetLevelSet(CellCenter[0, 0], CellCenter[0, 1]);

            if (Math.Abs(phi) <= 1.0e-5)
                return Math.Sqrt(dx*dx + dy*dy);
            else
                return 0.0;
        }


        public bool BoundaryPlusLevelsetTestSupported {
            get { return false; }
        }

        public double CellBoundaryPlusLevelsetArea(int jCell, string species, GridData g) {
            throw new NotImplementedException();
        }


    }
}
