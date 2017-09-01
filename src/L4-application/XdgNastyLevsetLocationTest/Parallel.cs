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
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XdgNastyLevsetLocationTest {
    class Parallel : ITest {

        public Parallel(double slope = 0.0, double intercept = 0.0)
            : this(new double[] { slope }, new double[] { intercept })
        { }

        public Parallel(IEnumerable<double> slopes, IEnumerable<double> intercepts) {

            Testcases = new List<Tuple<double, double>>(slopes.Count()*intercepts.Count());

            for (int i = 0; i < slopes.Count(); i++) {
                for (int k = 0; k < intercepts.Count(); k++) {
                    Testcases.Add(new Tuple<double, double>(slopes.ElementAt(i), intercepts.ElementAt(k)));
                }
            }

        }

        List<Tuple<double, double>> Testcases;
        int TestCaseCount = -1;

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


        double slope = double.NaN;
        double intercept = double.NaN;


        public GridCommons GetGrid() {
            double[] xNodes = new double[] { -1, 1 };
            double[] yNodes = new double[] { -2, 0, 2 };
            return Grid2D.Cartesian2DGrid(xNodes, yNodes);
        }

        double dx = 2;
        double dy = 2;

        public bool VolumeTestSupported {
            get { return true; }
        }

        public bool EdgeTestSupported {
            get { return false; }
        }

        public bool LevelsetTestSupported {
            get { return false; }
        }

        public double GetLevelSet(double x, double y) {
            return (y - intercept) + x*slope;
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
            MultidimensionalArray CellCenter = MultidimensionalArray.Create(1, 2);
            g.TransformLocal2Global(g.Cells.GetRefElement(jCell).Center, CellCenter, jCell);

            double phi = this.GetLevelSet(CellCenter[0, 0], CellCenter[0, 1]);
            double sign = phi*SpeciesSign(species);

            if (sign > 0)
                return dx*dy;
            else
                return 0.0;
        }

        public double EdgeArea(int iEdge, string species, GridData g) {
            throw new NotImplementedException();
        }

        public double LevelsetArea(int jCell, GridData g) {
            throw new NotImplementedException();
        }

        public bool BoundaryPlusLevelsetTestSupported {
            get {
                return false;
            }
        }

        public double CellBoundaryPlusLevelsetArea(int jCell, string species, GridData g) {
            throw new NotImplementedException();
            //Debug.Assert(g.SpatialDimension == 2);
            //MultidimensionalArray CellCenter = MultidimensionalArray.Create(1, 2);
            //g.TransformLocal2Global(MultidimensionalArray.Create(1, 2), CellCenter, jCell);

            //double phi = this.GetLevelSet(CellCenter[0, 0], CellCenter[0, 1]);
            //double sign = phi*SpeciesSign(species);

            //if (sign > 0)
            //    return 2.0*(dx + dy);
            //else
            //    return 0.0;
        }
    }
}
