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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils.Formula;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;

namespace BoSSS.Application.XNSE_Solver.Tests {


    class PolynomialTestForConvection : ITest {


         public Func<double[], double, double> GetPhi() {
            return ((_2D)(delegate (double x, double y) {
                return -y - (2.0 / 5.0) * x + (1.0 / 10.0) * x * x * x;
            })).Convert_xy2X().Convert_X2Xt();
        }

        public Func<double[], double, double> GetU(string species, int d) {
            if (d == 0) {
                return ((_2D)((x, y) => x * x)).Convert_xy2X().Convert_X2Xt();
            } else if (d == 1) {
                return ((_2D)((x, y) => -2*x*y)).Convert_xy2X().Convert_X2Xt();
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        public double dt {
            get {
                return 1.0;
            }
        }
        
        public GridCommons CreateGrid() {

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 6), GenericBlas.Linspace(-2, 2, 6));
            //var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 2), GenericBlas.Linspace(-2, 2, 2));

            grd.EdgeTagNames.Add(1, "Velocity_Inlet");
            grd.DefineEdgeTags(delegate(double[] _X) {
                return 1;
            });

            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("Velocity_Inlet", new AppControl.BoundaryValueCollection());
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => X[0] * X[0]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#A",
                (X, t) => -2.0 * X[0] * X[1]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => X[0] * X[0]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#B",
                (X, t) => -2.0 * X[0] * X[1]);

            return config;
        }

        public Func<double[], double, double> GetPress(string species) {
            return (X,t) => 0.0;
        }
        
        public double rho_A {
            get { return 100.0; }
        }

        public double rho_B {
            get { return 0.033; }
        }

        public double mu_A {
            get { return 0.0; }
        }

        public double mu_B {
            get { return 0.0; }
        }

        public ScalarFunction GetS(double time) {
            return ((_2D)((x, y) => 0)).Vectorize();
        }

        public Func<double[], double> GetF(string species, int d) {
            if (species == "A") {
                if (d == 0) {
                    return ((_2D)((x, y) => 2.0 * x * x * x)).Convert_xy2X();
                } else if (d == 1) {
                    return ((_2D)((x, y) => 2.0 * x * x * y)).Convert_xy2X();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            } else if (species == "B") {
                if (d == 0) {
                    return ((_2D)((x, y) => 2.0 * x * x * x)).Convert_xy2X();
                } else if (d == 1) {
                    return ((_2D)((x, y) => 2.0 * x * x * y)).Convert_xy2X();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        public double Sigma {
            get { return 0.0; }
        }

        public bool Material {
            get { return true; }
        }

        public bool steady {
            get { return true; }
        }

        public bool IncludeConvection {
            get { return true; }
        }

        public int LevelsetPolynomialDegree {
            get { return 4; }
        }

        public double[] AcceptableL2Error {
            get { return new double[] { 1.0e-5, 1.0e-5, 1.0e-5 }; }
        }

        public double[] AcceptableResidual {
            get { return new double[] { 1.0e-6, 1.0e-6, 1.0e-6 }; }
        }

        public int SpatialDimension {
            get {
                return 2;
            }
        }

    }
}
