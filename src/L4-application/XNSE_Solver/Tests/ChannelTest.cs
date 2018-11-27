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
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
//using BoSSS.Solution.Utils.Formula;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// a periodic channel flow, with water on bottom and air on top
    /// </summary>
    class ChannelTest : ITest {

        public bool Material {
            get {
                return true;
            }
        }

        public bool steady {
            get {
                return true;
            }
        }

        public bool IncludeConvection {
            get {
                return true;
            }
        }

        /// <summary>
        /// the zero-level-set is identical to the x-axis
        /// </summary>
        /// <param name="time"></param>
        /// <returns></returns>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double t) {
                var Coord = ROTinv.Transform(X[0], X[1]);
                double x = Coord[0];
                double y = Coord[1];

                return -y;
            };
        }

        
        public int LevelsetPolynomialDegree {
            get {
                return 1;
            }
        }

        public ChannelTest(double angle) {
            //double angle = 0.0;
            //double angle = 60.0 * Math.PI / 180.0;
            ROT = AffineTrafo.Some2DRotation(angle);
            ROTinv = ROT.Invert();
        }

        AffineTrafo ROT;

        AffineTrafo ROTinv;

        bool periodic = true;


        public Func<double[], double> GetF(string species, int d) {

            double rho = double.NaN;
            switch(species) {
                case "A": rho = rho_A; break;
                case "B": rho = rho_B; break;
                throw new ArgumentException();
            }

            if (mu_A == 0.0 && mu_B == 0) {
                return (X => 0.0);
            } else {

                double sc = Math.Min(this.mu_A, this.mu_B);
                double[] Fvec = new double[] { (1.0) * sc, 0 };
                var FvecT = ROT.Transform(Fvec);

                return (X => FvecT[d]/rho);
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            double a2, a1, a0;

            if (species == "A") {
                ParabolaCoeffs_A(out a2, out a1, out a0);
            } else if (species == "B") {
                ParabolaCoeffs_B(out a2, out a1, out a0);
            } else
                throw new ArgumentException();



            return ((_2D)(delegate (double _x, double _y) {

                var Coord = ROTinv.Transform(_x, _y);
                double y = Coord[1];

                Debug.Assert(Coord[0] >= -2);
                Debug.Assert(Coord[0] <= +2);
                Debug.Assert(Coord[1] >= -1);
                Debug.Assert(Coord[1] <= +1);

                double u = (a0 + a1 * y + a2 * y * y);
                //double u = 1 - y*y;
                var UT = ROT.Transform(u, 0.0);

                return UT[d];

            })).Convert_xy2X().Convert_X2Xt();
        }

        private void ParabolaCoeffs_B(out double a2, out double a1, out double a0) {

            double muA = this.mu_A;
            double muB = this.mu_B;

            if (muA <= 0.0 && muB <= 0.0) {
                muA = 0.001;
                muB = 0.001;
            }
            if (muA <= 0.0 != muB <= 0.0)
                throw new NotSupportedException();

            double sc = Math.Min(muA, muB);

            a0 = 1.0 / (muA + muB);
            a1 = (muB - muA) / (2.0 * (muB + muA) * muB);
            a2 = (-1.0) / (2 * muB);

            a0 *= sc;
            a1 *= sc;
            a2 *= sc;
        }

        private void ParabolaCoeffs_A(out double a2, out double a1, out double a0) {

            double muA = this.mu_A;
            double muB = this.mu_B;

            if (muA <= 0.0 && muB <= 0.0) {
                muA = 0.001;
                muB = 0.001;
            }
            if (muA <= 0.0 != muB <= 0.0)
                throw new NotSupportedException();

            double sc = Math.Min(muA, muB);

            a0 = 1.0 / (muA + muB);
            a1 = (muB - muA) / (2.0 * (muB + muA) * muA);
            a2 = (-1.0) / (2 * muA);

            a0 *= sc;
            a1 *= sc;
            a2 *= sc;
        }

        

        public double dt {
            get {
                return 1.0;
            }
        }

        public GridCommons CreateGrid() {

            //var yNodes = GenericBlas.Linspace(-1, 1, 9);

            //var yNodes1 = yNodes.GetSubVector(0, 4);
            //var yNodes2 = yNodes.GetSubVector(5, 4);
            //var _yNodes = ArrayTools.Cat(yNodes1, yNodes2);
            //var _yNodes = yNodes;

            var _yNodes = new double[] { -1, -0.8, -0.6, 0.6, 1.0 };
            //var _yNodes = GenericBlas.Linspace(-1, 1, 6);

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 4), _yNodes, periodicX: periodic);
            if (periodic) {
                grd.EdgeTagNames.Add(1, "wall_top");
                grd.EdgeTagNames.Add(2, "wall_bottom");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom wall
                        return 2;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top wall
                        return 1;


                    throw new ArgumentOutOfRangeException();
                    //return 1;
                });

                Console.WriteLine("ChannelTest, periodic.");

            } else {

                grd.EdgeTagNames.Add(1, "wall_top");
                grd.EdgeTagNames.Add(2, "wall_bottom");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");

                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom wall
                        return 2;

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top wall
                        return 1;

                    if (Math.Abs(x - (-2)) < 1.0e-6)
                        // inlet
                        return 3;

                    if (Math.Abs(x - (2)) < 1.0e-6)
                        // outlet
                        return 4;

                    throw new ArgumentOutOfRangeException();
                    //return 1;
                });
            }
            var grdT = grd.Transform(ROT);

            return grdT;
        }


        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("wall_top", new AppControl.BoundaryValueCollection());
            config.Add("wall_bottom", new AppControl.BoundaryValueCollection());

            if (!periodic) {
                if (!this.ROT.ApproximateEquals(AffineTrafo.Some2DRotation(0.0)))
                    throw new NotSupportedException();

                double A_a0, A_a1, A_a2, B_a0, B_a1, B_a2;
                this.ParabolaCoeffs_A(out A_a2, out A_a1, out A_a0);
                this.ParabolaCoeffs_B(out B_a2, out B_a1, out B_a0);

                config.Add("velocity_inlet", new AppControl.BoundaryValueCollection());
                config["velocity_inlet"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#A",
                    (X, t) => A_a0 + A_a1 * X[1] + A_a2 * X[1] * X[1]);
                config["velocity_inlet"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#B",
                    (X, t) => B_a0 + B_a1 * X[1] + B_a2 * X[1] * X[1]);

                config.Add("Pressure_Outlet", new AppControl.BoundaryValueCollection());
            }

            return config;
        }

        public Func<double[], double, double> GetPress(string species) {
            return (X, t) => 0.0;
        }


        /// <summary>
        /// specific weight, air
        /// </summary>
        public double rho_B {
            get {
                return 1.2;
            }
        }

        /// <summary>
        /// specific weight, water
        /// </summary>
        public double rho_A {
            get {
                return 1000;
            }
        }

        /// <summary>
        /// dynamic viscosity, air
        /// </summary>
        public double mu_B {
            get {
                return 17.1e-3;
            }
        }

        /// <summary>
        /// dynamic viscosity, water
        /// </summary>
        public double mu_A {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// surface tension of water (surface tension has no effect due to the planar interface)
        /// </summary>
        public double Sigma {
            get {
                //return 0.0;
                return 72.75e-3;
            }
        }
                       

        public double[] AcceptableL2Error {
            get {
                return new double[] { 1.0e-7, 1.0e-7, 1.0e-7 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return new double[] { 1.0e-7, 1.0e-7, 1.0e-7 };
            }
        }

        public int SpatialDimension {
            get {
                return 2;
            }
        }

        
    }




}
