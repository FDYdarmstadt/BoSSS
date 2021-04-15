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
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using System.Globalization;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;


namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// Basic test case for rotating the leve-set field
    /// in a constant velocity field
    /// </summary>
    class LevelSetRotationTest : IXNSElsTest {

        public bool TestImmersedBoundary => false;

        /// <summary>
        /// nix
        /// </summary>
        public Func<double[], double, double> GetPhi2() {
            throw new NotImplementedException(); // will never be called, as long as 'TestImmersedBoundary' == false;
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            throw new NotImplementedException();
        }

        double L = 1.0;
        double Radius = 0.4;
        double Uscale = 1;

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetRotationTest(int spatDim, int LevelSetDegree) {
            this.SpatialDimension = spatDim;
            this.LevelsetPolynomialDegree = LevelSetDegree;
        }


        public double dt {
            get {
                return -1.0;    // will be set in LevelSetTest() according to level set cfl 
            }
        }

        /// <summary>
        /// computes the timestep size according to the level-set CFL condition
        /// </summary>
        /// <param name="Resolution"></param>
        /// <param name="LSdegree"></param>
        /// <returns></returns>
        public double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel) {
            int gridCells1D = (9 * Resolution) * (AMRlevel + 1);
            double h = 1.0 * L / (double)gridCells1D;
            double dt = h / Uscale;
            dt /= (double)(LSdegree * LSdegree);

            return dt;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double getEndTime() {
            return 1.0; //TODO
        }

        /// <summary>
        /// creates a square in 2D and a cube in 3D
        /// </summary>
        /// <param name="Resolution"></param>
        /// <returns></returns>
        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();
            
            //double L = 1.0;
            var xNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);
            var yNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);
            var zNodes = GenericBlas.Linspace(-L, L, 9 * Resolution + 1);

            GridCommons grd;
            switch (SpatialDimension) {
                case 2:
                grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                break;
                case 3:
                grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                break;
                default:
                throw new ArgumentOutOfRangeException();
            }

            grd.EdgeTagNames.Add(1, "velocity_inlet");

            grd.DefineEdgeTags(delegate (double[] X) {
                double x = X[0], y = X[1];
                if (Math.Abs(x - (-L)) <= 1.0e-7)
                    // velocity inlet
                    return (byte)1;
                if (Math.Abs(x - (L)) <= 1.0e-7)
                    //  velocity inlet (
                    return (byte)1;
                if (Math.Abs(y - (-L)) <= 1.0e-7)
                    // bottom wall
                    return (byte)1;
                if (Math.Abs(y - (L)) <= 1.0e-7)
                    // top wall
                    return (byte)1;
                if (SpatialDimension == 3) {
                    double z = X[2];
                    if (Math.Abs(z - (-L)) <= 1.0e-7)
                        // back wall
                        return (byte)1;
                    if (Math.Abs(z - (L)) <= 1.0e-7)
                        // front wall
                        return (byte)1;
                }

                throw new ArgumentOutOfRangeException();
            });

            return grd;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("velocity_inlet", new AppControl.BoundaryValueCollection());
            config["velocity_inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => Uscale * X[1]);
            config["velocity_inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Uscale * X[1]);

            config["velocity_inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#A",
                (X, t) => -Uscale * X[0]);
            config["velocity_inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#B",
                (X, t) => -Uscale * X[0]);

            if (SpatialDimension == 3) {

                // todo

                //config.Add("velocity_inlet_back", new AppControl.BoundaryValueCollection());
                //config["velocity_inlet_back"].Evaluators.Add(
                //    VariableNames.Velocity_d(0) + "#A",
                //    (X, t) => Ux);
                //config["velocity_inlet_back"].Evaluators.Add(
                //    VariableNames.Velocity_d(0) + "#B",
                //    (X, t) => Ux);

                //config.Add("velocity_inlet_front", new AppControl.BoundaryValueCollection());
                //config["velocity_inlet_front"].Evaluators.Add(
                //    VariableNames.Velocity_d(0) + "#A",
                //    (X, t) => Ux);
                //config["velocity_inlet_front"].Evaluators.Add(
                //    VariableNames.Velocity_d(0) + "#B",
                //    (X, t) => Ux);
            }

            return config;
        }

        bool quadratic = false;

        /// <summary>
        /// Level-Set: at t=0 circle with radius R, at the center (0,0); non-moving.
        /// </summary>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double time) {


                double x = X[0], y = X[1];
                //x -= time * Uscale;

                double dist;
                switch (SpatialDimension) {
                    case 2: { dist = Math.Sqrt(x * x + y * y); break; }
                    case 3: { double z = X[2]; dist = Math.Sqrt(x * x + y * y + z * z); break; }
                    default:
                    throw new ArgumentOutOfRangeException();
                }

                return quadratic ? (Radius * Radius - dist * dist) : Radius - dist;

            };
        }

        /// <summary>
        /// velocity: constant rotational velocity field.
        /// </summary>
        public Func<double[], double, double> GetU(string species, int d) {
            switch (SpatialDimension) {
                case 2:
                if (d == 0)
                    return ((_2D)((x, y) => Uscale * y)).Convert_xy2X().Convert_X2Xt();
                else if (d == 1)
                    return ((_2D)((x, y) => -Uscale * x)).Convert_xy2X().Convert_X2Xt();
                else
                    throw new ArgumentOutOfRangeException();
                case 3:
                if (d == 0)
                    return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                else if (d == 1)
                    return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                else if (d == 2)
                    return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                else
                    throw new ArgumentOutOfRangeException();
                default:
                throw new ArgumentOutOfRangeException();
            }

        }

        /// <summary>
        /// pressure: 
        /// </summary>
        public Func<double[], double, double> GetPress(string species) {
            switch (SpatialDimension) {
                case 2: {
                    switch (species) {
                        case "A": return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                        case "B": return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                        default: throw new ArgumentException();
                    }
                }
                case 3: {
                    switch (species) {
                        case "A": return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                        case "B": return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                        default: throw new ArgumentException();
                    }
                }
                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// No bulk force.
        /// </summary>
        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }


        /// <summary>
        /// no surface tension effects during level set motion 
        /// </summary>
        public double Sigma {
            get {
                return 0.0;
            }
        }

        /// <summary>
        /// considering level set movement in single fluid (A=B)
        /// </summary>
        public double rho_A {
            get { return 1.0; }
        }

        /// <summary>
        ///  considering level set movement in single fluid (A=B)
        /// </summary>
        public double rho_B {
            get { return 1.0; }
        }

        /// <summary>
        ///  considering level set movement in single fluid (A=B)
        /// </summary>
        public double mu_A {
            get { return 1.0; }
        }

        /// <summary>
        ///  considering level set movement in single fluid (A=B)
        /// </summary>
        public double mu_B {
            get { return 1.0; }
        }


        public bool Material {
            get { return true; }
        }

        public bool steady {
            get { return false; }
        }

        public bool IncludeConvection {
            get { return false; }
        }

        public int LevelsetPolynomialDegree {
            get;
            private set;
        }

        public double[] AcceptableL2Error {
            get {
                  return new double[] { 1.0e-6, 1.0e-6, 1.0e-1 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return  new double[] { };
            }
        }


        /// <summary>
        /// returns spatial dimension
        /// </summary>
        public int SpatialDimension {
            get;
            private set;
        }

    }

}
