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
    /// Basic test case for advecting the level-set field
    /// in a constant velocity field
    /// </summary>
    public class LevelSetAdvectionTest : IXNSElsTest {

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

        double Ux = 1;
        double Uy = 0.5;
        double Uz = 0.7;

        double T = 0.4; // (Radius / Ux)

        double sign;

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetAdvectionTest(int spatDim, int LevelSetDegree, bool reversed) {
            this.SpatialDimension = spatDim;
            this.LevelsetPolynomialDegree = LevelSetDegree;
            sign = (reversed) ? -1 : 1;
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
            double dt = h / Ux;
            dt /= (double)(LSdegree * LSdegree);

            int timesteps = (int)Math.Ceiling(T / dt); // make sure that the singular point in time is exactly on one timestep

            return (T/(double)timesteps);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double getEndTime() {
            return (Radius / Ux);
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

            grd.EdgeTagNames.Add(1, "velocity_inlet_prescribed");


            grd.DefineEdgeTags(delegate (double[] X) {
                double x = X[0], y = X[1];
                if (Math.Abs(x - (-L)) <= 1.0e-8)
                    // velocity inlet left
                    return (byte)1;
                if (Math.Abs(x - (L)) <= 1.0e-8)
                    // velocity inlet right
                    return (byte)1;
                if (Math.Abs(y - (-L)) <= 1.0e-8)
                    // velocity inlet bottom
                    return (byte)1;
                if (Math.Abs(y - (L)) <= 1.0e-8)
                    // velocity inlet top
                    return (byte)1;
                if (SpatialDimension == 3) {
                    double z = X[2];
                    if (Math.Abs(z - (-L)) <= 1.0e-8)
                        // velocity inlet back
                        return (byte)1;
                    if (Math.Abs(z - (L)) <= 1.0e-8)
                        // velocity inlet front
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

            config.Add("velocity_inlet_prescribed", new AppControl.BoundaryValueCollection());
            config["velocity_inlet_prescribed"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => Math.Sign((Radius / Ux) - t) * sign * Ux);
            config["velocity_inlet_prescribed"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Math.Sign((Radius / Ux) - t) * sign * Ux);
            config["velocity_inlet_prescribed"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#A",
                (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uy);
            config["velocity_inlet_prescribed"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#B",
                (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uy);
            if (this.SpatialDimension == 3) {
                config["velocity_inlet_prescribed"].Evaluators.Add(
                    VariableNames.Velocity_d(2) + "#A",
                    (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uz);
                config["velocity_inlet_prescribed"].Evaluators.Add(
                    VariableNames.Velocity_d(2) + "#B",
                    (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uz);
            }


            return config;
        }

        bool quadratic = false;

        /// <summary>
        /// Level-Set movement.
        /// </summary>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double time) {

                double x0 = -sign * (Radius / 2.0);
                double y0 = -sign * (Radius / 4.0);
                double x = X[0] - x0, y = X[1] - y0;
                x -= Math.Sign((Radius / Ux) - time) * time * sign * Ux;
                y -= Math.Sign((Radius / Ux) - time) * time * sign * Uy;

                double dist;
                switch (SpatialDimension) {
                    case 2: { dist = Math.Sqrt(x * x + y * y); break; }
                    case 3: { double z0 = -sign * (Radius / 3.0); 
                        double z = X[2]; z -= Math.Sign((Radius / Ux) - time) * time * sign * Uz;
                        dist = Math.Sqrt(x * x + y * y + z * z); break; }
                    default:
                    throw new ArgumentOutOfRangeException();
                }

                return quadratic ? (Radius * Radius - dist * dist) : Radius - dist;

            };
        }

        /// <summary>
        /// velocity: constant velocity field with Ux in x-direction.
        /// </summary>
        public Func<double[], double, double> GetU(string species, int d) {
            switch (SpatialDimension) {
                case 2:
                    if (d == 0)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * sign * Ux; 
                    else if (d == 1)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uy;
                    else
                        throw new ArgumentOutOfRangeException();
                case 3:
                    if (d == 0)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * sign * Ux;
                    else if (d == 1)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uy;
                    else if (d == 2)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * sign * Uz;
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
                return new double[] { 1.0e-3, 1.0e-3 }; //, 1.0e-1, 1.0e-1, 1.0e-2, 1.0e-2 };
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


    /// <summary>
    /// Basic test case for advecting the level-set field
    /// in a constant velocity field
    /// </summary>
    class LevelSetAdvectionOnWallTest : IXNSElsTest {

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

        double Ux = 1;
        //double Uy = 0.5;
        double Uz = 0.7;

        double T = 0.8; // 2 * (Radius / Ux)

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetAdvectionOnWallTest(int spatDim, int LevelSetDegree) {
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
            double dt = h / Ux;
            dt /= (double)(LSdegree * LSdegree);

            int timesteps = (int)Math.Ceiling((T / 2.0) / dt); // make sure that the singular point in time is exactly on one timestep

            return ((T / 2.0) / (double)timesteps);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double getEndTime() {
            return 2.0 * (Radius / Ux);
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
            var yNodes = GenericBlas.Linspace(0, L, 4 * Resolution + 1);
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

            grd.EdgeTagNames.Add(1, "velocity_inlet_prescribed");
            grd.EdgeTagNames.Add(2, "navierslip_linear");

            grd.DefineEdgeTags(delegate (double[] X) {
                double x = X[0], y = X[1];
                if (Math.Abs(x - (-L)) <= 1.0e-8)
                    // velocity inlet left
                    return (byte)1;
                if (Math.Abs(x - (L)) <= 1.0e-8)
                    // velocity inlet right
                    return (byte)1;
                if (Math.Abs(y - (0)) <= 1.0e-8)
                    // velocity inlet bottom
                    return (byte)2;
                if (Math.Abs(y - (L)) <= 1.0e-8)
                    // velocity inlet top
                    return (byte)1;
                if (SpatialDimension == 3) {
                    double z = X[2];
                    if (Math.Abs(z - (-L)) <= 1.0e-8)
                        // velocity inlet back
                        return (byte)1;
                    if (Math.Abs(z - (L)) <= 1.0e-8)
                        // velocity inlet front
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

            config.Add("velocity_inlet_prescribed", new AppControl.BoundaryValueCollection());
            config["velocity_inlet_prescribed"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => Math.Sign((Radius / Ux) - t) * Ux);
            config["velocity_inlet_prescribed"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Math.Sign((Radius / Ux) - t) * Ux);

            if (this.SpatialDimension == 3) {
                config["velocity_inlet_prescribed"].Evaluators.Add(
                    VariableNames.Velocity_d(2) + "#A",
                    (X, t) => Math.Sign((Radius / Ux) - t) * Uz);
                config["velocity_inlet_prescribed"].Evaluators.Add(
                    VariableNames.Velocity_d(2) + "#B",
                    (X, t) => Math.Sign((Radius / Ux) - t) * Uz);
            }


            return config;
        }

        bool quadratic = false;

        /// <summary>
        /// Level-Set movement.
        /// </summary>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double time) {

                double x0 = -Radius / 2.0;
                double x = X[0] - x0, y = X[1];
                x -= Math.Sign((Radius / Ux) - time) * time * Ux;

                double dist;
                switch (SpatialDimension) {
                    case 2: { dist = Math.Sqrt(x * x + y * y); break; }
                    case 3: {
                        double z0 = -Radius / 3.0;
                        double z = X[2]; z -= Math.Sign((Radius / Ux) - time) * time * Uz;
                        dist = Math.Sqrt(x * x + y * y + z * z); break;
                    }
                    default:
                        throw new ArgumentOutOfRangeException();
                }

                return quadratic ? (Radius * Radius - dist * dist) : Radius - dist;

            };
        }

        /// <summary>
        /// velocity: constant velocity field with Ux in x-direction.
        /// </summary>
        public Func<double[], double, double> GetU(string species, int d) {
            switch (SpatialDimension) {
                case 2:
                    if (d == 0)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * Ux;
                    else if (d == 1)
                        return (X, t) => 0.0;
                    else
                        throw new ArgumentOutOfRangeException();
                case 3:
                    if (d == 0)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * Ux;
                    else if (d == 1)
                        return (X, t) => 0.0;
                    else if (d == 2)
                        return (X, t) => Math.Sign((Radius / Ux) - t) * Uz;
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
                        case "A": return (X, t) => 0.0;
                        case "B": return (X, t) => 0.0;
                        default: throw new ArgumentException();
                    }
                }
                case 3: {
                    switch (species) {
                        case "A": return (X, t) => 0.0;
                        case "B": return (X, t) => 0.0;
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
                return new double[] { };
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
