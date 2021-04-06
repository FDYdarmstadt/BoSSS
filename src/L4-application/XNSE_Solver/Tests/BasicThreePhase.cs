using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Tests {
    
    /// <summary>
    /// Elementary test for the three-phase (Fluid, Fluid, Solid) capabilities of <see cref="XNSE"/> with two level-sets.
    /// </summary>
    class BasicThreePhase : IXNSETest {

        public bool TestImmersedBoundary => true;

        public Func<double[], double, double> GetPhi2U(int d) {
            if(d==0) {
                return delegate (double[] X, double t) {
                    return 1.0;
                };
            } else {
                return (X, t) => 0.0;
            }
        }

        public BasicThreePhase(double R = 0.8, bool bConvection = true, bool bSteady = true, int spatDim = 2) {
            this.Radius = R;
            this.IncludeConvection = bConvection;
            this.steady = bSteady;
            this.spatialDimension = spatDim;
        }

        int spatialDimension;

        double Radius = 0.8;
        double Ux = 1.0;

        /// <summary>
        /// Level-Set: at t=0 circle with radius R, at the center (0,0); moving with Ux in x-direction.
        /// </summary>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double time) {

                double x = X[0], y = X[1];
                x -= time * Ux;

                switch (spatialDimension) {
                    case 2:
                        return (this.Radius * this.Radius - (x * x + y * y));
                    case 3:
                        double z = X[2];
                        return (this.Radius * this.Radius - (x * x + y * y + z * z));
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }; 
        }



        /// <summary>
        /// Second Level-Set
        /// </summary>
        public Func<double[], double, double> GetPhi2() {
            return delegate (double[] X, double time) {

                double x = X[0], y = X[1];
                double x0 = -1.45;
                x0 += time * Ux;

                return -x + x0;
            }; 
        }



        public bool IncludeConvection {
            get;
            private set;
        }

        /// <summary>
        /// velocity: not known, except for t=0;
        /// </summary>
        public Func<double[], double, double> GetU(string species, int d) {
            switch (spatialDimension) {
                case 2:
                    if (d == 0)
                        return ((_2D)((x, y) => Ux)).Convert_xy2X().Convert_X2Xt();
                    else if (d == 1)
                        return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                    else
                        throw new ArgumentOutOfRangeException();
                case 3:
                    if (d == 0)
                        return ((_3D)((x, y, z) => Ux)).Convert_xyz2X().Convert_X2Xt();
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

      
        public double dt {
            get {
                return 0.1;
            }
        }


        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();


            var xNodes = GenericBlas.Linspace(-1.5, 1.5, 9 * Resolution + 1);
            var yNodes = GenericBlas.Linspace(-1.5, 1.5, 9 * Resolution + 1);
            var zNodes = GenericBlas.Linspace(-1.5, 1.5, 9 * Resolution + 1);

            GridCommons grd;
            switch (spatialDimension) {
                case 2:
                    grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                    break;
                case 3:
                    grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
            } 

            grd.EdgeTagNames.Add(1, "wall_top");
            grd.EdgeTagNames.Add(2, "wall_bottom");
            grd.EdgeTagNames.Add(3, "velocity_inlet_pos");
            grd.EdgeTagNames.Add(4, "velocity_inlet_neg");
            if (spatialDimension == 3) {
                grd.EdgeTagNames.Add(5, "wall_front");
                grd.EdgeTagNames.Add(6, "wall_back");
            }

            grd.DefineEdgeTags(delegate(double[] X) {
                double x = X[0], y = X[1];
                if (Math.Abs(x - (-1.5)) <= 1.0e-7)
                    // velocity inlet
                    return (byte)3;
                if (Math.Abs(x - (1.5)) <= 1.0e-7)
                    //  velocity inlet (
                    return (byte)4;
                if (Math.Abs(y - (-1.5)) <= 1.0e-7)
                    // bottom wall
                    return (byte)2;
                if (Math.Abs(y - (+1.5)) <= 1.0e-7)
                    // top wall
                    return (byte)1;
                if (spatialDimension == 3) {
                    double z = X[2];
                    if (Math.Abs(z - (-1.5)) <= 1.0e-7)
                        // back wall
                        return (byte)6;
                    if (Math.Abs(z - (+1.5)) <= 1.0e-7)
                        // front wall
                        return (byte)5;
                }

                throw new ArgumentOutOfRangeException();
            });

            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("velocity_inlet_pos", new AppControl.BoundaryValueCollection());
            config["velocity_inlet_pos"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => Ux);
            config["velocity_inlet_pos"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Ux);

            config.Add("velocity_inlet_neg", new AppControl.BoundaryValueCollection());
            config["velocity_inlet_neg"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => Ux);
            config["velocity_inlet_neg"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Ux);

            config.Add("wall_bottom", new AppControl.BoundaryValueCollection());
            config["wall_bottom"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => 1.0);
            config["wall_bottom"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Ux);

            config.Add("wall_top", new AppControl.BoundaryValueCollection());
            config["wall_top"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => Ux);
            config["wall_top"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => Ux);

            if (spatialDimension == 3) {

                config.Add("wall_back", new AppControl.BoundaryValueCollection());
                config["wall_back"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#A",
                    (X, t) => 1.0);
                config["wall_back"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#B",
                    (X, t) => Ux);

                config.Add("wall_front", new AppControl.BoundaryValueCollection());
                config["wall_front"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#A",
                    (X, t) => Ux);
                config["wall_front"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#B",
                    (X, t) => Ux);
            }

            return config;
        }


        /// <summary>
        /// pressure: not known, except for t = 0;
        /// </summary>
        public Func<double[], double, double> GetPress(string species) {
            switch (spatialDimension) {
                case 2: {
                        switch (species) {
                            case "A": return ((_2D)((x, y) => (-1.0 / this.Radius) * Sigma)).Convert_xy2X().Convert_X2Xt();
                            case "B": return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                            default: throw new ArgumentException();
                        }
                    }
                case 3: {
                        switch (species) {
                            case "A": return ((_3D)((x, y, z) => (-2.0 / this.Radius) * Sigma)).Convert_xyz2X().Convert_X2Xt();
                            case "B": return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
                            default: throw new ArgumentException();
                        }
                    }
                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

        


        /// <summary>
        /// specific weight, air
        /// </summary>
        public double rho_A {
            get { return 1.2; }
        }

        /// <summary>
        /// specific weight, water
        /// </summary>
        public double rho_B {
            get { return 1000; }
        }

        /// <summary>
        /// dynamic viscosity, air * 10 (must be sufficiently large for stable steady-state solution)
        /// </summary>
        public double mu_A {
            get { return 17.1e-3*10; }
        }

        /// <summary>
        /// dynamic viscosity, water * 10 (must be sufficiently large for stable steady-state solution)
        /// </summary>
        public double mu_B {
            get { return 1.0*10; }
        }

        /// <summary>
        /// surface tension of water
        /// </summary>
        public double Sigma {
            get {
                return 72.75e-3;
                //return 0.0;
            }
        }

      

        /// <summary>
        /// No bulk force.
        /// </summary>
        public Func<double[],double> GetF(string species, int d) {
            return (X => 0.0);
        }

        

        public bool Material {
            get { return true; }
        }

        public bool steady {
            get;
            private set;
        }

        public int LevelsetPolynomialDegree {
            get { return 2; }
        }

        /// <summary>
        ///  tuned for grid resolution 1 (<see cref="CreateGrid(int)"/>), DG degree 2
        /// </summary>
        public double[] AcceptableL2Error {
            get { 
                var r = (spatialDimension == 2) ? new double[] { 1.0e-6, 1.0e-6, 5.0e-3 } : new double[] { 1.0e-6, 1.0e-6, 1.0e-6, 5.0e-3 };
                return r;
            }
        }

        /// <summary>
        /// tuned for grid resolution 1 (<see cref="CreateGrid(int)"/>), DG degree 2
        /// </summary>
        public double[] AcceptableResidual {
            get { 
                var r = (spatialDimension == 2) ? new double[] { 2.0e-5, 2.0e-5, 2.0e-5 } : new double[] { 2.0e-5, 2.0e-5, 2.0e-5, 2.0e-5 };
                if(!steady)
                    r.ScaleV(4);
                return r;
            }
        }

        public int SpatialDimension {
            get {
                return spatialDimension;
            }
        }

    }
}
