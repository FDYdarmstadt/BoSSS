using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Tests {
    
    /// <summary>
    /// Tests if the solver also works with no level-set in the domain
    /// </summary>
    class SinglePhaseChannel : ITest {

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
        /// level-set in Nirvana
        /// </summary>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double t) {
                return -1;
            };
        }

        
        public int LevelsetPolynomialDegree {
            get {
                return 1;
            }
        }

        public SinglePhaseChannel(double angle) {
            //double angle = 0.0;
            //double angle = 60.0 * Math.PI / 180.0;
            ROT = AffineTrafo.Some2DRotation(angle);
            ROTinv = ROT.Invert();
        }

        AffineTrafo ROT;

        AffineTrafo ROTinv;

      

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
            

            return ((_2D)(delegate (double _x, double _y) {

                var Coord = ROTinv.Transform(_x, _y);
                double y = Coord[1];

                double u = (1.0 - y * y);
                var UT = ROT.Transform(u, 0.0);

                return UT[d];

            })).Convert_xy2X().Convert_X2Xt();
        }

       


        public double dt {
            get {
                return 1.0;
            }
        }

        public bool periodic = false;

        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();

            // For Resolution = 1, 2, 3, ...

            int NoOfXcells = (int)Math.Round(Math.Pow(2, Resolution + 3)); // 16, 32, 64, ...
            int NoOfYcells = (int)Math.Round(Math.Pow(2, Resolution + 1)); // 4, 8, 16, ...
            
            

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(0, 10, NoOfXcells + 1), GenericBlas.Linspace(-1, +1, NoOfYcells + 1), periodicX: periodic);
            if (periodic) {
                 grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom wall
                        return "wall_bottom";

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top wall
                        return  "wall_top";

                    throw new ArgumentOutOfRangeException();
                });

                Console.WriteLine("ChannelTest, periodic.");

            } else {
                grd.DefineEdgeTags(delegate (double[] _X) {
                    var X = _X;
                    double x = X[0];
                    double y = X[1];

                    if (Math.Abs(y - (-1)) < 1.0e-6)
                        // bottom wall
                        return "wall_bottom";

                    if (Math.Abs(y - (+1)) < 1.0e-6)
                        // top wall
                        return  "wall_top";

                    if(Math.Abs(x - (0)) < 1.0e-6)
                        // inlet
                        return "Velocity_Inlet";

                    if (Math.Abs(x - (10)) < 1.0e-6)
                        // outlet
                        return "Pressure_Outlet";

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

                config.Add("velocity_inlet", new AppControl.BoundaryValueCollection());
                config["velocity_inlet"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#A",
                    (X, t) => 1.0 - X[1] * X[1]);
                config["velocity_inlet"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#B",
                    (X, t) => 1.0 - X[1] * X[1]);

                config.Add("Pressure_Outlet", new AppControl.BoundaryValueCollection());
            }

            return config;
        }

        public Func<double[], double, double> GetPress(string species) {
            return (X, t) => 0.0;
        }


        /// <summary>
        /// specific weight, set to 1.0
        /// </summary>
        public double rho_B {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// specific weight, set to 1.0
        /// </summary>
        public double rho_A {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// dynamic viscosity, set to 1.0
        /// </summary>
        public double mu_B {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// dynamic viscosity, set to 1.0
        /// </summary>
        public double mu_A {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// 0.0
        /// </summary>
        public double Sigma {
            get {
                return 0.0;
                
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
