using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// Test for the verification of the <see cref="Logging.SphericalHarmonicsLogging"/> post-processing:
    /// An initial level-set is defined using a specific set of modes for the spherical harmonics (<see cref="modes"/>) representation.
    /// It is verified that <see cref="Logging.SphericalHarmonicsLogging"/> is capable of recovering these values for the modes
    /// with a certain accuracy.
    /// </summary>
    /// <remarks>
    /// - The physics of this test-case is irrelevant
    /// - Fk, Oct. 2021; developed for DACH-cooperation project with TU Graz, Prof. Brenn.
    /// </remarks>
    class SphericalHarmonicsTest : IXNSETest {
        public double mu_A => 1e-2;

        public double mu_B => 1e-2;

        public double Sigma => 0.0;

        public bool TestImmersedBoundary => false;

        public int SpatialDimension => 3;

        public double dt => 1e-5;

        public double rho_A => 1000;

        public double rho_B => 1;

        public bool Material => true;

        public bool steady => true;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 3;

        public double[] AcceptableL2Error => new double[] { 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8 };

        public double[] AcceptableResidual => new double[] { 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8 };


        public bool ComputeOnQuarterDomain = false;

        public bool IsRotationalSymmetric = false;


        public GridCommons CreateGrid(int Resolution) {
            //double sz = 2;
            double sz = 0.65;
            
            var nodes = GenericBlas.Linspace(-sz, +sz, 20 * Resolution + 1);
            var grd = Grid3D.Cartesian3DGrid(nodes, nodes, nodes);
            if (ComputeOnQuarterDomain) {
                var nodesPos = GenericBlas.Linspace(0, +sz, 10 * Resolution + 1);
                grd = Grid3D.Cartesian3DGrid(nodesPos, nodesPos, nodes);
            }

            grd.DefineEdgeTags((Vector X) => "wall");


            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("wall", new AppControl.BoundaryValueCollection());
            //config["velocity_inlet_pos"].Evaluators.Add(
            //    VariableNames.Velocity_d(0) + "#A",
            //    (X, t) => Ux);
            //config["velocity_inlet_pos"].Evaluators.Add(
            //    VariableNames.Velocity_d(0) + "#B",
            //    (X, t) => Ux);

            return config;
        }

        public Func<double[], double> GetF(string species, int d) {
            return ((double[] X) => 0.0);
        }

        /*
        // fail on 20^3, (-0.5..0.5)^3
        internal (int l, int m, double Ylm)[] modes = new[] {
            (+0, +0, 1.0),
            (+1, -1, -0.1),
            (+1, +0, -0.2),
            (+1, +1, -0.3),
            (+2,  0, 0.4)
        };
        */

        internal (int l, int m, double Ylm)[] modesNonRotSym = new[] {
            (+0, +0, 1.0*SphericalHarmonics.Get_Nlm(0,0)),
            (+1, -1, -0.1*SphericalHarmonics.Get_Nlm(1,-1)),
            (+1, +0, -0.2*SphericalHarmonics.Get_Nlm(1,0)),
            (+1, +1, -0.3*SphericalHarmonics.Get_Nlm(1,1)),
            (+2,  0, 0.4*SphericalHarmonics.Get_Nlm(2,0))
        };

        internal (int l, int m, double Ylm)[] modesRotSym = new[] {
            (+0, +0, 1.0*SphericalHarmonics.Get_Nlm(0,0)),
            (+1, +0, -0.2*SphericalHarmonics.Get_Nlm(1,0)),
            (+2,  0, 0.4*SphericalHarmonics.Get_Nlm(2,0)),
            //(+3,  0, -0.3*SphericalHarmonics.Get_Nlm(3,0))
        };

        internal (int l, int m, double Ylm)[] modes; 


        /// <summary>
        /// Computes the error in given spherical modes against the values set in the test.
        /// </summary>
        public double ComputeModeError(double[] LinearModes) {
            double ret = 0.0;

            modes = IsRotationalSymmetric ? modesRotSym : modesNonRotSym;
            double[] Shallmodes = new double[Logging.SphericalHarmonicsLogging.SH_dim(modes.Max(tt => tt.l), IsRotationalSymmetric)];
            int LL = LinearModes.Length;
            int I = Math.Max( modes.Max(tt => tt.l), LL);

            for(int i = 0; i < Math.Max(LinearModes.Length, Shallmodes.Length); i++) {
                var (l, m) = Logging.SphericalHarmonicsLogging.SH_mappingInv(i, IsRotationalSymmetric);

                double ma = i < LL ? LinearModes[i] : 0.0;
                double mb = modes.SingleOrDefault(tt => tt.l == l && tt.m == m).Ylm;
                double err = ma - mb;

                Console.WriteLine($"Error on mode {(l, m)}: {err}");
                ret += err.Pow2();
            }

            return ret.Sqrt();
        }



        public Func<double[], double, double> GetPhi() {
            double DistSphere(double[] X, double time) {
                //double x = X[0];
                //double y = X[1];
                //double z = X[2];

                double R = X.L2Norm();
                var (u, v) = SphericalHarmonics.GetAngular(X);

                modes = IsRotationalSymmetric ? modesRotSym : modesNonRotSym;
                double r = 0;
                foreach(var tt in modes)
                    r += SphericalHarmonics.MyRealSpherical(tt.l, tt.m, u, v)*tt.Item3;

                return R - r;
            }

            return DistSphere;
        }

        public Func<double[], double, double> GetPhi2() {
            throw new NotImplementedException();
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            throw new NotImplementedException();
        }

        public Func<double[], double, double> GetPress(string species) {
            return ((double[] X, double t) => 0.0);
        }

        public Func<double[], double, double> GetU(string species, int d) {
            return ((double[] X, double t) => 0.0);
        }
    }
}
