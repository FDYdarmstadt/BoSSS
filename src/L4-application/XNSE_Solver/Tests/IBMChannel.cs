using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Tests {
    
    /// <summary>
    /// Laminar channel flow for the immersed boundary method:
    /// The walls are represented using the immersed boundary and the channel can be rotated.
    /// - For moving walls, a polynomial order of 1 should produce the exact solution under all angles
    /// - For static walls, a polynomial order of 2 should produce the exact solution under all angles
    /// </summary>
    public class IBMChannel : IXNSETest {

        public IBMChannel(double angle, bool movingWalls) {
            m_angle = angle;
            m_movingWalls = movingWalls;
        }

        double m_angle;
        bool m_movingWalls;

        /// <summary>
        /// transformation from the physical coordinates to the rotated coordinate system
        /// </summary>
        double rot_x(double phys_x, double phys_y) {
            return Math.Cos(m_angle) * phys_x + Math.Sin(m_angle) * phys_y;
        }

        /// <summary>
        /// transformation from the physical coordinates to the rotated coordinate system
        /// </summary>
        double rot_y(double phys_x, double phys_y) {
            return -Math.Sin(m_angle) * phys_x + Math.Cos(m_angle) * phys_y;
        }

       

        double rot_x(double[] phys_X) {
            return rot_x(phys_X[0], phys_X[1]);
        }
        double rot_y(double[] phys_X) {
            return rot_y(phys_X[0], phys_X[1]);
        }

        public double mu_A => 0.1;

        public double mu_B => 1e-20;

        public double Sigma => 0.0;

        public bool TestImmersedBoundary => true;

        public int SpatialDimension => 2;

        public double dt => 1e10;

        public double rho_A => 2;

        public double rho_B => 1e-20;

        public bool Material => true;

        public bool steady => true;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 2;

        public double[] AcceptableL2Error => new[] { 1e-8, 1e-8, 1e-8 };

        public double[] AcceptableResidual => new[] { 1e-8, 1e-8, 1e-8 };

        string Dirichlet = IncompressibleBcType.Velocity_Inlet.ToString() + "_Dirichlet";
        

        public GridCommons CreateGrid(int Resolution) {
            double[] Xnodes = GenericBlas.Linspace(-2, 2, 6 * Resolution + 1);
            double[] Ynodes = GenericBlas.Linspace(-2, 2, 6 * Resolution + 1);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);


            grd.DefineEdgeTags(delegate (double[] X) {
                return Dirichlet;
            });

            return grd;
        }


        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add(Dirichlet, new AppControl.BoundaryValueCollection());
            config[Dirichlet].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", (X, t) => this.Vel(X)[0]);
            config[Dirichlet].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", (X, t) => this.Vel(X)[1]);

            return config;
        }

        public Func<double[], double> GetF(string species, int d) {
            return (X) => 0.0;
        }

        public Func<double[], double, double> GetPhi() {
            // If one does not use the first level-set (i.e. Phi = 0),
            // the code should detect this and assume it to be negative (i.e. set species A everywhere).
            // here, we test this functionality
            return (X, t) => 0.0; 
        }


        double IBMlevSet(double[] X, double t) {
            double y = rot_y(X);
            return -1 + y * y;
        }

        public Func<double[], double, double> GetPhi2() {
            return IBMlevSet;
        }

        public double[] Vel(double[] X) {
            double y = rot_y(X);
            double vyy;
            if(m_movingWalls) {
                vyy = y;
            } else {
                vyy = 1 - y*y;
            }

            double[] V = new[] { Math.Cos(m_angle), Math.Sin(m_angle) };
            V.ScaleV(vyy);
            return V;
        }


        public Func<double[], double, double> GetPhi2U(int d) {
            return (X, t) => Vel(X)[d];
        }

        double Pressure_StaticWalls(double[] X, double t) {
            double x = rot_x(X);
            return -2 * mu_A * x;
        }

        public Func<double[], double, double> GetPress(string species) {
            if(m_movingWalls) {
                return (X, t) => 0.0;
            } else {
                return Pressure_StaticWalls;
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            return (X, t) => Vel(X)[d];
        }
    }
}
