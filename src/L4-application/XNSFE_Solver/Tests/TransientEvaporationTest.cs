using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSFE_Solver.Tests {
    class TransientEvaporationTest : IXNSFETest {

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

        double angle;
        AffineTrafo ROT;

        /// <summary>
        /// ctor..
        /// </summary>
        /// <param name="angle"></param>
        public TransientEvaporationTest(double angle = 0.0) {
            this.angle = angle;
            this.ROT = AffineTrafo.Some2DRotation(angle);
        }

        public double mu_A => 1.0;

        public double mu_B => 1.0;

        public double Sigma => 0.0;

        public double c_A => 1.0;

        // has to be zero, because we do not want to include a convective contribution to heatflux and temperature profile
        public double c_B => 0.0;

        public double k_A => 1.0;

        public double k_B => 0.1;

        public double T_sat => 100.0;

        public double h_vap => 100.0;

        public bool CheckT => true;

        public bool CheckE => false;

        public int SpatialDimension => 2;

        public double dt => 1e-2;

        public double rho_A => 1.0;

        public double rho_B => 0.1;

        public bool Material => false;

        public bool steady => false;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 2;

        public double[] AcceptableL2Error => new double[] { 1e-7, 1e-7, 1e-7, 1e-7 };

        public double[] AcceptableResidual => new double[] { 1e-7, 1e-7, 1e-7, 1e-7 };

        public double theta_e => Math.PI / 6;

        double L = 1.0;
        public GridCommons CreateGrid(int Resolution) {
            double[] Xnodes = GenericBlas.Linspace(0, L, Resolution + 1);
            double[] Ynodes = GenericBlas.Linspace(0, L, Resolution + 1);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

            if (!steady) {
                grd.EdgeTagNames.Add(2, "wall_ConstantHeatFlux_upper");
                grd.EdgeTagNames.Add(1, "pressure_Dirichlet_ZeroGradient_lower");
            } else {
                grd.EdgeTagNames.Add(2, "pressure_dirichlet_ConstantTemperature_upper");
                grd.EdgeTagNames.Add(1, "velocity_inlet_ZeroGradient_lower");
            }
            grd.EdgeTagNames.Add(3, "slipsymmetry_ZeroGradient");


            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if (Math.Abs(X[1]) <= 1.0e-8)
                    et = 1;
                if (Math.Abs(X[1] - L) <= 1.0e-8)
                    et = 2;
                if (Math.Abs(X[0]) <= 1.0e-8 || Math.Abs(X[0] - L) <= 1.0e-8)
                    et = 3;

                return et;
            });

            return grd.Transform(ROT);
        }


        double zi0 = 0.1;
        double qv = 10.0;
        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            double qv = 10.0;
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("slipsymmetry_ZeroGradient", new AppControl.BoundaryValueCollection());

            config.Add("wall_ConstantHeatFlux_upper", new AppControl.BoundaryValueCollection());
            config.Add("pressure_Dirichlet_ZeroGradient_lower", new AppControl.BoundaryValueCollection());

            config["wall_ConstantHeatFlux_upper"].Evaluators.Add("HeatFluxX#B", (X, t) => qv * Math.Sin(angle));
            config["wall_ConstantHeatFlux_upper"].Evaluators.Add("HeatFluxY#B", (X, t) => qv * Math.Cos(angle));
            config["pressure_Dirichlet_ZeroGradient_lower"].Evaluators.Add("Pressure#A", (X, t) => 0.0);


            return config;
        }

        public Func<double, double> GetE() => null;

        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }

        public Func<double[], double> GetQ(string species) {
            return (X => 0.0);
        }

        public double xSemiAxis0 => 10 * L;

        public double ySemiAxis0 => 1e-6;

        public double yCenter0 => L - zi0;
        public Func<double[], double, double> GetPhi() {
            return (X, t) =>  X[1] - yCenter0 + qv / (rho_B * h_vap) * t;
        }

        //public Func<double[], double, double> GetPhi() {
        //    return (X, t) => (zi0 + qv / (rho_B * h_vap) * t) - Math.Cos(angle) * X[1] + Math.Sin(angle) * X[0];
        //}

        public Func<double[], double, double> GetPress(string species) {
            double dp = -(qv / h_vap).Pow2() * (1 / this.rho_A - 1 / this.rho_B);
            switch (species) {
                case "A": { return (X, t) => 0.0; }
                case "B": { return (X, t) => 0.0 - dp; }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetT(string species) {
            switch (species) {
                case "A": { return (X, t) => this.T_sat; }
                case "B": { return (X, t) => this.T_sat + (qv / this.k_B) * ( X[1] - L +  (zi0 + qv / (rho_B * h_vap) * t)); }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            switch (species) {
                case "A": { return (X, t) => d == 0 ? -0.9 * Math.Sin(angle) : d == 1 ? -0.9 * Math.Cos(angle) : throw new ArgumentException(); }
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }
    }
}

