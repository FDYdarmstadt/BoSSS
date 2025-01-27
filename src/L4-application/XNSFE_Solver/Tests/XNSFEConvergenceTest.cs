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
    class XNSFEConvergenceTest : IXNSFETest {
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
        AffineTrafo ROTinv;

        /// <summary>
        /// ctor..
        /// </summary>
        /// <param name="angle"></param>
        public XNSFEConvergenceTest() {
            // tilt, so that both velocity components are non zero
            this.angle = Math.PI / 6;
            this.ROT = AffineTrafo.Some2DRotation(angle);
            this.ROTinv = this.ROT.Invert();
        }

        public double mu_A => 1.0;

        public double mu_B => 1.0;

        public double Sigma => 1.0;

        public double c_A => 1.0;

        public double c_B => 1.0;

        public double k_A => 1.0;

        public double k_B => 1.0 ;

        public double T_sat => 0.0;//100.0;

        public double h_vap => 1.0;

        public bool CheckT => true;

        public bool CheckE => false;

        public int SpatialDimension => 2;

        public double dt => 1e-1;

        public double rho_A => 1.0;

        public double rho_B => 0.1;

        public bool Material => false;

        public bool steady => true;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 2;

        public double[] AcceptableL2Error => new double[] { 1e-7, 1e-7, 1e-7, 1e-7 };

        public double[] AcceptableResidual => new double[] { 1e-7, 1e-7, 1e-7, 1e-7 };

        public GridCommons CreateGrid(int Resolution) {
            double[] Xnodes = GenericBlas.Linspace(-1, 1, Resolution + 1);
            double[] Ynodes = GenericBlas.Linspace(-1, 1, Resolution + 1);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicY: true);

            grd.EdgeTagNames.Add(1, "velocity_inlet_ConstantTemperature");
            grd.EdgeTagNames.Add(2, "pressure_outlet_ConstantTemperature");

            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if (Math.Abs(X[0] - Xnodes.First()) <= 1.0e-8)
                    et = 1;
                if (Math.Abs(X[0] - Xnodes.Last()) <= 1.0e-8)
                    et = 2;

                return et;
            });

            return grd.Transform(ROT);
        }

        double U_In => 1.0;
        double Massflux => this.rho_A * this.U_In;
        int n => 1;
        double p0 => 1.0;
        double T0 => 1.0;
        double alpha => this.Massflux * this.h_vap / (this.k_B * this.T0);

        public double theta_e => throw new NotImplementedException();

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {

            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("velocity_inlet_ConstantTemperature", new AppControl.BoundaryValueCollection());
            config.Add("pressure_outlet_ConstantTemperature", new AppControl.BoundaryValueCollection());

            config["velocity_inlet_ConstantTemperature"].Evaluators.Add("VelocityX#A", (X, t) => this.U_In * Math.Cos(this.angle));
            config["velocity_inlet_ConstantTemperature"].Evaluators.Add("VelocityY#A", (X, t) => this.U_In * Math.Sin(this.angle));
            config["velocity_inlet_ConstantTemperature"].Evaluators.Add("Temperature#A", (X, t) => this.T_sat);
            config["pressure_outlet_ConstantTemperature"].Evaluators.Add("Pressure#B", (X, t) => 0.0);
            config["pressure_outlet_ConstantTemperature"].Evaluators.Add("Temperature#B", (X,t) => this.T_sat + this.T0 * (1.0-Math.Exp(-alpha)));
            return config;
        }

        public Func<double, double> GetE() => null;

        public Func<double[], double> GetF(string species, int d) {
            switch (species) {
                case "A": { return X => ROT.Transform(new double[] { 1, 0 })[d] * this.p0 / this.rho_A * this.n * 2 * Math.PI * Math.Cos(this.n * 2 * Math.PI * ROTinv.Transform(X)[0]); }
                case "B": { return X => ROT.Transform(new double[] { 1, 0 })[d] * this.p0 / this.rho_B * this.n * 2 * Math.PI * Math.Cos(this.n * 2 * Math.PI * ROTinv.Transform(X)[0]); }
                default: { throw new ArgumentException(); }
            }
        }
        public Func<double[], double> GetQ(string species) {
            switch (species) {
                case "A": { return X => this.rho_A * this.c_A * this.U_In * this.n * 2 * Math.PI * Math.Sin(this.n * 2 * Math.PI * ROTinv.Transform(X)[0]) - this.k_A * (this.n * 2 * Math.PI).Pow2() * Math.Cos(this.n * 2 * Math.PI * ROTinv.Transform(X)[0]); }
                case "B": { return X => (this.rho_A * this.c_B * this.U_In * this.alpha + this.alpha.Pow2()) * this.T0 * Math.Exp(-this.alpha * ROTinv.Transform(X)[0]); }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetPhi() {
            return (X, t) => ROTinv.Transform(X)[0];
        }

        public Func<double[], double, double> GetPress(string species) {
            switch (species) {
                case "A": { return (X, t) => -this.Massflux.Pow2() * (1.0/this.rho_A - 1.0/this.rho_B) + this.p0 * Math.Sin(this.n * 2 * Math.PI * ROTinv.Transform(X)[0]); }
                case "B": { return (X, t) => this.p0 * Math.Sin(this.n * 2 * Math.PI * ROTinv.Transform(X)[0]); }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetT(string species) {
            switch (species) {
                case "A": { return (X, t) => this.T_sat + this.T0 * (1.0 - Math.Cos(n * 2 * Math.PI * ROTinv.Transform(X)[0])); }
                case "B": { return (X, t) => this.T_sat + this.T0 * (1.0 - Math.Exp(-alpha * ROTinv.Transform(X)[0])); }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            switch (species) {
                case "A": { return (X, t) => ROT.Transform(new double[] { this.U_In, 0 })[d]; }
                case "B": { return (X, t) => ROT.Transform(new double[] { this.rho_A / this.rho_B * this.U_In, 0 })[d]; }
                default: { throw new ArgumentException(); }
            }
        }
    }

}
