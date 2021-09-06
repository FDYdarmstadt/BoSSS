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

namespace BoSSS.Application.XNSFE_Solver.Tests {

    /// <summary>
    /// a simple Test for HeatTransport in a box shaped domain
    /// </summary>
    class HeatConductivityTest : IXHeatTest {

        double v;
        public HeatConductivityTest(double v = 0.0) {
            this.v = v;
        }

        public int SpatialDimension => 2;

        public double dt => 0.1;

        public double rho_A => 1.0;

        public double rho_B => 1.0;

        public bool Material => true;

        public bool steady => true;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 1;

        public double[] AcceptableL2Error =>  new double[] { 1.0e-7 };

        public double[] AcceptableResidual => new double[] { 1.0e-7 };

        public double c_A => 1.0;

        public double c_B => 1.0;

        public double k_A => 2.0;

        public double k_B => 1.0;

        public double T_sat => Double.MaxValue;

        public double h_vap => 0.0;

        public bool CheckT => true;

        public bool CheckE => false;

        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, 1, Resolution + 1), GenericBlas.Linspace(-1, 1, Resolution + 1), periodicX: false);

            grd.EdgeTagNames.Add(1, "ZeroGradient_top");
            grd.EdgeTagNames.Add(2, "ZeroGradient_bottom");
            grd.EdgeTagNames.Add(3, "ConstantTemperature_left");
            grd.EdgeTagNames.Add(4, "ConstantTemperature_right");

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

                if (Math.Abs(x - (-1)) < 1.0e-6)
                    // left
                    return 3;

                if (Math.Abs(x - (1)) < 1.0e-6)
                    // right
                    return 4;

                throw new ArgumentOutOfRangeException();
                //return 1;
            });
            return grd;
        }

        double T_l = 0.0, T_r = 1.0;
        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {

            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("ZeroGradient_top", new AppControl.BoundaryValueCollection());
            config.Add("ZeroGradient_bottom", new AppControl.BoundaryValueCollection());

            config.Add("ConstantTemperature_left", new AppControl.BoundaryValueCollection());
            config["ConstantTemperature_left"].Evaluators.Add(
                VariableNames.Temperature + "#A",
                (X, t) => T_l);
            config["ConstantTemperature_left"].Evaluators.Add(
                VariableNames.Temperature + "#B",
                (X, t) => T_l);

            config.Add("ConstantTemperature_right", new AppControl.BoundaryValueCollection());
            config["ConstantTemperature_right"].Evaluators.Add(
                VariableNames.Temperature + "#A",
                (X, t) => T_r);
            config["ConstantTemperature_right"].Evaluators.Add(
                VariableNames.Temperature + "#B",
                (X, t) => T_r);

            return config;
        }

        double dx = 0.0;
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double t) {
                double x = X[0];
                double y = X[1];

                return x - dx;
            };
        }

        public Func<double[], double, double> GetT(string species) {
            double k_A = this.k_A;
            double k_B = this.k_B;
            double dx = this.dx;
            double T_l = this.T_l;
            double T_r = this.T_r;
            double q = (T_r - T_l) / ((1 + dx) / k_A + (1 - dx) / k_B);
            switch (species) {
                case "A":
                    return (X, t) => q / k_A * (X[0] + 1.0) ;                    
                case "B":
                    return (X, t) => q / k_B * (X[0] - dx) + q / k_A * (dx + 1.0);                    
                default:
                    throw new ArgumentException();                    
            }           
        }

        public Func<double, double> GetE() => null;
    }

    /// <summary>
    /// Circular preheated domain in a larger quadratic area. No Heatflux over the boundaries.
    /// Stationary temperature should be constant over whole domain.
    /// </summary>
    class HeatDecayTest : IXHeatTest {

        public HeatDecayTest(double R = 0.8, double Q = 0.0) {
            this.r = R;
            this.q = Q;
        }

        public double c_A => 32.0;

        public double c_B => 4.8;

        public double k_A => 10.0;

        public double k_B => 0.04;

        public double T_sat => double.MaxValue;

        public double h_vap => 0.0;

        public int SpatialDimension => 2;

        public double dt => 0.1;

        public double rho_A => 3.4;

        public double rho_B => 10.0;

        public bool Material => true;

        public bool steady => false;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 4;

        public double[] AcceptableL2Error => new double[] { 1.0e-7 };

        public double[] AcceptableResidual => new double[] { 1.0e-7 };

        public bool CheckT => false;

        public bool CheckE => true;

        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 2)
                throw new ArgumentException();

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-1, 1, Resolution + 1), GenericBlas.Linspace(-1, 1, Resolution + 1), periodicX: false);

            grd.EdgeTagNames.Add(1, "ConstantHeatFlux_upper");
            grd.EdgeTagNames.Add(2, "ConstantHeatFlux_right");
            grd.EdgeTagNames.Add(3, "ConstantHeatFlux_lower");
            grd.EdgeTagNames.Add(4, "ConstantHeatFlux_left");

            grd.DefineEdgeTags(delegate (double[] _X) {
                var X = _X;
                double x = X[0];
                double y = X[1];

                if (Math.Abs(y - (+1)) < 1.0e-6)
                    // upper 
                    return 1;

                if (Math.Abs(y - (-1)) < 1.0e-6)
                    // lower 
                    return 3;

                if (Math.Abs(x - (+1)) < 1.0e-6)
                    // right
                    return 2;

                if (Math.Abs(x - (-1)) < 1.0e-6)
                    // left
                    return 4;

                throw new ArgumentOutOfRangeException();
                //return 1;
            });
            return grd;
        }

        double q;
        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("ConstantHeatFlux_upper", new AppControl.BoundaryValueCollection());
            config["ConstantHeatFlux_upper"].Evaluators.Add(
                VariableNames.HeatFluxY + "#B",
                (X, t) => q);

            config.Add("ConstantHeatFlux_right", new AppControl.BoundaryValueCollection());
            config["ConstantHeatFlux_right"].Evaluators.Add(
                VariableNames.HeatFluxX + "#B",
                (X, t) => q);

            config.Add("ConstantHeatFlux_lower", new AppControl.BoundaryValueCollection());
            config["ConstantHeatFlux_lower"].Evaluators.Add(
                VariableNames.HeatFluxY + "#B",
                (X, t) => -q);

            config.Add("ConstantHeatFlux_left", new AppControl.BoundaryValueCollection());
            config["ConstantHeatFlux_left"].Evaluators.Add(
                VariableNames.HeatFluxX + "#B",
                (X, t) => -q);

            return config;
        }

        double r;
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double t) {
                double x = X[0];
                double y = X[1];
                // circular signed distance
                return Math.Sqrt(Math.Pow(x,2) + Math.Pow(y,2)) - r;
            };
        }

        double T0 = 10.0;
        public Func<double[], double, double> GetT(string species) {
            double capacity_a = rho_A * c_A * Math.PI * Math.Pow(r, 2);
            double capacity_b = rho_B * c_B * (4.0 - Math.PI * Math.Pow(r, 2));
            // initial thermal energy
            double e0 = capacity_a * T0;
            // also the t -> inf solution is written here (assuming q = 0)
            switch (species) {
                case "A": { return (X, t) => t == 0.0 ? T0 : e0 / (capacity_a + capacity_b); }
                case "B": { return (X, t) => t == 0.0 ? 0.0 : e0 / (capacity_a + capacity_b); }
                default: { throw new ArgumentException(); }
            }
            
        }

        public Func<double, double> GetE() {
            double capacity_a = rho_A * c_A * Math.PI * Math.Pow(r, 2);
            double capacity_b = rho_B * c_B * (4.0 - Math.PI * Math.Pow(r, 2));
            // initial thermal energy
            double e0 = capacity_a * T0;
            double dE = 8 * q;
            return t => e0 + t * dE;
        }
    }
}
