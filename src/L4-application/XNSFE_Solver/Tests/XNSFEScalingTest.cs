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
    class XNSFEScalingTest : IXNSFETest {

        public bool TestImmersedBoundary => false;

        /// <summary>
        /// nix
        /// </summary>
        public Func<double[], double, double> GetPhi2() {
            throw new NotImplementedException();//return (X, t) => -(X[1] + 1.0); // will never be called, as long as 'TestImmersedBoundary' == false;
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            throw new NotImplementedException(); //return (X, t) => 0.0;
        }

        double angle;
        AffineTrafo ROT;

        int setup;
        bool equal;

        /// <summary>
        /// ctor..
        /// </summary>
        /// <param name="angle"></param>
        public XNSFEScalingTest(int Setup, bool EqualFluids) {
            this.angle = 0.0;
            this.ROT = AffineTrafo.Some2DRotation(angle);

            setup = Setup;
            equal = EqualFluids;
        }

        public double mu_A => 1.0;

        public double mu_B => equal ? 1.0 : 0.1;

        public double Sigma => setup == 1 ? 0.0 : 1.0;

        public double c_A => 1.0;

        public double c_B => equal ? 1.0 : 0.1;

        public double c_C => 1.0;

        public double k_A => 1.0;

        public double k_B => equal ? 1.0 : 0.1;

        public double k_C => 1.0;

        public double T_sat => 0.0;//100.0;

        public double h_vap => setup == 3 ? 100.0 : 0.0;

        public bool CheckT => true;
        
        public bool CheckE => false;

        public int SpatialDimension => 2;

        public double dt => 1e-1;

        public double rho_A => 1.0;

        public double rho_B => equal ? 1.0 : 0.1;

        public double rho_C => 1.0;

        public bool Material => false;

        public bool steady => true;

        public bool IncludeConvection => false;

        public int LevelsetPolynomialDegree => 2;

        public double[] AcceptableL2Error => new double[] { 1e-7, 1e-7, 1e-7, 1e-7 }; // dont rly care about these

        public double[] AcceptableResidual => new double[] { 1e-7, 1e-7, 1e-7, 1e-7 };

        public GridCommons CreateGrid(int Resolution) {
            double[] Xnodes = GenericBlas.Linspace(-1, 1, Resolution + 1);
            double[] Ynodes = GenericBlas.Linspace(-1, 1, Resolution + 2);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

            grd.EdgeTagNames.Add(1, "wall_ConstantHeatFlux");
            grd.EdgeTagNames.Add(2, "pressure_outlet_ConstantTemperature");

            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if (Math.Abs(X[1] - Ynodes.First()) <= 1.0e-8)
                    et = 1;
                if (Math.Abs(X[1] - Ynodes.Last()) <= 1.0e-8)
                    et = 2;

                return et;
            });

            return grd.Transform(ROT);
        }


        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {

            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("wall_ConstantHeatFlux", new AppControl.BoundaryValueCollection());
            config.Add("pressure_outlet_ConstantTemperature", new AppControl.BoundaryValueCollection());

            return config;
        }

        public Func<double, double> GetE() => null;

        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }
        public Func<double[], double> GetQ(string species) {
            return (X => 0.0);
        }

        double R => 0.9;

        public double theta_e => throw new NotImplementedException();

        public Func<double[], double, double> GetPhi() {
            return (X, t) => -X[1];
        }

        public Func<double[], double, double> GetPress(string species) {            
            switch (species) {
                case "A": { return (X, t) => 0.0; }
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetT(string species) {
            switch (species) {
                case "A": { return (X, t) => this.T_sat; }
                case "B": { return (X, t) => this.T_sat; }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            switch (species) {
                case "A": { return (X, t) => 0.0; }
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }
    }
}

