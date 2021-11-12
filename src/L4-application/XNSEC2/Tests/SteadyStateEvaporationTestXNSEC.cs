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

using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.XNSEC {

    internal class SteadyStateEvaporationTestXNSEC : IXNSECTest_Heat {
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

        private double angle;
        private AffineTrafo ROT;

        /// <summary>
        /// ctor..
        /// </summary>
        /// <param name="angle"></param>
        public SteadyStateEvaporationTestXNSEC(double angle = 0.0) {
            this.angle = angle;
            this.ROT = AffineTrafo.Some2DRotation(angle);
        }

        public double mu_A => 1.0;

        public double mu_B => 1.0;

        public double Sigma => 1.0 * 0;

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

        public double dt => 5e-4 * 1e100;

        public double rho_A => 1.0;

        public double rho_B => 0.1;

        public bool Material => false;

        public bool steady => true;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 2;

        public double[] AcceptableL2Error => new double[] { 1e-7, 1e-7, 1e-7, 1e-7/*, 1e-7, 1e-7*/ };

        public double[] AcceptableResidual => new double[] { 1e-7, 1e-7, 1e-7, 1e-7/*, 1e-7, 1e-7 */};

        public int NumberOfChemicalComponents => /*2*/1;

        public bool ChemicalReactionTermsActive => false;

        public bool EnableMassFractions => false;

        public bool EnableTemperature => true;

        public double[] GravityDirection => new double[] { 0, 0, 0 };

        private double L = 1.0;

        public GridCommons CreateGrid(int Resolution) {
            double[] Xnodes = GenericBlas.Linspace(0, L, Resolution + 1);
            double[] Ynodes = GenericBlas.Linspace(0, L, Resolution + 1);
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

            grd.EdgeTagNames.Add(1, "ScalarDirichlet_PressureOutlet");
            grd.EdgeTagNames.Add(2, "velocity_inlet_upper");



            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if (Math.Abs(X[1]) <= 1.0e-8)
                    et = 1;
                if (Math.Abs(X[1] - L) <= 1.0e-8)
                    et = 2;  

                return et;
            });

            return grd.Transform(ROT);
        }

        private double zi0 = 0.2;
        private double qv = 10.0;

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            double qv = 10.0;
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            //config.Add("slipsymmetry_ZeroGradient", new AppControl.BoundaryValueCollection());

            config.Add("ScalarDirichlet_PressureOutlet", new AppControl.BoundaryValueCollection());
            config.Add("velocity_inlet_upper", new AppControl.BoundaryValueCollection());

            //config["ScalarDirichlet_PressureOutlet"].Evaluators.Add("Temperature#B", (X, t) => this.T_sat + 10 * (qv / this.k_B) * (zi0 - Math.Cos(angle) * X[1] + Math.Sin(angle) * X[0]));
            config["ScalarDirichlet_PressureOutlet"].Evaluators.Add("Temperature", (X, t) => this.T_sat + (qv / this.k_B) * (zi0 - Math.Cos(angle) * X[1] + Math.Sin(angle) * X[0]));

            config["velocity_inlet_upper"].Evaluators.Add("VelocityX#A", (X, t) => 0.1 * Math.Sin(angle));
            config["velocity_inlet_upper"].Evaluators.Add("VelocityY#A", (X, t) => -0.1 * Math.Cos(angle));
            config["velocity_inlet_upper"].Evaluators.Add("Temperature#A", (X, t) => this.T_sat);

            if (EnableMassFractions) { 
            config["ScalarDirichlet_PressureOutlet"].Evaluators.Add("MassFraction0", (X, t) => 1.0);
            config["velocity_inlet_upper"].Evaluators.Add("MassFraction1", (X, t) => 1.0);
            }
            return config;
        }

        public Func<double, double> GetE() => null;

        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }

        public Func<double[], double> GetQ(string species) {
            return (X => 0.0);
        }

        public Func<double[], double, double> GetPhi() {
            return (X, t) => zi0 - Math.Cos(angle) * X[1] + Math.Sin(angle) * X[0];
        }

        public Func<double[], double, double> GetPress(string species) {
            double dp = -(qv / h_vap).Pow2() * (1 / this.rho_A - 1 / this.rho_B);
            switch (species) {
                case "A": { return (X, t) => 0.0 + dp; }
                case "B": { return (X, t) => 0.0; }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetT(string species) {
            throw new ArgumentException();
        }

        public Func<double[], double, double> GetTemperature(string species) {
            switch (species) {
                case "A": { return (X, t) => this.T_sat; }
                case "B": { return (X, t) => this.T_sat + (qv / this.k_B) * (zi0 - Math.Cos(angle) * X[1] + Math.Sin(angle) * X[0]); }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            switch (species) {
                case "A": { return (X, t) => d == 0 ? 0.1 * Math.Sin(angle) : d == 1 ? -0.1 * Math.Cos(angle) : throw new ArgumentException(); }
                case "B": { return (X, t) => d == 0 ? 1.0 * Math.Sin(angle) : d == 1 ? -1.0 * Math.Cos(angle) : throw new ArgumentException(); }
                default: { throw new ArgumentException(); }
            }
        }

        public Func<double[], double, double> GetMassFractions(string species, int q) {
            //if (q == 0) {
            //    return ((_3D)((t, x, y) => species == "A" ? 0.0 : 1.0)).Convert_txy2Xt();
            //} else if (q == 1) {
            //    return ((_3D)((t, x, y) => species == "B" ? 1.0 : 0.0)).Convert_txy2Xt();
            //} else {
            //    throw new ArgumentOutOfRangeException();
            //}
            return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();
        }
    }
}