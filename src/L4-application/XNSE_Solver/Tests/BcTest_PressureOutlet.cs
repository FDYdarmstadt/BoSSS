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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// The main issue of this example is to verify the correctness of the pressure-outlet boundary
    /// condition (<see cref="IncompressibleBcType.Pressure_Outlet"/>):
    /// For this example, the boundary condition
    /// \f[ 
    ///     \frac{1}{\mathrm{Rey}} \nabla \vec{u} \cdot \vec{n}_{\partial \Omega} - p \vec{n}_{\partial \Omega} = 0 \textrm{ on } \Gamma_{\mathrm{POlt}}
    /// \f]
    /// holds at \f$ x=1\f$ , but
    /// \f[ 
    ///     \frac{1}{\mathrm{Rey}} \nabla \vec{u} \cdot \vec{n}_{\partial \Omega} \neq 0, \ \ p \neq 0 \textrm{ on } \Gamma_{\mathrm{POlt}}.
    /// \f]
    /// (In this example the Reynolds number is just density over viscosity.)
    /// Since the terms of the boundary condition do not cancel out individually, but only in sum,
    /// this example is very suitable to verify the correctness of the <see cref="IncompressibleBcType.Pressure_Outlet"/>-implementation.
    /// </summary>
    class BcTest_PressureOutlet : ITest {


        public Func<double[], double, double> GetPhi() {
            return ((_2D)((x, y) => (x - 0.5))).Convert_xy2X().Convert_X2Xt();
        }

        public Func<double[], double, double> GetU(string species, int d) {
            if (d == 0) {
                return ((_2D)((x, y) => x)).Convert_xy2X().Convert_X2Xt();
            } else if (d == 1) {
                return ((_2D)((x, y) => -y)).Convert_xy2X().Convert_X2Xt();
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        public double dt {
            get {
                return 1.0;
            }
        }

        IncompressibleBcType outcond = IncompressibleBcType.Pressure_Outlet;

        public GridCommons CreateGrid() {
            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(0, 1, 4), GenericBlas.Linspace(0, 1, 4));

            grd.EdgeTagNames.Add(1, IncompressibleBcType.Velocity_Inlet.ToString());
            grd.EdgeTagNames.Add(2, outcond.ToString());
            grd.DefineEdgeTags(delegate (double[] _X) {
                double x = _X[0];
                if (Math.Abs(x - 1.0) < 1.0e-6)
                    return 2;

                return 1;
            });

            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("Velocity_Inlet", new AppControl.BoundaryValueCollection());
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => X[0]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#A",
                (X, t) => -X[1]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => X[0]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#B",
                (X, t) => -X[1]);

            config.Add(outcond.ToString(), new AppControl.BoundaryValueCollection());

            return config;
        }

        double Rey {
            get {
                Debug.Assert(this.mu_A == this.mu_B);
                Debug.Assert(this.rho_A == this.rho_B);

                return (this.rho_B / this.mu_A);
            }
        }

        public Func<double[], double, double> GetPress(string species) {
            return ((_2D)((x, y) => 1.0 / this.Rey)).Convert_xy2X().Convert_X2Xt();
        }

        public double rho_A {
            get {
                return 1.0;
            }
        }

        public double rho_B {
            get {
                return 1.0;
            }
        }

        public double mu_A {
            get {
                return 0.1;
            }
        }

        public double mu_B {
            get {
                return 0.1;
            }
        }
        
        public Func<double[], double> GetF(string species, int d) {
            return X => 0.0;
        }
        
        public double Sigma {
            get {
                return 0.0;
            }
        }

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
                return false;
            }
        }

        public int LevelsetPolynomialDegree {
            get {
                return 1;
            }
        }

        public double[] AcceptableL2Error {
            get {
                return new double[] { 1.0e-8, 1.0e-8, 1.0e-8 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return new double[] { 1.0e-8, 1.0e-8, 1.0e-8 };
            }
        }

        public int SpatialDimension {
            get {
                return 2;
            }
        }

        
    }
}
