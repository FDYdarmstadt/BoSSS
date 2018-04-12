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

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// This tests the correctness of the strain jump condition:
    /// \f[ 
    ///   \left\llbracket \mu \left( \vec{u} + \vec{u}^T \right) \vec{n}_{\mathfrak{I}} \right\rrbracket = 0
    /// \f]
    /// Important for this testcase is that here,
    /// the individual terms of the jump condition are nonzero and only the sum of them cancels out, i.e.
    /// \f[ 
    ///   \left\llbracket \mu \vec{u} \cdot \vec{n}_{\mathfrak{I}} \right\rrbracket \neq 0
    ///   \quad \textrm{and} \quad
    ///   \left\llbracket \mu \vec{u}^T \cdot \vec{n}_{\mathfrak{I}} \right\rrbracket \neq 0,
    /// \f]
    /// which is not the case in many other test-cases.
    /// </summary>
    class ViscosityJumpTest : ITest {

        public int SpatialDimension {
            get {
                return 2;
            }
        }

        public Func<double[],double,double> GetPhi() {
            return ((_3D)((time, x, y) => x + y)).Convert_txy2Xt();
        }

        public int LevelsetPolynomialDegree {
            get {
                return 1;
            }
        }

        public Func<double[],double,double> GetU(string species, int d) {
            if (d == 0) {
                return ((_3D)((t, x, y) => y)).Convert_txy2Xt();
            } else if (d == 1) {
                return ((_3D)((t, x, y) => -x)).Convert_txy2Xt();
            } else {
                throw new ArgumentOutOfRangeException();
            }
        }

        public double dt {
            get {
                return 0.0;
            }
        }

        public GridCommons CreateGrid() {

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 5), GenericBlas.Linspace(-2, 2, 5));
            //var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 3), GenericBlas.Linspace(-2, 2, 3));

            grd.EdgeTagNames.Add(1, "Velocity_Inlet");
            grd.DefineEdgeTags(delegate (double[] _X) {
                return 1;
            });

            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("Velocity_Inlet", new AppControl.BoundaryValueCollection());
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => X[1]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#A",
                (X, t) => -X[0]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => X[1]);
            config["Velocity_Inlet"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#B",
                (X, t) => -X[0]);

            return config;
        }

        public Func<double[],double,double> GetPress(string species) {
            return ((_3D)((t, x, y) => 0)).Convert_txy2Xt();
        }

        /// <summary>
        /// the density has no effect in this test (steady-state, no convection terms => density does not appear in the eq.
        /// and should have no effect on the matrix)
        /// </summary>
        public double rho_A {
            get {
                return 1.2;
            }
        }

        /// <summary>
        /// the density has no effect in this test (steady-state, no convection terms => density does not appear in the eq.
        /// and should have no effect on the matrix)
        /// </summary>
        public double rho_B {
            get {
                return 1.0e-3;
            }
        }

        /// <summary>
        /// this test will work for all combinations of viscosities
        /// </summary>
        public double mu_A {
            get {
                return 17.1e-3;
                //return 1.0;
            }
        }

        /// <summary>
        /// this test will work for all combinations of viscosities
        /// </summary>
        public double mu_B {
            get {
                return 1.0;
            }
        }

        public Foundation.ScalarFunction GetS(double time) {
            return ((_2D)((x, y) => 0)).Vectorize();
        }

        public Func<double[],double> GetF(string species, int d) {
            return (X => 0.0);
        }

        public Foundation.ScalarFunction GetSF(double time, int d) {
            return ((_2D)((x, y) => 0)).Vectorize();
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


        public double[] AcceptableL2Error {
            get {
                return new double[] { 5.0e-10, 5.0e-10, 1.0e-7 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return new double[] { 1.0e-7, 1.0e-7, 1.0e-7 };
            }
        }

    }
}
