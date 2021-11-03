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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;

//using BoSSS.Solution.Utils.Formula;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Application.XNSEC {

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
    class BcTest_PressureOutlet : IXNSECTest {

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

        private int m_SpatialDimension;

        public BcTest_PressureOutlet(int _SpatialDimension) {
            m_SpatialDimension = _SpatialDimension;
        }

        public Func<double[], double, double> GetPhi() {
            switch(m_SpatialDimension) {
                case 2:
                return ((_2D)((x, y) => (x - 0.5))).Convert_xy2X().Convert_X2Xt();
                case 3:
                return ((_3D)((x, y, z) => (x - 0.5))).Convert_xyz2X().Convert_X2Xt();
                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            switch(m_SpatialDimension) {
                case 2:
                if(d == 0) {
                    return ((_2D)((x, y) => x)).Convert_xy2X().Convert_X2Xt();
                } else if(d == 1) {
                    return ((_2D)((x, y) => -y)).Convert_xy2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
                case 3:
                if(d == 0) {
                    return ((_3D)((x, y, z) => x)).Convert_xyz2X().Convert_X2Xt();
                } else if(d == 1) {
                    return ((_3D)((x, y, z) => -y)).Convert_xyz2X().Convert_X2Xt();
                } else if(d == 2) {
                    return ((_3D)((x, y, z) => 0)).Convert_xyz2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        public double dt {
            get {
                return 1.0;
            }
        }

        IncompressibleBcType outcond = IncompressibleBcType.Pressure_Outlet;

        public GridCommons CreateGrid(int Resolution) {
            if(Resolution < 1)
                throw new ArgumentException();

            GridCommons grd;
            switch(m_SpatialDimension) {
                case 2:
                grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(0, 1, 3 * Resolution + 1), GenericBlas.Linspace(0, 1, 3 * Resolution + 1));
                break;
                case 3:
                grd = Grid3D.Cartesian3DGrid(GenericBlas.Linspace(0, 1, 3 * Resolution + 1),
                    GenericBlas.Linspace(0, 1, 3 * Resolution + 1), GenericBlas.Linspace(0, 1, 3 * Resolution + 1));
                break;
                default:
                throw new ArgumentOutOfRangeException();
            }

            grd.EdgeTagNames.Add(1, IncompressibleBcType.Velocity_Inlet.ToString());
            grd.EdgeTagNames.Add(2, outcond.ToString());
            grd.DefineEdgeTags(delegate (double[] _X) {
                double x = _X[0];
                if(Math.Abs(x - 1.0) < 1.0e-6)
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
            if(m_SpatialDimension == 3) {
                config["Velocity_Inlet"].Evaluators.Add(
                    VariableNames.Velocity_d(2) + "#A",
                    (X, t) => 0);
                config["Velocity_Inlet"].Evaluators.Add(
                    VariableNames.Velocity_d(2) + "#B",
                    (X, t) => 0);
            }

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
            switch(m_SpatialDimension) {
                case 2:
                return ((_2D)((x, y) => 1.0 / this.Rey)).Convert_xy2X().Convert_X2Xt();
                case 3:
                return ((_3D)((x, y, z) => 1.0 / this.Rey)).Convert_xyz2X().Convert_X2Xt();
                default:
                throw new ArgumentOutOfRangeException();
            }
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

  

        public Func<double[], double, double> GetMassFractions(string species, int q) {
            switch(m_SpatialDimension) {
                case 2:
                if(q == 0) {
                    return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
                default:
                throw new ArgumentOutOfRangeException();
            }
        }

        public Func<double[], double, double> GetTemperature(string species) {
            switch(m_SpatialDimension) {
                case 2:
                return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();

                default:
                throw new ArgumentOutOfRangeException();
            }
        }


        /// <summary>
        /// zero
        /// </summary>
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
                switch(m_SpatialDimension) {
                    case 2:
                    return new double[] { 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8 };
                    case 3:
                    return new double[] { 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8 };
                    default:
                    throw new ArgumentOutOfRangeException();
                }
            }
        }

        public double[] AcceptableResidual {
            get {
                switch(m_SpatialDimension) {
                    case 2:
                    return new double[] { 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8 };
                    case 3:
                    return new double[] { 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8 };
                    default:
                    throw new ArgumentOutOfRangeException();
                }
            }
        }

        public int SpatialDimension {
            get {
                return m_SpatialDimension;
            }
        }

        public int NumberOfChemicalComponents => 1;

        public bool ChemicalReactionTermsActive => false;

        public double[] GravityDirection => new double[] { 0, 0, 0 };
    }
}
