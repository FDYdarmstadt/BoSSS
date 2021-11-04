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

namespace BoSSS.Application.XNSEC {

    static public partial class FullNSEControlExamples {

        public class PseudoTwoDimensional_TwoPhaseFlow : IXNSECTest {
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

            private int m_SpatialDimension = 2;
            private bool m_DifferentFluids;
            private bool m_ViscosityActive;
            private bool m_RightPressureOutlet;

            public PseudoTwoDimensional_TwoPhaseFlow(bool differentFluids, bool ViscosityActive, bool RightPressureOutlet) {
                m_DifferentFluids = differentFluids;
                m_ViscosityActive = ViscosityActive;
                m_RightPressureOutlet = RightPressureOutlet;
            }

            public Func<double[], double, double> GetPhi() {
                return ((_2D)(delegate (double x, double y) {
                    return (x - 2.3);
                })).Convert_xy2X().Convert_X2Xt();
            }

            public Func<double[], double, double> GetU(string species, int d) {
                if (d == 0) {
                    if (species == "A") {
                        return ((_2D)((x, y) => 4.0)).Convert_xy2X().Convert_X2Xt();
                    } else if (species == "B") {
                        return ((_2D)((x, y) => m_DifferentFluids ? 2.0 : 4.0)).Convert_xy2X().Convert_X2Xt();
                    } else {
                        throw new ArgumentOutOfRangeException();
                    }
                } else if (d == 1) {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }

            public double dt {
                get {
                    return 1.0;
                }
            }

            public GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                double L = 5;
                double H = 2;
                double[] Xnodes = GenericBlas.Linspace(0, L, Resolution + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 4);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicY: true);

                grd.EdgeTagNames.Add(1, "velocity_inlet_left");
                if (m_RightPressureOutlet) {
                    grd.EdgeTagNames.Add(2, "Pressure_outlet");
                } else {
                    grd.EdgeTagNames.Add(2, "velocity_inlet_right");
                }
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 2;
                    return et;
                });

                return grd;
            }

            public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

                config.Add("velocity_inlet_left", new AppControl.BoundaryValueCollection());
                config["velocity_inlet_left"].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", (X, t) => 4.0);
                config["velocity_inlet_left"].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
                config["velocity_inlet_left"].Evaluators.Add(VariableNames.Temperature + "#A", (X, t) => 1.0);
                config["velocity_inlet_left"].Evaluators.Add(VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                if (!m_RightPressureOutlet) {
                    config.Add("velocity_inlet_right", new AppControl.BoundaryValueCollection());
                    config["velocity_inlet_right"].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", (X, t) => m_DifferentFluids ? 2.0 : 4.0);
                    config["velocity_inlet_right"].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
                    config["velocity_inlet_right"].Evaluators.Add(VariableNames.Temperature + "#A", (X, t) => 1.0);
                    config["velocity_inlet_right"].Evaluators.Add(VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                }

                return config;
            }

            public Func<double[], double, double> GetPress(string species) {
                if (species == "A") {
                    return ((_2D)((x, y) => m_DifferentFluids ? -2.0 : 0.0)).Convert_xy2X().Convert_X2Xt();
                } else if (species == "B") {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }

            public double rho_A {
                get { return 0.25; }
            }

            public double rho_B {
                get { return m_DifferentFluids ? 0.5 : 0.25; }
            }

            public double mu_A {
                get {
                    return 10.0;
                }
            }

            public double mu_B {
                get {
                    return 1.0;
                }
            }

            public ScalarFunction GetS(double time) {
                return (m_SpatialDimension == 2) ? ((_2D)((x, y) => 0)).Vectorize() : ((_3D)((x, y, z) => 0)).Vectorize();
            }

            public Func<double[], double> GetF(string species, int d) {
                if (d == 0) {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X();
                } else if (d == 1) {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetTemperature(string species) {
                switch (m_SpatialDimension) {
                    case 2:
                        return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();

                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetMassFractions(string species, int q) {
                switch (m_SpatialDimension) {
                    case 2:
                        if (q == 0) {
                            return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public double Sigma {
                get { return 0.0; }
            }

            public bool Material {
                get { return true; }
            }

            public bool steady {
                get { return true; }
            }

            public bool IncludeConvection {
                get { return true; }
            }

            public int LevelsetPolynomialDegree {
                get { return 4; }
            }

            public double[] AcceptableL2Error {
                get {
                    return (m_SpatialDimension == 2) ?
                      new double[] { 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5 } : new double[] { 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5 };
                }
            }

            public double[] AcceptableResidual {
                get {
                    return (m_SpatialDimension == 2) ?
                      new double[] { 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6 } : new double[] { 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6 };
                }
            }

            public int SpatialDimension {
                get {
                    return m_SpatialDimension;
                }
            }

            public int NumberOfChemicalComponents => 1;

            public bool ChemicalReactionTermsActive => false;

            public double[] GravityDirection => new double[] { 0.0, 0.0, 0.0 };
            public bool EnableMassFractions => false;
            public bool EnableTemperature => false;
        }

        public class PolynomialTestForConvection : IXNSECTest {
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

            public PolynomialTestForConvection(int _SpatialDimension) {
                m_SpatialDimension = _SpatialDimension;
            }

            public Func<double[], double, double> GetPhi() {
                switch (m_SpatialDimension) {
                    case 2:
                        return ((_2D)(delegate (double x, double y) {
                            //return 1;
                            return -y - (2.0 / 5.0) * x + (1.0 / 10.0) * x * x * x;
                        })).Convert_xy2X().Convert_X2Xt();
                    case 3:
                        return ((_3D)(delegate (double x, double y, double z) {
                            return -y - (2.0 / 5.0) * x + (1.0 / 10.0) * x * x * x;
                        })).Convert_xyz2X().Convert_X2Xt();
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetU(string species, int d) {
                switch (m_SpatialDimension) {
                    case 2:
                        if (d == 0) {
                            return ((_2D)((x, y) => x * x)).Convert_xy2X().Convert_X2Xt();
                        } else if (d == 1) {
                            return ((_2D)((x, y) => -2 * x * y)).Convert_xy2X().Convert_X2Xt();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    case 3:
                        if (d == 0) {
                            return ((_3D)((x, y, z) => x * x)).Convert_xyz2X().Convert_X2Xt();
                        } else if (d == 1) {
                            return ((_3D)((x, y, z) => -2 * x * y)).Convert_xyz2X().Convert_X2Xt();
                        } else if (d == 2) {
                            return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X().Convert_X2Xt();
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

            public GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                GridCommons grd;

                switch (m_SpatialDimension) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 5 * Resolution + 1),
                            GenericBlas.Linspace(-2, 2, 5 * Resolution + 1));
                        //grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 2), GenericBlas.Linspace(-2, 2, 2));
                        break;

                    case 3:
                        grd = Grid3D.Cartesian3DGrid(GenericBlas.Linspace(-2, 2, 5 * Resolution + 1),
                          GenericBlas.Linspace(-2, 2, 5 * Resolution + 1), GenericBlas.Linspace(-2, 2, 5 * Resolution + 1));
                        break;

                    default:
                        throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "Velocity_Inlet");
                grd.DefineEdgeTags(delegate (double[] _X) {
                    return 1;
                });

                return grd;
            }

            public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

                config.Add("Velocity_Inlet", new AppControl.BoundaryValueCollection());
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", (X, t) => X[0] * X[0]);
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", (X, t) => -2.0 * X[0] * X[1]);
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.Temperature + "#A", (X, t) => 1.0);
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.MassFraction0 + "#A", (X, t) => 1.0);

                config["Velocity_Inlet"].Evaluators.Add(VariableNames.Velocity_d(0) + "#B", (X, t) => X[0] * X[0]);
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.Velocity_d(1) + "#B", (X, t) => -2.0 * X[0] * X[1]);
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.Temperature + "#B", (X, t) => 1.0);
                config["Velocity_Inlet"].Evaluators.Add(VariableNames.MassFraction0 + "#B", (X, t) => 1.0);
                if (m_SpatialDimension == 3) {
                    config["Velocity_Inlet"].Evaluators.Add(VariableNames.Velocity_d(2) + "#A", (X, t) => 0.0);
                    config["Velocity_Inlet"].Evaluators.Add(VariableNames.Velocity_d(2) + "#B", (X, t) => 0.0);
                }

                return config;
            }

            public Func<double[], double, double> GetPress(string species) {
                return (X, t) => 0.0;
            }

            public double rho_A {
                get { return 100.0; }
            }

            public double rho_B {
                get { return 0.033; }
            }

            public double mu_A {
                get { return 0.0; }
            }

            public double mu_B {
                get { return 0.0; }
            }

            public ScalarFunction GetS(double time) {
                return (m_SpatialDimension == 2) ? ((_2D)((x, y) => 0)).Vectorize() : ((_3D)((x, y, z) => 0)).Vectorize();
            }

            public Func<double[], double> GetF(string species, int d) {
                switch (m_SpatialDimension) {
                    case 2:
                        if (d == 0) {
                            return ((_2D)((x, y) => 2.0 * x * x * x)).Convert_xy2X();
                        } else if (d == 1) {
                            return ((_2D)((x, y) => 2.0 * x * x * y)).Convert_xy2X();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    case 3:
                        if (d == 0) {
                            return ((_3D)((x, y, z) => 2.0 * x * x * x)).Convert_xyz2X();
                        } else if (d == 1) {
                            return ((_3D)((x, y, z) => 2.0 * x * x * y)).Convert_xyz2X();
                        } else if (d == 2) {
                            return ((_3D)((x, y, z) => 0.0)).Convert_xyz2X();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetTemperature(string species) {
                switch (m_SpatialDimension) {
                    case 2:
                        return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();

                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetMassFractions(string species, int q) {
                switch (m_SpatialDimension) {
                    case 2:
                        if (q == 0) {
                            return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public double Sigma {
                get { return 0.0; }
            }

            public bool Material {
                get { return true; }
            }

            public bool steady {
                get { return true; }
            }

            public bool IncludeConvection {
                get { return true; }
            }

            public int LevelsetPolynomialDegree {
                get { return 4; }
            }

            public double[] AcceptableL2Error {
                get {
                    return (m_SpatialDimension == 2) ?
                      new double[] { 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5 } : new double[] { 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5 };
                }
            }

            public double[] AcceptableResidual {
                get {
                    return (m_SpatialDimension == 2) ?
                      new double[] { 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6 } : new double[] { 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6 };
                }
            }

            public int SpatialDimension {
                get {
                    return m_SpatialDimension;
                }
            }

            public int NumberOfChemicalComponents => 1;

            public bool ChemicalReactionTermsActive => false;

            public double[] GravityDirection => new double[] { 0.0, 0.0, 0.0 };

            public bool EnableMassFractions => false;

            public bool EnableTemperature => false;
        }
    }
}