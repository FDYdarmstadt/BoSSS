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
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Application.XNSFE_Solver.Tests {

    interface ISlipTest {
        double slipI { get; }
    }

    /// <summary>
    /// a periodic shear flow, with two phases.
    /// Interfacial slip (noslip, slip, freeslip) with equal and unequal viscosities can be tested
    /// With a nonmaterial interface there will be a jump in tangential stresses.
    /// This means that for larger massfluxes the sign of the velocity gradient can switch!
    /// Convective terms are not regarded.
    /// More details in Annual Report 2023 "Rieckmann" (add arxiv when available), Case 2
    /// </summary>
    class InterfaceSlipTestLin : IXNSFETest, ISlipTest {

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

        public bool Material {
            get {
                return false;
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

        /// <summary>
        /// the zero-level-set is identical to the x-axis
        /// lower Phase is "A"
        /// </summary>
        /// <param name="time"></param>
        /// <returns></returns>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double t) {
                var Coord = ROTinv.Transform(X[0], X[1]);
                double x = Coord[0];
                double y = Coord[1];

                return y;
            };
        }


        public int LevelsetPolynomialDegree {
            get {
                return 1;
            }
        }

        public InterfaceSlipTestLin(double angle, double slipI, double viscosityratio, double massflux) {
            //double angle = 0.0;
            //double angle = 60.0 * Math.PI / 180.0;
            ROT = AffineTrafo.Some2DRotation(angle);
            ROTinv = ROT.Invert();

            this.mu_B = viscosityratio;
            this.slipI = slipI;
            this.massflux = massflux;
        }

        AffineTrafo ROT;

        AffineTrafo ROTinv;

        bool periodic = true;

        public Func<double[], double> GetF(string species, int d) {
            return X => 0.0;
        }

        private double UpperShear { get{ return 0.5 * Math.Exp(-slipI); } }

        double[] Coeffs(string spc) {
            double[] c = new double[2];
            if(spc == "A") {
                if(slipI == double.PositiveInfinity) {
                    c[0] = 0.0;
                    c[1] = 0.0;
                } else {
                    c[0] = (UpperShear * mu_B * (mu_A + mu_B - massflux * slipI)) / (mu_A * mu_B + mu_A * mu_A + massflux * mu_A * slipI);
                    c[1] = (UpperShear * mu_B * (mu_A + mu_B - massflux * slipI)) / (mu_A * mu_B + mu_A * mu_A + massflux * mu_A * slipI);
                }
            } else if (spc == "B") {               
                if (slipI == double.PositiveInfinity) {
                    c[0] = 0.0;
                    c[1] = 0.0; // owing to the boundary conditions available we cann not test a shear flow with a permeable no-slip wall on top
                } else {
                    c[0] = UpperShear;
                    c[1] = (UpperShear * mu_B * (mu_A + mu_B - massflux * slipI + 2 * mu_A * slipI)) / (mu_A * mu_B + mu_A * mu_A + massflux * mu_A * slipI);
                }
            } else {
                throw new ArgumentException();
            }
            return c;
        }

        public Func<double[], double, double> GetU(string species, int d) {
            double[] c = Coeffs(species);

            return ((_2D)(delegate (double _x, double _y) {

                var Coord = ROTinv.Transform(_x, _y);
                double y = Coord[1];

                Debug.Assert(Coord[0] >= -2);
                Debug.Assert(Coord[0] <= +2);
                Debug.Assert(Coord[1] >= -1);
                Debug.Assert(Coord[1] <= +1);

                double u = (c[0] * y + c[1]);
                double rho = 0.0;
                if(species == "A") {
                    rho = rho_A;
                }else if(species == "B") {
                    rho = rho_B;
                } else {
                    throw new ArgumentException();
                }
                var UT = ROT.Transform(u, massflux/rho);

                return UT[d];

            })).Convert_xy2X().Convert_X2Xt();
        }

        public double dt {
            get {
                return 1.0;
            }
        }

        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();

            //var yNodes = GenericBlas.Linspace(-1, 1, 9);

            //var yNodes1 = yNodes.GetSubVector(0, 4);
            //var yNodes2 = yNodes.GetSubVector(5, 4);
            //var _yNodes = ArrayTools.Cat(yNodes1, yNodes2);
            //var _yNodes = yNodes;

            var __yNodes = new double[] { -1, -0.8, -0.6, 0.6, 1.0 };
            var _yNodes = new double[(__yNodes.Length - 1) * Resolution + 1];
            for(int i = 0; i < __yNodes.Length -1; i++) {
                double[] part = GenericBlas.Linspace(__yNodes[i], __yNodes[i + 1], Resolution + 1);
                _yNodes.SetSubVector(part, i * Resolution, part.Length);
            }

            //var _yNodes = GenericBlas.Linspace(-1, 1, 6);

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 3*Resolution + 1), _yNodes, periodicX: periodic);
            if (periodic) {
                grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");
                grd.EdgeTagNames.Add(2, "velocity_inlet_ConstantHeatFlux");

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


                    throw new ArgumentOutOfRangeException();
                    //return 1;
                });

                Console.WriteLine("Interface slip, shearflow, nonmaterial, periodic.");

            } else {
                throw new NotImplementedException();
            }
            var grdT = grd.Transform(ROT);

            return grdT;
        }


        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("pressure_outlet_ConstantTemperature", new AppControl.BoundaryValueCollection());
            config["pressure_outlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.VelocityX_GradientY + "#B",
                    (X, t) => UpperShear);
            config["pressure_outlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.Temperature + "#B",
                    (X, t) => 0.0);
            config.Add("velocity_inlet_ConstantHeatFlux", new AppControl.BoundaryValueCollection());
            config["velocity_inlet_ConstantHeatFlux"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#A",
                    (X, t) => 0.0);
            config["velocity_inlet_ConstantHeatFlux"].Evaluators.Add(
                    VariableNames.Velocity_d(1) + "#A",
                    (X, t) => massflux/rho_A);
            config["velocity_inlet_ConstantHeatFlux"].Evaluators.Add(
                    VariableNames.HeatFluxVectorComponent(1) + "#A",
                    (X, t) => -massflux * h_vap);

            if (!periodic) {
                throw new NotImplementedException();               
            }

            return config;
        }

        public Func<double[], double, double> GetPress(string species) {
            return (X, t) => species == "A" ? -massflux*(massflux/rho_A - massflux/rho_B) : 0.0;
        }

        public Func<double[], double> GetQ(string species) {
            return X => 0.0;
        }

        public Func<double[], double, double> GetT(string species) {
            return (X, t) => species == "A" ? -(h_vap*massflux)/k_A * X[1] : 0.0;
        }

        public Func<double, double> GetE() => null;


        /// <summary>
        /// interfacial slip length
        /// </summary>
        public double slipI {
            private set;
            get;
        }

        /// <summary>
        /// massflux
        /// </summary>
        public double massflux {
            private set;
            get;
        }

        /// <summary>
        /// specific weight, air
        /// </summary>
        public double rho_B {
            get {
                return 1.2;
            }
        }

        /// <summary>
        /// specific weight, water
        /// </summary>
        public double rho_A {
            get {
                return 1000;
            }
        }


        /// <summary>
        /// dynamic viscosity, air
        /// </summary>
        public double mu_B {
            private set;
            get;
        }

        /// <summary>
        /// dynamic viscosity, water
        /// </summary>
        public double mu_A {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// surface tension of water (surface tension has no effect due to the planar interface)
        /// </summary>
        public double Sigma {
            get {
                //return 0.0;
                return 72.75e-3;
            }
        }
                       

        public double[] AcceptableL2Error {
            get {
                return new double[] { 5.0e-6, 1.0e-7, 1.0e-7, 1.0e-7 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return new double[] { 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10 };
            }
        }

        public int SpatialDimension {
            get {
                return 2;
            }
        }

        public double c_A => 2.5; // rescaled, otherwise the values "explode" in this testcase

        public double c_B => 1.85;

        public double k_A => 0.6;

        public double k_B => 0.02;

        public double T_sat => 0.0;

        public double h_vap => 12000;

        public bool CheckT => true;

        public bool CheckE => false;
    }

    /// <summary>
    /// a periodic shear flow, with two phases.
    /// Interfacial slip (noslip, slip, freeslip) with equal and unequal viscosities can be tested
    /// With a nonmaterial interface there will be a jump in tangential stresses.
    /// This means that for larger massfluxes the sign of the velocity gradient can switch!
    /// Convective terms are regarded.
    /// More details in Annual Report 2023 "Rieckmann" (add arxiv when available), Case 3
    /// </summary>
    class InterfaceSlipTestNonLin : IXNSFETest, ISlipTest {

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

        public bool Material {
            get {
                return false;
            }
        }

        public bool steady {
            get {
                return true;
            }
        }

        public bool IncludeConvection {
            get {
                return true;
            }
        }

        /// <summary>
        /// the zero-level-set is identical to the x-axis
        /// lower Phase is "A"
        /// </summary>
        /// <param name="time"></param>
        /// <returns></returns>
        public Func<double[], double, double> GetPhi() {
            return delegate (double[] X, double t) {
                var Coord = ROTinv.Transform(X[0], X[1]);
                double x = Coord[0];
                double y = Coord[1];

                return y;
            };
        }


        public int LevelsetPolynomialDegree {
            get {
                return 1;
            }
        }

        public InterfaceSlipTestNonLin(double angle, double slipI, double viscosityratio, double massflux) {
            //double angle = 0.0;
            //double angle = 60.0 * Math.PI / 180.0;
            ROT = AffineTrafo.Some2DRotation(angle);
            ROTinv = ROT.Invert();

            this.mu_B = viscosityratio;
            this.slipI = slipI;
            this.massflux = massflux;
            if(massflux == 0.0) {
                throw new ArgumentException("Testcase only defined for massflux unequal zero");
            }
        }

        AffineTrafo ROT;

        AffineTrafo ROTinv;

        bool periodic = true;

        public Func<double[], double> GetF(string species, int d) {
            return X => 0.0;
        }

        private double UpperShear { get { return 0.5 * Math.Exp(-slipI); } }

        double[] Coeffs(string spc) {
            double[] c = new double[2];
            if (spc == "A") {
                if (slipI == double.PositiveInfinity) {
                    c[0] = 0.0;
                    c[1] = 0.0;
                } else {
                    c[0] = (UpperShear * mu_B * Math.Exp(-massflux / mu_B) * (mu_A + mu_B - massflux * slipI)) / (massflux * (mu_A + mu_B + massflux * slipI));
                    c[1] = -(UpperShear * mu_B * Math.Exp(-massflux / mu_A) * Math.Exp(-massflux / mu_B) * (mu_A + mu_B - massflux * slipI)) / (massflux * (mu_A + mu_B + massflux * slipI));
                }
            } else if (spc == "B") {
                if (slipI == double.PositiveInfinity) {
                    c[0] = 0.0;
                    c[1] = 0.0; // owing to the boundary conditions available we cann not test a shear flow with a permeable no-slip wall on top
                } else {
                    c[0] = (UpperShear * mu_B * Math.Exp(-massflux / mu_B)) / massflux;
                    c[1] = -(UpperShear * mu_B * Math.Exp(-massflux / mu_A) * Math.Exp(-massflux / mu_B) * (mu_A + mu_B - massflux * slipI)) / (massflux * (mu_A + mu_B + massflux * slipI));
                }
            } else {
                throw new ArgumentException();
            }
            return c;
        }

        public Func<double[], double, double> GetU(string species, int d) {
            double[] c = Coeffs(species);

            return ((_2D)(delegate (double _x, double _y) {

                var Coord = ROTinv.Transform(_x, _y);
                double y = Coord[1];

                double rho = 0.0;
                double mu = 0.0;
                if (species == "A") {
                    rho = rho_A;
                    mu = mu_A;
                } else if (species == "B") {
                    rho = rho_B;
                    mu = mu_B;
                } else {
                    throw new ArgumentException();
                }
                double V0 = massflux / rho;

                Debug.Assert(Coord[0] >= -2);
                Debug.Assert(Coord[0] <= +2);
                Debug.Assert(Coord[1] >= -1);
                Debug.Assert(Coord[1] <= +1);

                double u = (c[0] * Math.Exp(rho*V0/mu * y) + c[1]);
                
                var UT = ROT.Transform(u, V0);

                return UT[d];

            })).Convert_xy2X().Convert_X2Xt();
        }

        public double dt {
            get {
                return 1.0;
            }
        }

        public GridCommons CreateGrid(int Resolution) {
            if (Resolution < 1)
                throw new ArgumentException();

            //var yNodes = GenericBlas.Linspace(-1, 1, 9);

            //var yNodes1 = yNodes.GetSubVector(0, 4);
            //var yNodes2 = yNodes.GetSubVector(5, 4);
            //var _yNodes = ArrayTools.Cat(yNodes1, yNodes2);
            //var _yNodes = yNodes;

            var __yNodes = GenericBlas.Linspace(-1.0, 1.0, 10);
            var _yNodes = new double[(__yNodes.Length - 1) * Resolution + 1];
            for (int i = 0; i < __yNodes.Length - 1; i++) {
                double[] part = GenericBlas.Linspace(__yNodes[i], __yNodes[i + 1], Resolution + 1);
                _yNodes.SetSubVector(part, i * Resolution, part.Length);
            }

            //var _yNodes = GenericBlas.Linspace(-1, 1, 6);

            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2, 3 * Resolution + 1), _yNodes, periodicX: periodic);
            if (periodic) {
                grd.EdgeTagNames.Add(1, "pressure_outlet_ConstantTemperature");
                grd.EdgeTagNames.Add(2, "velocity_inlet_ConstantTemperature");

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


                    throw new ArgumentOutOfRangeException();
                    //return 1;
                });

                Console.WriteLine("Interface slip, shearflow, nonmaterial, periodic.");

            } else {
                throw new NotImplementedException();
            }
            var grdT = grd.Transform(ROT);

            return grdT;
        }


        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
            config.Add("pressure_outlet_ConstantTemperature", new AppControl.BoundaryValueCollection());
            config["pressure_outlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.VelocityX_GradientY + "#B",
                    (X, t) => UpperShear);
            config["pressure_outlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.Temperature + "#B",
                    (X, t) => 0.0);
            config.Add("velocity_inlet_ConstantTemperature", new AppControl.BoundaryValueCollection());
            config["velocity_inlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.Velocity_d(0) + "#A",
                    (X, t) => 0.0);
            config["velocity_inlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.Velocity_d(1) + "#A",
                    (X, t) => massflux / rho_A);
            var TA = GetT("A");
            config["velocity_inlet_ConstantTemperature"].Evaluators.Add(
                    VariableNames.Temperature + "#A",
                    (X, t) => TA(new double[] { 0.0, -1.0 }, t));

            if (!periodic) {
                throw new NotImplementedException();
            }

            return config;
        }

        public Func<double[], double, double> GetPress(string species) {
            return (X, t) => species == "A" ? -massflux * (massflux / rho_A - massflux / rho_B) : 0.0;
        }

        public Func<double[], double> GetQ(string species) {
            return X => 0.0;
        }

        public Func<double[], double, double> GetT(string species) {
            double c0 = -(h_vap) / c_A;
            double c1 = (h_vap) / c_A;
            double V0 = massflux / rho_A;
            return (X, t) => species == "A" ? c0 * Math.Exp(rho_A * c_A * V0 /k_A * X[1]) + c1 : 0.0;
        }

        public Func<double, double> GetE() => null;


        /// <summary>
        /// interfacial slip length
        /// </summary>
        public double slipI {
            private set;
            get;
        }

        /// <summary>
        /// massflux
        /// </summary>
        public double massflux {
            private set;
            get;
        }

        /// <summary>
        /// specific weight, air
        /// </summary>
        public double rho_B {
            get {
                return 1.2;
            }
        }

        /// <summary>
        /// specific weight, water
        /// </summary>
        public double rho_A {
            get {
                return 1000;
            }
        }


        /// <summary>
        /// dynamic viscosity, air
        /// </summary>
        public double mu_B {
            private set;
            get;
        }

        /// <summary>
        /// dynamic viscosity, water
        /// </summary>
        public double mu_A {
            get {
                return 1.0;
            }
        }

        /// <summary>
        /// surface tension of water (surface tension has no effect due to the planar interface)
        /// </summary>
        public double Sigma {
            get {
                //return 0.0;
                return 72.75e-3;
            }
        }


        public double[] AcceptableL2Error {
            get {
                return new double[] { 5.0e-5, 1.0e-3, 1.0e-3, 1.0e-1 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return new double[] { 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10 };
            }
        }

        public int SpatialDimension {
            get {
                return 2;
            }
        }

        public double c_A => 2.5; // rescaled, otherwise the values "explode" in this testcase

        public double c_B => 1.85;

        public double k_A => 0.6;

        public double k_B => 0.02;

        public double T_sat => 0.0;

        public double h_vap => 12000;

        public bool CheckT => true;

        public bool CheckE => false;
    }


}
