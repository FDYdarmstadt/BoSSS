using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSE_Solver.Tests {
    class ParameterizedLevelSetTest : IXNSETest {

        private double m_DomainLength;
        private double m_Twall;
        public ParameterizedLevelSetTest(double _DomainLength, double _Twall) {
            m_DomainLength = _DomainLength;
            m_Twall = _Twall;
        }

        public bool TestImmersedBoundary => false;

        /// <summary>
        /// 
        /// </summary>
        public Func<double[], double, double> GetPhi2() {
            throw new NotImplementedException(); // will never be called, as long as 'TestImmersedBoundary' == false;
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            throw new NotImplementedException();
        }

        public int SpatialDimension {
            get {
                return 2;
            }
        }

        public Func<double[], double, double> GetPhi() {
            double xSemiAxis = 2.0;
            double ySemiAxis = 2.0;
            double yCenter = 2.0;
            return PhiFunc(xSemiAxis, ySemiAxis, yCenter);
        }

        private static Func<double[], double, double> PhiFunc(double xSemiAxis, double ySemiAxis, double yCenter) {
            //C.ParameterizedLevelSetControl = new ParameterizedLevelSetControl(2.0, 2.0, 2.0); // xSemiAxis0, ySemiAxis0, yCenter0
            return ((_3D)((time, x, y) => y - ySemiAxis - ySemiAxis * (1 - x * x / (xSemiAxis * xSemiAxis) ).Sqrt() )).Convert_txy2Xt();
        }

         public Func<double[], double, double> GetU(string species, int d) {
            //double L = 0.00001;
             if (d == 0) {
                return ((_3D)((t, x, y) => 0.0)).Convert_txy2Xt();
             } else if (d == 1) {
                 if (species == "A") {
                    return ((_3D)((t, x, y) => (0.1 - 0.1 * x.Pow2() / m_DomainLength.Pow2()))).Convert_txy2Xt();
                 } else if (species == "B") {
                    return ((_3D)((t, x, y) => 0.0)).Convert_txy2Xt();
                 } else {
                   throw new ArgumentOutOfRangeException();
                 }
             } else {
                throw new ArgumentOutOfRangeException();
             }
         }

        Func<double[], double, double> GetTemperature(string species) {
            //double T_wall = 163.75;
            if (species == "A") {
                return ((_3D)((t, x, y) => m_Twall)).Convert_txy2Xt();
            } else if (species == "B") {
                return ((_3D)((t, x, y) => 0.0)).Convert_txy2Xt();

            } else {
                throw new ArgumentOutOfRangeException();
            }

        }


        public double dt {
            get {
                return 0.05;
            }
        }

        public GridCommons CreateGrid(int Resolution) {
            //double L = 0.00001;
            if (Resolution < 1)
                throw new ArgumentException();

            double[] Xnodes = GenericBlas.Linspace(-m_DomainLength, m_DomainLength, 2 * Resolution + 1); 
            double[] Ynodes = GenericBlas.Linspace(0, 5 * m_DomainLength, 5 * Resolution + 1); 

            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes); // creation of 2D grid 


            grd.EdgeTagNames.Add(1, "NavierSlip_Linear_ZeroGradient");
            grd.EdgeTagNames.Add(2, "Pressure_Outlet_ZeroGradient");
            grd.EdgeTagNames.Add(3, "Velocity_Inlet_ConstantTemperature");

            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if ((Math.Abs(X[0] - Xnodes.First()) < 1e-8) || (Math.Abs(X[0] - Xnodes.Last()) < 1e-8))
                    return 1; // walls
                else if ((Math.Abs(X[1] - Ynodes.Last()) < 1e-8))
                    return 2; // upper border
                else if ((Math.Abs(X[1] - Ynodes.First()) < 1e-8))
                    return 3; // bottom

                return et;
            });

            return grd;
        }

public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

            config.Add("NavierSlip_Linear_ZeroGradient", new AppControl.BoundaryValueCollection());
            config.Add("Pressure_Outlet_ZeroGradient", new AppControl.BoundaryValueCollection());
            config.Add("Velocity_Inlet_ConstantTemperature", new AppControl.BoundaryValueCollection());

            //config.Add("freeslip", new AppControl.BoundaryValueCollection());
            //config.Add("navierslip_linear", new AppControl.BoundaryValueCollection());
            //config.Add("pressure_outlet", new AppControl.BoundaryValueCollection());

            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#A",
                (X, t) => 0.0);
            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add(
                VariableNames.Velocity_d(1) + "#A",
                (X, t) => (0.1 - 0.1 * X[0].Pow2() / m_DomainLength.Pow2()));
            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add(
                VariableNames.Velocity_d(0) + "#B",
                (X, t) => 0.0);
            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add(
                VariableNames.Temperature + "#A",
                (X, t) => m_Twall);
            config["Velocity_Inlet_ConstantTemperature"].Evaluators.Add(
                VariableNames.Temperature + "#B",
                (X, t) => 0.0);

            return config;
        }


        public bool isGravity {
            get {
                return false;
            }
        }

        public Func<double[], double> GetF(string species, int d) {
            return (X => 0.0);
        }

        public Func<double[], double, double> GetPress(string species) {
            return ((_3D)((t, x, y) => 0)).Convert_txy2Xt();
        }

        public bool setExtSol {
            get {
                return false;
            }
        }

        /// <summary>
        ///
        /// </summary>
        public double rho_A {
            get {
                return  600.0;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public double rho_B {
            get {
                return  9.0;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public double mu_A {
            get {
                return  0.013;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public double mu_B {
            get {
                return  0.0001019;
            }
        }

        public double Sigma {
            get {
                return 0.0;
            }
        }

        public double c_A {
            get {
                return 100;
            }
        }

        public double c_B {
            get {
                return 37;
            }
        }
        public double k_A {
            get {
                return 4.8;
            }
        }

        public double k_B {
            get {
                return 0.0251;
            }
        }

        public double hVap {
            get {
                return 10000;
            }
        }
        public double T_sat {
            get {
                return 0;
            }
        }

        public double theta_e {
            get {
                return Math.PI / 6.0;
            }
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

        public bool solveCoupledHeatEquation {
            get {
                return true;
            }
        }

        public bool IncludeRecoilPressure {
            get {
                return true;
            }
        }

        public double[] AcceptableL2Error {
            get {
                return new double[] { 1.0e-6, 1.0e-6, 0.01 };
            }
        }

        public double[] AcceptableResidual {
            get {
                return new double[] { 1.0e-6, 1.0e-6, 1.0e-6 };
            }
        }

        public int LevelsetPolynomialDegree => throw new NotImplementedException();
    }
}