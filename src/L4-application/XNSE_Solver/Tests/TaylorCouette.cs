using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Tests {


    /// <summary>
    /// Two-phase Taylor Couette flow, used also as a testscase in publication;
    /// 
    /// </summary>
    public class TaylorCouette : IXNSETest {

        public enum Mode {

            /// <summary>
            /// Entire domain filled with species A
            /// </summary>
            TestIBM = 1,

            /// <summary>
            /// Phase boundary within the domain, <see cref="TestImmersedBoundary"/> is false;
            /// Note: at the time of writing, we are missing the implementation of the required viscosity mode, 
            /// so we cannot test IBM and 2-phase simultaneously.
            /// </summary>
            Test2Phase = 2
        }


        public TaylorCouette(Mode m) {

            double Ra = double.NaN, Ri = double.NaN, Rm = double.NaN;
            switch(m) {
                case Mode.Test2Phase:
                Ra = 2; // größer Inkreis
                Ri = Math.Sqrt(2) / 2; // kleinster Inkreis im Gebiet
                Rm = 0.5 * (Ra + Ri); ; // in the middle
                this.TestImmersedBoundary = false;
                break;

                case Mode.TestIBM:
                Ri = Math.Sqrt(2) / 2; // kleinster Inkreis im Gebiet
                Rm = Math.Sqrt(2 * 2 + 2 * 2) + 0.5; // ausserhalb
                Ra = Rm + 1; // noch weiter draußen
                this.TestImmersedBoundary = true;
                break;
            }
            exS = new ExactSol(Ri, Rm, Ra);


        }


        public double mu_A => ExactSol.muA;

        public double mu_B => ExactSol.muB;

        public double Sigma => ExactSol.sigma;

        public bool TestImmersedBoundary {
            get;
            private set;
        }



        public int SpatialDimension => 2;

        public double dt => throw new NotImplementedException();

        public double rho_A => ExactSol.rhoA;

        public double rho_B => ExactSol.rhoB;

        public bool Material => true;

        public bool steady => true;

        public bool IncludeConvection => true;

        public int LevelsetPolynomialDegree => 2;


        ExactSol exS;

        /// <summary>
        /// 
        /// </summary>
        class ExactSol {
            public const double Ui = 2; // inner radius tangential speed
            public const double Ua = 1; // outer radius tangential speed
            public const double rhoA = 0.1;
            public const double rhoB = 1.3;
            public const double muA = 0.1;
            public const double muB = 0.2;
            public const double sigma = 0.9;

            public double Ri = Math.Sqrt(2) / 2;
            
            public double Ra = 2;

            public double Rm; //  = (Ri + Ra) / 2;

            public ExactSol(double _Ri, double _Rm, double _Ra) {
                this.Rm = _Rm;
                this.Ra = _Ra;
                this.Ri = _Ri;

                _C1A = (Ra.Pow2() * Ri * Ui * muA - Ra.Pow2() * Ri * Ui * muB + Ra * Rm.Pow2() * Ua * muB - Ri * Rm.Pow2() * Ui * muA) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA);
                _C1B = (Ra * Ri.Pow2() * Ua * muA - Ra * Ri.Pow2() * Ua * muB + Ra * Rm.Pow2() * Ua * muB - Ri * Rm.Pow2() * Ui * muA) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA);
                _C2A = Ri * Ra * Rm.Pow2() * muB * (Ra * Ui - Ri * Ua) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA);
                _C2B = Ra * Ri * Rm.Pow2() * muA * (Ra * Ui - Ri * Ua) / (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA);
                _C3A = (1.0 / 2.0) * (-Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA.Pow2() * rhoA 
                    - Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA.Pow2() * rhoB 
                    + Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muB.Pow2() * rhoA 
                    + Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muB.Pow2() * rhoB 
                    - 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muB.Pow2() * rhoB 
                    + 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA.Pow2() * rhoA 
                    + 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow2() * muA * muB * sigma 
                    + 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow2() * muA * muB * sigma 
                    - 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(4) * muA * muB * sigma 
                    - Ra.Pow2() * Rm.Pow(7) * Ua.Pow2() * muB.Pow2() * rhoA + Ra.Pow2() * Rm.Pow(7) * Ua.Pow2() * muB.Pow2() * rhoB 
                    - Ri.Pow2() * Rm.Pow(7) * Ui.Pow2() * muA.Pow2() * rhoA + Ri.Pow2() * Rm.Pow(7) * Ui.Pow2() * muA.Pow2() * rhoB 
                    - 4 * Ra.Pow(4) * Ri.Pow(4) * muA * muB * sigma 
                    - 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow2() * muB.Pow2() * sigma 
                    - 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow2() * muA.Pow2() * sigma + 2 * Ra.Pow(4) * Ri.Pow(4) * muA.Pow2() * sigma 
                    + 2 * Ra.Pow(4) * Ri.Pow(4) * muB.Pow2() * sigma + 2 * Ra.Pow(4) * Rm.Pow(4) * muB.Pow2() * sigma 
                    + 2 * Ri.Pow(4) * Rm.Pow(4) * muA.Pow2() * sigma + 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA.Pow2() * rhoB * Math.Log(Rm) 
                    - 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muB.Pow2() * rhoA * Math.Log(Rm) 
                    - 4 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muB.Pow2() * rhoA * Math.Log(Rm) 
                    + 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muA * muB * rhoB * Math.Log(Rm) 
                    - 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muA * muB * rhoB * Math.Log(Rm) 
                    + 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA * muB * rhoA * Math.Log(Rm) 
                    + 4 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA.Pow2() * rhoB * Math.Log(Rm) 
                    - 2 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muA * muB * rhoA + 2 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA * muB * rhoB 
                    + 2 * Ra * Ri * Rm.Pow(7) * Ua * Ui * muA * muB * rhoA - 2 * Ra * Ri * Rm.Pow(7) * Ua * Ui * muA * muB * rhoB 
                    - 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA * muB * rhoA * Math.Log(Rm) 
                    + 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA * muB * rhoA * Math.Log(Rm) 
                    - 4 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA * muB * rhoB * Math.Log(Rm) 
                    + 4 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muA * muB * rhoB * Math.Log(Rm) 
                    - 4 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA * muB * rhoA * Math.Log(Rm) 
                    + 4 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muB.Pow2() * rhoA * Math.Log(Rm) 
                    - 4 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muA.Pow2() * rhoB * Math.Log(Rm) 
                    + 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muB.Pow2() * rhoA * Math.Log(Rm) 
                    - 4 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA.Pow2() * rhoB * Math.Log(Rm) 
                    + 2 * Ra.Pow(4) * Ri.Pow2() * Rm.Pow(3) * Ui.Pow2() * muA * muB * rhoA 
                    + 2 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muA.Pow2() * rhoB 
                    - 2 * Ra.Pow(3) * Ri.Pow(3) * Rm.Pow(3) * Ua * Ui * muB.Pow2() * rhoA 
                    + 2 * Ra.Pow(3) * Ri * Rm.Pow(5) * Ua * Ui * muB.Pow2() * rhoA 
                    - 2 * Ra.Pow2() * Ri.Pow(4) * Rm.Pow(3) * Ua.Pow2() * muA * muB * rhoB 
                    + 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ua.Pow2() * muA * muB * rhoB 
                    - 2 * Ra.Pow2() * Ri.Pow2() * Rm.Pow(5) * Ui.Pow2() * muA * muB * rhoA 
                    - 2 * Ra * Ri.Pow(3) * Rm.Pow(5) * Ua * Ui * muA.Pow2() * rhoB
                    ) / (Rm * (Ra.Pow2() * Ri.Pow2() * muA - Ra.Pow2() * Ri.Pow2() * muB + Ra.Pow2() * Rm.Pow2() * muB - Ri.Pow2() * Rm.Pow2() * muA).Pow2());
                _C3B = 0;

                if(_C1A.IsNaNorInf())
                    throw new ArithmeticException();
                if(_C2A.IsNaNorInf())
                    throw new ArithmeticException();
                if(_C3A.IsNaNorInf())
                    throw new ArithmeticException();
                if(_C1B.IsNaNorInf())
                    throw new ArithmeticException();
                if(_C2B.IsNaNorInf())
                    throw new ArithmeticException();
                if(_C3B.IsNaNorInf())
                    throw new ArithmeticException();
            }

            double _C1A, _C1B, _C2A, _C2B, _C3A, _C3B;

            /// <summary>
            /// velocity, phase A, radial coordinates
            /// </summary>
            public double vA(double r) {
                return _C1A * r + _C2A / r;
            }

            /// <summary>
            /// velocity, phase B, radial coordinates
            /// </summary>
            public double vB(double r) {
                return _C1B * r + _C2B / r;
            }

            /// <summary>
            /// pressure, phase A, radial coordinates
            /// </summary>
            public double pA(double r) {
                return 0.5 * rhoA * _C1A.Pow2() * r.Pow2() - 0.5 * rhoA * _C2A.Pow2() / (r.Pow2()) + 2 * rhoA * _C1A * _C2A * Math.Log(r) + _C3A;
            }

            /// <summary>
            /// pressure, phase B, radial coordinates
            /// </summary>
            public double pB(double r) {
                return 0.5 * rhoB * _C1B.Pow2() * r.Pow2() - 0.5 * rhoB * _C2B.Pow2() / (r.Pow2()) + 2 * rhoB * _C1B * _C2B * Math.Log(r) + _C3B;
            }


            public double UA1(double[] X, double t) { 
                return (-X[1] / X.L2Norm()) * vA(X.L2Norm()); 
            }
            
            public double UA2(double[] X, double t) { 
                return (+X[0] / X.L2Norm()) * vA(X.L2Norm()); 
            }
            
            public double UB1(double[] X, double t) { 
                return (-X[1] / X.L2Norm()) * vB(X.L2Norm()); 
            }
            
            public double UB2(double[] X, double t) { 
                return (+X[0] / X.L2Norm()) * vB(X.L2Norm()); 
            }

            public double PA(double[] X, double t) { 
                double p = pA(X.L2Norm());
                if(p.IsNaNorInf())
                    throw new ArithmeticException();
                return p;
            }
            
            public double PB(double[] X, double t) { 
                double p = pB(X.L2Norm());
                if(p.IsNaNorInf())
                    throw new ArithmeticException();
                return p;
            }

        }




        public double[] AcceptableL2Error => new double[] { 1, 1, 1 };

        public double[] AcceptableResidual => new double[] { 1, 1, 1 };

        string innerWallTag = IncompressibleBcType.Velocity_Inlet.ToString() + "_inner";
        
        string outerWallTag = IncompressibleBcType.Velocity_Inlet.ToString() + "_outer";
        
        public GridCommons CreateGrid(int Resolution) {


            double[] Xnodes = GenericBlas.Linspace(-2, 2, 8 * Resolution + 1);
            double[] Ynodes = GenericBlas.Linspace(-2, 2, 8 * Resolution + 1);
            var cutOut = new BoundingBox(new double[] { -0.5, -0.5 }, new double[] { +0.5, +0.5 });
            var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, CutOuts: cutOut);
            grd.EdgeTagNames.Add(1, innerWallTag);
            grd.EdgeTagNames.Add(2, outerWallTag);

            grd.DefineEdgeTags(delegate (double[] X) {
                byte et = 0;
                if(Math.Abs(X[0] - (-0.5)) <= 1.0e-8 || Math.Abs(X[0] - (+0.5)) <= 1.0e-8
                    || Math.Abs(X[1] - (-0.5)) <= 1.0e-8 || Math.Abs(X[1] - (+0.5)) <= 1.0e-8)
                    return innerWallTag;
                if(Math.Abs(X[0] - (-2)) <= 1.0e-8 || Math.Abs(X[0] - (+2)) <= 1.0e-8
                    || Math.Abs(X[1] - (-2)) <= 1.0e-8 || Math.Abs(X[1] - (+2)) <= 1.0e-8)
                    return outerWallTag;

                throw new ArgumentOutOfRangeException("error in DefineEdgeTags");
            });

            return grd;
        }

        public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
            var config = new Dictionary<string, AppControl.BoundaryValueCollection>();


            config.Add(innerWallTag, new AppControl.BoundaryValueCollection());
            config[innerWallTag].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", exS.UA1);
            config[innerWallTag].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", exS.UA2);


            config.Add(outerWallTag, new AppControl.BoundaryValueCollection());
            if(TestImmersedBoundary) {
                config[outerWallTag].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", exS.UA1);
                config[outerWallTag].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", exS.UA2);
            } else {
                config[outerWallTag].Evaluators.Add(VariableNames.Velocity_d(0) + "#A", exS.UA1);
                config[outerWallTag].Evaluators.Add(VariableNames.Velocity_d(1) + "#A", exS.UA2);
            }

            return config;
        }

        public Func<double[], double> GetF(string species, int d) {
            return (double[] X) => 0.0;
        }

        /// <summary>
        /// Circle at radius <see cref="Rm"/>;
        /// - interior is negative, i.e. species A
        /// - exterior is positive, i.e. species B
        /// </summary>
        /// <returns></returns>
        public Func<double[], double, double> GetPhi() {
            
            Func<double[], double, double> phiFunc =  (X, t) => X.L2NormPow2() - exS.Rm.Pow2();  // quadratic form

            double[] Xtest = new[] { Math.Cos(0.77) * exS.Rm, Math.Sin(0.77) * exS.Rm };
            var valTest = phiFunc(Xtest, 0.0);


            return phiFunc;
        }

        public Func<double[], double, double> GetPhi2() {
            return (X, t) => - X.L2NormPow2() + (exS.Ri + 0.1).Pow2();  // quadratic form; fluid domain is negative w.r.t. second level-set
        }

        public Func<double[], double, double> GetPhi2U(int d) {
            switch(d) {
                case 0: return exS.UA1;
                case 1: return exS.UA2;
                default: throw new ArgumentOutOfRangeException();
            }
        }

        public Func<double[], double, double> GetPress(string species) {
            switch(species) {
                case "A": return exS.PA;
                case "B": return exS.PB;
                default: throw new ArgumentOutOfRangeException();
            }
        }

        public Func<double[], double, double> GetU(string species, int d) {
            new ValueTuple<string, int>("A", 0);

            switch(d) {
                case 0:
                switch(species) {
                    case "A": return exS.UA1;
                    case "B": return exS.UB1;
                    default: throw new ArgumentOutOfRangeException();
                }

                case 1:
                switch(species) {
                    case "A": return exS.UA2;
                    case "B": return exS.UB2;
                    default: throw new ArgumentOutOfRangeException();
                }

                default: throw new ArgumentOutOfRangeException();
            }
        }
    }
}
