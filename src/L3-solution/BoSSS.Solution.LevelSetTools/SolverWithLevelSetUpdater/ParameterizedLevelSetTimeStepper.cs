using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.LevelSetTools.ParameterizedLevelSet;

namespace BoSSS.Solution.LevelSetTools.ParameterizedLevelSet {

    /// <summary>
    /// options for the time-discretization of the Parameterized level-set evolution
    /// </summary>
    public enum Parameterized_Timestepper {
        /// <summary>
        /// Explicit Euler
        /// </summary>
        ExplicitEuler,
    }

    public class ParameterizedLevelSetTimeStepper {


        //internal double[] current_PLSproperty;

        //internal double[] previous_PLSproperty;

        protected Parameterized_Timestepper Timestepper;

        internal RungeKuttaScheme RKscheme;

        /// 
        /// </summary>
        /// <param name="Control"></param>
        /// <param name="currentState"></param>
        public ParameterizedLevelSetTimeStepper(ParameterizedLevelSetControl Control, RungeKuttaScheme _RKscheme) {

            this.Timestepper = Control.Timestepper;
            this.RKscheme = _RKscheme;
        }

        public (double xSemi1, double ySemi1, double yCenter1) MoveLevelSet(double dt, double forceX, double xSemi0, double ySemi0, double yCenter0) {

            double[] CoeffValues0 = new double[3] { 0, 0, 0 };
            double[] ParamValues0 = new double[3] { xSemi0, ySemi0, yCenter0 };
            double[] LowerBoundaryValues = new double[3] { -10, -10, -10 };
            double[] UpperBoundaryValues = new double[3] { 10, 10, 10 };

            var CoeffValues1 = GradientMinimizationMethod.Minimi(CoeffValues0, ParamValues0, MyRealF.F, NumericalGradient.Differentiate(MyRealF.F), LowerBoundaryValues, UpperBoundaryValues, dt, forceX);

            var xSemi1 = CoeffValues1[0] * dt + ParamValues0[0];
            var ySemi1 = CoeffValues1[1] * dt + ParamValues0[1];
            var yCenter1 = CoeffValues1[2] * dt + ParamValues0[2];

            return (xSemi1, ySemi1, yCenter1);
        }

    }
    /// <summary>
    /// 
    /// </summary>
    public static class GradientMinimizationMethod {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="InitialValueForCoeff"></param>
        /// <param name="EllipsePar0"></param>
        /// <param name="F"></param>
        /// <param name="Differentiate"></param>
        /// <param name="LowerSafeguard"></param>
        /// <param name="UpperSafeguard"></param>
        /// <param name="dt"></param>
        /// <param name="forceX"></param>
        /// <returns></returns>
        public static double[] Minimi(double[] InitialValueForCoeff, double[] EllipsePar0, Func<double[], double[], double, double, double> F, Func<double[], double[], double, double, double[]> Differentiate, double[] LowerSafeguard, double[] UpperSafeguard, double dt, double forceX) {

            int L = InitialValueForCoeff.Length;
            double[] CoeffOldIter = new double[L];
            double[] CoeffNewIter = new double[L];
            for (int i = 0; i < L; ++i) {
                CoeffOldIter[i] = InitialValueForCoeff[i];
                CoeffNewIter[i] = InitialValueForCoeff[i];
            }
            double[] derivOldIter = new double[L];
            double[] derivNewIter = Differentiate(CoeffOldIter, EllipsePar0, dt, forceX);
            double eps = 1e-16;
            double lambda = 0.001;
            Console.WriteLine("iter {0}: [{1}]", 0, string.Join(", ", CoeffOldIter));

            double funcNewIter = F(CoeffNewIter, EllipsePar0, dt, forceX);
            double funcOldIter = 0.0;
            double numerator = 0.0;
            double denumerator = 0.0;
            double exitCriteria = 0.0;

            for (int i = 1; i < 100; ++i) {
                //Update Schema
                funcOldIter = funcNewIter;
                for (int d = 0; d < derivOldIter.Length; ++d) {
                    CoeffOldIter[d] = CoeffNewIter[d];
                    derivOldIter[d] = derivNewIter[d];
                    CoeffNewIter[d] -= lambda * derivOldIter[d];
                }
                funcNewIter = F(CoeffNewIter, EllipsePar0, dt, forceX);

                //Lambda and exit criteria calculation
                derivNewIter = Differentiate(CoeffNewIter, EllipsePar0, dt, forceX);
                numerator = 0.0;
                denumerator = 0.0;
                exitCriteria = 0.0;
                for (int j = 0; j < CoeffOldIter.Length; ++j) {
                    exitCriteria += derivNewIter[j] * derivNewIter[j];
                    numerator += (CoeffNewIter[j] - CoeffOldIter[j]) * (derivNewIter[j] - derivOldIter[j]);
                    denumerator += (derivNewIter[j] - derivOldIter[j]) * (derivNewIter[j] - derivOldIter[j]);
                }
                numerator = Math.Abs(numerator);
                lambda = numerator / denumerator;

                //Console Output
                Console.WriteLine("iter {0}: [{1}]", i, string.Join(", ", CoeffNewIter));

                //Exit loop condition
                if (exitCriteria < eps) {
                    break;
                }
            }

            Console.WriteLine($"Minimum of function: {funcNewIter} at  X value: [{string.Join(", ", CoeffNewIter)}]");
            return CoeffNewIter;

        }

    }

    public static class NumericalGradient {

        static public Func<double[], double[], double, double, double[]> Differentiate(Func<double[], double[], double, double, double> F) {

            return delegate (double[] X, double[] EllPar0, double dt, double forceX) {
                double gridStepSize = 1.0e-10;
                int L = X.Length;
                double[] dF_dX = new double[L];
                double[] XLeft = new double[L];
                double[] XRight = new double[L];

                for (int i = 0; i < L; ++i) {
                    XLeft[i] = X[i];
                    XRight[i] = X[i];
                }

                for (int i = 0; i < L; ++i) {
                    if (i > 0) {
                        XLeft[i - 1] = X[i - 1];
                        XRight[i - 1] = X[i - 1];
                    }

                    XLeft[i] -= gridStepSize;
                    XRight[i] += gridStepSize;

                    //Console.WriteLine($" {X[i ]}   {XRight[i ]} {XLeft[i]}");
                    dF_dX[i] = 0.5 * (F(XRight, EllPar0, dt, forceX) - F(XLeft, EllPar0, dt, forceX)) / gridStepSize;
                }
                return dF_dX;
            };
        }


    }

    public static class IntegrationOverEllipse {

        static double stepNumber = 100000;


        public static double IntegralCalculation(double[] EllipseCoeff, double[] EllipsePar0, Func<double[], double[], double, double, double, double, double> f, double alphaMin, double alphaMax, double dt, double forceX) {

            Console.WriteLine($"IntegralCalculation {EllipseCoeff[0]},{EllipseCoeff[1]},{EllipseCoeff[2]} {alphaMin} {alphaMax} {stepNumber}");
            double a = EllipsePar0[0];
            double b = EllipsePar0[1];
            double yc = EllipsePar0[2];

            Func<double, double> x = (alpha) => a * Math.Cos(alpha);
            Func<double, double> y = (alpha) => b * Math.Sin(alpha) + yc;
            Func<double, double> DerivativeX = (alpha) => -a * Math.Sin(alpha);
            Func<double, double> DerivativeY = (alpha) => b * Math.Cos(alpha);


            double gramDeter(double alpha) {
                return Math.Sqrt(DerivativeX(alpha) * DerivativeX(alpha) + DerivativeY(alpha) * DerivativeY(alpha));
            }


            double funcForIntegration(double alpha) {
                double gramDeter_alpha = gramDeter(alpha);
                double q_alpha = f(EllipseCoeff, EllipsePar0, x(alpha), y(alpha), dt, forceX);

                if (double.IsNaN(gramDeter_alpha) || double.IsInfinity(gramDeter_alpha) || double.IsInfinity(q_alpha) || double.IsInfinity(q_alpha))
                    throw new ArithmeticException($" {EllipseCoeff} {gramDeter_alpha}, {q_alpha}");
                return gramDeter_alpha * q_alpha;
            }

            double stepSize = 0.5 * (alphaMax - alphaMin) / stepNumber;
            double integralResult = (funcForIntegration(alphaMin) + funcForIntegration(alphaMax)) / 3;

            for (int i = 1; i <= stepNumber; i++) {

                integralResult += (4 * funcForIntegration(alphaMin + (2 * i - 1) * stepSize) + 2 * funcForIntegration(alphaMin + 2 * i * stepSize)) / 3;
            }

            //var integral = integralResult * stepSize;
            Console.WriteLine("integral = " + integralResult * stepSize);
            return integralResult * stepSize;
        }

    }

    public static class MyRealF {
        public static double f(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double force) {

            // double timeStep = 0.2;
            double time0 = 0.2;
            double time1 = time0 + timeStep;

            double SquareRoot = (EllipsePar0[0].Pow2() - x.Pow2()).Sqrt(); //sqrt(a^2 - x^2)
                                                                           //Console.WriteLine("SquareRoot = " + SquareRoot);

            //double[] Velocity = new double[Dim] {x, y};
            //double[] Grad = new double[Dim] {-x * EllipsePar0[1]  / ( EllipsePar0[0] * SquareRoot), 1.0};

            double term1 = -EllipseCoeff[2] * time1 / timeStep;
            double term2 = EllipseCoeff[1] * time1 / (timeStep * EllipsePar0[0]) * SquareRoot;
            double term3 = EllipsePar0[1] * x.Pow2() * EllipseCoeff[0] * time1 / (timeStep * EllipsePar0[0].Pow2() * SquareRoot);

            //Console.WriteLine($"functerm { term1},{-0.2 * EllipseCoeff[2]},{term2} {0.2 * EllipseCoeff[1] * SquareRoot / EllipsePar0[0] } ");
            //return (EllipseCoeff[0] - 1.0) * x.Pow2() + EllipseCoeff[2];
            //return -2 * EllipseCoeff[2] + 2 * EllipseCoeff[1] * SquareRoot / EllipsePar0[0] + 1.0 + 2 * EllipsePar0[1] * x.Pow2() * EllipseCoeff[0] / (EllipsePar0[0].Pow2() * SquareRoot);
            //return term1 + term2 + term3 + Velocity[0] * Grad[0] + Velocity[1] * Grad[1];
            return term1 + term2 + term3 + force;
        }

        public static double fsqr(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double force) {
            var fval = f(EllipseCoeff, EllipsePar0, x, y, timeStep, force);
            return fval * fval;
        }

        public static double F(double[] EllipseCoeff, double[] EllipsePar0, double dt, double forceX) {
            const double alphaMin = -Math.PI * 5 / 6;
            const double alphaMax = -Math.PI / 6;

            return IntegrationOverEllipse.IntegralCalculation(EllipseCoeff, EllipsePar0, fsqr, alphaMin, alphaMax, dt, forceX);
        }

    }
}

