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
using BoSSS.Foundation.Quadrature;

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

        protected Parameterized_Timestepper Timestepper;


        /// 
        /// </summary>
        /// <param name="Control"></param>
        /// <param name="currentState"></param>
        public ParameterizedLevelSetTimeStepper(ParameterizedLevelSetControl Control) {

            this.Timestepper = Control.Timestepper;
        }

        public double[] MoveLevelSet(double dt, double time, double forceX, double xSemi0, double ySemi0, double yCenter0, IGridData gdat, CellQuadratureScheme levSetQuadScheme, int quadOrder) {

            double[] CoeffValues0 = new double[3] { 0, 0, 0 };
            double[] ParamValues0 = new double[3] { xSemi0, ySemi0, yCenter0 };
            double[] LowerBoundaryValues = new double[3] { -10, -10, -10 };
            double[] UpperBoundaryValues = new double[3] { 10, 10, 10 };

            var myf = new MyRealF(levSetQuadScheme, gdat, quadOrder);


            var CoeffValues1 = GradientMinimizationMethod.Minimi(CoeffValues0, ParamValues0, myf.F, NumericalGradient.Differentiate(myf.F), LowerBoundaryValues, UpperBoundaryValues, dt, time, forceX);

            double[] ParamValues1 = new double[CoeffValues1.Length];
            for (int i = 0; i < CoeffValues1.Length; i++) {
                ParamValues1[i] = CoeffValues1[i] * dt + ParamValues0[i]; // Multiply each element by the scalar
            }


            return ParamValues1;
        }

    }

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
        public static double[] Minimi(double[] InitialValueForCoeff, double[] EllipsePar0, Func<double[], double[], double, double, double, double> F, Func<double[], double[], double, double, double, double[]> Differentiate, double[] LowerSafeguard, double[] UpperSafeguard, double dt, double time, double forceX) {

            int L = InitialValueForCoeff.Length;
            double[] CoeffOldIter = new double[L];
            double[] CoeffNewIter = new double[L];
            for (int i = 0; i < L; ++i) {
                CoeffOldIter[i] = InitialValueForCoeff[i];
                CoeffNewIter[i] = InitialValueForCoeff[i];
            }
            double[] derivOldIter = new double[L];
            double[] derivNewIter = Differentiate(CoeffOldIter, EllipsePar0, dt, time, forceX);
            double eps = 1e-20;
            double lambda = 0.001;
            Console.WriteLine("iter {0}: [{1}]", 0, string.Join(", ", CoeffOldIter));

            double funcNewIter = F(CoeffNewIter, EllipsePar0, dt, time, forceX);
            double funcOldIter = 0.0;
            double numerator = 0.0;
            double denumerator = 0.0;
            double exitCriteria = 0.0;
            int itermax = 5000;

            for (int i = 1; i < itermax; ++i) {
                //Update Schema
                funcOldIter = funcNewIter;
                for (int d = 0; d < derivOldIter.Length; ++d) {
                    CoeffOldIter[d] = CoeffNewIter[d];
                    derivOldIter[d] = derivNewIter[d];
                    CoeffNewIter[d] -= lambda * derivOldIter[d];
                }
                funcNewIter = F(CoeffNewIter, EllipsePar0, dt, time, forceX);

                //Lambda and exit criteria calculation
                derivNewIter = Differentiate(CoeffNewIter, EllipsePar0, dt, time, forceX);
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
                if (i == itermax) {
                    throw new ApplicationException("Increase the number of iterations for the gradient method!");
                }
            }

            Console.WriteLine($"Minimum of function: {funcNewIter} at  X value: [{string.Join(", ", CoeffNewIter)}]");
            return CoeffNewIter;

        }

    }
 
    public static class NumericalGradient {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="F"></param>
        /// <returns></returns>
        static public Func<double[], double[], double, double, double, double[]> Differentiate(Func<double[], double[], double, double, double, double> F) {

            return delegate (double[] X, double[] EllPar0, double dt, double time, double forceX) {
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
                    dF_dX[i] = 0.5 * (F(XRight, EllPar0, dt, time, forceX) - F(XLeft, EllPar0, dt, time, forceX)) / gridStepSize;
                }
                return dF_dX;
            };
        }


    }

    public static class IntegrationOverEllipse {

        static double stepNumber = 100000;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="EllipseCoeff"></param>
        /// <param name="EllipsePar0"></param>
        /// <param name="f"></param>
        /// <param name="alphaMin"></param>
        /// <param name="alphaMax"></param>
        /// <param name="dt"></param>
        /// <param name="forceX"></param>
        /// <returns></returns>
        /// <exception cref="ArithmeticException"></exception>
        public static double IntegralCalculation(double[] EllipseCoeff, double[] EllipsePar0, Func<double[], double[], double, double, double, double, double, double> f, double alphaMin, double alphaMax, double dt, double time, double forceX) {

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
                double q_alpha = f(EllipseCoeff, EllipsePar0, x(alpha), y(alpha), dt, time, forceX);

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

    public class MyRealF {

        CellQuadratureScheme levSetQuadScheme;
        IGridData gdat;
        int quadOrder;

        public MyRealF(CellQuadratureScheme _levSetQuadScheme, IGridData _gdat, int _quadOrder) {
            levSetQuadScheme = _levSetQuadScheme;
            gdat = _gdat;
            quadOrder = _quadOrder;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="EllipseCoeff"></param>
        /// <param name="EllipsePar0"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="timeStep"></param>
        /// <param name="force"></param>
        /// <returns></returns>
        static double f(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double time, double force) {

            // double timeStep = 0.2;
            //double time0 = 0;
            //double time1 = time0 + timeStep;

            double SquareRoot = (EllipsePar0[0].Pow2() - x.Pow2()).Sqrt(); //sqrt(a^2 - x^2)
                                                                           //Console.WriteLine("SquareRoot = " + SquareRoot);

            //double[] Velocity = new double[Dim] {x, y};
            //double[] Grad = new double[Dim] {-x * EllipsePar0[1]  / ( EllipsePar0[0] * SquareRoot), 1.0};

            double term1 = -EllipseCoeff[2] * timeStep / timeStep;
            double term2 = EllipseCoeff[1] * timeStep / (timeStep * EllipsePar0[0]) * SquareRoot;
            double term3 = EllipsePar0[1] * x.Pow2() * EllipseCoeff[0] * timeStep / (timeStep * EllipsePar0[0].Pow2() * SquareRoot);

            //Console.WriteLine($"functerm { term1},{-0.2 * EllipseCoeff[2]},{term2} {0.2 * EllipseCoeff[1] * SquareRoot / EllipsePar0[0] } ");
            //return (EllipseCoeff[0] - 1.0) * x.Pow2() + EllipseCoeff[2];
            //return -2 * EllipseCoeff[2] + 2 * EllipseCoeff[1] * SquareRoot / EllipsePar0[0] + 1.0 + 2 * EllipsePar0[1] * x.Pow2() * EllipseCoeff[0] / (EllipsePar0[0].Pow2() * SquareRoot);
            //return term1 + term2 + term3 + Velocity[0] * Grad[0] + Velocity[1] * Grad[1];
            return term1 + term2 + term3 + force;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="EllipseCoeff"></param>
        /// <param name="EllipsePar0"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="timeStep"></param>
        /// <param name="force"></param>
        /// <returns></returns>
        static double fsqr(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double time, double force) {
            var fval = f(EllipseCoeff, EllipsePar0, x, y, timeStep, time, force);
            return fval * fval;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="EllipseCoeff"></param>
        /// <param name="EllipsePar0"></param>
        /// <param name="dt"></param>
        /// <param name="forceX"></param>
        /// <returns></returns>
        public double F(double[] EllipseCoeff, double[] EllipsePar0, double dt, double time, double forceX) {
            




            void VectorizedEvaluate(int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) {
                var PhysicalNodes = gdat.GlobalNodes.GetValue_Cell(rule.Nodes, i0, Length);
                for (int i = 0; i < Length; i++) {
                    for (int k = 0; k < rule.NoOfNodes; k++) {
                        double x = PhysicalNodes[i, k, 0];
                        double y = PhysicalNodes[i, k, 1];
                        EvalResult[i, k, 0] = fsqr(EllipseCoeff, EllipsePar0, x, y, dt, time, forceX);
                    }
                }
            }

            double TotalIntegral = 0;
            void SafeResult(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    //int jCell = i + i0;
                    TotalIntegral += ResultsOfIntegration[i, 0];
                }
            }


            var q = CellQuadrature.GetQuadrature(new int[] { 1 }, gdat, levSetQuadScheme.Compile(gdat, quadOrder),
                VectorizedEvaluate, SafeResult);
            q.Execute();

            {
                const double alphaMin = -Math.PI * 5 / 6;
                const double alphaMax = -Math.PI / 6;
                double TotalIntegral_RefVal = IntegrationOverEllipse.IntegralCalculation(EllipseCoeff, EllipsePar0, fsqr, alphaMin, alphaMax, dt, time, forceX);
                Console.WriteLine(TotalIntegral_RefVal - TotalIntegral);
            }

            return TotalIntegral;


        }

    }
}

