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
using System.Diagnostics;

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

        public double[] MoveLevelSet(double dt, double time, SinglePhaseField[] meanVelocity, double xSemi0, double ySemi0, double yCenter0, IGridData gdat, CellQuadratureScheme levSetQuadScheme, int quadOrder) {

            double[] CoeffValues0 = new double[3] { 0, 0, 0 };
            double[] ParamValues0 = new double[3] { xSemi0, ySemi0, yCenter0 };
            double[] LowerBoundaryValues = new double[3] { -1, -5, -5 };
            double[] UpperBoundaryValues = new double[3] { 1, 5, 5 };

            var myf = new MyRealF(levSetQuadScheme, gdat, quadOrder);


            var CoeffValues1 = GradientMinimizationMethod.Minimi(CoeffValues0, ParamValues0, myf.F, NumericalGradient.Differentiate(myf.F), LowerBoundaryValues, UpperBoundaryValues, dt, time, meanVelocity);

            double[] ParamValues1 = new double[CoeffValues1.Length];
            for (int i = 0; i < CoeffValues1.Length; i++) {
                 ParamValues1[i] = CoeffValues1[i] * dt + ParamValues0[i]; // Multiply each element by the scalar
                 if ((ParamValues1[i] > UpperBoundaryValues[i]) || (ParamValues1[i] < LowerBoundaryValues[i])) {
                    ParamValues1[i] = ParamValues0[i];
                    CoeffValues1[i] = 0;
                }
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
        public static double[] Minimi(double[] InitialValueForCoeff, double[] EllipsePar0, Func<double[], double[], double, double, SinglePhaseField[], double> F, Func<double[], double[], double, double, SinglePhaseField[], double[]> Differentiate, double[] LowerSafeguard, double[] UpperSafeguard, double dt, double time, SinglePhaseField[]  meanVelocity) {

            int L = InitialValueForCoeff.Length;
            double[] CoeffOldIter = new double[L];
            double[] CoeffNewIter = new double[L];
            for (int i = 0; i < L; ++i) {
                CoeffOldIter[i] = InitialValueForCoeff[i];
                CoeffNewIter[i] = InitialValueForCoeff[i];
            }
            double[] derivOldIter = new double[L];
            double[] derivNewIter = Differentiate(CoeffOldIter, EllipsePar0, dt, time, meanVelocity);
            double eps = 1e-16;
            double lambda = 0.001;
            Console.WriteLine("iter {0}: [{1}]", 0, string.Join(", ", CoeffOldIter));

            double funcNewIter = F(CoeffNewIter, EllipsePar0, dt, time, meanVelocity);
            double funcOldIter = 0.0;
            double numerator = 0.0;
            double denumerator = 0.0;
            double exitCriteria = 0.0;
            int itermax = 150000;

            for (int i = 1; i < itermax; ++i) {
                //Update Schema
                funcOldIter = funcNewIter;
                for (int d = 0; d < derivOldIter.Length; ++d) {
                    CoeffOldIter[d] = CoeffNewIter[d];
                    derivOldIter[d] = derivNewIter[d];
                    CoeffNewIter[d] -= lambda * derivOldIter[d];
                }
                funcNewIter = F(CoeffNewIter, EllipsePar0, dt, time, meanVelocity);

                //Lambda and exit criteria calculation
                derivNewIter = Differentiate(CoeffNewIter, EllipsePar0, dt, time, meanVelocity);
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
                //if (CoeffNewIter[0] < 0) {
                //    CoeffNewIter[0]=0;
                //}
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
        static public Func<double[], double[], double, double, SinglePhaseField[], double[]> Differentiate(Func<double[], double[], double, double, SinglePhaseField[], double> F) {

            return delegate (double[] X, double[] EllPar0, double dt, double time, SinglePhaseField[] meanVelocity) {
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
                    dF_dX[i] = 0.5 * (F(XRight, EllPar0, dt, time, meanVelocity) - F(XLeft, EllPar0, dt, time, meanVelocity)) / gridStepSize;
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
        public static double IntegralCalculation(double[] EllipseCoeff, double[] EllipsePar0, Func<double[], double[], double, double, double, double, double[], double> f, double alphaMin, double alphaMax, double dt, double time, double[] velocity) {

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
                double q_alpha = f(EllipseCoeff, EllipsePar0, x(alpha), y(alpha), dt, time, velocity);

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
            Console.WriteLine("integral_simpson = " + integralResult * stepSize);
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
        static double f(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double time, double[] velocity) {

            // double timeStep = 0.2;
            //double time0 = 0;
            //double time1 = time0 + timeStep;

            //double SquareRoot = (EllipsePar0[0].Pow2() - x.Pow2()).Sqrt(); //sqrt(a^2 - x^2)
            //                                                               //Console.WriteLine("SquareRoot = " + SquareRoot);

            ////double[] Velocity = new double[Dim] {x, y};
            //double Gradx =  -x * EllipsePar0[1] / (EllipsePar0[0] * SquareRoot);
            //double Grady = 1.0;
            ////Debugger.Launch();
            //double Grad = (Gradx.Pow2() + Grady.Pow2()).Sqrt();

            //double term1 = -EllipseCoeff[2];
            //double term2 = EllipseCoeff[1] / EllipsePar0[0] * SquareRoot;
            //double term3 = EllipsePar0[1] * x.Pow2() * EllipseCoeff[0]  / ( EllipsePar0[0].Pow2() * SquareRoot);

            double SquareRoot = (EllipsePar0[1].Pow2() * (1 - x.Pow2() / EllipsePar0[0].Pow2())).Sqrt();
            double numer =  EllipsePar0[1] * EllipseCoeff[1] * EllipsePar0[0] * (EllipsePar0[0].Pow2() - x.Pow2()) + x.Pow2() * EllipsePar0[1].Pow2() * EllipseCoeff[0];
            double DerivFuncEllipse = -EllipseCoeff[2] + numer / (SquareRoot * EllipsePar0[0].Pow2() * EllipsePar0[0]);
            double Gradx =  -EllipsePar0[1].Pow2() * x / (EllipsePar0[0].Pow2() * SquareRoot);
            double Grady = 1;
            //double Grad = (Gradx.Pow2() + Grady.Pow2()).Sqrt();
            double res  = velocity[0] * Gradx + velocity[1] * Grady;
            //Console.WriteLine($"forceX: {res} ");
            return DerivFuncEllipse + res;

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
        static double fsqr(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double time, double[] velocity) {
            var fval = f(EllipseCoeff, EllipsePar0, x, y, timeStep, time, velocity);
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
        public double F(double[] EllipseCoeff, double[] EllipsePar0, double dt, double time, SinglePhaseField[] meanVelocity) {
            




            void VectorizedEvaluate(int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) {
                var PhysicalNodes = gdat.GlobalNodes.GetValue_Cell(rule.Nodes, i0, Length);
                var VelocityValues = MultidimensionalArray.Create(Length, rule.NoOfNodes, 2);
                for (int d = 0; d < 2; d++) {
                    meanVelocity[d].Evaluate(i0, Length, rule.Nodes, VelocityValues.ExtractSubArrayShallow(-1, -1, d));
                }
                for (int i = 0; i < Length; i++) {
                    for (int k = 0; k < rule.NoOfNodes; k++) {
                        double x = PhysicalNodes[i, k, 0];
                        double y = PhysicalNodes[i, k, 1];
                        double[] velocity = new double[] { VelocityValues[i, k, 0], VelocityValues[i, k, 1] };
                        EvalResult[i, k, 0] = fsqr(EllipseCoeff, EllipsePar0, x, y, dt, time, velocity);
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
                double angle = Math.Atan((3.0 / 16).Sqrt());
                double alphaMin = -Math.PI + angle ;
                double alphaMax = -angle;
                //double TotalIntegral_RefVal = IntegrationOverEllipse.IntegralCalculation(EllipseCoeff, EllipsePar0, fsqr, alphaMin, alphaMax, dt, time, velocity);
                Console.WriteLine($"integral_bosss: {TotalIntegral}");
               // Console.WriteLine($"integral_diff: {TotalIntegral_RefVal - TotalIntegral}");
            }

            return TotalIntegral;


        }

    }
}

