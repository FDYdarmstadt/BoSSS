using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using ilPSP.LinSolvers;
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
        ExplicitEuler,
    }
   
    public class ParameterizedLevelSetTimeStepper {

        protected Parameterized_Timestepper Timestepper;

   
        public ParameterizedLevelSetTimeStepper(ParameterizedLevelSetControl Control) {

            this.Timestepper = Control.Timestepper;
        }


        /// <summary>
        /// calculation of parameters for the current time step
        /// </summary>
        /// <param name="dt"></param> time step size
        /// <param name="time"></param> 
        /// <param name="meanVelocity"></param> velocity from a flow solver
        /// <param name="xSemi0"></param> x-SemiAxis of ellipse
        /// <param name="ySemi0"></param> y-SemiAxis of ellipse
        /// <param name="yCenter0"></param> y-coordinate of center of ellipse
        /// <param name="gdat"></param>
        /// <param name="levSetQuadScheme"></param>
        /// <param name="quadOrder"></param>
        /// <returns></returns>
        public double[] MoveLevelSet(double dt, double time, SinglePhaseField[] meanVelocity, double xSemi0, double ySemi0, double yCenter0, IGridData gdat, CellQuadratureScheme levSetQuadScheme, int quadOrder) {

            double[] CoeffValues0 = new double[3] { 0.0, 0.0, 0.0};
            double[] ParamValues0 = new double[3] { xSemi0, ySemi0, yCenter0 };

            var myf = new MyRealF(levSetQuadScheme, gdat, quadOrder);


            var CoeffValues1 = GradientMinimizationMethod.Minimi(CoeffValues0, ParamValues0, myf.F, NumericalGradient.Differentiate(myf.F), dt, time, meanVelocity);
            //var CoeffValues1 = GradientMinimizationMethod.Minimi_Newton(CoeffValues0, ParamValues0, myf.F, NumericalGradient.CalculateHessian(myf.F), NumericalGradient.Differentiate(myf.F), dt, time, meanVelocity);

            double[] ParamValues1 = new double[CoeffValues1.Length];
            for (int i = 0; i < CoeffValues1.Length; i++) {
                 ParamValues1[i] = CoeffValues1[i] * dt + ParamValues0[i]; 
            }
            return ParamValues1;
        }

    }

    public static class GradientMinimizationMethod {

       /// <summary>
       /// minimization method 
       /// </summary>
        public static double[] Minimi(double[] InitialValueForCoeff, double[] EllipsePar0, Func<double[], double[], double, double, SinglePhaseField[], double> F, Func<double[], double[], double, double, SinglePhaseField[], double[]> Differentiate, double dt, double time, SinglePhaseField[]  meanVelocity) {

            double[] CoeffNewIter = (double[])InitialValueForCoeff.Clone();
            //double[] CoeffLookAhead = (double[])InitialValueForCoeff.Clone();
            double lambda = 20.0;// 0.165;//20

            double[] deriv = new double[InitialValueForCoeff.Length]; // derivative of function
            double eps = 1e-06;//GenericBlas.MachineEps.Sqrt(); 
            double gamma = 0.99; // learning rate
            double[] grad_acc = new double[3] { 0.0, 0.0, 0.0}; // accumulated gradient values

            double[] G = new double[3] { 0.0, 0.0, 0.0 }; //sum of squares of gradients
            double[] correctedG = new double[3] { 0.0, 0.0, 0.0 }; //corrected sum of squares of gradients for zero-bias correction
            double[] lambda_ad = new double[3] { 0.0, 0.0, 0.0 }; //adaptive learning rate
            double beta = 0.999; // decay rate for zero-bias correction


            int itermax = 10000000; // maximum number of iterations
            Stopwatch stw = Stopwatch.StartNew(); // to check time of execution

            for (int i = 1; i < itermax; ++i) {
                //Update Schema
                deriv = Differentiate(CoeffNewIter, EllipsePar0, dt, time, meanVelocity); 
                for (int d = 0; d < deriv.Length; ++d) {
                    G[d] += deriv[d] * deriv[d];
                    correctedG[d] = G[d] / (1 - Math.Pow(beta, i));
                    lambda_ad[d] = lambda / (Math.Sqrt(correctedG[d] + eps));
                    
                    grad_acc[d] = grad_acc[d] * gamma + lambda_ad[d] * deriv[d];
                    CoeffNewIter[d] -= grad_acc[d]; // to update parameters
                }

                //Exit criteria calculation
                double exitCriteria = 0.0;
                for (int j = 0; j < CoeffNewIter.Length; ++j) {
                    exitCriteria += deriv[j] * deriv[j];
                }
                if (CoeffNewIter[0].IsNaNorInf() || CoeffNewIter[1].IsNaNorInf() || CoeffNewIter[2].IsNaNorInf())
                    throw new ArithmeticException($"breakdown in minimizer: {CoeffNewIter[0]}, {CoeffNewIter[1]}, {CoeffNewIter[2]}");


                Console.WriteLine("iter {0}: \t [{1}]", i, string.Join(", ", CoeffNewIter));
                Console.WriteLine("            \t [{1}] \t {2} ", i, string.Join(", ", deriv), deriv.L2Norm());
                Console.WriteLine("            \t[{1}]", i, string.Join(", ", correctedG));
                Console.WriteLine("            \t[{1}]", i, string.Join(", ", lambda_ad));

                //Exit loop condition
                if (Math.Sqrt(exitCriteria) < eps) {
                    break;
                }
                if (i == itermax) {
                    throw new ApplicationException("Increase the number of iterations for the gradient method!");
                }
               
            }
            stw.Stop();
            Console.WriteLine($"Elapsed Time minimization: {stw.Elapsed}");
            stw.Reset();

            Console.WriteLine($"Minimum of function: {F(CoeffNewIter, EllipsePar0, dt, time, meanVelocity)} at  X value: [{string.Join(", ", CoeffNewIter)}]");
            return CoeffNewIter;

        }

        // Minimi_Newton is not finished yet
        public static double[] Minimi_Newton(double[] InitialValueForCoeff, double[] EllipsePar0, Func<double[], double[], double, double, SinglePhaseField[], double> F, Func<double[], double[], double, double, SinglePhaseField[], MultidimensionalArray> CalculateHessian, Func<double[], double[], double, double, SinglePhaseField[], double[]> Differentiate, double dt, double time, SinglePhaseField[] meanVelocity) {

            double[] CoeffOldIter = (double[])InitialValueForCoeff.Clone();
            double[] CoeffNewIter = (double[])InitialValueForCoeff.Clone();
            double eps = 1e-12;
            int L = InitialValueForCoeff.Length;

            double[] deriv = Differentiate(CoeffOldIter, EllipsePar0, dt, time, meanVelocity);
            MultidimensionalArray Hess = CalculateHessian(CoeffOldIter, EllipsePar0, dt, time, meanVelocity);
            MultidimensionalArray Hess_Invert = MultidimensionalArray.Create(L, L);
            Hess.InvertTo(Hess_Invert);

            double[] Res_Multipl = new double[L];

            for (int d = 0; d < L; ++d) {
                Res_Multipl[d] = 0;
                for (int j = 0; j < L; j++) {
                    Res_Multipl[d] += Hess_Invert[d, j] * deriv[j];
                }
            }

            //BlockMsrMatrix.Multiply(res, Hess, derivOldIter);
            int itermax = 150000;
            for (int i = 1; i < itermax; ++i) {

                double exitCriteria = 0.0;

                for (int j = 0; j < L; ++j) {
                    CoeffNewIter[j] = CoeffOldIter[j] - Res_Multipl[j];
                    exitCriteria += (CoeffNewIter[j] - CoeffOldIter[j]).Pow2();
                }

                Console.WriteLine("iter {0}: [{1}]", i, string.Join(", ", CoeffNewIter));
                CoeffOldIter = CoeffNewIter;
                deriv = Differentiate(CoeffOldIter, EllipsePar0, dt, time, meanVelocity);
                Hess = CalculateHessian(CoeffOldIter, EllipsePar0, dt, time, meanVelocity);

                //Exit loop condition
                if (exitCriteria.Sqrt() < eps) {
                    break;
                }
            }
            double funcNewIter = F(CoeffNewIter, EllipsePar0, dt, time, meanVelocity);
            Console.WriteLine($"Minimum of function: {funcNewIter} at  X value: [{string.Join(", ", CoeffNewIter)}]");
            return CoeffNewIter;

        }

    }


    /// <summary>
    /// calculation of the numerical derivative, which is used in the minimization method
    /// </summary>
    public static class NumericalGradient {

    static public Func<double[], double[], double, double, SinglePhaseField[], double[]> Differentiate(Func<double[], double[], double, double, SinglePhaseField[], double> F) {

        return delegate (double[] X, double[] EllPar0, double dt, double time, SinglePhaseField[] meanVelocity) {
            double gridStepSize = GenericBlas.MachineEps.Sqrt();  //1e-14; 
            int L = X.Length;
            double[] dF_dX = new double[L];


            double[] XLeft = (double[])X.Clone();
            double[] XRight = (double[])X.Clone();

            for (int i = 0; i < L; ++i) {
                if (i > 0) {
                    XLeft[i - 1] = X[i - 1];
                    XRight[i - 1] = X[i - 1];
                }

                XLeft[i] -= gridStepSize;
                XRight[i] += gridStepSize;

                dF_dX[i] = 0.5 * (F(XRight, EllPar0, dt, time, meanVelocity) - F(XLeft, EllPar0, dt, time, meanVelocity)) / gridStepSize;
            }

            return dF_dX;
 
        };
    }

      // CalculateHessian is not finished yet
    static public Func<double[], double[], double, double, SinglePhaseField[], MultidimensionalArray> CalculateHessian(Func<double[], double[], double, double, SinglePhaseField[], double> F) {

        return delegate (double[] X, double[] EllPar0, double dt, double time, SinglePhaseField[] meanVelocity) {
            double gridStepSize = 1e-12;
            int L = X.Length;
            MultidimensionalArray hessian = MultidimensionalArray.Create(L, L);
            double[] X_plus_h = new double[L];
            double[] X_minus_h = new double[L];

            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {

                    X_plus_h = (double[])X.Clone();
                    X_minus_h = (double[])X.Clone();

                    X_plus_h[i] += gridStepSize;
                    X_plus_h[j] += gridStepSize;
                    X_minus_h[i] -= gridStepSize;
                    X_minus_h[j] -= gridStepSize;

                    double f_plus_plus = F(X_plus_h, EllPar0, dt, time, meanVelocity);
                    X_plus_h[j] -= 2 * gridStepSize;  // adjust jth component only
                    double f_plus_minus = F(X_plus_h, EllPar0, dt, time, meanVelocity);
                    X_minus_h[j] += 2 * gridStepSize;  // revert jth component to minus and add 2*h
                    double f_minus_plus = F(X_minus_h, EllPar0, dt, time, meanVelocity);
                    X_minus_h[j] -= 2 * gridStepSize; // adjust jth component only
                    double f_minus_minus = F(X_minus_h, EllPar0, dt, time, meanVelocity);

                    if (i == j) {
                        hessian[i, j] = (f_plus_plus - 2 * F(X, EllPar0, dt, time, meanVelocity) + f_minus_minus) / gridStepSize.Pow2();
                    } else {
                        hessian[i, j] = (f_plus_plus - f_plus_minus - f_minus_plus + f_minus_minus) / (4 * gridStepSize.Pow2());
                    }

                }
            }

            return hessian;
        };
    }

    }
    //public static class IntegrationOverEllipse {

    //    static double stepNumber = 100000;

    //    public static double IntegralCalculation(double[] EllipseCoeff, double[] EllipsePar0, Func<double[], double[], double, double, double, double, double[], double> f, double alphaMin, double alphaMax, double dt, double time, double[] velocity) {

    //        Console.WriteLine($"IntegralCalculation {EllipseCoeff[0]},{EllipseCoeff[1]},{EllipseCoeff[2]} {alphaMin} {alphaMax} {stepNumber}");
    //        double a = EllipsePar0[0];
    //        double b = EllipsePar0[1];
    //        double yc = EllipsePar0[2];

    //        Func<double, double> x = (alpha) => a * Math.Cos(alpha);
    //        Func<double, double> y = (alpha) => b * Math.Sin(alpha) + yc;
    //        Func<double, double> DerivativeX = (alpha) => -a * Math.Sin(alpha);
    //        Func<double, double> DerivativeY = (alpha) => b * Math.Cos(alpha);


    //        double gramDeter(double alpha) {
    //            return Math.Sqrt(DerivativeX(alpha) * DerivativeX(alpha) + DerivativeY(alpha) * DerivativeY(alpha));
    //        }


    //        double funcForIntegration(double alpha) {
    //            double gramDeter_alpha = gramDeter(alpha);
    //            double q_alpha = f(EllipseCoeff, EllipsePar0, x(alpha), y(alpha), dt, time, velocity);

    //            if (double.IsNaN(gramDeter_alpha) || double.IsInfinity(gramDeter_alpha) || double.IsInfinity(q_alpha) || double.IsInfinity(q_alpha))
    //                throw new ArithmeticException($" {EllipseCoeff} {gramDeter_alpha}, {q_alpha}");
    //            return gramDeter_alpha * q_alpha;
    //        }

    //        double stepSize = 0.5 * (alphaMax - alphaMin) / stepNumber;
    //        double integralResult = (funcForIntegration(alphaMin) + funcForIntegration(alphaMax)) / 3;

    //        for (int i = 1; i <= stepNumber; i++) {

    //            integralResult += (4 * funcForIntegration(alphaMin + (2 * i - 1) * stepSize) + 2 * funcForIntegration(alphaMin + 2 * i * stepSize)) / 3;
    //        }

    //        //var integral = integralResult * stepSize;
    //        Console.WriteLine("integral_simpson = " + integralResult * stepSize);
    //        return integralResult * stepSize;
    //    }

    //}


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
        /// calculation of the function f obtained after substitution the level-set function in the level-set evolution equation
        /// </summary>
        static double f(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double time, double[] velocity) {

            double SquareRoot = (EllipsePar0[1].Pow2() * (1 - x.Pow2() / EllipsePar0[0].Pow2())).Sqrt();
            double numer = EllipsePar0[1] * EllipseCoeff[1] * EllipsePar0[0] * (EllipsePar0[0].Pow2() - x.Pow2()) + x.Pow2() * EllipsePar0[1].Pow2() * EllipseCoeff[0];
            double DerivFuncEllipse = -EllipseCoeff[2] + numer / (SquareRoot * EllipsePar0[0].Pow2() * EllipsePar0[0]);
            double Gradx = -EllipsePar0[1].Pow2() * x / (EllipsePar0[0].Pow2() * SquareRoot);
            double Grady = 1.0;
            double res = velocity[0] * Gradx + velocity[1] * Grady;
            return DerivFuncEllipse + res;

        }

        /// <summary>
        /// squaring of the function f
        /// </summary>
        static double fsqr(double[] EllipseCoeff, double[] EllipsePar0, double x, double y, double timeStep, double time, double[] velocity) {
            var fval = f(EllipseCoeff, EllipsePar0, x, y, timeStep, time, velocity);
            return fval * fval;
        }


        /// <summary>
        /// calculation of the numerical integral from the square of the function f over the interface
        /// </summary>
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
               
                //Console.WriteLine($"integral_bosss: {TotalIntegral}");
                // Console.WriteLine($"integral_diff: {TotalIntegral_RefVal - TotalIntegral}");
            }

            return TotalIntegral;


        }

    }
}

