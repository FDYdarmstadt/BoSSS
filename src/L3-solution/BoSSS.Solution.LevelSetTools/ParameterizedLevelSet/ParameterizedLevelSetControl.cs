using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Solution.Timestepping;

namespace BoSSS.Solution.LevelSetTools.ParameterizedLevelSet {

    /// <summary>
    /// Encapsulation of Options for FourierLevSet
    /// </summary>
    [Serializable]
    abstract public class ParameterizedLevelSetControl : ILevSetControl {

        public abstract double[] CurveParametes { get; }

        /// <summary>
        /// See <see cref="Parameterized_Timestepper"/>
        /// </summary>
        public Parameterized_Timestepper Timestepper = Parameterized_Timestepper.ExplicitEuler;
        
        abstract public double PhiFunc(double[] X);
    }

    [Serializable]
    public class ParameterizedLevelSetControlEllipse : ParameterizedLevelSetControl {

        /// <summary>
        /// x-semiAxis of elliptic interface
        /// </summary>
        public double xSemiAxis;
        /// <summary>
        /// y-semiAxis of elliptic interface
        /// </summary>
        public double ySemiAxis;

        /// <summary>
        /// y-coordinate of center of ellipse
        /// </summary>
        public double yCenter;


        public ParameterizedLevelSetControlEllipse(double xSemiAxis, double ySemiAxis, double yCenter) {
            this.xSemiAxis = xSemiAxis;
            this.ySemiAxis = ySemiAxis;
            this.yCenter = yCenter;

        }

        public override double[] CurveParametes => new double[] { xSemiAxis, ySemiAxis, yCenter };

        public override double PhiFunc(double[] X) {
            return X[1] - yCenter + ySemiAxis * (1 - X[0].Pow2() / xSemiAxis.Pow2()).Sqrt();
        }
    }


    [Serializable]
    public class ParameterizedLevelSetControPolynomial : ParameterizedLevelSetControl {

        


        public ParameterizedLevelSetControPolynomial(double xSemiAxis, double ySemiAxis, double yCenter) {


        }

        public override double[] CurveParametes => throw new NotImplementedException();

        public override double PhiFunc(double[] X) {
            throw new NotImplementedException();
        }
    }





    public static class ParameterizedLevelSetFactory {

        public static ParameterizedLevelSetTimeStepper Build_Timestepper(ParameterizedLevelSetControl Control) {

            switch (Control.Timestepper) {
                case Parameterized_Timestepper.ExplicitEuler:
                    break;
                default:
                    throw new ArgumentException();
            }
            return new ParameterizedLevelSetTimeStepper(Control);

        }

    }

}
