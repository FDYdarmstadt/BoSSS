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
    public class ParameterizedLevelSetControl : ILevSetControl {

        /// <summary>
        /// 
        /// </summary>
        public double  xSemiAxis;
        /// <summary>
        /// 
        /// </summary>
        public double ySemiAxis;

        /// <summary>
        /// 
        /// </summary>
        public double yCenter;

        /// <summary>
        /// See <see cref="Parameterized_Timestepper"/>
        /// </summary>
        public Parameterized_Timestepper Timestepper = Parameterized_Timestepper.ExplicitEuler;

        /// <summary>
        /// 
        /// </summary>
        public ParameterizedLevelSetControl(double xSemiAxis, double ySemiAxis, double yCenter) {
            this.xSemiAxis = xSemiAxis;
            this.ySemiAxis = ySemiAxis;
            this.yCenter = yCenter;

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
