using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSERO_Solver {
    class ParticleSuperEllipsoidFlat : ParticleSuperEllipsoid{

        double m_Length;
        double m_Exponent;
        double m_Thickness;

        public ParticleSuperEllipsoidFlat(
            IMotion motion, 
            double length, 
            double thickness, 
            int superEllipsoidExponent, 
            double[] startPos, 
            double startAngl = 0, 
            double activeStress = 0, 
            double[] startTransVelocity = null, 
            double startRotVelocity = 0)
            : base(motion, length, thickness, superEllipsoidExponent, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            if (startPos.Length != 2)
                throw new ArgumentOutOfRangeException("Spatial dimension does not fit particle definition");

            m_Length = length;
            m_Thickness = thickness;
            m_Exponent = superEllipsoidExponent;
        }

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {
            double alpha = -(Motion.GetAngle(0));
            double r;
            r = - Math.Pow(
                    Math.Pow(((X[0] - Postion[0]) * Math.Cos(alpha) - (X[1] - Postion[1]) * Math.Sin(alpha)) / m_Length, m_Exponent)
                    + Math.Pow(((X[0] - Postion[0]) * Math.Sin(alpha) + (X[1] - Postion[1]) * Math.Cos(alpha)) / m_Thickness, m_Exponent) 
                , 1 / m_Exponent)
                + 1;
            r *= 0.2;
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
            return r;
        }
    }
}
