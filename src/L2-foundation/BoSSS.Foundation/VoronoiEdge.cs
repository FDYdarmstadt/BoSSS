using BoSSS.Foundation.Grid.Voronoi;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System;

namespace BoSSS.Foundation.Grid.Voronoi
{
    static class VoronoiEdge
    {
        //Galiliean-invariant cosmological hydronamical simulations on a moving mesh, p17, equation 33
        public static double NormalVelocity(double[] posR, double[] velR, double[] posL, double[] velL, double[] x, Vector normal)
        {
            double meanVelocity = (velR[0] + velL[0]) / 2 * normal[0] + (velR[1] + velL[1]) / 2 * normal[1];
            meanVelocity -= RotationVelocity(posR, velR, posL, velL, x);
            return meanVelocity;
        }

        //Galiliean-invariant cosmological hydronamical simulations on a moving mesh, p17, equation 32
        //Maybe caching would improve performance
        static double RotationVelocity(double[] posR, double[] velR, double[] posL, double[] velL, double[] x)
        {
            //Eq 32
            Vector rR_minus_rL = new Vector(
                posR[0] - posL[0],
                posR[1] - posL[1]);
            Vector w_L_minus_w_R = new Vector(
                velL[0] - velR[0],
                velL[1] - velR[1]);
            Vector f_minus_r_l_plus_r_L = new Vector(
                x[0] - (posR[0] + posL[0]) / 2.0,
                x[1] - (posR[1] + posL[1]) / 2.0);
            double result = w_L_minus_w_R * f_minus_r_l_plus_r_L;
            result /= rR_minus_rL.L2Norm();
            return result;
        }

        public static void TestRotationVelocity()
        {
            double rotation = RotationVelocity(
                new double[] { 1, 0 },
                new double[] { 0, 1 },
                new double[] { -1, 0 },
                new double[] { 0, -1 },
                new double[] { 0, -1 });
            double expectedRotation = 1;
            //Assert.IsTrue(Math.Abs(rotation - expectedRotation) < 1e-12,
            //    "Speed of rotation is not correct.");
        }
    }
}
