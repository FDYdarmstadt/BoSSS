using System;
using System.Collections.Generic;
using System.Text;

namespace ilPSP.Utils {


    /// <summary>
    /// Spherical Harmonic Functions and Legendre Functions, see
    /// https://en.wikipedia.org/wiki/Spherical_harmonics,
    /// explicit representations taken from https://de.wikipedia.org/wiki/Kugelfl%C3%A4chenfunktionen#Darstellung
    /// </summary>
    public static class SphericalHarmonics {
        public static double factorial(int m) {
            if(m < 0)
                throw new ArgumentException();
            double r = 1.0;
            for(int i = 1; i <= m; i++)
                r *= i;
            return r;
        }

        /// <summary>
        /// some base to an integer power
        /// </summary>
        public static double Pow(double a, int m) {
            if(m < 0)
                throw new ArgumentException();
            double r = 1.0;
            for(int i = 0; i < m; i++)
                r *= a;
            return r;
        }

        /// <summary>
        /// Polar coordinates (e.g. input for <see cref="MyRealSpherical"/>) from Cartesian 
        /// </summary>
        /// <param name="X"></param>
        /// <returns>
        /// - theta: inclination angle, from 0 to $` \pi `$
        /// - phi: azimuth, from 0 to $` 2 \pi `$
        /// </returns>
        public static (double theta, double phi) GetAngular(Vector X) {
            Vector X0 = X.Normalize();

            double theta = Math.Acos(X0.z);
            double phi = Math.Atan2(X0.y, X0.x);

            if(phi < 0)
                phi += 2 * Math.PI;

            if(phi < 0 || phi >= Math.PI * 2 || phi.IsNaNorInf())
                throw new ApplicationException($"phi computation failed: {phi} for {X0} -- {X}");
            if(theta < 0 || theta > Math.PI || theta.IsNaNorInf())
                throw new ApplicationException($"theta computation failed: {theta} for {X0} -- {X}");

            Vector X0_recovered = new Vector(
                Math.Cos(phi) * Math.Sin(theta),
                Math.Sin(phi) * Math.Sin(theta),
                Math.Cos(theta));

            double errDist = X0.Dist(X0_recovered);
            if(errDist >= 1.0e-6)
                throw new ApplicationException($"angular coordinate computation failed: {X} -> {X0} -> ({phi},{theta}) -> {X0_recovered}");

            return (theta, phi);
        }


        static SphericalHarmonics() {
            int l_max = 10;
            Legendre_m = new double[l_max];
            for(int m = 0; m < l_max; m++) {
                Legendre_m[m] = Pow(-1, m) * factorial(2 * m) / (Pow(2, m) * factorial(m));
            }

            Spherical_Nlm = (double[,]) Array.CreateInstance(typeof(double), new int[] { l_max, 2 * l_max + 2 }, new int[] { 0, -l_max });
            for(int l = 0; l < l_max; l++) {
                for(int m = -l; m <= l; m++) {
                    Spherical_Nlm[l,m] = (1 / Math.Sqrt(2 * Math.PI)) * Math.Sqrt(((2.0 * l + 1.0) / 2.0) * factorial(l - m) / factorial(l + m)); 
                }
            }
        }

        static double[] Legendre_m;
       

        /// <summary>
        /// Legendre functions
        /// </summary>
        /// <param name="l"></param>
        /// <param name="m">
        /// For <paramref name="m"/>=0, one obtains Legendre Polynomials
        /// </param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double MyLegendre(int l, int m, double x) {
            double f, a;

            if(l < m) {
                //f = 0;
                return 0.0;
            } else if(m < 0) {
                a = MyLegendre(l, -m,  x);
                f = a / ((Pow(-1, m)) * factorial(l - m) / factorial(l + m));
                return f;
            } else if(l == m) {
                // compute (1 - x^2)^(m/2)
                if(m % 2 == 0)
                    a = Pow(1 - x * x, m / 2);
                else
                    a = Pow(Math.Sqrt(1 - x * x), m);

                f = Legendre_m[m] * a;  //Pow(-1, m) * factorial(2 * m) / (Pow(2,m) * factorial(m)) * a;
                return f;
            } else {
                f = x * (2 * l - 1) * MyLegendre(l - 1, m, x) - (l + m - 1) * MyLegendre(l - 2, m, x);
                f = f / (l - m);
                return f;
            }
        }

        static double[,] Spherical_Nlm;

        /// <summary>
        /// 
        /// </summary>
        static public double Get_Nlm(int l, int m) {
            return Spherical_Nlm[l, Math.Abs(m)];
        }


        /// <summary>
        /// Real value (real as in complex vs. real numbers) spherical harmonics.
        /// </summary>
        /// <param name="l"></param>
        /// <param name="m"></param>
        /// <param name="theta">
        /// <see cref="GetAngular"/>
        /// </param>
        /// <param name="phi">
        /// <see cref="GetAngular"/>
        /// </param>
        /// <returns></returns>
        public static double MyRealSpherical(int l, int m, double theta, double phi) {
            if(theta < 0 || theta > Math.PI) {
                throw new ArgumentException();
            }
            if(phi < 0 || phi > 2*Math.PI) {
                throw new ArgumentException();
            }

            // avoid the scaling to conform with the cooperation project with group Prof. Brenn (Oscillating Droplet)
            double Nlm = 1.0; 
            //Nlm = Spherical_Nlm[l, Math.Abs(m)]; // Math.Sqrt((2 * l + 1) / 2 * factorial(l - m) / factorial(l + m));

            double Ylm;
            if(m >= 0) {
                Ylm =  Nlm * MyLegendre(l, m, Math.Cos(theta)) * Math.Cos(m * phi);
            } else {
                Ylm =  Nlm * MyLegendre(l, -m, Math.Cos(theta)) * Math.Sin(-m * phi);
            }
            return Ylm;
        }



    }
}
