using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.LevelSetTools.TestCases {
    public class Popcorn {

        double radius;
        double sigma;
        double A;

        public Popcorn(double radius, double sigma, double A = 2.0) {
            this.radius = radius;
            this.sigma = sigma;
            this.A = A;
        }

        /// <summary>
        /// Signed-Distance Level-Set Function for a Torus
        /// </summary>
        /// <param name="X">Cartesian Coordinates</param>
        /// <returns>Level-Set Value at X</returns>
        public double SignedDistance(double[] X) {
            double SumOfSeries = 0;
            for (int i = 0; i <= 11; i++) {
                double[] xk = Xk(i);
                SumOfSeries += A * Math.Exp(-((X[0] - xk[0]).Pow2() + (X[1] - xk[1]).Pow2() + (X[2] - xk[2]).Pow2()) / sigma * sigma);
            }
            return Math.Sqrt(X[0].Pow2() + X[1].Pow2() + X[2].Pow2() - radius - SumOfSeries);
        }


        ///// <summary>
        ///// Signed-Distance 2D Level-Set Function for a Torus
        ///// </summary>
        ///// <param name="X">Cartesian Coordinates</param>
        ///// <returns>Level-Set Value at X</returns>
        //public double SignedDistance2D(double[] X) {
        //    return Math.Sqrt((Math.Sqrt(X[0].Pow2() + X[1].Pow2()) - bigR).Pow2()) - smallR;
        //}

        /// <summary>
        /// Get pre-determined values for sum series
        /// </summary>
        /// <param name="X">Cartesian Coordinates</param>
        /// <returns>Level-Set Value at X</returns>
        public double[] Xk(int k) {
            double[] Xk = new double[3];
            if (k <= 4) {
                Xk[0] = 2 * Math.Cos(2 * k * Math.PI / 5) * radius / Math.Sqrt(5);
                Xk[1] = 2 * Math.Sin(2 * k * Math.PI / 5) * radius / Math.Sqrt(5);
                Xk[2] = 1 * radius / Math.Sqrt(5);
            } else if (k <= 9) {
                Xk[0] = 2 * Math.Cos((2 * (k - 5) - 1) * Math.PI / 5) * radius / Math.Sqrt(5);
                Xk[1] = 2 * Math.Sin((2 * (k - 5) - 1) * Math.PI / 5) * radius / Math.Sqrt(5);
                Xk[2] = -1 * radius / Math.Sqrt(5);
            } else if (k == 10) {
                Xk[2] = 1 * radius;
            } else if (k == 11) {
                Xk[2] = -1 * radius;
            } else {
                throw new NotImplementedException();
            }
            return Xk;

        }
    }
}
