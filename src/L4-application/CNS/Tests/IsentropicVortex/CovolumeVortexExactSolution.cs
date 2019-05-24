/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using System;

namespace CNS.Tests.IsentropicVortex {

    /// <summary>
    /// Exact solution for isentropic vortex in a covolume gas; see Mueller2014
    /// </summary>
    static class CovolumeVortexExactSolution {

        /// <summary>
        /// Computes the solution at (x, t).
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="u"></param>
        /// <param name="t"></param>
        /// <param name="c"></param>
        /// <param name="c1"></param>
        /// <returns></returns>
        public static StateVector GetSolution(double x, double y, double u, double t, CNSControl c, double c1) {
            double xBar = x - u * t;
            double r = Math.Sqrt(xBar * xBar + y * y);
            double theta = Math.Atan2(y, xBar);
            CovolumeGas eos = c.EquationOfState as CovolumeGas;
            double gamma = eos.HeatCapacityRatio;
            double covolume = eos.Covolume;

            double rho = GetDensity(r, eos, c1);
            double p = Math.Pow(rho / (1 - covolume * rho), gamma);
            double normU = r * Math.Exp(0.5 * (1.0 - r * r));
            Vector U = new Vector(
                u - Math.Sin(theta) * normU,
                Math.Cos(theta) * normU,
                0.0);

            return new StateVector(
                c.GetMaterial(),
                rho,
                rho * U,
                p * (1.0 - covolume * rho) / (gamma - 1.0) + 0.5 * rho * U * U);
        }

        /// <summary>
        /// Computes the density from a non-linear equation using a simple
        /// bisection method
        /// </summary>
        /// <param name="r"></param>
        /// <param name="eos"></param>
        /// <param name="c1"></param>
        /// <returns></returns>
        private static double GetDensity(double r, CovolumeGas eos, double c1) {
            double gamma = eos.HeatCapacityRatio;
            double covolume = eos.Covolume;

            double densityLeft = 0.0 + 1e-14;
            double densityRight = 1.0 / covolume - 1e-14;
            double fOld = double.PositiveInfinity;
            int maxIterations = 1000;

            double density = double.NaN;
            double f = double.NaN;

            int n = 0;
            while (true) {
                if (n >= maxIterations) {
                    throw new Exception("No convergence");
                }
                
                density = 0.5 * (densityRight + densityLeft);
                f = F(density, r, eos, c1);

                if (Math.Abs(fOld - f) < 1e-13) {
                    break;
                } else if (f > 0.0) {
                    densityRight = density;
                } else {
                    densityLeft = density;
                }

                fOld = f;
                n++;
            }

            return density;
        }

        /// <summary>
        /// The non-linear equation F=0 to be solved
        /// </summary>
        /// <param name="density"></param>
        /// <param name="r"></param>
        /// <param name="eos"></param>
        /// <param name="c1"></param>
        /// <returns></returns>
        private static double F(double density, double r, CovolumeGas eos, double c1) {
            double gamma = eos.HeatCapacityRatio;
            double covolume = eos.Covolume;

            return (gamma - covolume * density) / (gamma - 1.0) / density *
                    Math.Pow(density / (1.0 - covolume * density), gamma)
                - c1 + 0.5 * Math.Exp(1.0 - r * r);
        }
    }
}
