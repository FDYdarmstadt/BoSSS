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

using System;
using System.Collections.Generic;

namespace BoSSS.Foundation {

    /// <summary>
    /// Provides a Function to Calculate Quadrature Orders
    /// </summary>
    static public class QuadOrderFunc {

        /// <summary>
        /// Twice the maximum degree of any of the domain, parameter or codomain fields.
        /// </summary>
        /// <returns></returns>
        public static Func<int[], int[], int[], int> MaxDegTimesTwo() {
            return (DomDegs, ParamDegs, CoDomDegs) => 2 * Math.Max(Math.Max(MaxOr0(DomDegs), MaxOr0(ParamDegs)), MaxOr0(CoDomDegs));
        }

        /// <summary>
        /// A fixed quadrature order, which completely ignored the DG polynomial degrees of any DG field.
        /// </summary>
        /// <param name="order">fixed quadrature order.</param>
        /// <returns>A function which always returns the value <paramref name="order"/>.</returns>
        public static Func<int[], int[], int[], int> FixedOrder(int order) {
            return ((DomDegs, ParamDegs, CoDomDegs) => order);
        }

        /// <summary>
        /// Assumes a linear integrand, resp. flux.
        /// </summary>
        /// <returns>
        /// A function that computes an integration degree according to
        /// 1 + $codomainDegree + $parameterDegree
        /// and rounding the result to the next highest even number
        /// </returns>
        public static Func<int[], int[], int[], int> Linear() {
            return SumOfMaxDegrees(1, false);
        }

        /// <summary>
        /// Assumes a non-linear integrand where the highest term behaves like
        /// \f[h^{<paramref name="NonLinDeg"/>}\f]
        /// </summary>
        /// <param name="NonLinDeg"></param>
        /// <returns>
        /// A function that computes an integration degree according to
        /// <paramref name="NonLinDeg"/> + $codomainDegree + $parameterDegree
        /// and rounding the result to the next highest even number
        /// </returns>
        public static Func<int[], int[], int[], int> NonLinear(int NonLinDeg) {
            return SumOfMaxDegrees(NonLinDeg, false);
        }

        /// <summary>
        /// Just as <see cref="NonLinear(int)"/>, but ignoring the degree of
        /// the parameter fields.
        /// </summary>
        /// <param name="NonLinDeg"></param>
        /// <param name="RoundUp"></param>
        /// <returns></returns>
        public static Func<int[], int[], int[], int> NonLinearWithoutParameters(int NonLinDeg, bool RoundUp = false) {
            return ((DomDegs, ParamDegs, CoDomDegs) => (QuadOrderFunc.SumOfMaxDegreesWithRounding(DomDegs, new int[0], CoDomDegs, NonLinDeg, RoundUp)));
        }

        /// <summary>
        /// Version of <see cref="NonLinear(int)"/> where rounding can be
        /// controlled
        /// </summary>
        /// <param name="NonLinDeg"></param>
        /// <param name="RoundUp"></param>
        /// <returns></returns>
        public static Func<int[], int[], int[], int> SumOfMaxDegrees(int NonLinDeg = 1, bool RoundUp = true) {
            return ((DomDegs, ParamDegs, CoDomDegs) => (SumOfMaxDegreesWithRounding(DomDegs, ParamDegs, CoDomDegs, NonLinDeg, RoundUp)));
        }
        
        /// <summary>
        /// used by <see cref="SumOfMaxDegrees(int, bool)"/>.
        /// </summary>
        private static int SumOfMaxDegreesWithRounding(int[] DomainDegrees, int[] ParamDegrees, int[] CoDomainDegrees, int NonlinDegree, bool RoundUp) {
            int order = MaxOr0(DomainDegrees) * NonlinDegree + MaxOr0(ParamDegrees) + MaxOr0(CoDomainDegrees);
            if (RoundUp) {
                if ((order % 2) != 0) {
                    order++;
                }
            }
            return order;
        }

        /// <summary>
        /// Returns 0 for an empty list
        /// </summary>
        static int MaxOr0<T>(T l) where T : IEnumerable<int> {
            int r = 0;
            foreach (int i in l) {
                r = Math.Max(i, r);
            }
            return r;
        }
    }
}
