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

namespace ilPSP {

    /// <summary>
    /// Extension methods for <see cref="int"/> values
    /// </summary>
    public static class IntExtensions {

        /// <summary>
        /// Calculates the factorial of <paramref name="n"/>
        /// </summary>
        /// <param name="n">
        /// A non-negative number.
        /// </param>
        /// <returns>\f$ n!\f$ </returns>
        public static int Factorial(this int n) {
            if (n < 0) {
                throw new ArgumentException("Number must be non-negative");
            }

            int result = 1;
            while (n > 1) {
                result *= n;
                n--;
            }
            return result;
        }

        /// <summary>
        /// Calculates the binomial coefficient
        /// \f$ \binom{n}{k}\f$ , often also referred to as
        /// "<paramref name="n"/> choose <paramref name="k"/>"
        /// </summary>
        /// <param name="n">A non-negative number</param>
        /// <param name="k">A non-negative number</param>
        /// <returns>
        /// \f$ \binom{n}{k} = \frac{n!}{k!(n-k)!}\f$ 
        /// </returns>
        public static int Choose(this int n, int k) {
            if (n < 0) {
                throw new ArgumentException("Number must be non-negative", "n");
            }

            if (k < 0) {
                throw new ArgumentException("Number must be non-negative", "k");
            }

            if (k > n) {
                throw new ArgumentException("Argument must be bigger than n=" + n, "k");
            }

            if (k > n - k) {
                k = n - k;
            }

            int coefficient = 1;
            for (int i = 0; i < k; i++) {
                coefficient *= (n - i) / (i + 1);
            }
            return coefficient;
        }

        /// <summary>
        /// executes the <paramref name="LoopResult"/> for all values from 0
        /// (including) to <paramref name="L"/> (excluding);
        /// </summary>
        /// <returns>
        /// Array of length <paramref name="L"/>,
        /// the i-th entry is equal to <paramref name="LoopResult"/>(i).
        /// </returns>
        public static T[] ForLoop<T>(this int L, Func<int, T> LoopResult) {
            T[] ret = new T[L];
            for (int i = 0; i < L; i++) {
                ret[i] = LoopResult(i);
            }
            return ret;
        }
    }
}
