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
using ilPSP;

namespace BoSSS.Foundation {

    /// <summary>
    /// delegate for steady-state vector functions in D-dimensional space,
    /// vectorized definition.
    /// </summary>
    /// <param name="input">positions in space at which the function should be evaluated;
    /// - 1st index: point index;
    /// - 2nd index: spatial coordinate vector (from 0 to D-1);
    /// </param>
    /// <param name="output">result of function evaluation;
    /// - 1st index: point index, corresponds with 1st index of <paramref name="input"/>
    /// - 2nd index: spatial direction, corresponds with 2nd index of <paramref name="input"/>
    /// </param>
    public delegate void VectorFunction(MultidimensionalArray input, MultidimensionalArray output);

    /// <summary>
    /// delegate for time-dependent vector functions in D-dimensional space,
    /// vectorized definition.
    /// </summary>
    /// <param name="input">positions in space at which the function should be evaluated;
    /// - 1st index: point index;
    /// - 2nd index: spatial coordinate vector (from 0 to D-1);
    /// </param>
    /// <param name="output">result of function evaluation;
    /// - 1st index: point index, corresponds with 1st index of <paramref name="input"/>
    /// - 2nd index: spatial direction, corresponds with 2nd index of <paramref name="input"/>
    /// </param>
    /// <param name="time">
    /// time
    /// </param>
    public delegate void VectorFunctionTimeDep(MultidimensionalArray input, double time, MultidimensionalArray output);

    /// <summary>
    /// Extensions functions for <see cref="ScalarFunction"/>
    /// </summary>
    public static class VectorFunctionExt {

        /// <summary>
        /// fixing the time for some time-depended function <paramref name="f"/>.
        /// </summary>
        public static VectorFunction SetTime(this VectorFunctionTimeDep f, double t) {
            return delegate (MultidimensionalArray input, MultidimensionalArray output) {
                f(input, t, output);
            };
        }

    }
}
