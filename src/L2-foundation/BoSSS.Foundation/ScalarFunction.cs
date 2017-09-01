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
    /// delegate for scalar functions in D-dimensional space;
    /// vectorized definition;
    /// </summary>
    /// <param name="input">positions in space at which the function should be evaluated;
    /// 1st index: point index;
    /// 2nd index: spatial coordinate vector (from 0 to D-1);
    /// </param>
    /// <param name="output">result of function evaluation;
    /// 1st index: point index, corresponds with 1st index of <paramref name="input"/>
    /// </param>
    public delegate void ScalarFunction(MultidimensionalArray input, MultidimensionalArray output);
}
