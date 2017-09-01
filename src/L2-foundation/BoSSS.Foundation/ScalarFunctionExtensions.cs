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
using System.Diagnostics;
using System.Linq;
using System.Text;
using ilPSP;

namespace BoSSS.Foundation {
    /// <summary>
    /// extend! extend! extend!
    /// </summary>
    public static class ScalarFunctionExtensions {

        private class Map1 {
            internal Func<double, double> T;
            internal ScalarFunctionEx f;

            internal void Tf(int j0, int Len, NodeSet N, MultidimensionalArray result) {
                f(j0, Len, N, result);
                Debug.Assert(Len == result.GetLength(0));
                Debug.Assert(result.Dimension == 2);
                result.ApplyAll(this.T);
            }
        }

        /// <summary>
        /// creates a <see cref="ScalarFunctionEx"/> which represents the comosition
        /// \f$ T \circ f\f$ 
        /// </summary>
        /// <param name="f">original function \f$ f\f$ </param>
        /// <param name="T">transformation \f$ T\f$ </param>
        /// <returns></returns>
        static public ScalarFunctionEx Map(this ScalarFunctionEx f, Func<double, double> T) {
            return (new Map1() { T = T, f = f }).Tf;
        }

    }
}
