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
using System.Linq;
using System.Text;

namespace ilPSP {
    public partial class MultidimensionalArray {


        /// <summary>
        /// Matrix/Matrix multiplication (works only on 2D arrays).
        /// </summary>
        /// <returns><paramref name="A"/>*<paramref name="B"/> </returns>
        static public MultidimensionalArray operator *(MultidimensionalArray A, MultidimensionalArray B) {
            MultidimensionalArray R = MultidimensionalArray.Create(A.NoOfRows, B.NoOfCols);
            R.Multiply(1.0, A, B, 0.0, "ij", "ik", "kj");
            return R;
        }


        /// <summary>
        /// addition of two arrays
        /// </summary>
        static public MultidimensionalArray operator +(MultidimensionalArray A, MultidimensionalArray B) {
            var R = MultidimensionalArray.Create(A.Lengths);
            R.Acc(1.0, A);
            R.Acc(1.0, B);
            return R;
        }

        /// <summary>
        /// subtraction of two arrays
        /// </summary>
        static public MultidimensionalArray operator -(MultidimensionalArray A, MultidimensionalArray B) {
            var R = MultidimensionalArray.Create(A.Lengths);
            R.Acc(1.0, A);
            R.Acc(-1.0, B);
            return R;
        }

    }
}
