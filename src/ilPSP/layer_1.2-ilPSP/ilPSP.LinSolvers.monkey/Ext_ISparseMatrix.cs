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

namespace ilPSP.LinSolvers {
    
    /// <summary>
    /// extension methods for <see cref="ISparseMatrix"/>.
    /// </summary>
    static public class Ext_ISparseMatrix {

        /// <summary>
        /// Sparse Matrix-Vector product in parallel
        /// performs <paramref name="acc"/> = <paramref name="acc"/>*<paramref name="beta"/> + <paramref name="alpha"/>*this*<paramref name="a"/>;
        /// </summary>
        static public void SpMVpara<VectorType1, VectorType2>(this MsrMatrix M, double alpha, VectorType1 a, double beta, VectorType2 acc) 
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> 
        {
            Debug.Assert(M.RowPartitioning.IsMutable == false);
            Debug.Assert(M.ColPartition.IsMutable == false);

            if (M.RowPartitioning.MpiSize > 1) {
                // parallel

                monkey.CPU.RefMatrix rm = new monkey.CPU.RefMatrix(M);

                rm.SpMV(alpha, a, beta, acc);
                rm.Dispose();
            } else {
                M.SpMV(alpha, a, beta, acc);
            }
        }

    }
}
