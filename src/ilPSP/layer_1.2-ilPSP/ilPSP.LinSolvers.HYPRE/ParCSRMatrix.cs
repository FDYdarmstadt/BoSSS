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
using System.Text;
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {

    /// <summary>
    /// 
    /// </summary>
    internal class ParCSRMatrix {

        /// <summary>
        /// y = alpha*Matrix*x + beta*y
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRMatrixMatvec")]
        public static extern int ParCSRMatrixMatvec(double alpha, T_ParCSR_matrix ParCSRMatrix, T_ParCRS_vector ParVector_x, double beta, T_ParCRS_vector ParVector_y);
    }
}
