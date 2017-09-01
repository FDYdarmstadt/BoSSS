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

namespace ilPSP.LinSolvers.HYPRE.Wrappers {

    /// <summary>
    /// delegate for setup function of HYPRE solvers and preconditioners
    /// </summary>
    delegate int Setup(T_Solver solver, T_ParCSR_matrix Matrix_A, T_ParCRS_vector Vector_b, T_ParCRS_vector Vector_x);

    /// <summary>
    /// delegate for solver function of HYPRE solvers and preconditioners
    /// </summary>
    delegate int Solve(T_Solver solver, T_ParCSR_matrix Matrix_A, T_ParCRS_vector Vector_b, T_ParCRS_vector Vector_x);
    
    /// <summary>
    /// pointer to solver object
    /// </summary>
    struct T_Solver {
        public IntPtr p;
    }

    /// <summary>
    /// Pointer to Matrix object
    /// </summary>
    struct T_IJMatrix {
        public IntPtr p;
    }

    /// <summary>
    /// Pointer to vector object
    /// </summary>
    struct T_IJVector {
        public IntPtr p;
    }

    /// <summary>
    /// assembled matrix
    /// </summary>
    struct T_ParCSR_matrix {
        public IntPtr p;
    }

    /// <summary>
    /// assembled vector
    /// </summary>
    struct T_ParCRS_vector {
        public IntPtr p;
    }


}
