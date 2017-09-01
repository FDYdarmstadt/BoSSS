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
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {
    
    
    class Utilities {
        [DllImport("HYPRE")]
        extern public static int HYPRE_GetError();

        /* Check if the given error flag contains the given error code */
        [DllImport("HYPRE")]
        extern public static int HYPRE_CheckError(int hypre_ierr, int hypre_error_code);

        /* Return the index of the argument (counting from 1) where
           argument error (HYPRE_ERROR_ARG) has occured */
        [DllImport("HYPRE")]
        public static extern int HYPRE_GetErrorArg();

        /* Describe the given error flag in the given string */
        [DllImport("HYPRE")]
        unsafe public static extern void HYPRE_DescribeError(int hypre_ierr, byte* descr);

        /* Clears the hypre error flag */
        [DllImport("HYPRE")]
        public static extern int HYPRE_ClearAllErrors();

        /* Clears the given error code from the hypre error flag */
        [DllImport("HYPRE")]
        public static extern int HYPRE_ClearError(int hypre_error_code);
    }
}
