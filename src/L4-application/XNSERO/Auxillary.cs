/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using ilPSP;
using System;

namespace BoSSS.Application.XNSERO_Solver {
    /// <summary>
    /// Helper class for the FSISolver. Contains additional methods for testing.
    /// </summary>
    [Serializable]
    public class Auxillary {
        internal void TestArithmeticException(double[] variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(Vector variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void TestArithmeticException(double variable, string variableName) {
            ThrowIsNaNException(variable, variableName);
            ThrowIsInfinityException(variable, variableName);
        }

        internal void ThrowIsNaNException(double[] variable, string variableName) {
            for (int i = 0; i < variable.Length; i++) {
                if (double.IsNaN(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
            }
        }

        internal void ThrowIsNaNException(Vector variable, string variableName) {
            for (int i = 0; i < variable.Dim; i++) {
                if (double.IsNaN(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
            }
        }

        internal void ThrowIsInfinityException(Vector variable, string variableName) {
            for (int i = 0; i < variable.Dim; i++) {
                if (double.IsInfinity(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
            }
        }

        internal void ThrowIsNaNException(double variable, string variableName) {
            if (double.IsNaN(variable))
                throw new ArithmeticException("Error during update of " + variableName + ", value is NaN.");
        }


        internal void ThrowIsInfinityException(double[] variable, string variableName) {
            for (int i = 0; i < variable.Length; i++) {
                if (double.IsInfinity(variable[i]))
                    throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
            }
        }

        internal void ThrowIsInfinityException(double variable, string variableName) {
            if (double.IsInfinity(variable))
                throw new ArithmeticException("Error during update of " + variableName + ", value is infinity.");
        }

        
    }
}
