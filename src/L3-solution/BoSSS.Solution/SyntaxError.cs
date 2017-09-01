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

namespace BoSSS.Solution.Utils.Formula {
    /**
     * Thrown by {@link FormulaParser} when error occurs.
     * @author udav
     */
    public class SyntaxError : ParseException {
        /// <summary>  </summary>
        /// <param name="errorCode"></param>
        /// <param name="parsedString"></param>
        /// <param name="position"></param>
        public SyntaxError(int errorCode, String parsedString, int position) :
            base(errorCode, parsedString, position) {
        }

        /// <summary> </summary>
        /// <param name="errorCode"></param>
        /// <param name="expectedToken"></param>
        /// <param name="parsedString"></param>
        /// <param name="position"></param>
        public SyntaxError(int errorCode, String[] expectedToken, String parsedString, int position) :
            base(errorCode, parsedString, position) {
        }
    }
}
