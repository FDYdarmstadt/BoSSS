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
     * Thrown when {@link Tokenizer} or {@link FormulaParser} occurs error.
     * @author udav
     */
    public class ParseException : Exception {
        /// <summary> </summary>
        public const int SUCCESSFUL = 0; // successful
        /// <summary> </summary>
        public const int UNEXPECTED_TOKEN = 1;
        /// <summary> </summary>
        public const int EXPECTED_IDENTIFIER = 2;
        /// <summary> </summary>
        public const int EXPECTED_NUMBER = 4;
        /// <summary> </summary>
        public const int EXPECTED_OPERATOR = 8;
        /// <summary> </summary>
        public const int EXPECTED_SEPARATOR = 16;
        /// <summary> </summary>
        public const int EXPECTED_TOKEN = 32;
        /// <summary> </summary>
        public const int EXPECTED_NOTHING = 64;

        private int errorCode;
        private String[] expectedTokens;
        private String parsedString;
        private int position;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="errorCode"></param>
        /// <param name="parsedString"></param>
        /// <param name="position"></param>
        public ParseException(int errorCode, String parsedString, int position) :
            this(errorCode, new String[0], parsedString, position) {
            
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="errorCode"></param>
        /// <param name="expectedTokens"></param>
        /// <param name="parsedString"></param>
        /// <param name="position"></param>
        public ParseException(int errorCode, String[] expectedTokens, String parsedString, int position) {
            this.errorCode = errorCode;
            this.expectedTokens = expectedTokens;
            this.parsedString = parsedString;
            this.position = position;
        }

        /// <summary> </summary>
        /// <returns>Error code</returns>
        public int getErrorCode() {
            return errorCode;
        }

        /**
         * @return String to be parsed
         */

        public String getParsedString() {
            return parsedString;
        }

        /**
         * @return Position in string to be parsed, in which error occurs.
         */
        public int getPosition() {
            return position;
        }

        /**
         * @return Copy of expected tokens list.
         */
        public String[] getExpectedTokens() {
            return (string[]) expectedTokens.Clone();
        }
    }
}
