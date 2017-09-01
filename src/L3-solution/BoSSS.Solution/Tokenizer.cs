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
using System.Text;

namespace BoSSS.Solution.Utils.Formula {

    /**
     * Tokenizer is a lexical analizer.
     * It breakes input string into tokens. If tokenizer occurs error,
     * it throws {@link LexicalError} exception. Tokenizer can read
     * identifier, separator, operator or number.
     * @author udav
     */
    class Tokenizer {
        public const int UNKNOWN = -1;
        public const int IDENTIFIER = 0;
        public const int NUMBER = 1;  // number literal
        public const int OPERATOR = 2;
        public const int SEPARATOR = 3;

        private const String operators = "+-*/^=%";
        private const String separator = ",()";
        private /*final*/ char[] _string;

        /** Buffer for token store */
        private StringBuilder buffer  = new StringBuilder(32000);
        /** Position at string */
        private int curPos = 0;
        /** Type of current token */
        private int tokenType = UNKNOWN;

        /***
         * Creates new tokenizer for parsing of the specified string.
         * @param string string for parsing
         */
        public Tokenizer(string __string) {
            this._string = __string.ToCharArray();
        }

        /**
         * @return Parsed string
         */
        public String getParsedString() {
            return new string(_string);
        }

        /***
         * @return Current position
         */
        public int getCurrentPosition() {
            return curPos;
        }

        /**
         * @return Returns true if the parsed string contains more tokens; otherwise returns false.
         */
        public bool hasNext() {
            whiteSpaces();
            return curPos < _string.Length;
        }

        /**
         * @return Current token type
         */
        public int getTokenType() {
            return tokenType;
        }

        /**
         * Reads next token. The token's type can be taken from method {@link #getTokenType()}.
         * If error occurs or no more tokens, throws {@link LexicalError}.
         * @return token
         * @throws LexicalError throws wher error occurs or no more tokens.
         */
        public string token() {
            whiteSpaces();
            buffer.Length = 0;

            if (operators.IndexOf(_string[curPos]) >= 0) {
                return _operator();
            } else if (separator.IndexOf(_string[curPos]) >= 0) {
                return _separator();
            } else if (char.IsLetter(_string[curPos])) {
                return identifier();
            } else if (_string[curPos] == '.' || char.IsDigit(_string[curPos])) {
                return number();
            } else {
                throw new LexicalError(LexicalError.UNEXPECTED_TOKEN, new String(_string), curPos);
            }
        }

        /**
         * Reads separator. Method {@link #getTokenType()} must returns type {@link #SEPARATOR}.
         * If error occurs or no more tokens, throws {@link LexicalError}.
         * @return separator
         * @throws LexicalError throws wher error occurs or no more tokens.
         */
        public String _separator() {
            whiteSpaces();
            buffer.Length = 0;
            if (separator.IndexOf(_string[curPos]) < 0) {
                throw new LexicalError(LexicalError.EXPECTED_SEPARATOR, new String(_string), curPos);
            }
            buffer.Append(_string[curPos]);
            ++curPos;
            tokenType = SEPARATOR;
            return buffer.ToString();
        }

        /**
         * Reads operator. Method {@link #getTokenType()} must returns type {@link #OPERATOR}.
         * If error occurs or no more tokens, throws {@link LexicalError}.
         * @return operator
         * @throws LexicalError throws wher error occurs or no more tokens.
         */
        public String _operator() {
            whiteSpaces();
            buffer.Length = 0;
            if (operators.IndexOf(_string[curPos]) < 0) {
                throw new LexicalError(LexicalError.EXPECTED_OPERATOR, new String(_string), curPos);
            }
            buffer.Append(_string[curPos]);
            ++curPos;
            tokenType = OPERATOR;
            return buffer.ToString();
        }

        /**
         * Reads identifier. Identifier is a char sequence, starting with symbol which follows
         * symbols or numbers. Method {@link #getTokenType()} must returns type {@link #IDENTIFIER}.
         * If error occurs or no more tokens, throws {@link LexicalError}.
         * @return identifier
         * @throws LexicalError throws wher error occurs or no more tokens.
         */
        public String identifier()  {
		    buffer.Length = 0;
		    whiteSpaces();
		    if (!char.IsLetter(_string[curPos])) {
			    throw new LexicalError(LexicalError.EXPECTED_IDENTIFIER, new String(_string), curPos);
		    }
		    buffer.Append(_string[curPos]);
		    ++curPos;
		    while (curPos < _string.Length &&
				    char.IsLetterOrDigit(_string[curPos])) {
			    buffer.Append(_string[curPos]);
			    ++curPos;
		    }
		    if (curPos < _string.Length && char.IsLetterOrDigit(_string[curPos])) {
			    throw new LexicalError(LexicalError.EXPECTED_IDENTIFIER, new String(_string), curPos);
		    }
		    tokenType = IDENTIFIER;
		    return buffer.ToString();
	    }

        /**
         * Reads number. Method {@link #getTokenType()} must returns type {@link #NUMBER}.
         * If error occurs or no more tokens, throws {@link LexicalError}.
         * @return number
         * @throws LexicalError throws wher error occurs or no more tokens.
         */
        public String number() {
            buffer.Length = 0;
            whiteSpaces();
            if (_string[curPos] != '.' && !char.IsDigit(_string[curPos])) {
                throw new LexicalError(LexicalError.EXPECTED_NUMBER, new String(_string), curPos);
            }
            if (_string[curPos] == '.') {
                buffer.Append('.');
                ++curPos;
                intNumber();
            } else {
                intNumber();
                if (curPos < _string.Length && _string[curPos] == '.') {
                    buffer.Append('.');
                    ++curPos;
                    if (curPos < _string.Length && char.IsDigit(_string[curPos])) {
                        intNumber();
                    }
                }
            }
            tokenType = NUMBER;
            return buffer.ToString();
        }

        /**
         * Eat white spaces
         */
        private void whiteSpaces() {
            while (curPos < _string.Length &&
                    char.IsWhiteSpace(_string[curPos])) {
                ++curPos;
            }
        }

        private void intNumber() {
            if (!char.IsDigit(_string[curPos])) {
                throw new LexicalError(LexicalError.EXPECTED_NUMBER, new String(_string), curPos);
            }
            buffer.Append(_string[curPos]);
            ++curPos;
            while (curPos < _string.Length && char.IsDigit(_string[curPos])) {
                buffer.Append(_string[curPos]);
                ++curPos;
            }

            if (curPos < _string.Length && char.IsLetter(_string[curPos])) {
                throw new LexicalError(LexicalError.EXPECTED_NUMBER, new String(_string), curPos);
            }
        }
    }
}