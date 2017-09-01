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
     * Interface contains instructions generated
     * by formula parser for formula interpreter.
     * @author udav
     */
    enum Instruction {
        PUSH = 0,
        RSLV = 1,
        PLUS = 3,
        MINUS = 4,
        MULT = 5,
        DIV = 6,
        POW = 7,
        UNARY_PLUS = 8,
        UNARY_MINUS = 9,
        RND = 10,
        MOD = 11,
        FUNC = 12,

        SIN = 100,
        COS = 101,
        TAN = 102,
        ASIN = 103,
        ACOS = 104,
        ATAN = 105,
        EXP = 106,
        LOG = 107,
        LOG10 = 108,
        SQRT = 109,
        CEIL = 110,
        FLOOR = 111,
        ABS = 112,
        MAX = 113,
        MIN = 114,
        SINH = 115,
        COSH = 116,
        TANH = 117,
        JUMP = 118,
        ATAN2 = 119,
        SIGN = 120,
        DILOG = 121,
        JACOBI_SN = 122,
        JACOBI_CN = 123,
        JACOBI_DN = 124,
    }
}