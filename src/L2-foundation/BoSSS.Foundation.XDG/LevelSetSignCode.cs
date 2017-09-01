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

namespace BoSSS.Foundation.XDG {
    /// <summary>
    /// encodes the sign of up to four level sets in four bits,
    /// (at one specific point, while <see cref="LevelsetSign"/> and <see cref="LevelsetCellSignCode"/>
    /// define all possibilities in <em>one cell</em>.);
    /// </summary>
    public struct LevelSetSignCode {
        /// <summary>
        /// see <see cref="Val"/>;
        /// </summary>
        internal int val;

        /// <summary>
        /// Bit 0 corresponds with sign of level set 0, Bit 1 corresponds ...
        /// </summary>
        public int Val {
            get {
                return val;
            }
            set {
                if (val < 0 || val >= 16)
                    throw new ArgumentOutOfRangeException();
                val = value;
            }
        }

        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            if (obj is LevelSetSignCode) {
                return (((LevelSetSignCode)obj).val == this.val);
            } else {
                return false;
            }
        }

        /// <summary>
        /// hasch code
        /// </summary>
        public override int GetHashCode() {
            return val;
        }

        /// <summary>
        /// equality
        /// </summary>
        public static bool operator ==(LevelSetSignCode a, LevelSetSignCode b) {
            return (a.val == b.val);
        }

        /// <summary>
        /// inequality
        /// </summary>
        public static bool operator !=(LevelSetSignCode a, LevelSetSignCode b) {
            return (a.val != b.val);
        }



        /// <summary>
        /// Encodes the signs of the level set values into a byte code
        /// </summary>
        /// <param name="levelSetValues">
        /// The values of the level sets at a given point
        /// </param>
        /// <returns>
        /// A byte code where the n-th bit (starting from the least significant
        /// bit) is 0 if <paramref name="levelSetValues"/>[n] is less than zero
        /// and 1 otherwise.
        /// </returns>
        public static LevelSetSignCode ComputeLevelSetBytecode(params double[] levelSetValues) {
            if (levelSetValues.Length > 4) {
                throw new ArgumentException("Currently, a maximum of 4 level sets is supported", "levelSetValues");
            }

            int result = 0;
            if (levelSetValues.Length > 0 && levelSetValues[0] >= 0) {
                result += 1;
            }
            if (levelSetValues.Length > 1 && levelSetValues[1] >= 0) {
                result += 2;
            }
            if (levelSetValues.Length > 2 && levelSetValues[2] >= 0) {
                result += 4;
            }
            if (levelSetValues.Length > 3 && levelSetValues[3] >= 0) {
                result += 8;
            }

            LevelSetSignCode res;
            res.val = result;
            return res;
        }


        /// <summary>
        /// extracts the sign of one level set
        /// </summary>
        /// <param name="LevSetIdx"></param>
        /// <returns>
        /// true: positive sign
        /// false: negative sign
        /// </returns>
        public bool GetSign(int LevSetIdx) {
            if (LevSetIdx < 0 || LevSetIdx >= 4)
                throw new IndexOutOfRangeException();
            return (val & (0x1 << LevSetIdx)) != 0;
        }

        /// <summary>
        /// manipulates the sign of one level set;
        /// </summary>
        /// <param name="LevSetIdx"></param>
        /// <param name="sign">
        /// true: positive sign
        /// false: negative sign
        /// </param>
        public void SetSign(int LevSetIdx, bool sign) {
            if (LevSetIdx < 0 || LevSetIdx >= 4)
                throw new IndexOutOfRangeException();

            int mask = (0x1 << LevSetIdx);
            int NotMask = ~mask;
            mask = sign ? mask : 0;

            val = (val & NotMask) | mask;
        }
    }
}
