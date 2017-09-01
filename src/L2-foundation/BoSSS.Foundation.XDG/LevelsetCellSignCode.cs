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
    /// Encodes, for one cell, the sign of all four level sets;
    /// </summary>
    public struct LevelsetCellSignCode {

        /// <summary>
        /// 3adic representation of the <see cref="LevelsetSign"/> for all four level sets;
        /// </summary>
        internal int lsSig;

        /// <summary>
        /// extracts the level set sign for level set no. <paramref name="LevSetIdx"/>
        /// </summary>
        public LevelsetSign GetSign(int LevSetIdx) {
            int _s3 = lsSig / 27;
            int rem = _s3 * 3;
            if (LevSetIdx == 3)
                return ((LevelsetSign)(_s3 - 1));

            int _s2 = lsSig / 9 - rem;
            rem += _s2;
            rem *= 3;
            if (LevSetIdx == 2)
                return ((LevelsetSign)(_s2 - 1));

            int _s1 = lsSig / 3 - rem;
            rem += _s1;
            rem *= 3;
            if (LevSetIdx == 1)
                return ((LevelsetSign)(_s1 - 1));

            int _s0 = lsSig - rem;
            if (LevSetIdx == 0)
                return ((LevelsetSign)(_s0 - 1));

            throw new ArgumentOutOfRangeException();
        }

        /// <summary>
        /// manipulates the level set sign for level set no. <paramref name="LevSetIdx"/>
        /// </summary>
        public void SetSign(int LevSetIdx, LevelsetSign sign) {

            int oldSign = (int)GetSign(LevSetIdx) + 1;

            int bs = 1;
            for (int i = 0; i < LevSetIdx; i++)
                bs *= 3;

            int newSign = (int)sign + 1;

            lsSig += (newSign - oldSign) * bs;
        }


        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            if (obj is SpeciesId) {
                return (((LevelsetCellSignCode)obj).lsSig == this.lsSig);
            } else {
                return false;
            }
        }

        /// <summary>
        /// hasch code
        /// </summary>
        public override int GetHashCode() {
            return lsSig;
        }

        /// <summary>
        /// equality
        /// </summary>
        public static bool operator ==(LevelsetCellSignCode a, LevelsetCellSignCode b) {
            return (a.lsSig == b.lsSig);
        }

        /// <summary>
        /// inequality
        /// </summary>
        public static bool operator !=(LevelsetCellSignCode a, LevelsetCellSignCode b) {
            return (a.lsSig != b.lsSig);
        }

        /// <summary>
        /// true, if <paramref name="cd"/> is contained in this 
        /// </summary>
        public bool IsContained(LevelSetSignCode cd, int NoOfSignificant_LevSets) {
            for (int i = 0; i < NoOfSignificant_LevSets; i++) {
                LevelsetSign ls_i = GetSign(i);
                if (ls_i == LevelsetSign.Both)
                    continue;

                LevelsetSign sg = cd.GetSign(i) ? LevelsetSign.Positive : LevelsetSign.Negative;
                if (sg != ls_i)
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Extracts the level set cell sign code
        /// </summary>
        /// <param name="RegionCode"></param>
        /// <returns>
        /// the reduced region code for the "full" region code <paramref name="RegionCode"/>;
        /// this is a number between 0 (including) and 81 (excluding).
        /// </returns>
        static public LevelsetCellSignCode Extract(ushort RegionCode) {
            int r = 0;

            if ((RegionCode & 0xf000) > 0x8000)
                r += 2; // level set 4 purly positive
            if ((RegionCode & 0xf000) == 0x8000)
                r += 1; // level set 4 crosses the cell (positive and negative)
            r *= 3;

            if ((RegionCode & 0x0f00) > 0x0800)
                r += 2; // level set 3 purly positive
            if ((RegionCode & 0x0f00) == 0x0800)
                r += 1; // level set 3 crosses the cell (positive and negative)
            r *= 3;

            if ((RegionCode & 0x00f0) > 0x0080)
                r += 2; // level set 2 purly positive
            if ((RegionCode & 0x00f0) == 0x0080)
                r += 1; // level set 2 crosses the cell (positive and negative)
            r *= 3;

            if ((RegionCode & 0x000f) > 0x0008)
                r += 2; // level set 1 purly positive
            if ((RegionCode & 0x000f) == 0x0008)
                r += 1; // level set 1 crosses the cell (positive and negative)

            LevelsetCellSignCode ret;
            ret.lsSig = r;
            return ret;
        }
    }
}
