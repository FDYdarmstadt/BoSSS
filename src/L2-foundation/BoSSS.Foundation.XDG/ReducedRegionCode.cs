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
    /// Encodes the cut-cell - state of a cell, for all four level sets. 
    /// </summary>
    public struct ReducedRegionCode {

        /// <summary>
        /// 3adic representation of the <em>reduced region</em> for all four level sets;
        /// </summary>
        internal int rrc;

        /// <summary>
        /// For one level set, the <em>reduced region</em> is one of the numbers 0,1 or 2:
        /// <list type="bullet">
        ///   <item><b>2:</b> the cell if in the positive FAR region: <i>distance</i> == 7, i.e. the <i>code</i> is 0xf</item>
        ///   <item><b>1:</b>the cell if in the negative FAR region: <i>distance</i> == -1, i.e the <i>code</i> is 0x1</item>
        ///   <item><b>0:</b>the cell is cuttet or in the near reagion: -6 &lt; <i>distance</i> &lt; 6</item>
        /// </list>
        /// The reduced region code is build from the reduced regions of each level set 
        /// by encoding them into a 3-adic representation.
        /// </summary>
        public int _Val {
            get {
                return rrc;
            }
            set {
                if (value < 0 || value >= 81)
                    throw new ArgumentOutOfRangeException();
                rrc = value;
            }
        }

        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            if (obj is SpeciesId) {
                return (((ReducedRegionCode)obj).rrc == this.rrc);
            } else {
                return false;
            }
        }

        /// <summary>
        /// hasch code
        /// </summary>
        public override int GetHashCode() {
            return rrc;
        }

        /// <summary>
        /// equality
        /// </summary>
        public static bool operator ==(ReducedRegionCode a, ReducedRegionCode b) {
            return (a.rrc == b.rrc);
        }

        /// <summary>
        /// inequality
        /// </summary>
        public static bool operator !=(ReducedRegionCode a, ReducedRegionCode b) {
            return (a.rrc != b.rrc);
        }

        /// <summary>
        /// Extracts the reduced region code.
        /// </summary>
        /// <param name="RegionCode"></param>
        /// <returns>
        /// the reduced region code for the "full" region code <paramref name="RegionCode"/>;
        /// this is a number between 0 (including) and 81 (excluding).
        /// </returns>
        static public ReducedRegionCode Extract(ushort RegionCode) {
            int r = 0;

            if ((RegionCode & 0xf000) == 0xf000)
                r += 2; // code for level set 4 is +FAR
            if ((RegionCode & 0xf000) == 0x1000)
                r += 1; // code for level set 4 is -FAR
            r *= 3;

            if ((RegionCode & 0x0f00) == 0x0f00)
                r += 2; // code for level set 3 is +FAR
            if ((RegionCode & 0x0f00) == 0x0100)
                r += 1; // code for level set 3 is -FAR
            r *= 3;

            if ((RegionCode & 0x00f0) == 0x00f0)
                r += 2; // code for level set 2 is +FAR
            if ((RegionCode & 0x00f0) == 0x0010)
                r += 1; // code for level set 2 is -FAR
            r *= 3;

            if ((RegionCode & 0x000f) == 0x000f)
                r += 2; // code for level set 1 is +FAR
            if ((RegionCode & 0x000f) == 0x0001)
                r += 1; // code for level set 1 is -FAR

            ReducedRegionCode ret;
            ret.rrc = r;
            return ret;
        }
    }


}
