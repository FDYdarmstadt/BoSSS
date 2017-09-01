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
using System.IO;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// see <see cref="LevelSetTracker.UpdateTracker(int, int[])"/>;
    /// </summary>
    public class LevelSetCFLException : Exception {

        static string ComposeMessage(int[] failCnt) {
            StringWriter stw = new StringWriter();
            stw.Write("failed on Level Set(s): ");
            bool beistrich = false;
            for (int i = 0; i < failCnt.Length; i++) {
                if (failCnt[i] > 0) {
                    if (beistrich)
                        stw.Write(", ");
                    stw.Write(i + " (" + failCnt[i] + " cells)");
                }
                stw.Write(";");
            }

            return stw.ToString();
        }

        /// <summary>
        /// ctor
        /// </summary>
        internal LevelSetCFLException(int[] __fail)
            : base(ComposeMessage(__fail)) {
            m_IndicesOfViolatingLevSets = __fail;
        }

        int[] m_IndicesOfViolatingLevSets;

        /// <summary>
        /// index: Level Set index <br/>
        /// content: Number of cells in which the CFL condition is violated.
        /// </summary>
        public int[] NoOfViolatedCellsPerLevSet {
            get {
                return m_IndicesOfViolatingLevSets;
            }
        }
    }
}
