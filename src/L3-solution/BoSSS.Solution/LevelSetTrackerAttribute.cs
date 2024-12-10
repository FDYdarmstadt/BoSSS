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

namespace BoSSS.Solution {
    /// <summary>
    /// Attribute for decorating the <see cref="BoSSS.Foundation.XDG.LevelSetTracker"/>;
    /// </summary>
    [AttributeUsage(AttributeTargets.Field, AllowMultiple = false, Inherited = false)]
    public class LevelSetTrackerAttribute : Attribute {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="speciesTable">
        /// specified which level set signs map to which species;
        /// e.g. for three species, one may set<c>---:A --+:B -+-:A -++:C +**:D</c>.
        /// </param>
        /// <param name="NearCellWidth"></param>
        public LevelSetTrackerAttribute(string speciesTable, int NearCellWidth) {
            m_SpeciesTable = speciesTable;
            m_NearCellWidth = NearCellWidth;
        }

        string m_SpeciesTable;

        /// <summary>
        /// desired width of near region for the level-set tracker
        /// </summary>
        internal int m_NearCellWidth;

        internal Array GetSpeciesTable(int NoOfLevelSets) {
            return BoSSS.Foundation.XDG.LevelSetTracker.GetSpeciesTable(m_SpeciesTable, NoOfLevelSets);
        }


    }

}
