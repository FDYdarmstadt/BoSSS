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

using BoSSS.Foundation.XDG;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    public abstract class AMRLevelIndicatorWithLevelset : AMRLevelIndicator {


        /// <summary>
        /// control object
        /// </summary>
        new protected SolverWithLevelSetUpdaterControl Control {
            get {
                return (SolverWithLevelSetUpdaterControl)(base.Control);
            }
        }


        /// <summary>
        /// the level-set which represents the fluid interface
        /// </summary>
        protected LevelSet LevSet {
            get {
                return (LevelSet)(this.LsTrk.LevelSets[0]);
            }
        }

        public override bool Equals(object obj) {
            if (!base.Equals(obj))
                return false;
            
            return true;
        }

        public override int GetHashCode() {
            return base.GetHashCode();
        }


    }
}
