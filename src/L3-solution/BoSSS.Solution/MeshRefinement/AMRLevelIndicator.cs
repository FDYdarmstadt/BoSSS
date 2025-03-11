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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {

    /// <summary>
    /// Base-class for adaptive mesh refinement level indication.
    /// Objects of this type can be added to the list `activeAMRlevelIndicators`,
    /// and the respective post-processing (implemented in the derived class) will be executed.
    /// </summary>
    [Serializable]
    public abstract class AMRLevelIndicator {
        [DataMember]
        public int maxRefinementLevel = 1;


        /// <summary>
        /// reference to solver application class
        /// </summary>
        [JsonIgnore]
        protected IApplication SolverMain {
            get;
            private set;
        }

        /// <summary>
        /// reference to application control
        /// </summary>
        [JsonIgnore]
        protected Control.AppControl Control {
            get {
                return SolverMain.ControlBase;
            }
        }


        /// <summary>
        /// Access to the level set tracker 
        /// </summary>
        [JsonIgnore]
        protected LevelSetTracker LsTrk {
            get {
                return SolverMain.LsTrk;
            }
        }


        /// <summary>
        /// 
        /// </summary>
        [JsonIgnore]
        protected BoSSS.Foundation.Grid.Classic.GridData GridData {
            get {
                return (Foundation.Grid.Classic.GridData)(this.SolverMain.GridData);
            }
        }


        [JsonIgnore]
        [NonSerialized]
        bool m_Setupdone = false;

        /// <summary>
        /// 
        /// </summary>
        virtual public void Setup(IApplication solverMain) {
            if (m_Setupdone)
                throw new NotSupportedException("Setup Routine may only be called once.");
            m_Setupdone = true;

            this.SolverMain = solverMain;
        }


        /// <summary>
        /// checks whether a cell needs to be refined or may be coarsen.
        /// </summary>
        /// <returns> on return the values define:
        ///     -1: cell may be coarsen
        ///      0: no changes in cell 
        ///     +1: cell needs to be refined
        /// the array index corresponds to the local cell index    
        /// </returns>
        public abstract int[] DesiredCellChanges();

        /// <summary>
        /// comparison
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as AMRLevelIndicator;
            if (other == null)
                return false;

            if (this.GetType() != other.GetType())
                return false;

            return this.maxRefinementLevel == other.maxRefinementLevel;
        }

        /// <summary>
        /// 
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }
    }

}
