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
using BoSSS.Foundation;

namespace BoSSS.Solution.XNSECommon {
    
    
    
    public class OperatorConfiguration {
        public PhysicalParameters physParams;

        /// <summary>
        /// advanced operator configuration
        /// </summary>
        public DoNotTouchParameters dntParams;

        /// <summary>
        /// Controls the domain variables that the operator should contain. <br/>
        /// This controls only the formal operator shape, not the actual components.
        /// index: 0 -- Velocity components, 1 -- pressure component; 
        /// </summary>
        public bool[] DomBlocks = new bool[2];

        /// <summary>
        /// Controls the codomain variables that the operator should contain. <br/>
        ///This controls only the formal operator shape, not the actual components.
        /// index: 0 -- momentum equation, 1 -- continuity equation; 
        /// </summary>
        public bool[] CodBlocks = new bool[2];

        /// <summary>
        /// include transport operator
        /// </summary>
        public bool Transport;

        /// <summary>
        /// include viscous operator
        /// </summary>
        public bool Viscous;

        /// <summary>
        /// include pressure gradient
        /// </summary>
        public bool PressureGradient;

        /// <summary>
        /// Use the surface force term -- mostly only useful for manufactured solutions.
        /// </summary>
        public bool ArtificialsurfaceForce = false;

        /// <summary>
        /// include continuity equation
        /// </summary>
        public bool continuity;

        /// <summary>
        /// Switch to turn velocity extension on/off.
        /// </summary>
        public bool UseXDG4Velocity = true;
    }
}
