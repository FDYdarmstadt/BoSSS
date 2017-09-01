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
    /// 
    /// </summary>
    public enum IOListOption {
        /// <summary>
        /// field will always be added to the <see cref="Application{T}.IOFields"/> set,
        /// no matter what the control file says
        /// </summary>
        Always,

        /// <summary>
        /// field will be added to the <see cref="Application{T}.IOFields"/> set,
        /// if specified in the control file
        /// </summary>
        ControlFileDetermined,

        /// <summary>
        /// field will never be added to the <see cref="Application{T}.IOFields"/> set,
        /// no matter what the control file says
        /// </summary>
        Never
    }
}
