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
using System.Diagnostics;

namespace BoSSS.Foundation.Quadrature {
    
    /// <summary>
    /// common methods for alll kinds of quadrature objects
    /// </summary>
    public interface IQuadrature {

        /// <summary>
        /// execute the quadrature
        /// </summary>
        void Execute();

        /// <summary>
        /// Additional timers, which may be derived in derived classes, that are added to the performance analysis.
        /// </summary>
        Stopwatch[] CustomTimers {
            get;
            set;
        }

        /// <summary>
        /// Names of the custom timers (see <see cref="CustomTimers"/>), for reference reasons in the log files.
        /// Index correlates with <see cref="CustomTimers"/>;
        /// </summary>
        string[] CustomTimers_Names {
            get;
            set;
        }

        /// <summary>
        /// If one wants to establish a call tree among the custom timers, one can use this.
        /// </summary>
        int[] CustomTimers_RootPointer {
            get;
            set;
        }
    }
}
