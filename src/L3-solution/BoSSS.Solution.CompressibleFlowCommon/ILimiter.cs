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

using BoSSS.Foundation;
using System.Collections.Generic;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    /// <summary>
    /// A generic sloper limiter
    /// </summary>
    public interface ILimiter {

        /// <summary>
        /// The sensor used by this slope limiter to detected troubled cells
        /// </summary>
        IShockSensor Sensor {
            get;
        }

        /// <summary>
        /// Limits the values of the primal variables
        /// </summary>
        /// <param name="program"></param>
        void LimitFieldValues(IEnumerable<DGField> ConservativeVariables, IEnumerable<DGField> DerivedFields);
    }
}
