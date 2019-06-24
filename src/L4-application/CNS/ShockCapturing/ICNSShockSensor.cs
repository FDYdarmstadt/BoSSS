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
using BoSSS.Foundation.Grid;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using System;
using System.Collections;
using System.Collections.Generic;

namespace CNS.ShockCapturing {

    /// <summary>
    /// Defines a sensor that yields large positive values in regions with
    /// strong oscillations and small values otherwise
    /// </summary>
    public interface ICNSShockSensor : IShockSensor {

        /// <summary>
        /// Updates the sensor values in all cells
        /// </summary>
        void UpdateSensorValues(IEnumerable<DGField> fieldSet, ISpeciesMap speciesMap, CellMask cellMask);
    }
}
