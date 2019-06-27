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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    public interface IShockSensor {

        /// <summary>
        /// Returns the current value of the sensor in cell <paramref name="cell"/>
        /// </summary>
        /// <param name="cell"></param>
        /// <returns>The sensor value as a double</returns>
        double GetSensorValue(int cell);

        /// <summary>
        /// Returns the sensor values stored in a DG field
        /// </summary>
        /// <returns>Sensor values in a DG field</returns>
        DGField GetSensorField();
    }

    /// <summary>
    /// Extension methods for <see cref="IShockSensor"/>
    /// </summary>
    public static class IShockSensorExtensions {

        /// <summary>
        /// Returns a cell mask containing all cells that are considered shock
        /// if the sensor limit is given <paramref name="sensorLimit"/>
        /// </summary>
        /// <returns></returns>
        public static CellMask GetShockedCellMask(this IShockSensor sensor, IGridData gridData, double sensorLimit, double cellSize, int dgDegree) {
            BitArray shockedCellArray = new BitArray(gridData.iLogicalCells.NoOfLocalUpdatedCells);
            for (int i = 0; i < gridData.iLogicalCells.NoOfLocalUpdatedCells; i++) {
                double sensorValue = sensor.GetSensorValue(i);
                shockedCellArray[i] = sensorValue > (sensorLimit * cellSize / dgDegree);
            }
            return new CellMask(gridData, shockedCellArray);
        }
    }
}
