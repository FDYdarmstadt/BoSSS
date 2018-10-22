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
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using System;
using System.Collections.Generic;

namespace CNS.ShockCapturing {

    public class ConservativeFVMLimiter : ILimiter {

        private double sensorLimit;

        private double cellSize;

        private int dgDegree;

        public ConservativeFVMLimiter(IShockSensor sensor, double sensorLimit, double cellSize, int dgDegree) {
            this.Sensor = sensor;
            this.sensorLimit = sensorLimit;
            this.cellSize = cellSize;
            this.dgDegree = dgDegree;
        }

        public IShockSensor Sensor {
            get;
            private set;
        }

        public void LimitFieldValues(IEnumerable<DGField> fielsSet) {

            IProgram<CNSControl> program;
            CellMask shockedCells = Sensor.GetShockedCellMask(program.GridData, sensorLimit, cellSize, dgDegree);



            foreach (DGField field in program.WorkingSet.ConservativeVariables) {
                for (int i = 0; i < field.GridDat.iLogicalCells.NoOfLocalUpdatedCells; i++) { // loop over all cells
                    foreach (Chunk chunk in shockedCells) { // loop over all cells in mask *inside* loop over cells => quadratic runtime
                        foreach (int cell in chunk.Elements) { // where is the 'cell' index actually used?? 
                            for (int j = 1; j < field.Coordinates.NoOfCols; j++) {
                                field.Coordinates[i, j] = 0.0;
                            }
                        }
                    }
                }
            }

            throw new NotImplementedException("Implementation is seriously flawed, quadratic runtime w.r.t. number of cells: check implementation!");

        }
    }
}
