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
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using System.Collections.Generic;
using System.Linq;

namespace CNS.ShockCapturing {

    public class PrimitiveFVMLimiter : ILimiter {

        private double sensorLimit;

        private double cellSize;

        private int dgDegree;

        private IProgram<CNSControl> program;

        public PrimitiveFVMLimiter(IShockSensor sensor, double sensorLimit, double cellSize, int dgDegree, IProgram<CNSControl> program) {
            this.Sensor = sensor;
            this.sensorLimit = sensorLimit;
            this.cellSize = cellSize;
            this.dgDegree = dgDegree;
            this.program = program;
        }

        public IShockSensor Sensor {
            get;
            private set;
        }

        public void LimitFieldValues(IEnumerable<DGField> ConservativeVariables, IEnumerable<DGField> DerivedFields) {

            var GridData = ConservativeVariables.First().GridDat;

            // Make sure primitive fields are up-to-date
            for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                CNSVariables.Velocity[d].UpdateFunction(
                    DerivedFields.Single(f => f.Identification == CNSVariables.Velocity[d].Name),
                    CellMask.GetFullMask(GridData),
                    program);
            }
            CNSVariables.Pressure.UpdateFunction(
                DerivedFields.Single(f => f.Identification == CNSVariables.Pressure.Name),
                CellMask.GetFullMask(GridData),
                program);

            // Limit and store primitive fields
            string[] primitiveFieldNames = { CompressibleVariables.Density, CNSVariables.Velocity.xComponent, CNSVariables.Velocity.yComponent, CNSVariables.Pressure };
            DGField[] primitiveFields = new DGField[primitiveFieldNames.Length];
            CellMask shockedCells = Sensor.GetShockedCellMask(GridData, sensorLimit, cellSize, dgDegree);
            int k = 0;
            foreach (string name in primitiveFieldNames) {
                DGField field = ConservativeVariables.
                    Concat(DerivedFields).
                    Where(f => f.Identification.Equals(name)).
                    Single();

                foreach (Chunk chunk in shockedCells) {
                    foreach (int cell in chunk.Elements) {
                        for (int j = 1; j < field.Coordinates.NoOfCols; j++) {
                            field.Coordinates[cell, j] = 0.0;
                        }
                    }
                }

                primitiveFields[k] = field;
                k++;
            }

            // Update conservative variables by using limited primitive variables only in shocked cells
            int D = CompressibleEnvironment.NumberOfDimensions;

            for (int d = 0; d < D; d++) {
                DGField mom_d = ConservativeVariables.Single(f => f.Identification == CompressibleVariables.Momentum[d].Name);
                mom_d.Clear(shockedCells);
                mom_d.ProjectFunction(
                    1.0,
                    (X, U, j) => U[0] * U[1 + d],
                    new CellQuadratureScheme(true, shockedCells),
                    primitiveFields);
            }

            // Update total energy
            DGField Energy = ConservativeVariables.Single(f => f.Identification == CompressibleVariables.Energy.Name);
            Energy.Clear(shockedCells);
            Energy.ProjectFunction(
                1.0,
                delegate (double[] X, double[] U, int jCell) {
                    double K = 0.0;
                    for (int d = 0; d < D; d++) {
                        K += U[d + 1] * U[d + 1];
                    }
                    return U[D + 1] / (program.Control.EquationOfState.HeatCapacityRatio - 1.0) + 0.5 * U[0] * K;
                },
                new CellQuadratureScheme(true, shockedCells),
                primitiveFields);
        }
    }
}
