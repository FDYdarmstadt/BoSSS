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
using CNS.IBM;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CNS.ShockCapturing {

    public class PerssonSensor : IShockSensor {

        private Variable sensorVariable;

        private double sensorLimit;

        private double[] sensorValues;

        public PerssonSensor(Variable sensorVariable, double sensorLimit) {
            this.sensorVariable = sensorVariable;
            this.sensorLimit = sensorLimit;
        }

        public void UpdateSensorValues(CNSFieldSet fieldSet, ISpeciesMap speciesMap, CellMask cellMask) {
            DGField fieldToTest = fieldSet.ConservativeVariables.
                Concat(fieldSet.DerivedFields.Values).
                Where(f => f.Identification == sensorVariable.Name).
                Single();
            int degree = fieldToTest.Basis.Degree;
            int noOfCells = fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if (sensorValues == null || sensorValues.Length != noOfCells) {
                sensorValues = new double[noOfCells];
            }

            IMatrix coordinatesTimesMassMatrix;
            IMatrix coordinatesTruncatedTimesMassMatrix;
            if (speciesMap is ImmersedSpeciesMap ibmMap) {
                BlockMsrMatrix massMatrix = ibmMap.GetMassMatrixFactory(fieldToTest.Mapping).MassMatrix;

                // Old
                DGField temp = fieldToTest.CloneAs();
                massMatrix.SpMV(1.0, fieldToTest.CoordinateVector, 0.0, temp.CoordinateVector);
                coordinatesTimesMassMatrix = temp.Coordinates;

                // Neu
                DGField uTruncated = fieldToTest.CloneAs();

                // Set all coordinates to zero
                for (int cell = 0; cell < uTruncated.Coordinates.NoOfRows; cell++) {
                    for (int coordinate = 0; coordinate < uTruncated.Coordinates.NoOfCols; coordinate++) {
                        uTruncated.Coordinates[cell, coordinate] = 0;
                    }
                }

                // Copy only the coordiantes that belong to the highest modes 
                foreach (int cell in cellMask.ItemEnum) {
                    foreach (int coordinate in fieldToTest.Basis.GetPolynomialIndicesForDegree(cell, degree)) {
                        uTruncated.Coordinates[cell, coordinate] = fieldToTest.Coordinates[cell, coordinate];
                    }
                }

                // Calculate M times u
                DGField vecF_Field = fieldToTest.CloneAs();
                massMatrix.SpMV(1.0, uTruncated.CoordinateVector, 0.0, vecF_Field.CoordinateVector);
                coordinatesTruncatedTimesMassMatrix = vecF_Field.Coordinates;
            } else {
                // Mass matrix is identity
                coordinatesTimesMassMatrix = fieldToTest.Coordinates;
                coordinatesTruncatedTimesMassMatrix = fieldToTest.Coordinates;
            }
            //IMatrix coordinatesTimesMassMatrix = fieldToTest.Coordinates;

            // This is equivalent to norm(restrictedField) / norm(originalField)
            // Note: THIS WILL FAIL IN CUT CELLS
            foreach (int cell in cellMask.ItemEnum) {
                double numerator = 0.0;
                foreach (int coordinate in fieldToTest.Basis.GetPolynomialIndicesForDegree(cell, degree)) {
                    //numerator += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                    numerator += fieldToTest.Coordinates[cell, coordinate] * coordinatesTruncatedTimesMassMatrix[cell, coordinate];
                }

                double denominator = 0.0;
                for (int coordinate = 0; coordinate < fieldToTest.Basis.Length; coordinate++) {
                    //denominator += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                    denominator += fieldToTest.Coordinates[cell, coordinate] * coordinatesTimesMassMatrix[cell, coordinate];
                }

                Debug.Assert(denominator != 0, "Persson sensor: Denominator is zero!");

                Debug.Assert(!(numerator / denominator).IsNaN(), "Persson sensor: Sensor value is NaN!");

                Debug.Assert(numerator / denominator >= 0, "Persson sensor: Sensor value is negative!");

                sensorValues[cell] = numerator / denominator;
            }
        }

        public double GetSensorValue(int cellIndex) {
            return sensorValues[cellIndex];
        }

    }
}