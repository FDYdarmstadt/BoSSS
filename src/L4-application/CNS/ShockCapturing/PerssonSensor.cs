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
using CNS.IBM;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Collections.Generic;
using System.Linq;

namespace CNS.ShockCapturing {

    public class PerssonSensor : IShockSensor {

        private string m_sensorVariableName;

        private double sensorLimit;

        private double[] sensorValues;

        public PerssonSensor(string sensorVariableName, double sensorLimit) {
            this.m_sensorVariableName = sensorVariableName;
            this.sensorLimit = sensorLimit;
        }

        public void UpdateSensorValues(IEnumerable<DGField> fieldSet, ISpeciesMap speciesMap, CellMask cellMask) {
            using (new FuncTrace()) {
                DGField fieldToTest = fieldSet.
                    Where(f => f.Identification == m_sensorVariableName).
                    Single();
                int degree = fieldToTest.Basis.Degree;
                int noOfCells = fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                if (sensorValues == null || sensorValues.Length != noOfCells) {
                    sensorValues = new double[noOfCells];
                }

                IMatrix coordinatesTimesMassMatrix;
                IMatrix coordinatesTruncatedTimesMassMatrix;
                if (speciesMap is ImmersedSpeciesMap ibmMap) {
                    // Note: This has to be the _non_-agglomerated mass matrix
                    // because we still live on the non-agglomerated mesh at this
                    // point
                    BlockMsrMatrix massMatrix = ibmMap.GetMassMatrixFactory(fieldToTest.Mapping).NonAgglomeratedMassMatrix;

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

                //cellMask.SaveToTextFile("fluidCells.txt");

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

                    double result;
                    if (denominator == 0.0) {
                        result = 0.0;
                    } else {
                        result = numerator / denominator;
                    }

                    //Debug.Assert(denominator != 0, "Persson sensor: Denominator is zero!");

                    //Debug.Assert(!(numerator / denominator).IsNaN(), "Persson sensor: Sensor value is NaN!");

                    //Debug.Assert(numerator / denominator >= 0, "Persson sensor: Sensor value is negative!");

                    sensorValues[cell] = result;
                }
            }
        }

        /// <summary>
        /// Version for XDG
        /// </summary>
        /// <param name="fieldToTest"></param>
        /// <param name="massMatrix"></param>
        /// <param name="cellMask"></param>
        public void UpdateSensorValues(SinglePhaseField fieldToTest, BlockMsrMatrix massMatrix, CellMask cellMask) {
            int degree = fieldToTest.Basis.Degree;
            int noOfCells = fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if (sensorValues == null || sensorValues.Length != noOfCells) {
                sensorValues = new double[noOfCells];
            }

            IMatrix coordinatesTimesMassMatrix;
            IMatrix coordinatesTruncatedTimesMassMatrix;
            // Note: This has to be the _non_-agglomerated mass matrix
            // because we still live on the non-agglomerated mesh at this
            // point
            //MassMatrixFactory massMatrixFactory = levelSetTracker.GetXDGSpaceMetrics(new SpeciesId[] { a_speciesID, b_speciesID }, nonLinearQuadratureDegree).MassMatrixFactory;

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
            //IMatrix coordinatesTimesMassMatrix = fieldToTest.Coordinates;

            //cellMask.SaveToTextFile("fluidCells.txt");

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

                double result;
                if (denominator == 0.0) {
                    result = 0.0;
                } else {
                    result = numerator / denominator;
                }

                //Debug.Assert(denominator != 0, "Persson sensor: Denominator is zero!");

                //Debug.Assert(!(numerator / denominator).IsNaN(), "Persson sensor: Sensor value is NaN!");

                //Debug.Assert(numerator / denominator >= 0, "Persson sensor: Sensor value is negative!");

                sensorValues[cell] = result;
            }
        }

        public double GetSensorValue(int cellIndex) {
            return sensorValues[cellIndex];
        }

    }
}