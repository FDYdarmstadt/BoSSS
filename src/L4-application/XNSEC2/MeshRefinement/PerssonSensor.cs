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
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.XNSEC {


    public class PerssonSensor2 {

        private string sensorVariable;

        private ConventionalDGField sensorField;

        /// <summary>
        /// The field which should be tested for high energy modes
        /// </summary>
        public ConventionalDGField fieldToTestRestricted;

        /// <summary>
        /// Sensor field
        /// </summary>
        public PerssonSensor2(DGField fieldToTest) {
            sensorVariable = fieldToTest.Identification;
            sensorField = new SinglePhaseField(
                new Basis(fieldToTest.GridDat, 0), "sensor");
            fieldToTestRestricted = new SinglePhaseField(
                new Basis(fieldToTest.GridDat, fieldToTest.Basis.Degree - 1));
        }

        /// <summary>
        /// Get the mean value of sensor filed
        /// </summary>
        public double GetValue(int cellIndex) {
            return sensorField.GetMeanValue(cellIndex);
        }

        /// <summary>
        /// Get the actual sensor field
        /// </summary>
        public DGField GetField() {
            return sensorField;
        }

        /// <summary>
        /// Update the actual sensor field
        /// </summary>
        public void Update(DGField fieldToTest) {
            sensorField.Clear();
            fieldToTestRestricted.Clear();

            fieldToTestRestricted.AccLaidBack(1.0, fieldToTest);

            DGField difference = fieldToTestRestricted - fieldToTest;

            foreach (Chunk chunk in CellMask.GetFullMask(fieldToTest.GridDat)) {
                for (int i = 0; i < chunk.Len; i++) {
                    int cell = i + chunk.i0;

                    CellMask singleCellMask = new CellMask(fieldToTest.GridDat, Chunk.GetSingleElementChunk(cell));

                    CellQuadratureScheme scheme = new CellQuadratureScheme(domain: singleCellMask);
                    var rule = scheme.SaveCompile(fieldToTest.GridDat, 2 * fieldToTest.Basis.Degree);

                    double[] a = difference.LocalLxError((ScalarFunction)null, null, rule);
                    double[] b = fieldToTest.LocalLxError((ScalarFunction)null, null, rule);
                    double mySensorValue = a[0] / b[0];

                    // not using L2-norm, but rather only the scalar product
                    mySensorValue = mySensorValue * mySensorValue;

                    if (double.IsNaN(mySensorValue)) {
                        mySensorValue = 0.0;
                    }

                    sensorField.SetMeanValue(cell, mySensorValue);

                }
            }
        }
    }

    /*
        public class PerssonSensor {

            private double[] sensorValues;

            private LevelSetTracker levelSetTracker;

            private SinglePhaseField sensorField;

            public PerssonSensor(LevelSetTracker levelSetTracker) {
                this.levelSetTracker = levelSetTracker;
                this.sensorField = new SinglePhaseField(new Basis(levelSetTracker.GridDat, degree: 0), "AAA");
            }

            public void UpdateSensorValues(DGField field, CellMask cellMask) {
                using (new FuncTrace()) {
                    int degree = field.Basis.Degree;
                    int noOfCells = field.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                    if (sensorValues == null || sensorValues.Length != noOfCells) {
                        sensorValues = new double[noOfCells];
                    }

                    // Clear sensor field just to be sure that there are not any old entries
                    this.sensorField.Clear();

                    IMatrix coordinatesTimesMassMatrix;
                    IMatrix coordinatesTruncatedTimesMassMatrix;

                    bool massMatrixisIdentity = false;

                    if (!massMatrixisIdentity) {
                        // Note: This has to be the _non_-agglomerated mass matrix
                        // because we still live on the non-agglomerated mesh at this point
                        SpeciesId speciesID = this.levelSetTracker.GetSpeciesId("A");
                        MassMatrixFactory massMatrixFactory = this.levelSetTracker.GetXDGSpaceMetrics(new SpeciesId[] { speciesID }, CutCellsQuadOrder: 6).MassMatrixFactory;
                        BlockMsrMatrix massMatrix = massMatrixFactory.GetMassMatrix(field.Mapping, inverse: false);

                        // Old
                        DGField temp = field.CloneAs();
                        massMatrix.SpMV(1.0, field.CoordinateVector, 0.0, temp.CoordinateVector);
                        coordinatesTimesMassMatrix = temp.Coordinates;

                        // Neu
                        DGField uTruncated = field.CloneAs();

                        // Set all coordinates to zero
                        for (int cell = 0; cell < uTruncated.Coordinates.NoOfRows; cell++) {
                            for (int coordinate = 0; coordinate < uTruncated.Coordinates.NoOfCols; coordinate++) {
                                uTruncated.Coordinates[cell, coordinate] = 0;
                            }
                        }

                        // Copy only the coordiantes that belong to the highest modes 
                        foreach (int cell in cellMask.ItemEnum) {
                            foreach (int coordinate in field.Basis.GetPolynomialIndicesForDegree(cell, degree)) {
                                uTruncated.Coordinates[cell, coordinate] = field.Coordinates[cell, coordinate];
                            }
                        }

                        // Calculate M times u
                        DGField vecF_Field = field.CloneAs();
                        massMatrix.SpMV(1.0, uTruncated.CoordinateVector, 0.0, vecF_Field.CoordinateVector);
                        coordinatesTruncatedTimesMassMatrix = vecF_Field.Coordinates;
                    } else {
                        // Mass matrix is identity
                        coordinatesTimesMassMatrix = field.Coordinates;
                        coordinatesTruncatedTimesMassMatrix = field.Coordinates;
                    }


                    // This is equivalent to norm(restrictedField) / norm(originalField)
                    // Note: THIS WILL FAIL IN TRUE XDG CUT CELLS WITH TWO SPECIES
                    foreach (int cell in cellMask.ItemEnum) {
                        double numerator = 0.0;
                        foreach (int coordinate in field.Basis.GetPolynomialIndicesForDegree(cell, degree)) {
                            //numerator += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                            numerator += field.Coordinates[cell, coordinate] * coordinatesTruncatedTimesMassMatrix[cell, coordinate];
                        }

                        double denominator = 0.0;
                        for (int coordinate = 0; coordinate < field.Basis.Length; coordinate++) {
                            //denominator += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                            denominator += field.Coordinates[cell, coordinate] * coordinatesTimesMassMatrix[cell, coordinate];
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

                        this.sensorField.SetMeanValue(cell, result);
                    }
                }
            }

            public double GetSensorValue(int cellIndex) {
                return sensorValues[cellIndex];
            }

            public SinglePhaseField GetSensorField() {
                return this.sensorField;
            }
        }

        */


}