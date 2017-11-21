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
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CNS.ShockCapturing {

    public class PerssonSensor : IShockSensor {

        private Variable sensorVariable;

        private double sensorLimit;

        private double[] sensorValues;

        private int fileCount = 0;

        private int lineCount = 0;

        //private StreamWriter sensorOnOffWriter = null;

        //private DGField fieldToTestRestricted;

        ICompositeQuadRule<QuadRule> compositeRule;

        public PerssonSensor(Variable sensorVariable, double sensorLimit) {
            this.sensorVariable = sensorVariable;
            this.sensorLimit = sensorLimit;
        }

        public void UpdateSensorValues(CNSFieldSet fieldSet) {
            DGField fieldToTest = fieldSet.ConservativeVariables.
                Concat(fieldSet.DerivedFields.Values).
                Where(f => f.Identification == sensorVariable.Name).
                Single();
            int degree = fieldToTest.Basis.Degree;
            int noOfCells = fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if (sensorValues == null || sensorValues.Length != noOfCells) {
                sensorValues = new double[noOfCells];
            }

            // This is equivalent to norm(restrictedField) / norm(originalField)
            // Note: THIS WILL FAIL IN CUT CELLS
            for (int cell = 0; cell < noOfCells; cell++) {
                double numerator = 0.0;
                foreach (int coordinate in fieldToTest.Basis.GetPolynomialIndicesForDegree(cell, degree)) {
                    numerator += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                }

                double denominator = 0.0;
                for (int coordinate = 0; coordinate < fieldToTest.Basis.Length; coordinate++) {
                    denominator += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                }

                Debug.Assert(denominator != 0, "Persson sensor: Denominator is zero!");

                sensorValues[cell] = numerator / denominator;
            }
            

            //if (fieldToTestRestricted == null || fieldToTestRestricted.Basis.Degree != degree || sensorValues.Length != noOfCells) {
            //    fieldToTestRestricted = new SinglePhaseField(
            //        new Basis(fieldToTest.GridDat, degree));
            //    CellQuadratureScheme scheme = new CellQuadratureScheme();
            //    compositeRule = scheme.SaveCompile(fieldToTest.GridDat, 2 * fieldToTest.Basis.Degree);
            //}

            //int i0 = 1;
            //int degreeMinusOne = degree - 1;
            //for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
            //    i0 *= (degreeMinusOne + d + 1);
            //}
            //for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
            //    i0 /= (d + 1);
            //}

            //for (int cell = 0; cell < noOfCells; cell++) {
            //    for (int coordinate = i0; coordinate < fieldToTest.Basis.Length; coordinate++) {
            //        fieldToTestRestricted.Coordinates[cell, coordinate] = fieldToTest.Coordinates[cell, coordinate];
            //    }
            //}

            //// Note: This is the _only_ expensive part!
            //double[] a = fieldToTestRestricted.LocalLxError((ScalarFunction)null, null, compositeRule);
            //double[] b = fieldToTest.LocalLxError((ScalarFunction)null, null, compositeRule);


            ////double[] cellCenters = ((GridData)fieldSet.Density.GridDat).Cells.CellCenter.ExtractSubArrayShallow(-1, 0).To1DArray();

            ////if (sensorOnOffWriter == null) {
            ////    sensorOnOffWriter = CreateNewStreamWriter();
            ////}

            //sensorLimit = sensorLimit / Math.Pow(fieldToTest.Basis.Degree, 4);
            //for (int cell = 0; cell < noOfCells; cell++) {
            //    double mySensorValue = a[cell] / b[cell];

            //    #region AV on in entire domain (checked for each time step)
            //    //double maxSensor = sensorValues.GetMeanValue(0);
            //    //if (mySensorValue > maxSensor) {
            //    //    for (int j = 0; j < fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells; j++) {
            //    //        sensorValues.SetMeanValue(j, mySensorValue);
            //    //    }
            //    //}
            //    #endregion

            //    // comment when using region AV on in entire domain (checked for each time step)
            //    sensorValues[cell] = mySensorValue;

            //    // Store AV On/Off-Values to extra file
            //    //if (mySensorValue > sensorLimit) {
            //    //    sensorOnOffWriter.Write(cellCenters[i + chunk.i0] + "\t" + "{0:0.000000000} \r\n", physTime);
            //    //    lineCount++;
            //    //    if (FileIsFull(5000)) {
            //    //        sensorOnOffWriter.Flush();
            //    //        sensorOnOffWriter.Close();
            //    //        sensorOnOffWriter = CreateNewStreamWriter();
            //    //    }
            //    //}
            //}
            ////sensorOnOffWriter.Flush();
        }

        private StreamWriter CreateNewStreamWriter() {
            fileCount++;
            string streamPath = "AV_OnOff_" + fileCount + ".txt";
            //string streamPath = @"\\fdyprime\userspace\geisenhofer\AV_OnOff_" + fileCount + ".txt";
            StreamWriter streamWriter = new StreamWriter(streamPath);
            lineCount = 1;
            streamWriter.WriteLine("x" + "\t" + "t");
            lineCount++;
            streamWriter.Flush();
            return streamWriter;
        }

        private bool FileIsFull(int maxNumOfLines) {
            if (lineCount % maxNumOfLines == 0) {
                return true;
            } else
                return false;
        }

        public double GetSensorValue(int cellIndex) {
            return sensorValues[cellIndex];
        }

    }
}