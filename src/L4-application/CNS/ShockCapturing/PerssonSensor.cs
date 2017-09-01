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
using BoSSS.Foundation.Quadrature;
using System;
using System.IO;
using System.Linq;

namespace CNS.ShockCapturing {

    public class PerssonSensor : IShockSensor {

        private Variable sensorVariable;

        private double sensorLimit;

        private DGField sensorValues;

        private int fileCount = 0;

        private int lineCount = 0;

        //private StreamWriter sensorOnOffWriter = null;

        public PerssonSensor(Variable sensorVariable, double sensorLimit) {
            this.sensorVariable = sensorVariable;
            this.sensorLimit = sensorLimit;
        }

        public void UpdateSensorValues(CNSFieldSet fieldSet) {
            DGField fieldToTest = fieldSet.ConservativeVariables.
                Concat(fieldSet.DerivedFields.Values).
                Where(f => f.Identification == sensorVariable.Name).
                Single();

            DGField fieldToTestRestricted = new SinglePhaseField(
                new Basis(fieldToTest.GridDat, fieldToTest.Basis.Degree - 1));
            fieldToTestRestricted.AccLaidBack(1.0, fieldToTest);

            DGField difference = fieldToTestRestricted - fieldToTest;

            sensorValues = new SinglePhaseField(
                new Basis(fieldToTest.GridDat, 0));

            CellQuadratureScheme scheme = new CellQuadratureScheme();
            var rule = scheme.SaveCompile(fieldToTest.GridDat, 2 * fieldToTest.Basis.Degree);
            double[] a = difference.LocalLxError((ScalarFunction)null, null, rule);


            double[] b = fieldToTest.LocalLxError((ScalarFunction)null, null, rule);
            //var clone = fieldToTest.CloneAs();
            //foreach (Chunk chunk in CellMask.GetFullMask(clone.GridDat)) {
            //    for (int i = 0; i < chunk.Len; i++) {
            //        int cell = i + chunk.i0;
            //        clone.SetMeanValue(cell, 1.0);
            //    }
            //}
            //double[] b = clone.LocalLxError((ScalarFunction)null, null, rule);

            int cellIndex = 0;

            sensorLimit = sensorLimit / Math.Pow(fieldToTest.Basis.Degree, 4);

            double[] cellCenters = ((GridData)fieldSet.Density.GridDat).Cells.CellCenter.ExtractSubArrayShallow(-1, 0).To1DArray();

            //if (sensorOnOffWriter == null) {
            //    sensorOnOffWriter = CreateNewStreamWriter();
            //}

            foreach (Chunk chunk in CellMask.GetFullMask(fieldToTest.GridDat)) {
                for (int i = 0; i < chunk.Len; i++) {
                    int cell = i + chunk.i0;

                    //CellMask singleCellMask = new CellMask(fieldToTest.GridDat, Chunk.GetSingleElementChunk(cell));
                    //double a = difference.L2Norm(singleCellMask);
                    //double b = fieldToTest.L2Norm(singleCellMask);
                    double mySensorValue = a[cellIndex] / b[cellIndex];

                    // not using L2-norm --> only the scalar product
                    //mySensorValue = mySensorValue * mySensorValue;

                    #region AV on in entire domain (checked for each time step)
                    //double maxSensor = sensorValues.GetMeanValue(0);
                    //if (mySensorValue > maxSensor) {
                    //    for (int j = 0; j < fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells; j++) {
                    //        sensorValues.SetMeanValue(j, mySensorValue);
                    //    }
                    //}
                    #endregion

                    // comment when using region AV on in entire domain (checked for each time step)
                    sensorValues.SetMeanValue(cell, mySensorValue);

                    // Store AV On/Off-Values to extra file
                    //if (mySensorValue > sensorLimit) {
                    //    sensorOnOffWriter.Write(cellCenters[i + chunk.i0] + "\t" + "{0:0.000000000} \r\n", physTime);
                    //    lineCount++;
                    //    if (FileIsFull(5000)) {
                    //        sensorOnOffWriter.Flush();
                    //        sensorOnOffWriter.Close();
                    //        sensorOnOffWriter = CreateNewStreamWriter();
                    //    }
                    //}

                    cellIndex++;
                }
                //sensorOnOffWriter.Flush();
            }
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
            return sensorValues.GetMeanValue(cellIndex);
        }

    }
}