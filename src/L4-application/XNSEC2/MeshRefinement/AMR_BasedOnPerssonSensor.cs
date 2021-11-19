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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// refinement of cells based on a sensor
    /// </summary>
    [Serializable]
    public class AMR_BasedOnPerssonSensor : AMRLevelIndicatorWithLevelset {

        [DataMember]
        private string FieldName;

        /// <summary>
        /// Empty constructor for serialization
        /// </summary>
        private AMR_BasedOnPerssonSensor() { }

        public AMR_BasedOnPerssonSensor(string fieldname, int _refLevel) {
            sensor = new PerssonSensor();
            FieldName = fieldname;
            maxRefinementLevel = _refLevel;
        }

        [DataMember]
        protected Sensor sensor;

        [DataMember]
        protected DGField FieldForSensor {
            get {
                return (XDGField)((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == FieldName).SingleOrDefault();
            }
        }

        public override int[] DesiredCellChanges() {
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            int D = GridData.SpatialDimension;
            var field = (XDGField)this.FieldForSensor;

            if(field == null) // this should take care of the startup AMR, where temperature is still not defined
                return levels;
            
            sensor.UpdateSensorValues((XDGField)this.FieldForSensor);

            double[] sensorValues = new double[J];

            for(int j = 0; j < J; j++) {
                sensorValues[j] = sensor.GetSensorValue(j);
            }

            var sensorField = new SinglePhaseField(new Basis(GridData, 2), "Sensor");
            var Sensor_aftermodifications = new SinglePhaseField(new Basis(GridData, 2), "Sensor_aftermodifications");
            foreach(Chunk chunk in CellMask.GetFullMask(sensorField.GridDat)) {
                for(int i = 0; i < chunk.Len; i++) {
                    int cell = chunk.i0 + i;
                    sensorField.SetMeanValue(cell, sensorValues[i]);
                }
            }

            //var meanval = sensorField.GetMeanValueTotal(null);
            //  sensorField.AccConstant(-1.0 * meanval);

            double minSensorValue; double maxSensorValue;
            sensorField.GetExtremalValues(out minSensorValue, out maxSensorValue);

            double denom = maxSensorValue;

            for(int j = 0; j < J; j++) {
                double sensorval = sensor.GetSensorValue(j);
                int currentLevel = GridData.Grid.Cells[j].RefinementLevel;
                if((sensorval / denom > 0.5 && sensorval > 1e-10) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                }
            }
            Console.WriteLine("Number of refined cells around areas where pearson sensor for variable" + FieldName+ "is: " + levels.Sum() + ".\n");
            //foreach(Chunk chunk in CellMask.GetFullMask(sensorField.GridDat)) {
            //    for(int i = 0; i < chunk.Len; i++) {
            //        int cell = chunk.i0 + i;
            //        Sensor_aftermodifications.SetMeanValue(cell, SensorAfterModificationsValues[i]);
            //    }
            //}

            //Tecplot.PlotFields(new DGField[] { sensorField, Sensor_aftermodifications }, "SensorDebugging", 0.0, 3);
            return levels;
        }
    }

    public class PerssonSensor : Sensor {
        private double[] sensorValues;

        public PerssonSensor() {
        }

        public override double GetSensorValue(int cellIndex) {
            return sensorValues[cellIndex];
        }

        /// <summary>
        /// Update the actual sensor field
        /// </summary>
        public override void UpdateSensorValues(DGField fieldToTest) {
            int noOfCells = fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            if(sensorValues == null || sensorValues.Length != noOfCells) {
                sensorValues = new double[noOfCells];
            }

            //ConventionalDGField fieldToTestRestricted = new SinglePhaseField(new Basis(fieldToTest.GridDat, fieldToTest.Basis.Degree - 1));
            XDGField fieldToTestRestricted = new XDGField(new XDGBasis(((XDGField)fieldToTest).Basis.Tracker, fieldToTest.Basis.Degree - 1));
            fieldToTestRestricted.Clear();
            fieldToTestRestricted.AccLaidBack(1.0, fieldToTest);

            DGField difference = fieldToTestRestricted - fieldToTest;

            foreach(Chunk chunk in CellMask.GetFullMask(fieldToTest.GridDat)) {
                for(int i = 0; i < chunk.Len; i++) {
                    int cell = i + chunk.i0;

                    CellMask singleCellMask = new CellMask(fieldToTest.GridDat, Chunk.GetSingleElementChunk(cell));

                    CellQuadratureScheme scheme = new CellQuadratureScheme(domain: singleCellMask);
                    var rule = scheme.SaveCompile(fieldToTest.GridDat, 2 * fieldToTest.Basis.Degree);

                    double[] a = difference.LocalLxError((ScalarFunction)null, null, rule);
                    double[] b = fieldToTest.LocalLxError((ScalarFunction)null, null, rule);
                    double mySensorValue = a[0] / b[0];

                    // not using L2-norm, but rather only the scalar product
                    mySensorValue = mySensorValue * mySensorValue;

                    if(double.IsNaN(mySensorValue)) {
                        mySensorValue = 0.0;
                    }
                    sensorValues[cell] = mySensorValue;
                }
            }
        }
    }
}