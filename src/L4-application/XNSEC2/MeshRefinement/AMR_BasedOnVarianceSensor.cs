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
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using ilPSP;
using System;
using System.Linq;

namespace BoSSS.Application.XNSEC {

 

    /// <summary>
    /// refinement of cells based on a gradient of a variable
    /// </summary>
    [Serializable]
    public class AMR_BasedOnVarianceSensor : AMRLevelIndicatorWithLevelset {
        public string FieldName;

        /// <summary>
        /// Empty constructor for serialization
        /// </summary>
        private AMR_BasedOnVarianceSensor() { }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mySensorType"></param>
        public AMR_BasedOnVarianceSensor(SensorType mySensorType) {
            switch(mySensorType) {
                case SensorType.PearsonSensor:
                sensor =  new PerssonSensor();
                break;
                case SensorType.VarianceSensor:
                sensor = new VarianceSensor();
                break;
                default:
                throw new NotImplementedException("Invalid sensor type");
            }
        }

        public enum SensorType {
            PearsonSensor,

            VarianceSensor,

            none
        }


        protected Sensor sensor;


        protected DGField FieldForSensor {
            get {
                return (XDGField)((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == FieldName).Single();
            }
        }

        public override int[] DesiredCellChanges() {
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];
            int D = GridData.SpatialDimension;

            //sensor.UpdateSensorValues(((XDGField)this.FieldForSensor).GetSpeciesShadowField("A")); // How should i do this?
            sensor.UpdateSensorValues((XDGField)this.FieldForSensor); // How should i do this?

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

            foreach(Chunk chunk in CellMask.GetFullMask(sensorField.GridDat)) {
                for(int i = 0; i < chunk.Len; i++) {
                    int cell = chunk.i0 + i;
                    var value = sensorField.GetMeanValue(cell);
                    sensorField.SetMeanValue(cell, value);
                }
            }

         

            double minSensorValue; double maxSensorValue;
            sensorField.GetExtremalValues(out minSensorValue, out maxSensorValue);


            double denom =  maxSensorValue - minSensorValue;

            double[] SensorAfterModificationsValues = new double[J];
            for(int j = 0; j < J; j++) {
                //double sensorval = Math.Abs(sensor.GetSensorValue(j) - minSensorValue);
                double sensorval = sensor.GetSensorValue(j) - minSensorValue;
                SensorAfterModificationsValues[j] = sensorval;
                int currentLevel = GridData.Grid.Cells[j].RefinementLevel;
                if((sensorval / denom > 0.7 || sensorval / denom < 0.3) && currentLevel < maxRefinementLevel) {
                    levels[j] = 0;
                }
       
            }

            foreach(Chunk chunk in CellMask.GetFullMask(sensorField.GridDat)) {
                for(int i = 0; i < chunk.Len; i++) {
                    int cell = chunk.i0 + i;
                    Sensor_aftermodifications.SetMeanValue(cell, SensorAfterModificationsValues[i]);
                }
            }

            Tecplot.PlotFields(new DGField[] { sensorField, Sensor_aftermodifications }, "SensorDebugging", 0.0, 3);
            return levels;
        }
    }
 
    public class VarianceSensor : Sensor {
        private double[] sensorValues;

        public VarianceSensor() {
        }

        public override double GetSensorValue(int cellIndex) {
            return sensorValues[cellIndex];
        }

        public override void UpdateSensorValues(DGField fieldToTest_xdg) {
            var fieldToTest = ((XDGField)fieldToTest_xdg).GetSpeciesShadowField("A");
            int noOfCells = fieldToTest.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if(sensorValues == null || sensorValues.Length != noOfCells) {
                sensorValues = new double[noOfCells];
            }

            CellMask cellMask = CellMask.GetFullMask(fieldToTest.GridDat);
            foreach(int cell in cellMask.ItemEnum) {
                double constantCoordinate = 0.0;

                foreach(int coordinate in fieldToTest.Basis.GetPolynomialIndicesForDegree(cell, 0)) {
                    constantCoordinate += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                }
                double cellCoordinates = 0.0;
                for(int coordinate = 0; coordinate < fieldToTest.Basis.Length; coordinate++) {
                    cellCoordinates += fieldToTest.Coordinates[cell, coordinate] * fieldToTest.Coordinates[cell, coordinate];
                }

                double result = Math.Sqrt(cellCoordinates - constantCoordinate) / (Math.Sqrt(fieldToTest.GridDat.iLogicalCells.GetCellVolume(cell)));

                sensorValues[cell] = result;
            }
        }
    }

}