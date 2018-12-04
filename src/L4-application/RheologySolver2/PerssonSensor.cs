using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Application.Rheology {

    public class PerssonSensor {

        private string sensorVariable;

        private ConventionalDGField sensorField;

        public ConventionalDGField fieldToTestRestricted;

        public PerssonSensor(DGField fieldToTest) {
            sensorVariable = fieldToTest.Identification;
            sensorField = new SinglePhaseField(
                new Basis(fieldToTest.GridDat, 0), "sensor");
            fieldToTestRestricted = new SinglePhaseField(
                new Basis(fieldToTest.GridDat, fieldToTest.Basis.Degree - 1));
        }

        public double GetValue(int cellIndex) {
            return sensorField.GetMeanValue(cellIndex);
        }

        public DGField GetField() {
            return sensorField;
        }

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
}