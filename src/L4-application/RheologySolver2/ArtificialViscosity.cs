using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System;
using BoSSS.Foundation.Quadrature;
using ilPSP;

namespace BoSSS.Application.Rheology {

    public static class ArtificialViscosity {

        public static double GetViscosity(int jCell, double cellSize, int dgDegree, double perssonSensor, double sensorLimit, double maxViscosity) {

            double sensorValue = Math.Log10(perssonSensor + 1e-15);
            double limit = Math.Log10(sensorLimit/Math.Pow(dgDegree,4));

            double epsilonE;
            double kappa = 1.0;
            double refMaxViscosity = maxViscosity;

            if (sensorValue < limit - kappa) {
                epsilonE = 0.0;
            } else if (sensorValue > limit + kappa) {
                epsilonE = refMaxViscosity;
            } else {
                epsilonE = 0.5 * refMaxViscosity * (1.0 + Math.Sin(0.5 * Math.PI * (sensorValue - limit) / kappa));
            }

            epsilonE = epsilonE * cellSize / dgDegree;

            return epsilonE;
        }


        public static void ProjectArtificalViscosityToDGField(SinglePhaseField avField, PerssonSensor sensor, double sensorLimit, double maxViscosity, CellMask cellMask = null) {

            MultidimensionalArray h_min = avField.GridDat.iGeomCells.h_min;
            int p = sensor.fieldToTestRestricted.Basis.Degree + 1;

            if (cellMask == null) {
                cellMask = CellMask.GetFullMask(avField.GridDat);
            }

            avField.Clear();
            foreach (int cell in cellMask.ItemEnum) {

                double localViscosity = GetViscosity(cell, h_min[cell], p, sensor.GetValue(cell), sensorLimit, maxViscosity);

                avField.SetMeanValue(cell, localViscosity);
            }
        }

    }
}