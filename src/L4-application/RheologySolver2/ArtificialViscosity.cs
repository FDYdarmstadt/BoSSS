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


using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System;
using BoSSS.Foundation.Quadrature;
using ilPSP;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Artificial viscosity for shock capturing
    /// </summary>
    public static class ArtificialViscosity {

        /// <summary>
        /// calculating the value of the artificial viscosity
        /// </summary>
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

        /// <summary>
        /// project the value of the artificial viscosity onto a DG field
        /// </summary>
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