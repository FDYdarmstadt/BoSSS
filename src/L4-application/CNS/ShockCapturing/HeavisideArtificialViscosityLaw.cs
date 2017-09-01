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

using System;
using System.Linq;

namespace CNS.ShockCapturing {

    public class HeavisideArtificialViscosityLaw : IArtificialViscosityLaw {

        private IShockSensor sensor;

        private double cellSize;

        private double dgDegree = 0;

        private double sensorLimit;

        private double epsilon0;

        public HeavisideArtificialViscosityLaw(IShockSensor sensor, double cellSize, int dgDegree, double sensorLimit, double epsilon0) {
            this.sensor = sensor;
            this.cellSize = cellSize;
            this.sensorLimit = sensorLimit;
            this.epsilon0 = epsilon0;
        }

        public double GetViscosity(int jCell, double cellSize, StateVector state) {
            double S0 = sensorLimit / Math.Pow(dgDegree, 4);
            double Se = sensor.GetSensorValue(jCell);

            if (Se > S0) {
                return epsilon0 * cellSize / dgDegree;
            } else {
                return 0.0;
            }
        }
    }
}
