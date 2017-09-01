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

namespace CNS.ShockCapturing {

    public class SmoothedHeavisideArtificialViscosityLaw : IArtificialViscosityLaw {

        private IShockSensor sensor;
        
        private double dgDegree;

        private double sensorLimit;

        private double epsilon0;

        private double kappa;

        public SmoothedHeavisideArtificialViscosityLaw(IShockSensor sensor, int dgDegree, double sensorLimit, double epsilon0, double kappa) {
            this.sensor = sensor;
            this.dgDegree = dgDegree;
            this.sensorLimit = sensorLimit;
            this.epsilon0 = epsilon0;
            this.kappa = kappa;
        }
      
        public double GetViscosity(int jCell, double cellSize, StateVector state) {
            double s0 = Math.Log10(sensorLimit / (double)Math.Pow(dgDegree, 4));
            double se = Math.Log10(sensor.GetSensorValue(jCell) + 1e-15);

            double epsilonE;
            if (se < s0 - kappa) {
                epsilonE = 0.0;
            } else if (se > s0 + kappa) {
                epsilonE = epsilon0;
            } else {
                epsilonE = 0.5 * epsilon0 * (1.0 + Math.Sin(0.5 * Math.PI * (se - s0) / kappa));
            }
            
            double lambdaMax = state.SpeedOfSound + state.Velocity.Abs();
            double fudgeFactor = 0.5;   // Kloeckner (2011)
            epsilonE = fudgeFactor * epsilonE * lambdaMax * cellSize / dgDegree;
            return epsilonE;
        }
    }
}
