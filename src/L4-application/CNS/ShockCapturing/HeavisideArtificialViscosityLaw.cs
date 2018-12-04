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

using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using System;

namespace CNS.ShockCapturing {

    /// <summary>
    /// Artificial viscosity law where artificial viscosity is either fully on
    /// or off, without a smooth transitioning
    /// </summary>
    public class HeavisideArtificialViscosityLaw : IArtificialViscosityLaw {

        private IShockSensor sensor;

        private double dgDegree;

        private double sensorLimit;

        private double refViscosity;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="sensor">
        /// A shock sensor
        /// </param>
        /// <param name="dgDegree">
        /// Approximation degree of the field to which the given
        /// <paramref name="sensor"/> is applied
        /// </param>
        /// <param name="refSensorLimit">
        /// Base value of the hard sensor limit (before scaling) that
        /// distinguishes between AV or no AV
        /// </param>
        /// <param name="refViscosity">
        /// Base value of the viscosity (before scaling) that should be used in
        /// shocked cells.
        /// </param>
        public HeavisideArtificialViscosityLaw(IShockSensor sensor, int dgDegree, double refSensorLimit, double refViscosity) {
            this.sensor = sensor;
            this.sensorLimit = refSensorLimit / (double)Math.Pow(dgDegree, 4);
            this.refViscosity = refViscosity;
            this.dgDegree = dgDegree;
        }

        /// <summary>
        /// The viscosity to be used in the given cell
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="cellSize"></param>
        /// <param name="state"></param>
        /// <returns></returns>
        public double GetViscosity(int jCell, double cellSize, StateVector state) {
            if (sensor.GetSensorValue(jCell) > sensorLimit) {
                return refViscosity * cellSize / dgDegree;
            } else {
                return 0.0;
            }
        }

        /// <summary>
        /// True, if <see cref="GetViscosity(int, double, StateVector)"/>
        /// returns non-zero value for the given cell
        /// </summary>
        /// <param name="jCell"></param>
        /// <returns></returns>
        public bool IsShocked(int jCell) {
            return (sensor.GetSensorValue(jCell) > sensorLimit);
        }
    }
}
