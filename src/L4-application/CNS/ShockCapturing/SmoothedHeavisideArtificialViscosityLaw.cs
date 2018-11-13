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
    /// An artificial viscosity law with a smooth transitioning between full
    /// AV (sensor above upper threshold) cells and AV (sensor below lower
    /// threshold).
    /// </summary>
    public class SmoothedHeavisideArtificialViscosityLaw : IArtificialViscosityLaw {

        private IShockSensor sensor;

        private double dgDegree;

        private double sensorLimit;

        private double refMaxViscosity;

        private double kappa;

        private double? lambdaMax;

        private double? fudgeFactor;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="sensor">
        /// Some shock sensor
        /// </param>
        /// <param name="dgDegree">
        /// Approximation degree of the field to which the given
        /// <paramref name="sensor"/> is applied
        /// </param>
        /// <param name="refSensorLimit">
        /// Base value of the sensor limit (before scaling) that defines the
        /// amount of AV to be used. More precisely, this values defines the
        /// center of smooth Heaviside approximation (i.e., exactly half of the
        /// specified maximum <paramref name="refMaxViscosity"/> will be
        /// applied for the given sensor value)
        /// </param>
        /// <param name="refMaxViscosity">
        /// Base value of the maximum viscosity (before scaling) that should be
        /// used in shocked cells.
        /// </param>
        /// <param name="kappa">
        /// Width (in terms of "range of sensor values") of the smoothed
        /// transition zone between full AV and no AV
        /// </param>
        /// <param name="lambdaMax">
        /// Optional: Additional scaling coefficient for the AV. If none is
        /// given, an estimate will be used according the local flow state
        /// </param>
        /// <param name="fudgeFactor">
        /// Correction factor, typically set to 0.5 for the compressible Euler equations (Kloeckner et al. 2011)
        /// </param>
        public SmoothedHeavisideArtificialViscosityLaw(IShockSensor sensor, int dgDegree, double refSensorLimit, double refMaxViscosity, double kappa, double? lambdaMax = null, double? fudgeFactor = null) {
            this.sensor = sensor;
            this.dgDegree = dgDegree;
            this.sensorLimit = Math.Log10(refSensorLimit / (double)Math.Pow(dgDegree, 4));
            this.refMaxViscosity = refMaxViscosity;
            this.kappa = kappa;
            this.lambdaMax = lambdaMax;
            this.fudgeFactor = fudgeFactor;
        }

        /// <summary>
        /// Returns AV acoording to the configured sensor using a smoothed
        /// (sine wave) AV activation profile.
        /// </summary>
        /// <param name="jCell">
        /// The considered cell
        /// </param>
        /// <param name="cellSize">
        /// A measure for the size of the cell
        /// </param>
        /// <param name="state">
        /// The local flow state
        /// </param>
        /// <returns>
        /// For some sensor value \f$ s \f$, the artificial viscosity
        /// \f$ \nu \f$ is calculated as
        /// \f[ 
        ///     \nu = \nu_0 \begin{cases}
        ///         0     & \text{ if } s \lt s_0 - \kappa\\
        ///         1 & \text{ if } s \gt s_0 + \kappa\\
        ///         \frac{1}{2} \left(1 + \sin(\frac{1}{2 \kappa} \pi (s - s_0) \right)
        ///     \end{cases}
        /// \f]
        /// </returns>
        public double GetViscosity(int jCell, double cellSize, StateVector state) {
            double sensorValue = Math.Log10(sensor.GetSensorValue(jCell) + 1e-15);

            double epsilonE;
            if (sensorValue < sensorLimit - kappa) {
                epsilonE = 0.0;
            } else if (sensorValue > sensorLimit + kappa) {
                epsilonE = refMaxViscosity;
            } else {
                epsilonE = 0.5 * refMaxViscosity * (1.0 + Math.Sin(0.5 * Math.PI * (sensorValue - sensorLimit) / kappa));
            }

            double lambdaMax = this.lambdaMax ?? state.SpeedOfSound + state.Velocity.Abs();
            //lambdaMax = 20; //DMR
            //lambdaMax = 2; //Shock Tube

            double fudgeFactor = this.fudgeFactor ?? 0.5;   // Kloeckner (2011)

            epsilonE = fudgeFactor * epsilonE * lambdaMax * cellSize / dgDegree;

            return epsilonE;
        }

        /// <summary>
        /// Returns true of AV is non-zero, i.e. if sensor value is lower than
        /// the lower AV threshold used in 
        /// <see cref="GetViscosity(int, double, StateVector)"/>
        /// </summary>
        /// <param name="jCell"></param>
        /// <returns></returns>
        public bool IsShocked(int jCell) {
            double se = Math.Log10(sensor.GetSensorValue(jCell) + 1e-15);
            return (se >= sensorLimit - kappa);
        }
    }
}
