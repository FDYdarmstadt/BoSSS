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
using BoSSS.Foundation;
using BoSSS.Solution;
using CNS.Exception;

namespace CNS.TimeStepping {

    /// <summary>
    /// Special time stepper for semi-implicit schemes. Aggregates an explicit
    /// and an implicit timestepper and advances both at once.
    /// </summary>
    public class SemiImplicitTimeStepper : ITimeStepper {

        /// <summary>
        /// The explicit part
        /// </summary>
        private ITimeStepper explicitTimeStepper;

        /// <summary>
        /// The implicit part
        /// </summary>
        private ITimeStepper implicitTimeStepper;

        /// <summary>
        /// Constructs a new semi-implicit timestepper by aggregating
        /// <paramref name="explicitTimeStepper"/> and
        /// <paramref name="implicitTimeStepper"/>. Note that both
        /// timesteppers are required to have the same
        /// <see cref="CoordinateMapping"/>.
        /// </summary>
        /// <param name="explicitTimeStepper">The explicit part</param>
        /// <param name="implicitTimeStepper">The implicit part</param>
        public SemiImplicitTimeStepper(ITimeStepper explicitTimeStepper, ITimeStepper implicitTimeStepper) {
            if (explicitTimeStepper.Mapping != implicitTimeStepper.Mapping) {
                throw new ArgumentException("Both time steppers currently need to have the same coordinate mapping", "implicitTimeStepper");
            }

            this.explicitTimeStepper = explicitTimeStepper;
            this.implicitTimeStepper = implicitTimeStepper;
        }

        #region ITimeStepper Members

        /// <summary>
        /// Returns the current time which has to be equal for both parts
        /// </summary>
        public double Time {
            get {
                double time = explicitTimeStepper.Time;
                if (implicitTimeStepper.Time != time) {
                    throw new InternalErrorException("Explicit and implicit part are at different time levels. This should not have happened");
                }
                return time;
            }
        }

        /// <summary>
        /// Resets both timesteppers
        /// </summary>
        /// <param name="NewTime"><see cref="ITimeStepper"/></param>
        public void ResetTime(double NewTime) {
            explicitTimeStepper.ResetTime(NewTime);
            implicitTimeStepper.ResetTime(NewTime);
        }

        /// <summary>
        /// First performs the explicit part, then then the implicit part.
        /// </summary>
        /// <param name="dt"></param>
        public double Perform(double dt) {
            explicitTimeStepper.Perform(dt);
            implicitTimeStepper.Perform(dt);
            return dt;
        }

        /// <summary>
        /// Returns the mapping which should be identical for both parts.
        /// </summary>
        public CoordinateMapping Mapping {
            get {
                // always equal to implicitTimeStepper.Mapping
                return explicitTimeStepper.Mapping;
            }
        }

        /// <summary>
        /// Returns the DG coordinates which should be identical for both
        /// parts.
        /// </summary>
        public CoordinateVector DGCoordinates {
            get {
                // always equal to implicitTimeStepper.DGCoordinates
                return implicitTimeStepper.DGCoordinates;
            }
        }

        #endregion
    }
}
