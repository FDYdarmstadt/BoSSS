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

namespace BoSSS.Solution.CompressibleFlowCommon.MaterialProperty {

    /// <summary>
    /// Tait equation of state for stiff fluids (e.g. water; for example, see
    /// FarhatEtAl2008). The generic form of the equation is
    /// \f$ 
    /// p = \alpha \rho^\kappa - \pi
    /// \f$ 
    /// with given constants \f$ \alpha\f$ ,
    /// \f$ \kappa\f$  and
    /// \f$ \pi\f$ .
    /// </summary>
    public class TaitLaw : IEquationOfState {

        /// <summary>
        /// Model parameter \f$ \alpha\f$ 
        /// </summary>
        private readonly double alpha;

        /// <summary>
        /// Model parameter \f$ \pi\f$ 
        /// </summary>
        private readonly double pi;

        /// <summary>
        /// Creates the equation of state for a stiff fluid with the given
        /// model parameters
        /// </summary>
        /// <param name="heatCapacityRatio">
        /// The heat capacity ratio
        /// </param>
        /// <param name="alpha">
        /// Pressure scaling
        /// </param>
        /// <param name="pi">
        /// Pressure offset
        /// </param>
        public TaitLaw(double heatCapacityRatio, double alpha, double pi) {
            this.HeatCapacityRatio = heatCapacityRatio;
            this.alpha = alpha;
            this.pi = pi;
        }

        /// <summary>
        /// Creates the equation of state for a stiff fluid with the given heat
        /// capacity ratio, the reference density, the reference pressure and a
        /// a pressure constant that has been determined at the reference state.
        /// </summary>
        /// <param name="heatCapacityRatio">
        /// The heat capacity ratio
        /// </param>
        /// <param name="referenceDensity">
        /// The reference density at which measurements have been taken
        /// </param>
        /// <param name="referencePressure">
        /// The reference pressure at which measurements have been taken
        /// </param>
        /// <param name="pressureConstant">
        /// A pressure constant that determines the compressibility of the
        /// fluid
        /// </param>
        /// <remarks>
        /// The input parameters typically strongly vary in size, so there is
        /// guarantee about the precision of the calculated constants.
        /// However, this should hardly matter because the measurement error is
        /// probably much higher anyway...
        /// </remarks>
        public TaitLaw(double heatCapacityRatio, double referenceDensity, double referencePressure, double pressureConstant) {
            this.HeatCapacityRatio = heatCapacityRatio;
            this.alpha = (referencePressure + pressureConstant / heatCapacityRatio) / Math.Pow(referenceDensity, heatCapacityRatio);
            this.pi = pressureConstant / heatCapacityRatio;
        }

        /// <summary>
        /// Instantiates the specific equation of state for water
        /// </summary>
        public static TaitLaw Water {
            get {
                return new TaitLaw(7.15, 1E+5, 3.31E+8, 1E+3);
            }
        }

        #region IEquationOfState Members

        /// <summary>
        /// The heat capacity ratio.
        /// </summary>
        public double HeatCapacityRatio {
            get;
            private set;
        }

        /// <summary>
        /// Compute the pressure from the given density.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetPressure"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// p = \alpha \rho^\kappa - \pi
        /// \f$ 
        /// </returns>
        public double GetPressure(StateVector state) {
            return alpha * Math.Pow(state.Density, HeatCapacityRatio) + pi;
        }

        /// <summary>
        /// Calculates the temperature \f$ T\f$  via
        /// \f$ T = (\kappa - 1) e\f$ .
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetTemperature"/>
        /// </param>
        /// <returns>
        /// \f$ (\kappa - 1) e\f$ 
        /// where
        /// \f$ e\f$  = <see name="StateVector.InnerEnergy"/>
        /// and \f$ \kappa\f$  is the heat capacity ratio.
        /// </returns>
        /// <remarks>
        /// While the inner energy of a Tait fluid is in principle not defined,
        /// we can make some assumptions that are explained in
        /// <see cref="GetInnerEnergy"/>
        /// </remarks>
        public double GetTemperature(StateVector state) {
            return (HeatCapacityRatio - 1.0) * state.SpecificInnerEnergy;
        }

        /// <summary>
        /// Calculates the inner energy for a Tait fluid under the assumption
        /// that it models an isentropic stiffened gas.
        /// </summary>
        /// <param name="density"></param>
        /// <param name="pressure"></param>
        /// <returns>
        /// \f$ \frac{p + \kappa \pi}{(\kappa - 1)}\f$  where
        /// \f$ \kappa\f$  is the heat capacity ratio.
        /// </returns>
        /// <remark>
        /// In principle, the inner energy cannot be calculated from the Tait
        /// law itself since it is, in some sense, incomplete. However, a Tait
        /// fluid can be interpreted as the special case of an
        /// <i>isentropic</i> stiffened gas, which implies that the inner
        /// energy can be computed from the stiffened gas law (cf.
        /// FarhatEtAl2008 and Mueller2014).
        /// </remark>
        public double GetInnerEnergy(double density, double pressure) {
            return (pressure + HeatCapacityRatio * pi) / (HeatCapacityRatio - 1.0);
        }

        /// <summary>
        /// Calculates the local speed of sound from the given density. For
        /// details, e.g. see FarhatEtAl2008.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// c = \sqrt{\gamma \alpha \rho^{gamma-1}}
        /// \f$ 
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            return Math.Sqrt(HeatCapacityRatio * alpha * Math.Pow(state.Density, HeatCapacityRatio - 1));
        }

        /// <summary>
        /// Always returns 1 since the equation of state dictates isentropic
        /// state changes, see FarhatEtAl2008
        /// </summary>
        /// <param name="state"></param>
        /// <returns></returns>
        public double GetEntropy(StateVector state) {
            return 1.0;
        }

        #endregion
    }
}
