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


namespace BoSSS.Solution.CompressibleFlowCommon.MaterialProperty {

    /// <summary>
    /// Interface for an equation of state which allows for the calculation of
    /// some flow characteristics that depend on the equation of state. In
    /// particular, we limit ourselves to <i>incomplete</i> equations of state
    /// (cf. MenikoffPlohr1989) where the pressure is given by a function of
    /// the density and the inner energy, i.e.
    /// \f$ 
    /// p = p(\rho, e),
    /// \f$ 
    /// that can <b>always</b> be reformulated as
    /// \f$ 
    /// e = e(\rho, p).
    /// \f$ 
    /// Note that this does not hold for the density in general, which means
    /// that the density is not necessarily uniquely defined by given values
    /// for pressure and inner energy.
    /// </summary>
    public interface IEquationOfState {

        /// <summary>
        /// Ratio of specific heat capacities
        /// </summary>
        double HeatCapacityRatio {
            get;
        }

        /// <summary>
        /// Calculates the pressure base on the given state
        /// </summary>
        /// <param name="state">
        /// The flow state in some point of the domain.
        /// </param>
        /// <returns>The pressure</returns>
        double GetPressure(StateVector state);

        /// <summary>
        /// Calculates the temperature based on the given state.
        /// </summary>
        /// <param name="state">
        /// The flow state in some point of the domain.
        /// </param>
        /// <returns>The temperature</returns>
        double GetTemperature(StateVector state);

        /// <summary>
        /// Calculates the inner energy for a given <paramref name="density"/>
        /// and <paramref name="pressure"/>.
        /// </summary>
        /// <param name="density">
        /// The local fluid density
        /// </param>
        /// <param name="pressure">
        /// The local fluid pressure
        /// </param>
        /// <returns>
        /// The inner energy \f$ \rho e = \rho \cdot e(\rho, p)\f$ 
        /// </returns>
        double GetInnerEnergy(double density, double pressure);

        /// <summary>
        /// Calculates the local speed of sound based on the given state.
        /// </summary>
        /// <param name="state">
        /// The flow state in some point of the domain.
        /// </param>
        /// <returns>The speed of sound</returns>
        double GetSpeedOfSound(StateVector state);

        /// <summary>
        /// Calculates the local entropy based on the given state.
        /// </summary>
        /// <param name="state">
        /// The flow state in some point of the domain.
        /// </param>
        /// <returns>The entropy</returns>
        double GetEntropy(StateVector state);
    }
}
