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


namespace CNS.MaterialProperty {

    /// <summary>
    /// Summarizes all material properties of a fluid. At the moment, these are
    /// the equation of state and the viscosity law. Other may follow later
    /// (Thermal conductivity? Stress tensor? Who knows...)
    /// </summary>
    public class Material {

        /// <summary>
        /// Constructs a new material.
        /// </summary>
        /// <param name="control">
        /// <see cref="CNSControl"/>
        /// </param>
        public Material(CNSControl control) {
            this.Control = control;
            this.EquationOfState = control.EquationOfState;
            this.ViscosityLaw = control.ViscosityLaw;
        }

        /// <summary>
        /// Associated CNScontrol file to generate the material
        /// </summary>
        public readonly CNSControl Control;

        /// <summary>
        /// The equation of state linking pressure, density and inner energy
        /// </summary>
        public IEquationOfState EquationOfState {
            get;
            private set;
        }

        /// <summary>
        /// The dimensionless viscosity as a function of the dimensionless
        /// temperature.
        /// </summary>
        public IViscosityLaw ViscosityLaw {
            get;
            private set;
        }
    }
}
