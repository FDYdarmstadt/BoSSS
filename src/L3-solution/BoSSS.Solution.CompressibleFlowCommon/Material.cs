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
    /// Summarizes all material properties of a fluid. At the moment, these are
    /// the equation of state and the viscosity law. Other may follow later
    /// (Thermal conductivity? Stress tensor? Who knows...)
    /// </summary>
    public class Material {

        /// <summary>
        /// Constructs a new material.
        /// </summary>
        public Material(IEquationOfState __EquationOfState, IViscosityLaw __ViscosityLaw, double __MachNumber, double __ReynoldsNumber, double __PrandtlNumber, double __FroudeNumber, double __ViscosityRatio) {
            //this.Control = control;
            this.EquationOfState = __EquationOfState;
            this.ViscosityLaw = __ViscosityLaw;
            MachNumber = __MachNumber;
            ReynoldsNumber = __ReynoldsNumber;
            FroudeNumber = __FroudeNumber;
            PrandtlNumber = __PrandtlNumber;
            ViscosityRatio = __ViscosityRatio;
        }


        /// <summary>
        /// The configured Mach Number in the far field.
        /// </summary>
        public double MachNumber { get; private set; }

        /// <summary>
        /// The configured Reynolds number in the far field, ignored if Euler equations are solved.
        /// </summary>
        public double ReynoldsNumber { get; private set; }

        /// <summary>
        /// The configured Prandtl number in the far field, ignored if Euler equations are solved;
        /// </summary>
        public double PrandtlNumber { get; private set; }

        /// <summary>
        /// The ratio of a characteristic flow velocity to the velocity of a gravitational wave.
        /// </summary>
        public double FroudeNumber { get; private set; }

        /// <summary>
        /// The ratio of the bulk viscosity to the shear viscosity.
        /// </summary>
        public double ViscosityRatio { get; private set; }


        ///// <summary>
        ///// Associated CNScontrol file to generate the material
        ///// </summary>
        //public readonly CNSControl Control;

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
