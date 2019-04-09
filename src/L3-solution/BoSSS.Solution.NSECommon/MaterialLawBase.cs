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

using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Base class for implementing different material laws.
    /// </summary>
    public abstract class MaterialLaw {

        /// <summary>
        /// Returns density dependent on some scalar,
        /// e.g. temperature or level-set.
        /// </summary>
        /// <param name="phi">scalar</param>
        /// <returns>density</returns>
        public abstract double GetDensity(params double[] phi);

        /// <summary>
        /// Returns scalar dependent density.
        /// Implements <see cref="BoSSS.Foundation.Func"/>.
        /// </summary>
        /// <param name="X">Not used.</param>
        /// <param name="phi">scalar</param>
        /// <param name="jCell">Not used.</param>
        /// <returns>Density</returns>
        public double GetDensity(double[] X, double[] phi, int jCell) {
            return this.GetDensity(phi);
        }

        /// <summary>
        /// Returns dynamic viscosity dependent on some scalar,
        /// e.g. temperature or level-set.
        /// </summary>
        /// <param name="phi">scalar</param>
        /// <returns>dynamic viscosity</returns>
        public abstract double GetViscosity(double phi);

        /// <summary>
        /// Returns scalar dependent dynamic viscosity.
        /// Implements <see cref="BoSSS.Foundation.Func"/>.
        /// </summary>
        /// <param name="X">Not used.</param>
        /// <param name="phi">scalar</param>
        /// <param name="jCell">Not used.</param>
        /// <returns>Dynamic viscosity</returns>
        public double GetViscosity(double[] X, double[] phi, int jCell) {
            if (phi.Length != 1)
                throw new ArgumentException();

            return this.GetViscosity(phi[0]);
        }


        /// <summary>
        /// Returns thermodynamic pressure as function of inital mass and temperature.
        /// </summary>
        /// <param name="InitialMass"></param>
        /// <param name="Temperature"></param>
        /// <returns></returns>
        public abstract double GetMassDeterminedThermodynamicPressure(double InitialMass, SinglePhaseField Temperature);


        #region CoupledLaxFriedrichs

        /// <summary>
        /// Returns lambda, i.e. penalty parameter for coupled Lax-Friedrichs flux.
        /// Not to be confused with heat conductivity!!!
        /// </summary>
        /// <param name="VelocityMean"></param>
        /// <param name="Normal"></param>
        /// <param name="ScalarMean"></param>
        /// <returns>Absolute value of lambda.</returns>
        public abstract double GetLambda(double[] VelocityMean, double[] Normal, double ScalarMean);

        /// <summary>
        /// Derivative of density with respect to scalar.
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public abstract double DiffRho_Temp(double phi);

        #endregion

        /// <summary>
        /// Paramaters for <see cref="BoSSS.Foundation.IEquationComponent.ParameterOrdering"/>
        /// </summary>
        public abstract IList<string> ParameterOrdering {
            get;
        }
    }
}
