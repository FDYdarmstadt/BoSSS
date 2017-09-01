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
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Material law for multiphase flows.
    /// </summary>
    public class MaterialLawMultiphase : MaterialLaw {

        double rho1;
        double rho2;
        double alpha;
        double mu1;
        double mu2;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="rho1"></param>
        /// <param name="rho2"></param>
        /// <param name="mu1"></param>
        /// <param name="mu2"></param>
        public MaterialLawMultiphase(double rho1, double rho2, double mu1, double mu2) {
            this.rho1 = rho1;
            this.rho2 = rho2;
            this.alpha = rho2 / rho1;
            this.mu1 = mu1;
            this.mu2 = mu2;
        }

        /// <summary>
        /// Returns the density as a function of the level-set.
        /// </summary>
        /// <param name="phi">Level-set</param>
        /// <returns></returns>
        public override double GetDensity(params double[] phi) {
            double density = rho1 * phi[0] + rho2 * (1 - phi[0]);
            return ClipProperty(phi[0], density, rho1, rho2);
        }

        /// <summary>
        /// Returns the viscosity as a function of the level-set.
        /// </summary>
        /// <param name="phi">Level-set</param>
        /// <returns></returns>
        public override double GetViscosity(double phi) {
            double mu = mu1 * phi + mu2 * (1 - phi);
            return ClipProperty(phi, mu, mu1, mu2);
        }

        /// <summary>
        /// Clips some phyisical property <paramref name="prop"/>
        /// for unphysical values of <paramref name="phi"/>,
        /// i.e. greater than 1.0 or less than 0.0
        /// </summary>
        /// <param name="phi">Level-Set value</param>
        /// <param name="prop">Value of the property, which should be clipped</param>
        /// <param name="prop1">Value of the property for <paramref name="phi"/> = 1.0</param>
        /// <param name="prop2">Value of the property for <paramref name="phi"/> = 0.0</param>
        /// <returns>
        /// The clipped value.
        /// </returns>
        private double ClipProperty(double phi, double prop, double prop1, double prop2) {
            if ((phi > 0.0) && (phi < 1.0)) {
                return prop;
            } else if (phi >= 1.0) {
                return prop1;
            } else if (phi <= 0.0) {
                return prop2;
            } else {
                throw new ArgumentException();
            }
        }

        public override double GetLambda(double[] VelocityMean, double[] Normal, double ScalarMean) {
            throw new NotImplementedException();

            //Debug.Assert(VelocityMean.Length == Normal.Length, "Mismatch in dimensions!");            

            //double V_n = 0.0;
            //for (int d = 0; d < VelocityMean.Length; d++)
            //    V_n += VelocityMean[d] * Normal[d];

            //double Lambda;            

            //switch (MixLaw) {
            //    case MixingLaw.Linear:
            //        Lambda = V_n * (3.0 * ScalarMean * (1.0 - alpha) - 2.0) / (2.0 * ScalarMean * (1.0 - alpha) - 1.0);
            //        break;
            //    case MixingLaw.NonLinear:
            //        Lambda = V_n * (2.0 * alpha + ScalarMean * (1.0 - alpha)) / alpha;
            //        break;
            //    default:
            //        throw new ArgumentException();
            //}

            //return Math.Abs(Lambda);
        }

        public override double DiffRho_Temp(double phi) {
            throw new NotImplementedException();

            //double DrhoDT;

            //switch (MixLaw) {
            //    case MixingLaw.Linear:
            //        return DrhoDT = rho2 - rho1;
            //    case MixingLaw.NonLinear:
            //        return DrhoDT = rho1 * rho2 * (rho2 - rho1) / Math.Pow(rho1 * phi + rho2 * (1.0 - phi), 2.0);
            //    default:
            //        throw new ArgumentException();
            //}
        }
    }
}
