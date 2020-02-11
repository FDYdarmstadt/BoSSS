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
using System.Diagnostics;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Lambda for Lax-Friedrichs flux-
    /// </summary>
    static public class LambdaConvection {

        /// <summary>
        /// Generic Lambda for convective operators of form \nabla(U Scalar)
        /// </summary>
        /// <param name="ScalarMean"></param>
        /// <param name="VelocityMean"></param>
        /// <param name="Normal"></param>
        /// <param name="FactorTwo">
        /// True: Lambda = 2.0 * rho * V_n (used for momentum equation).
        /// False: Lambda = rho * V_n (used for temperature equation).
        /// </param>
        /// <returns>        
        /// </returns>
        static public double GetLambda(double[] VelocityMean, double[] Normal, bool FactorTwo, double ScalarMean = 1) {
            Debug.Assert(VelocityMean.Length == Normal.Length, "Mismatch in dimensions!");

            double V_n = 0.0;
            for (int d = 0; d < VelocityMean.Length; d++)
                V_n += VelocityMean[d] * Normal[d];

            double Lambda = ScalarMean * V_n;
            if (FactorTwo)
                Lambda *= 2.0;

            if (double.IsNaN(Lambda) || double.IsInfinity(Lambda))
                throw new NotFiniteNumberException();

            return Math.Abs(Lambda);
        }

        /// <summary>
        /// Lambda for convective operators with variable density,
        /// i.e. momentum equation for multiphase and Low-Mach flows
        /// and temperature advection for Low-Mach flows.
        /// </summary>
        /// <param name="ScalarMean"></param>
        /// <param name="VelocityMean"></param>
        /// <param name="Normal"></param>
        /// <param name="EoS"></param>
        /// <param name="FactorTwo">
        /// True: Lambda = 2.0 * rho * V_n (used for momentum equation).
        /// False: Lambda = rho * V_n (used for temperature equation).
        /// </param>
        /// <returns>        
        /// </returns>
        static public double GetLambda(double[] VelocityMean, double[] Normal, MaterialLaw EoS, bool FactorTwo, params double[] ScalarMean) {
            double rhoMean = EoS.GetDensity(ScalarMean);
            return GetLambda(VelocityMean, Normal, FactorTwo, rhoMean);
        }

    }
}
