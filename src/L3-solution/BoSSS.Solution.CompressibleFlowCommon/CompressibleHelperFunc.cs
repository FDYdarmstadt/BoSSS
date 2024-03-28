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

using MathNet.Numerics.Distributions;
using System;
using static BoSSS.Solution.CompressibleFlowCommon.Variable;

namespace BoSSS.Solution.CompressibleFlowCommon {
    public static class CompressibleHelperFunc
    {
        public static double ComputeVelocity(double rho, double p,double M,double gamma)
        {
            return M* Math.Sqrt(gamma * p / rho);
        }
        /// <summary>
        /// computes all post shock values from values left of the shock
        /// </summary>
        /// all left values and gamma
        /// <param name="rhoL"></param>
        /// <param name="uL"></param>
        /// <param name="pL"></param>
        /// <param name="ML"></param>
        /// <param name="gamma"></param>
        /// <returns></returns>
        public static (double rhoR, double uR,double pR,double cR,double MR) ComputeNormalShockWaveRelations(double rhoL, double uL, double pL, double ML,double gamma)
        {
            double rhoR = (gamma + 1) * ML * ML / (2 + (gamma - 1) * ML * ML) * rhoL;
            double pR = (1 + 2 * gamma / (gamma + 1) * (ML * ML - 1)) * pL;
            double uR = (2 + (gamma - 1) * ML * ML) / ((gamma + 1) * ML * ML) * uL;    // (1)
            double cR= Math.Sqrt(gamma * pR / rhoR);
            double MR = uR/cR;
            return ( rhoR, uR, pR,cR, MR);
        }
    }
}
