/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer StRoeSTmungsdynamik (chair of fluid dynamics)

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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics.Mcmc;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using XESF.Fluxes;

namespace XESTSF.Fluxes {
    /* INFO
            Implementation Guide Luine => s. Toro p.376
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            See p.376:
            "   Before solving these equations we note that in the purely one–dimensional case
                v˜ = w˜ = 0, α˜3 = α˜4 = 0 (wavestrength), K˜ (3) = K˜ (4) = 0 (11.67)           "
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      */
    

    public class RoeSTDensityFlux : RoeSTBaseFlux { 

        public RoeSTDensityFlux(IBoundaryConditionMap boundaryMap, Material material,double s_alpha,DGField[] _previous_u, int m_D=2) : base(boundaryMap, material,s_alpha, _previous_u,m_D) {
        }

        // Method-Header copied from HLLCDensityFlux.cs
        public override double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {

            //seperate space-time normal
            var normal_x = new double[normal.Length - 1];
            for(int i = 0; i < normal_x.Length; i++) {
                normal_x[i] = normal[i];
            }
            double n_t = normal[normal.Length - 1];


            // +++++++++++ Compute Roe stabilization term k_j ++++++++++
            // obtain Roe Flux quantities
            (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) = Setup(normal_x, Uin, Uout, D-1,n_t);

            //[[U]] =Uin - Uout
            double[] Ujmp = Uin.CloneAs();
            Ujmp.AccV(-1, Uout);

            //V_0^{-1}*[[U]]
            var result = V0_inv.MatVecMul(1.0, Ujmp);

            // Assemble row of Roe matrix corresponding to Density
            double[] V0 = new double[D+1];
            V0[0] = 1;
            for(int i = 0; i < D; i++) {
                V0[i+1] = normal[i];
            }
            V0[D] = 1;

            // V0*EV*(V_0^{-1}*[[U]]) - EV = matrix with eigenvalues
            double k_j = 0;
            for(int i = 0; i < D +1; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }

            // +++++++++++ Compute left and right fluxes FL,FR ++++++++++
            Vector n = new Vector(normal);
            double FL = Flux(Uin, D) * n;
            double FR = Flux(Uout, D) * n;

            // ++++++++++ Add all together +++++++++++++++
            double densityFlux = 0.5 * (FL + FR) + 0.5 * k_j;

            //Debug.Assert(densityFlux > 0);
            //Debug.Assert(!densityFlux.IsNaN());
            return densityFlux;

        }

        /// <summary>
        /// Space-time density flux
        /// </summary>
        /// <param name="U"></param>
        /// <param name="D"></param>
        /// <returns></returns>
        public override Vector Flux(double[] U, int D) {
            // ++++++++++++ space component ++++++++++++
            Vector fluxvector = new Vector(D);
            for(int d = 0; d < D-1; d++) {
                // state.Momentum
                double flux = U[d + 1];
                Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                fluxvector[d] = flux;
            }
            // ++++++++++++ time component ++++++++++++
            fluxvector[D - 1] = U[0];
            return fluxvector;
        }


    }
}