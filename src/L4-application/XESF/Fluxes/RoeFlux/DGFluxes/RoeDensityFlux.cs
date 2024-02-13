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

namespace XESF.Fluxes {
    /* INFO
            Implementation Guide Line => s. Toro p.376
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            See p.376:
            "   Before solving these equations we note that in the purely one–dimensional case
                v˜ = w˜ = 0, α˜3 = α˜4 = 0 (wave strength), K˜ (3) = K˜ (4) = 0 (11.67)           "
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      */
    

    public class RoeDensityFlux : RoeBaseFlux { 

        public RoeDensityFlux(IBoundaryConditionMap boundaryMap, Material material,double s_alpha, int m_D=2) : base(boundaryMap, material,s_alpha,m_D) {
        }

        // Method-Header copied from HLLCDensityFlux.cs
        public double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {

            Vector n = new Vector(normal);
            double FL = Flux(Uin, D) * n;
            double FR = Flux(Uout, D) * n;

            (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) = Setup(normal, Uin, Uout, D);
            double vN = vAverage * (new Vector(normal));
            double gamma = this.material.EquationOfState.HeatCapacityRatio;

            double[] Udiff = Uin.CloneAs();
            Udiff.AccV(-1, Uout);
            var result = V0_inv.MatVecMul(1.0, Udiff);

            // +++++++++++++++ V0 +++++++++++++++
            double[] V0 = new double[D+2];
            V0[0] = 1;
            for(int i = 0; i < D; i++) {
                V0[i+1] = normal[i];
            }
            V0[D + 1] = 1;

            double k_j = 0;
            for(int i = 0; i < D + 2; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }
            double densityFlux = 0.5 * (FL + FR) + 0.5 * k_j;
            return densityFlux;

        }


        public Vector Flux(double[] U, int D) {
            Vector fluxvector = new Vector(D);
            for(int d = 0; d < D; d++) {
                // state.Momentum
                double flux = U[d + 1];
                Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                fluxvector[d] = flux;
            }
            return fluxvector;
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            // Copied from HLLCDensityFlux
            double[] Uout = GetExternalState(ref inp, Uin);
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (Vin);
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            // Copied from HLLCDensityFlux
            return InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D) * (_vIN - _vOUT);
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            // Kinda copied from HLLCDensityFlux
            Vector fluxvector = Flux(U, cpv.D);
            return -1.0 * fluxvector * GradV;
        }

    }
}