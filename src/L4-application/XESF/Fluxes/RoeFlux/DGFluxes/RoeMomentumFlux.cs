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
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XESF.Fluxes {

    public class RoeMomentumFlux : RoeBaseFlux {

        private int component;
        public RoeMomentumFlux(IBoundaryConditionMap boundaryMap, int component, Material material, double s_alpha,int D=2) : base(boundaryMap, material, s_alpha,D) {
            this.component = component;
        }

        public double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {

            //(double[] waveStrength, double[] eigenVal, double[,] eigenVec) = SetupXSplit(normal, Uin, Uout, D);
            
            Vector n = new Vector(normal);
            double FL = Flux(Uin, D)*n;
            double FR = Flux(Uout, D)*n;

            (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) = Setup(normal, Uin, Uout, D);

            double[] Udiff = Uin.CloneAs();
            Udiff.AccV(-1, Uout);
            var result = V0_inv.MatVecMul(1.0, Udiff);

            // +++++++++++++++ V0 +++++++++++++++
            double[] V0 = new double[D + 2];
            V0[0] = vAverage[component] - speedOfSoundAverage * normal[component];
            for(int i = 1; i < D+1; i++) {
                if(component == i-1) {
                    V0[i] = (vAverage[component] - normal[component]) * normal[i-1] + 1;
                } else {
                    V0[i] = (vAverage[component] - normal[component]) * normal[i-1];
                }
            }
            V0[D + 1] = vAverage[component] + speedOfSoundAverage * normal[component];
            // +++++++++++++++ end +++++++++++++++
   
            double k_j = 0;
            for(int i = 0; i < D + 2; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }


            double momentumFlux = 0.5 * (FL + FR) + 0.5 * k_j;

            return momentumFlux;

        }

        public Vector Flux(double[] U, int D) {
            Vector Output = new Vector(D);
            int L = U.Length;
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;

            double density = U[0];
            double energy = U[D + 1];

            double momentumSquared = 0.0;
            for(int d = 0; d < D; d++) {
                Output[d] += U[component + 1] * U[d + 1] / density;
                momentumSquared += U[d + 1] * U[d + 1];
            }

            double flux = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density) / gammaMachSquared;
            //Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
            Output[component] += flux;
            return Output;
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            // Copied from HLLCMomentumFlux
            double[] Uout = GetExternalState(ref inp, Uin);
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (Vin);
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            // Copied from HLLCMomentumFlux
            return InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D) * (_vIN - _vOUT);
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector fluxvector = Flux(U, cpv.D);
            return -1.0 * fluxvector * GradV;
        }
    }
}
