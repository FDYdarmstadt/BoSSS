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
using XESF.Fluxes;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.CodeAnalysis.CSharp.Syntax;

namespace XESF.Fluxes {

    public class RoeEnergyFlux : RoeBaseFlux{

        public RoeEnergyFlux(IBoundaryConditionMap boundaryMap, Material material, double s_alpha,int D) : base(boundaryMap, material, s_alpha,D) {
        }

        public double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {

            Vector n = new Vector(normal);
            double FL = Flux(Uin, D) * n;
            double FR = Flux(Uout, D) * n;

            (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) = Setup(normal, Uin, Uout, D);
            double vN = vAverage * (new Vector(normal));

            double[] Udiff = Uin.CloneAs();
            Udiff.AccV(-1, Uout);
            var result = V0_inv.MatVecMul(1.0, Udiff);

            // +++++++++++++++ V0 +++++++++++++++
            double[] V0 = new double[D + 2];
            V0[0] = enthalpyAverage - vN * speedOfSoundAverage;
            for(int i = 1; i < D + 1; i++) {
                V0[i] = vAverage[i - 1] + (vAverage.L2NormPow2() / 2 - vN) * normal[i - 1];
            }
            V0[D + 1] = enthalpyAverage + vN * speedOfSoundAverage;

            double k_j = 0;
            for(int i = 0; i < D + 2; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }

            //OLD stuff
            //(MultidimensionalArray V0_inv_old, MultidimensionalArray eigenVals_old,MultidimensionalArray V0_old)= OldSetup(normal, Uin, Uout, D); 
            //var prod = V0_old * eigenVals_old * V0_inv_old;
            //var k_4 = prod.MatVecMul(1.0, Udiff);
            //if(Math.Abs(k_4[D+1]-k_j)>1e-15) {
            //    Console.WriteLine("*************** Fehler im Code: Energy Flux ***************");
            //    var result_old=V0_inv_old.MatVecMul(1.0, Udiff);
            //}

            double energyFlux = 0.5 * (FL + FR) + 0.5 * k_j;

            //if(energyFlux.IsNaN()) {
            //    Console.WriteLine("*************** Fehler im Code: Energy Flux ***************"); 
            //}
            return energyFlux;


        }


        // Copied from HLLCEnergyFlux 
        //!!!!
        // xxx noch checken+überarbeiten !!!! ist noch im raw-Zustand
        //!!!!
        public Vector Flux(double[] U, int D) {
            Vector Output = new Vector(D);
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;


            double density = U[0];
            double energy = U[D + 1];

            double momentumSquared = 0.0;
            for(int d = 0; d < D; d++) {
                momentumSquared += U[d + 1] * U[d + 1];
            }
            double pressure = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density);

            //return state.Velocity * (state.Energy + state.Pressure);
            for(int d = 0; d < D; d++) {
                double flux = U[d + 1] / density * (energy + pressure);
                //Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                Output[d] += flux;
            }
            return Output;
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            // Copied from HLLCEnergyFlux
            double[] Uout = GetExternalState(ref inp, Uin);
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (Vin);
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            // Copied from HLLCEnergyFlux
            return InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D) * (_vIN - _vOUT);
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            // Copied from HLLCEnergyFlux
            // xxx am Ende auf Korrektheit überprüfen mit dem Code der noch implemntiert wird/wurde !!
            Vector fluxvector = Flux(U, cpv.D);
            return -1.0 * fluxvector * GradV;
        }
    }
}
