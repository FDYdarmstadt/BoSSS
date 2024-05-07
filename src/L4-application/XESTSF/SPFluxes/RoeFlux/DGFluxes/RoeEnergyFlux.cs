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
using XESTSF.Fluxes;
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
using XESF.Fluxes;

namespace XESTSF.Fluxes {

    public class RoeSTEnergyFlux : RoeSTBaseFlux{

        public RoeSTEnergyFlux(IBoundaryConditionMap boundaryMap, Material material, double s_alpha, DGField[] _previous_u, int m_D = 2) : base(boundaryMap, material, s_alpha, _previous_u, m_D) {
        }

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
            double vN = vAverage * (new Vector(normal_x));

            //[[U]] =Uin - Uout
            double[] Ujmp = Uin.CloneAs();
            Ujmp.AccV(-1, Uout);

            //V_0^{-1}*[[U]]
            var result = V0_inv.MatVecMul(1.0, Ujmp);

            // Assemble row of Roe matrix corresponding to Energy
            double[] V0 = new double[D + 1];
            V0[0] = enthalpyAverage - vN * speedOfSoundAverage;
            for(int i = 1; i < D ; i++) {
                V0[i] = vAverage[i - 1] + (vAverage.L2NormPow2() / 2 - vN) * normal[i - 1];
            }
            V0[D] = enthalpyAverage + vN * speedOfSoundAverage;

            // V0*EV*(V_0^{-1}*[[U]]) - EV = matrix with eigenvalues
            double k_j = 0;
            for(int i = 0; i < D + 1; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }

            // +++++++++++ Compute left and right fluxes FL,FR ++++++++++
            Vector n = new Vector(normal);
            double FL = Flux(Uin, D) * n;
            double FR = Flux(Uout, D) * n;

            // ++++++++++ Add all together +++++++++++++++
            double energyFlux = 0.5 * (FL + FR) + 0.5 * k_j;
            //Debug.Assert(energyFlux > 0);
            Debug.Assert(!double.IsNaN(energyFlux) || !double.IsInfinity(energyFlux));

            return energyFlux;


        }


        /// <summary>
        /// Space-time energy flux
        /// </summary>
        /// <param name="U"></param>
        /// <param name="D"></param>
        /// <returns></returns>
        public override Vector Flux(double[] U, int D) {

            // ++++++++++++ space component ++++++++++++
            Vector Output = new Vector(D);
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;
            double density = U[0];
            double energy = U[D];
            double momentumSquared = 0.0;
            for(int d = 0; d < D-1; d++) {
                momentumSquared += U[d + 1] * U[d + 1];
            }
            double pressure = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density);

            for(int d = 0; d < D-1; d++) {
                double flux = U[d + 1] / density * (energy + pressure);
                Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                Output[d] += flux;
            }

            // ++++++++++++ time component ++++++++++++
            Output[D - 1] = energy;
            return Output;
        }


    }
}
