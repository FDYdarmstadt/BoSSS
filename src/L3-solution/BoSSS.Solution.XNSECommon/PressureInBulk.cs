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
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Solution.XNSECommon.Operator.Pressure {
    
    /// <summary>
    /// pressure gradient, bulk phase, by central differences
    /// </summary>
    public class PressureInBulk : PressureGradientLin_d {

        public PressureInBulk(int _d, IncompressibleMultiphaseBoundaryCondMap bcMap, double rhoA, double rhoB)
            : base(_d, bcMap) {
            m_rhoA = rhoA;
            m_rhoB = rhoB;
            base.pressureFunction = null;
            this.m_bcMap = bcMap;
        }

        double m_rhoA;
        double m_rhoB;
        IncompressibleMultiphaseBoundaryCondMap m_bcMap;


        public void SetParameter(string speciesName, SpeciesId SpcId, MultidimensionalArray LenScales) {
            switch (speciesName) {
                case "A": oneOverRho = 1.0 / m_rhoA; SetBndfunction("A"); break;
                case "B": oneOverRho = 1.0 / m_rhoB; SetBndfunction("B"); break;
            default: throw new ArgumentException("Unknown species.");
            }
        }


        void SetBndfunction(string S) {
            base.pressureFunction = this.m_bcMap.bndFunction[VariableNames.Pressure + "#" + S];
        }

        double oneOverRho = double.NaN;


        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            //switch (m_varMode) {
            //    case EquationAndVarMode.u_p_2:
            //    case EquationAndVarMode.mom_p:
                return base.BorderEdgeFlux(ref inp, Uin);

            //    case EquationAndVarMode.u_p:
            //    return base.BorderEdgeFlux(ref inp, Uin)*oneOverRho;

            //    case EquationAndVarMode.u_Psi:
            //    Uin[0] /= oneOverRho;
            //    return base.BorderEdgeFlux(ref inp, Uin)*oneOverRho;

            //    default:
            //    throw new NotImplementedException();
            //}
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            base.Flux(ref inp, U, output);
            //switch (m_varMode) {
            //    case EquationAndVarMode.u_p:
            //    BLAS.dscal(output.Length, oneOverRho, output, 1);
            //    break;

            //    case EquationAndVarMode.u_p_2:
            //    case EquationAndVarMode.mom_p:
            //    case EquationAndVarMode.u_Psi:
            //    // nop
            //    break;

            //    default:
            //    throw new NotImplementedException();
            //}
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            //switch (m_varMode) {
            //    case EquationAndVarMode.u_p:
            //    return base.InnerEdgeFlux(ref inp, Uin, Uout)*oneOverRho;

            //    case EquationAndVarMode.u_p_2:
            //    case EquationAndVarMode.mom_p:
            //    case EquationAndVarMode.u_Psi:
                return base.InnerEdgeFlux(ref inp, Uin, Uout);

            //    default:
            //    throw new NotImplementedException();
            //}
        }
    }
}

