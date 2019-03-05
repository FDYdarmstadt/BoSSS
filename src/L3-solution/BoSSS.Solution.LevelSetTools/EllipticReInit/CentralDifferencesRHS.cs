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
using System.Diagnostics;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Solution.LevelSetTools.EllipticReInit {


    /// <summary>
    /// The Right Hand Side of the EllipticReinitEquation:
    /// Div(DiffusionRate(Abs(Grad Phi)) Grad(Phi)  )
    /// </summary>
    public class CentralDifferencesRHSForm : RHSForm {

        /// <summary>
        /// ctor
        /// </summary>
        public CentralDifferencesRHSForm(double PenaltyBase, LevelSetTracker LSTrck) : base(PenaltyBase, LSTrck) { }
                

        public new IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

        public override TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.GradUxV);
            }
        }


        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIn, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vA, double[] Grad_vB) {
            double Acc = 0.0;
           
            double AbsGrad_uIn = 0;
            double AbsGrad_uOut = 0;

            for (int d = 0; d < inp.D; d++) {
                AbsGrad_uIn += Grad_uIn[0, d] * Grad_uIn[0, d];
                AbsGrad_uOut += Grad_uOut[0, d] * Grad_uOut[0, d];
            }

            AbsGrad_uIn = Math.Sqrt(AbsGrad_uIn);
            AbsGrad_uOut = Math.Sqrt(AbsGrad_uOut);
            
            //Central Differences
            bool eval = NearFieldBitMask[inp.jCellIn] && NearFieldBitMask[inp.jCellOut];
            for (int d = 0; d < inp.D; d++) {
                // Central Differences
                Acc += 0.5 * (DiffusionRate(AbsGrad_uIn, eval) * Grad_uIn[0, d]* inp.Normale[d] + DiffusionRate(AbsGrad_uOut, eval) * Grad_uOut[0, d]* inp.Normale[d]) * (vIn - vOut);
                // Inner Values
                //Acc += DiffusionRate(AbsGrad_uA, eval) * Grad_uA[0, d] * inp.Normale[d] * vA - DiffusionRate(AbsGrad_uB, eval) * Grad_uB[0, d] * inp.Normale[d] * vB;

                // consistency term
                //The Volume Integrals are already non-symmetric, so theres nothing to symmetrize here
                //Acc += 0.5 * (DiffusionRate(AbsGrad_uA, eval) * Grad_vA[d] + DiffusionRate(AbsGrad_uB, eval) * Grad_vB[d]) * (uA[0] - uB[0]) * inp.Normale[d];  // symmetry term
            }
            return Acc;

        }
    }
}
