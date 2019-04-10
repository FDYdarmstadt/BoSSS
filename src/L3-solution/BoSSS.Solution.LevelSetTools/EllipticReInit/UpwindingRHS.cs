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
    /// \f[ \nabla \cdot ( \text{DiffusionRate}( | \nabla \varphi | )  \nabla \varphi ) \f]
    /// </summary>
    public class EllipticReInitUpwindForm_RHS : RHSForm {

        /// <summary>
        /// 
        /// </summary>
        public EllipticReInitUpwindForm_RHS(double PenaltyBase, LevelSetTracker LSTrck): base(PenaltyBase, LSTrck) {
            this.D = LSTrck.GridDat.SpatialDimension;
            LSTrck.Subscribe(this);
        }

        int D;
        
        /// <summary>
        /// %
        /// </summary>
        public override IList<string> ParameterOrdering
        {
            get
            {
                if (D == 2) {
                    return new string[] { "OldLevelSet", "MeanLevelSetGradient[0]", "MeanLevelSetGradient[1]" };
                }
                else ///3D
                {
                    return new string[] { "OldLevelSet", "MeanLevelSetGradient[0]", "MeanLevelSetGradient[1]", "MeanLevelSetGradient[2]" };
                }
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public override TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.GradUxV);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIn, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIn, double[] Grad_vOut) {
            double Acc = 0.0;

            double AbsGrad_uIn = 0;
            double AbsGrad_uOut = 0;

            // levelSet*gradPhi*n
            double DirectionSelector = 0;

            double DirectionIn = 0;
            double DirectionOut = 0;
            double AbsMeanGradUIn = 0;
            double AbsMeanGradUOut = 0;

            for (int d = 0; d < inp.D; d++) {
                AbsGrad_uIn += Grad_uIn[0, d] * Grad_uIn[0, d];
                AbsGrad_uOut += Grad_uOut[0, d] * Grad_uOut[0, d];

                DirectionIn += inp.Parameters_IN[d + 1] * inp.Normale[d] * Math.Sign(inp.Parameters_IN[0]);
                DirectionOut += inp.Parameters_OUT[d + 1] * inp.Normale[d] * Math.Sign(inp.Parameters_OUT[0]);

                AbsMeanGradUIn += inp.Parameters_IN[d + 1] * inp.Parameters_IN[d + 1];
                AbsMeanGradUOut += inp.Parameters_OUT[d + 1] * inp.Parameters_OUT[d + 1];
            }

            AbsGrad_uIn = Math.Sqrt(AbsGrad_uIn);
            AbsGrad_uOut = Math.Sqrt(AbsGrad_uOut);

            AbsMeanGradUIn = Math.Sqrt(AbsMeanGradUIn);
            AbsMeanGradUOut = Math.Sqrt(AbsMeanGradUOut);
            //DirectionIn /= AbsMeanGradUIn;
            //DirectionOut /= AbsMeanGradUOut;


            bool cut = CutCells[inp.jCellIn] && CutCells[inp.jCellOut];

            // The information is always advancing from the cut-cell into the outer cells
            //bool CellInCut = LSTrck._Regions.GetCutCellMask().GetBitMask()[inp.jCellIn];
            //bool CellOutCut = LSTrck._Regions.GetCutCellMask().GetBitMask()[inp.jCellOut];
            bool CellInCut = CutCells[inp.jCellIn];
            bool CellOutCut = CutCells[inp.jCellOut];

            if (CellInCut && CellOutCut) {
                //central differences
                for (int d = 0; d < inp.D; d++) {
                    //return 0;
                    Acc += 0.5 * (DiffusionRate(AbsGrad_uIn, cut) * Grad_uIn[0, d] * inp.Normale[d] + DiffusionRate(AbsGrad_uOut, cut) * Grad_uOut[0, d] * inp.Normale[d]) * (vIn-vOut);
                }
                return Acc;
            }
            else if(CellInCut && !CellOutCut) {
                DirectionSelector = 1;
            }
            else if (!CellInCut && CellOutCut) {
                DirectionSelector = -1;
            }
            else {
                DirectionSelector = DirectionIn  + DirectionOut ;
                //if (DirectionIn * DirectionOut < 0 && DirectionIn > 0) {
                //    return 0;
                //}
            }

            if (DirectionSelector > 0) {
                for (int d = 0; d < inp.D; d++) {
                    Acc += 0.5 * (DiffusionRate(AbsGrad_uIn, cut) * Grad_uIn[0, d] * inp.Normale[d] + DiffusionRate(AbsGrad_uOut, cut) * Grad_uOut[0, d] * inp.Normale[d]) * (- vOut);
                }
            }
            else // if (DirectionSelector <= 0) 
                {
                for (int d = 0; d < inp.D; d++) {
                    Acc += 0.5 * (DiffusionRate(AbsGrad_uIn, cut) * Grad_uIn[0, d] * inp.Normale[d] + DiffusionRate(AbsGrad_uOut, cut) * Grad_uOut[0, d] * inp.Normale[d]) * inp.Normale[d] * (vIn );
                }
            }
            //else {
            //    /// Do nothing
            //}
            return Acc;

        }

        

    }
}
