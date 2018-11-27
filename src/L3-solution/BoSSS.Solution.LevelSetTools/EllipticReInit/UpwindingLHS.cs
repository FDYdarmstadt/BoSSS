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
    /// The Laplace Operator for the Elliptic Reinitialization
    /// \f$ \operatorname{div}(\operatorname{grad} \varphi) \f$
    /// No Boundary Conditions are set -> Boundary Conditions are determined by Interface only
    /// </summary>
    class EllipticReInitUpwindForm_Laplace : SIPLaplace, IObserver<LevelSetTracker.LevelSetRegions> {
        public EllipticReInitUpwindForm_Laplace(double PenaltyBase, LevelSetTracker LSTrck) : base(PenaltyBase, LSTrck.GridDat.Cells.PenaltyLengthScales, VariableNames.LevelSet) {
            this.D = LSTrck.GridDat.SpatialDimension;
            this.m_CutCells = LSTrck.Regions.GetCutCellMask().GetBitMask();
            this.LSTrck = LSTrck;
            LSTrck.Subscribe(this);
        }

        LevelSetTracker LSTrck;
        System.Collections.BitArray m_CutCells;
        int D;

        public override IList<string> ParameterOrdering
        {
            get
            {   if (D == 2) {
                    return new string[] { "OldLevelSet", "MeanLevelSetGradient[0]", "MeanLevelSetGradient[1]" };
                }
                else ///3D
                {
                    return new string[] { "OldLevelSet", "MeanLevelSetGradient[0]", "MeanLevelSetGradient[1]", "MeanLevelSetGradient[2]" };
                }
            }
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIn, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIn, double[] Grad_vOut) {
            double Acc = 0.0;
            double pnlty = 2 * base.GetPenalty(inp.jCellIn, inp.jCellOut);

            // levelSet*gradPhi*n
            double DirectionSelector = 0;

            double DirectionIn = 0;
            double DirectionOut = 0;
            double AbsMeanGradUIn = 0;
            double AbsMeanGradUOut = 0;

            for (int d = 0; d < inp.D; d++) {
                DirectionIn += inp.Parameters_IN[d+1] * inp.Normale[d] * Math.Sign(inp.Parameters_IN[0]);
                DirectionOut += inp.Parameters_OUT[d+1] * inp.Normale[d] * Math.Sign(inp.Parameters_OUT[0]);
                AbsMeanGradUIn += inp.Parameters_IN[d + 1]* inp.Parameters_IN[d + 1];
                AbsMeanGradUOut += inp.Parameters_OUT[d + 1]* inp.Parameters_OUT[d + 1];
            }

            AbsMeanGradUIn = Math.Sqrt(AbsMeanGradUIn);
            AbsMeanGradUOut = Math.Sqrt(AbsMeanGradUOut);
            //DirectionIn /= AbsMeanGradUIn;
            //DirectionOut /= AbsMeanGradUOut;


            double penaltyfactor = 1;

            // levelSet*gradPhi*n

            // The information is always advancing from the cut-cell into the outer cells
            //bool CellInCut = LSTrck._Regions.GetCutCellMask().GetBitMask()[inp.jCellIn];
            //bool CellOutCut = LSTrck._Regions.GetCutCellMask().GetBitMask()[inp.jCellOut];
            bool CellInCut = m_CutCells[inp.jCellIn];
            bool CellOutCut = m_CutCells[inp.jCellOut];

            if (CellInCut && CellOutCut) {
                //central differences
                for (int d = 0; d < inp.D; d++) {
                    Acc += 0.5 * (Grad_vIn[d] + Grad_vOut[d]) * inp.Normale[d] * (uIn[0] - uOut[0]);
                    Acc += 0.5 * (Grad_uIn[0, d] + Grad_uOut[0, d]) * inp.Normale[d] * (vIn - vOut);
                }
                Acc -= (uIn[0] - uOut[0]) * (vIn - vOut) * penaltyfactor * pnlty;
                return Acc;
            }
            else if (CellInCut && !CellOutCut) {
                DirectionSelector = 1;
            }
            else if (!CellInCut && CellOutCut) {
                DirectionSelector = -1;
            }
            else {
                DirectionSelector = DirectionIn + DirectionOut;
                //if (DirectionIn * DirectionOut < 0 && DirectionIn >0 ) {
                //    return 0;
                //}
            }

            

            if (DirectionSelector > 0) {
                for (int d = 0; d < inp.D; d++) {
                    Acc += Grad_vOut[d] * inp.Normale[d] * (uIn[0] - uOut[0]);
                    Acc += 0.5 * (Grad_uIn[0, d] + Grad_uOut[0, d]) * inp.Normale[d] * (-vOut);
                }

                Acc -= (uIn[0] - uOut[0]) * (-vOut) * penaltyfactor * pnlty;
            }
            else //if (DirectionSelector <= 0)
            {

                for (int d = 0; d < inp.D; d++) {
                    Acc += Grad_vIn[d] * inp.Normale[d] * (uIn[0] - uOut[0]);
                    Acc += 0.5 * (Grad_uIn[0, d] + Grad_uOut[0, d]) * inp.Normale[d] * (vIn);
                }
                Acc -= (uIn[0] - uOut[0]) * (vIn) * penaltyfactor * pnlty;
            }
            //else {
            //    /// Do nothing
            //}
            return Acc;

        }

        public void OnNext(LevelSetTracker.LevelSetRegions value) {
            m_CutCells = LSTrck.Regions.GetCutCellMask().GetBitMask();
        }

        public void OnError(Exception error) {
            // Do Nothing
        }

        public void OnCompleted() {
            // Do nothing
        }

        /// <summary>
        /// Here is some more Code doing nothing for performance reasons
        /// like this, the boundary terms are not even evaluated
        /// </summary>
        public new TermActivationFlags BoundaryEdgeTerms
        {
            get
            {
                return (TermActivationFlags.None);
            }
        }



        public new double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] uA, double[,] Grad_uA, double vA, double[] Grad_vA) {
            return 0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return false;
        }
    }

}
