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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Pressure stabilization for incompressible flows using an equal-order formulation.
    /// D. A. D. Pietro and A. Ern, Mathematical Aspects of Discontinuous Galerkin Methods. Springer Berlin Heidelberg, 2012.
    /// (Chapter 6.2.4.2)
    /// </summary>
    public class PressureStabilizationAtLevelSet : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;
        double PressureStabilizationFactor;
        protected double ReynoldsA;
        protected double ReynoldsB;

        /// <summary>
        /// Ctor.
        /// </summary>
        public PressureStabilizationAtLevelSet(LevelSetTracker lstrk, double PressureStabilizationFactor, double _reynoldsA, double _reynoldsB) {
            this.m_LsTrk = lstrk;
            this.PressureStabilizationFactor = PressureStabilizationFactor;
            this.ReynoldsA = _reynoldsA;
            this.ReynoldsB = _reynoldsB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Pressure };
            }
        }

        MultidimensionalArray h_max_Edge;
        public double InnerEdgeForm(ref CommonParams inp,
            double[] UA, double[] UB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double Flx_InCell, Flx_OutCell;

            double PosCellLengthScale = PosLengthScaleS[inp.jCellOut];
            double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

            double h_max = Math.Max(NegCellLengthScale, PosCellLengthScale);
            double penalty = PressureStabilizationFactor * h_max;
            Flx_InCell = penalty * (UA[0] * ReynoldsA - UB[0] * ReynoldsB);
            Flx_OutCell = -Flx_InCell;

            double res = 0;

            res = Flx_InCell* vA +Flx_OutCell * vB;

            return res;

        }

        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }

     
    }
}
