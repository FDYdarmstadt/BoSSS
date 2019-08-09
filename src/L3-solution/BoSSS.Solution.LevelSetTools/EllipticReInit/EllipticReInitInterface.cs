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
    /// Penalty Operator for the Zero Boundary Condition at the Interface
    /// \f[ 
    /// a(u,v) = \alpha \int_{\Gamma} u v   \mathrm{dS}
    /// \f]
    /// </summary>
    public class EllipticReInitInterfaceForm : ILevelSetForm, ILevelSetEquationComponentCoefficient {
        readonly double PenaltyBase;

        readonly LevelSetTracker LSTrk;

        /// <summary>
        /// old ctor
        /// </summary>
        public EllipticReInitInterfaceForm(double PenaltyBase, LevelSetTracker LSTrk) {
            this.PenaltyBase = PenaltyBase;
            this.LSTrk = LSTrk;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return (TermActivationFlags.UxV);
            }
        }

        MultidimensionalArray NegCellLengthScaleS;
        MultidimensionalArray PosCellLengthScaleS;

        /// <summary>
        /// Called by 
        /// <see cref="XSpatialOperatorMk2.XEvaluatorNonlin.Evaluate{Tout}(double, double, Tout, double[])"/>
        /// resp.
        /// <see cref="XSpatialOperatorMk2.XEvaluatorLinear.ComputeMatrix{M, V}(M, V)"/>.
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegCellLengthScaleS = csA.CellLengthScales;
            if(csB != null)
                PosCellLengthScaleS = csB.CellLengthScales;
        }


        /// <summary>
        /// The penalty at the interface enforcing phi=0 at the old position
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double NegCellLengthScale = NegCellLengthScaleS[inp.jCell];
            double PosCellLengthScale = (PosCellLengthScaleS != null) ? PosCellLengthScaleS[inp.jCell] : NegCellLengthScaleS[inp.jCell];

            double hmin;
            if(NegCellLengthScale.IsNaN()) {
                hmin = PosCellLengthScale;
            } else if(PosCellLengthScale.IsNaN()) {
                hmin = NegCellLengthScale;
            } else {
                hmin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            }

            return -2 * PenaltyBase / hmin * (uA[0]) * (vA);

        }

        /// <summary>
        /// Only depends on the level-set
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.LevelSet };
            }
        }

        /// <summary>
        /// empty
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

        /// <summary>
        /// 0
        /// </summary>
        public int LevelSetIndex {
            get { return 0; }
        }

        /// <summary>
        /// B
        /// </summary>
        public SpeciesId PositiveSpecies {
            get { return this.LSTrk.GetSpeciesId("B"); }
        }

        /// <summary>
        /// A
        /// </summary>
        public SpeciesId NegativeSpecies {
            get { return this.LSTrk.GetSpeciesId("A"); }
        }
    }
}
