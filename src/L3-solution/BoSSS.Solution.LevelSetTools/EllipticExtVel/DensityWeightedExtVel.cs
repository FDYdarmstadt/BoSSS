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
/* =======================================================================
Copyright 2018 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;


namespace BoSSS.Solution.LevelSetTools.EllipticExtension {
    /// <summary>
    /// Penalty Operator for the Zero Boundary Condition at the Interface
    /// \f[ 
    /// a(u,v) = \alpha \int_{\Gamma} u v   \mathrm{dS}
    /// \f]
    /// </summary>
    public class DensityWeightedExtVel : ILevelSetForm {

        double PenaltyBase;
        LevelSetTracker LSTrk;
        double[] Weights;

        public DensityWeightedExtVel(double PenaltyBase, LevelSetTracker LSTrk, double[] rho) {
            this.PenaltyBase = PenaltyBase;
            this.LSTrk = LSTrk;
            double SumOfDensities = rho.Sum();
            Weights = new double[2];
            Weights[0] = rho[0] / SumOfDensities;
            Weights[1] = rho[1] / SumOfDensities;
        }


        /// <summary>
        /// Penalty Term enforcing the duirichlet value at the interface
        /// Note: this Form is written only in terms of uA, since there is no XDG-field involved
        /// </summary>
        /// <param name="inp">inp.ParamsNeg[0] is the Dirichlet value from the parameter-field</param>
        /// <param name="uA">the unknown</param>
        /// <param name="uB">not needed</param>
        /// <param name="Grad_uA">not needed</param>
        /// <param name="Grad_uB">not needed</param>
        /// <param name="vA">test function</param>
        /// <param name="vB">not needed</param>
        /// <param name="Grad_vA">not needed</param>
        /// <param name="Grad_vB">not needed</param>
        /// <returns>the evaluated penalty flux</returns>
        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double hmin;
            if(inp.NegCellLengthScale.IsNaN()) {
                hmin = inp.PosCellLengthScale;
            } else if(inp.PosCellLengthScale.IsNaN()) {
                hmin = inp.NegCellLengthScale;
            } else {
                hmin = Math.Min(inp.NegCellLengthScale, inp.PosCellLengthScale);
            }

            return PenaltyBase * 2 / hmin * (uA[0] - (inp.ParamsNeg[0] * Weights[0] + inp.ParamsPos[0] * Weights[1])) * (vA);

        }



        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "Extension" };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { "InterfaceValue" };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.LSTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.LSTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.V);
            }
        }

    }
}
