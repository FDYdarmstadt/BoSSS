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
    class EllipticReInitInterfaceForm : ILevelSetForm {
        double PenaltyBase;
        LevelSetTracker LSTrk;

        public EllipticReInitInterfaceForm(double PenaltyBase, LevelSetTracker LSTrk) {
            this.PenaltyBase = PenaltyBase;
            this.LSTrk = LSTrk;
        }

        public TermActivationFlags LevelSetTerms
        {
            get
            {
                return (TermActivationFlags.UxV);
            }
        }



        /// <summary>
        /// The penalty at the interface enforcing phi=0 at the old position
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="uA"></param>
        /// <param name="uB"></param>
        /// <param name="Grad_uA"></param>
        /// <param name="Grad_uB"></param>
        /// <param name="vA"></param>
        /// <param name="vB"></param>
        /// <param name="Grad_vA"></param>
        /// <param name="Grad_vB"></param>
        /// <returns></returns>
        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double hmin;
            if (inp.NegCellLengthScale.IsNaN()) {
                hmin = inp.PosCellLengthScale;
            }
            else if (inp.PosCellLengthScale.IsNaN()) {
                hmin = inp.NegCellLengthScale;
            }
            else {
                hmin = Math.Min(inp.NegCellLengthScale, inp.PosCellLengthScale);
            }

            return - 2* PenaltyBase /  hmin * (uA[0]) * (vA);

        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.LevelSet };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { };
            }
        }

        public int LevelSetIndex
        {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies
        {
            get { return this.LSTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies
        {
            get { return this.LSTrk.GetSpeciesId("A"); }
        }
    }
}
