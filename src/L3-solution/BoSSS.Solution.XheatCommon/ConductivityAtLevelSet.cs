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
using System.Threading.Tasks;

using ilPSP;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;


namespace BoSSS.Solution.XheatCommon {


    public class ConductivityAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ConductivityAtLevelSet(LevelSetTracker lstrk, double _kA, double _kB, double _penalty, double _Tsat) {
            this.m_LsTrk = lstrk;
            this.kA = _kA;
            this.kB = _kB;
            this.penalty = _penalty;
            this.Tsat = _Tsat;
            this.m_D = lstrk.GridDat.SpatialDimension;

        }

        double kA;
        double kB;

        double penalty;

        double Tsat;

        int m_D;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for(int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            Debug.Assert(!(double.IsInfinity(hCutCellMin) || double.IsNaN(hCutCellMin)));

            if(hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            double Ret = 0.0;

            Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - vB);                           // consistency term
            //Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (vA - 0);                           // consistency term
            //Ret -= 0.5 * (kA * Grad_uA_xN + kB * Grad_uB_xN) * (0 - vB);                           // consistency term

            Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - uB[0]);                     // symmetry term
            //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (uA[0] - Tsat);                     // symmetry term
            //Ret -= 0.5 * (kA * Grad_vA_xN + kB * Grad_vB_xN) * (Tsat - uB[0]);                     // symmetry term

            Ret += (penalty / hCutCellMin) * (uA[0] - uB[0]) * (vA - vB) * (Math.Abs(kA) > Math.Abs(kB) ? kA : kB); // penalty term
            //Ret += (penalty / hCutCellMin) * (uA[0] - Tsat) * (vA - vB) * (Math.Abs(kA) > Math.Abs(kB) ? kA : kB); // penalty term
            //Ret += (penalty / hCutCellMin) * (Tsat - uB[0]) * (vA - vB) * (Math.Abs(kA) > Math.Abs(kB) ? kA : kB); // penalty term


            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
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

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.Temperature }; }
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
            get {
                return null;
            }
        }


    }


}
