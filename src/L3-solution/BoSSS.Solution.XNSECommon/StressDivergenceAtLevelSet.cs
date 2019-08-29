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
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;
using System.Collections;

namespace BoSSS.Solution.XNSECommon.Operator.Viscosity {

    /// <summary>
    /// 
    /// </summary>
    public class StressDivergenceAtLevelSet : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public StressDivergenceAtLevelSet(LevelSetTracker lstrk, double _reynoldsA, double _reynoldsB, double[] _penalty1, double _penalty2, int _component,
            bool _staticInt = false) {
            this.m_LsTrk = lstrk;
            this.muA = 1 / _reynoldsA;
            this.muB = 1 / _reynoldsB;
            this.penalty1 = _penalty1;
            this.penalty2 = _penalty2;
            this.component = _component;
            this.m_D = lstrk.GridDat.SpatialDimension;
            this.staticInt = _staticInt;
        }

        double muA;
        double muB;
        double[] penalty1;
        double penalty2;
        int component;
        int m_D;

        bool staticInt;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
            double[] TA, double[] TB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            Debug.Assert(this.ArgumentOrdering.Count == 3);
            //Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            //Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            //Debug.Assert(Grad_uA.GetLength(1) == D);
            //Debug.Assert(Grad_uB.GetLength(1) == D);

            double[] Grad_uA_xN = new double[2], Grad_uB_xN = new double[2];
            double Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    Grad_uA_xN[dd] += Grad_uA[dd, d] * N[d];
                    Grad_uB_xN[dd] += Grad_uB[dd, d] * N[d];
                }
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            Debug.Assert(TA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(TB.Length == this.ArgumentOrdering.Count);


            double res = 0;


            res += 0.5 * (TA[0]* muA + TB[0]* muB) * N[0] + 0.5 * (TA[1]* muA + TB[1]* muB) * N[1]; // central difference for stress divergence

            switch (component) {
                case 0:
                    res += -penalty2 / hCellMin * (TA[2]* muA - TB[2]* muB);
                    break;
                case 1:
                    res += -penalty2 / hCellMin * (TA[2] * muA - TB[2]* muB);
                    break;
                default:
                    throw new NotImplementedException();


            }

            return  res;

            //switch (m_ViscosityImplementation) {
            //    // old Form (H-Implementation)
            //    case ViscosityImplementation.H: {

            //double wA;
            //double wB;
            //double wPenalty;
            //if (!weighted) {
            //    wA = 0.5;
            //    wB = 0.5;
            //    wPenalty = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //} else {
            //    wA = muB / (muA + muB);
            //    wB = muA / (muA + muB);
            //    wPenalty = muA * muB / (muA + muB);
            //}

            ////double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //if (!staticInt) {
            //    Ret -= (wA * muA * Grad_uA_xN[component] + wB * muB * Grad_uB_xN[component]) * (vA - vB);                           // consistency term
            //    Ret -= (wA * muA * Grad_vA_xN + wB * muB * Grad_vB_xN) * (uA[component] - uB[component]);     // symmetry term
            //    Ret += (penalty / hCutCellMin) * (uA[component] - uB[component]) * (vA - vB) * wPenalty; // penalty term
            //    // Transpose Term
            //    for (int i = 0; i < D; i++) {
            //        Ret -= (wA * muA * Grad_uA[i, component] + wB * muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
            //        Ret -= (wA * muA * Grad_vA[i] + wB * muB * Grad_vB[i]) * (uA[i] - uB[i]) * N[component];  // symmetry term
            //    }

            //} else {

            //    //wall
            //    Ret -= (wA * muA * Grad_uA_xN[component]) * (vA);                           // consistency term
            //    Ret -= (wA * muA * Grad_vA_xN) * (uA[component]);     // symmetry term
            //    Ret += (penalty / hCutCellMin) * (uA[component] - 0) * (vA) * muA; // penalty term

            //    Ret += (wB * muB * Grad_uB_xN[component]) * (vB);                           // consistency term
            //    Ret += (wB * muB * Grad_vB_xN) * (uB[component]);     // symmetry term
            //    Ret += (penalty / hCutCellMin) * (0 - uB[component]) * (0 - vB) * muB; // penalty term
            //    // Transpose Term
            //    for (int d = 0; d < D; d++) {
            //        Ret -= (wA * muA * Grad_uA[d, component]) * (vA) * N[d];  // consistency term
            //        Ret -= (wA * muA * Grad_vA[d]) * (uA[d]) * N[component];  // symmetry term
            //        Ret += (wB * muB * Grad_uB[d, component]) * (vB) * N[d];  // consistency term
            //        Ret += (wB * muB * Grad_vB[d]) * (uB[d]) * N[component];  // symmetry term
            //    }
            //}
            //return Ret;
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }

        //private static bool rem = true;

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get {
                switch (component) {
                    case 0:
                        return new string[] { VariableNames.StressXX, VariableNames.StressXY, VariableNames.VelocityX };
                    case 1:
                        return new string[] { VariableNames.StressXY, VariableNames.StressYY, VariableNames.VelocityY };
                    default:
                        throw new NotImplementedException();
                }
            }
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
