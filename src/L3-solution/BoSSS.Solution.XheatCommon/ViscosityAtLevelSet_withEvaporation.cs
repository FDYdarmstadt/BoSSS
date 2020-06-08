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

namespace BoSSS.Solution.XheatCommon {

    /// <summary>
    /// 
    /// </summary>
    public class ViscosityAtLevelSet_FullySymmetric_withEvap : EvaporationAtLevelSet {


        public ViscosityAtLevelSet_FullySymmetric_withEvap(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component, 
            ThermalParameters thermParams, double _sigma) 
            : base(lstrk.GridDat.SpatialDimension, lstrk, thermParams, _sigma) {

            this.muA = _muA;
            this.muB = _muB;
            this.penalty = _penalty;
            this.component = _component;

        }

        double muA;
        double muB;
        double penalty;
        int component;



        /// <summary>
        /// default-implementation
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];


            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < m_D; d++) {
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }
            double Ret = 0.0;

            double PosCellLengthScale = PosLengthScaleS[inp.jCellOut];
            double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;


            double M = ComputeEvaporationMass(inp.Parameters_IN, inp.Parameters_OUT, N, inp.jCellIn);
            if (M == 0.0)
                return 0.0;

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);


            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component];     // symmetry term
            Ret -= (penalty / hCutCellMin) * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[component] * (vA - vB) * muMax; // penalty term
                                                                                                               // Transpose Term
            for (int i = 0; i < m_D; i++) {
                //Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                Ret += 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * N[component] * M * ((1 / m_rhoA) - (1 / m_rhoB)) * N[i];
            }

            return -Ret;
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;


        public override void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            base.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);

            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

        }


        public override TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.GradV | TermActivationFlags.V;
            }
        }

        //public override IList<string> ParameterOrdering {
        //    get {
        //        return ArrayTools.Cat(VariableNames.HeatFlux0Vector(m_D), VariableNames.Temperature0, VariableNames.Curvature, VariableNames.DisjoiningPressure);
        //    }
        //}


    }


}
