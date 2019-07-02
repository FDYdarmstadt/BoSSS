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


        public ViscosityAtLevelSet_FullySymmetric_withEvap(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component, double _rhoA, double _rhoB,
            ThermalParameters thermParams, double _Rint, double _sigma) {
            //double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penalty = _penalty;
            this.component = _component;
            this.D = lstrk.GridDat.SpatialDimension;

            this.rhoA = _rhoA;
            this.rhoB = _rhoB;

            this.kA = thermParams.k_A;
            this.kB = thermParams.k_B;
            this.hVapA = thermParams.hVap_A;
            this.Rint = _Rint;

            this.Tsat = thermParams.T_sat;
            this.sigma = _sigma;
            this.pc = thermParams.pc;
        }

        double muA;
        double muB;
        double penalty;
        int component;

        double rhoA;
        double rhoB;


        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, int jCell) {

            double qEvap = ComputeHeatFlux(paramsNeg, paramsPos, N, jCell);

            if (qEvap == 0.0)
                return 0.0;

            double hVap = (hVapA > 0) ? hVapA : -hVapA;
            double M = qEvap / hVap;

            //Console.WriteLine("mEvap - GeneralizedViscosityAtLevelSet_FullySymmetric: {0}", M);

            return M;

        }


        /// <summary>
        /// default-implementation
        /// </summary>
        public override double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];


            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }
            double Ret = 0.0;

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;


            double M = ComputeEvaporationMass(inp.ParamsNeg, inp.ParamsPos, N, inp.jCell);
            if (M == 0.0)
                return 0.0;

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);


            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1 / rhoA) - (1 / rhoB)) * N[component];     // symmetry term
            Ret -= (penalty / hCutCellMin) * M * ((1 / rhoA) - (1 / rhoB)) * N[component] * (vA - vB) * muMax; // penalty term
                                                                                                               // Transpose Term
            for (int i = 0; i < D; i++) {
                //Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                Ret += 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * N[component] * M * ((1 / rhoA) - (1 / rhoB)) * N[i];
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

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, D), VariableNames.Temperature, "Curvature", "DisjoiningPressure"); //;
            }
        }


    }


}
