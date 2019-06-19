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
    public class ViscosityAtLevelSet_FullySymmetric : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ViscosityAtLevelSet_FullySymmetric(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component, 
            bool _staticInt = false, bool _weighted = false) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penalty = _penalty;
            this.component = _component;
            this.m_D = lstrk.GridDat.SpatialDimension;
            this.staticInt = _staticInt;
            this.weighted = _weighted;
        }

        double muA;
        double muB;
        double penalty;
        int component;
        int m_D;

        bool staticInt;
        bool weighted;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double[] Grad_uA_xN = new double[2], Grad_uB_xN = new double[2];
            double Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    Grad_uA_xN[dd] += Grad_uA[dd, d] * N[d];
                    Grad_uB_xN[dd] += Grad_uB[dd, d] * N[d];
                }
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

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);
            //switch (m_ViscosityImplementation) {
            //    // old Form (H-Implementation)
            //    case ViscosityImplementation.H: {

            double wA;
            double wB;
            double wPenalty;
            if(!weighted) {
                wA = 0.5;
                wB = 0.5;
                wPenalty = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            } else { 
                wA = muB / (muA + muB);
                wB = muA / (muA + muB);
                wPenalty = muA*muB / (muA + muB);
            } 

            //double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            if(!staticInt) {
                Ret -= (wA * muA * Grad_uA_xN[component] + wB * muB * Grad_uB_xN[component]) * (vA - vB);                           // consistency term
                Ret -= (wA * muA * Grad_vA_xN + wB * muB * Grad_vB_xN) * (uA[component] - uB[component]);     // symmetry term
                Ret += (penalty / hCutCellMin) * (uA[component] - uB[component]) * (vA - vB) * wPenalty; // penalty term
                // Transpose Term
                for(int i = 0; i < D; i++) {
                    Ret -= (wA * muA * Grad_uA[i, component] + wB * muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                    Ret -= (wA * muA * Grad_vA[i] + wB * muB * Grad_vB[i]) * (uA[i] - uB[i]) * N[component];  // symmetry term
                }

            } else {
                ////free slip
                //for(int d = 0; d < D; d++) {
                //    Ret -= N[d] * (wA * muA * Grad_uA_xN[d]) * (vA) * N[component];                           // consistency term
                //    Ret -= N[component] * (wA * muA * Grad_vA_xN) * (uA[d]) * N[d];     // symmetry term
                //    Ret += (penalty / hCutCellMin) * (uA[d] - 0) * N[d] * (vA) * N[component] * muA; // penalty term

                //    Ret += N[d] * (wB * muB * Grad_uB_xN[d]) * (vB) * N[component];                           // consistency term
                //    Ret += N[component] * (wB * muB * Grad_vB_xN) * (uB[d]) * N[d];     // symmetry term
                //    Ret += (penalty / hCutCellMin) * (0 - uB[d]) * N[d] * (0 - vB) * N[component] * muB; // penalty term
                //}
                //// Transpose Term
                //for(int dN = 0; dN < D; dN++) {
                //    for(int dD = 0; dD < D; dD++) {
                //        Ret -= N[dN] * (wA * muA * Grad_uA[dD, dN]) * N[dD] * (vA) * N[component];  // consistency term
                //        Ret -= N[dN] * (wA * muA * Grad_vA[dN]) * N[component] * (uA[dD]) * N[dD];  // symmetry term
                //        Ret += N[dN] * (wB * muB * Grad_uB[dD, dN]) * N[dD] * (vB) * N[component];  // consistency term
                //        Ret += N[dN] * (wB * muB * Grad_vB[dN]) * N[component] * (uB[dD]) * N[dD];  // symmetry term
                //    }
                //}

                //wall
                Ret -= (wA * muA * Grad_uA_xN[component]) * (vA);                           // consistency term
                Ret -= (wA * muA * Grad_vA_xN) * (uA[component]);     // symmetry term
                Ret += (penalty / hCutCellMin) * (uA[component] - 0) * (vA) * muA; // penalty term

                Ret += (wB * muB * Grad_uB_xN[component]) * (vB);                           // consistency term
                Ret += (wB * muB * Grad_vB_xN) * (uB[component]);     // symmetry term
                Ret += (penalty / hCutCellMin) * (0 - uB[component]) * (0 - vB) * muB; // penalty term
                // Transpose Term
                for(int d = 0; d < D; d++) {
                    Ret -= (wA * muA * Grad_uA[d, component]) * (vA) * N[d];  // consistency term
                    Ret -= (wA * muA * Grad_vA[d]) * (uA[d]) * N[component];  // symmetry term
                    Ret += (wB * muB * Grad_uB[d, component]) * (vB) * N[d];  // consistency term
                    Ret += (wB * muB * Grad_vB[d]) * (uB[d]) * N[component];  // symmetry term
                }
            }
                        


            //break;^^
            //    }
            //    // SWIP-form nach DiPietro/Ern:
            //    case ViscosityImplementation.SWIP:{
            //        Ret -= ( muB * muA * Grad_uA_xN + muA* muB* Grad_uB_xN) / (muA+muB) * (vA - vB);
            //        Ret -= (muB * muA * Grad_vA_xN + muA* muB* Grad_vB_xN)  / (muA + muB) * (uA[component] - uB[component]);
            //        Ret += (penalty / hCutCellMin) * (uA[component] - uB[component]) * (vA - vB) *(2.0*muA*muB/(muA + muB));
            //        // Transpose-Term
            //        for (int i = 0; i < D; i++) {
            //            Ret -= (muB * muA * Grad_uA[i, component] + muA * muB * Grad_uB[i, component]) / (muA + muB) * (vA - vB) * N[i];  // consistency term
            //            Ret -= (muB * muA * Grad_vA[i] + muA * muB * Grad_vB[i]) / (muA + muB) * (uA[i] - uB[i]) * N[component];  // symmetry term
            //        }
            //        break;
            //    }
            //    default: { throw new ArgumentException(); }
            //}
            /*
            {
                double Acc = 0.0;
                for (int i = 0; i < D; i++) {
                    Acc += 0.5 * (muA*Grad_uA[i, i] + muB*Grad_uB[i, i]) * (vA - vB) * N[component];  // consistency term
                    Acc += 0.5 * (muA*Grad_vA[component] + muB*Grad_vB[component]) * (uA[i] - uB[i]) * N[i];  // symmetry term
                }
                
                //double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
                //Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax; // penalty term

                Ret += Acc*(2.0/3.0);

            } // */


            return Ret;
        }

        

        //public static Stopwatch Optimized = new Stopwatch();
        //public static Stopwatch Classic = new Stopwatch();



       // public override void EdgeForm(LevSetIntParams inp, MultidimensionalArray Koeff_UxV, MultidimensionalArray Koeff_NablaUxV, MultidimensionalArray Koeff_UxNablaV, MultidimensionalArray Koeff_NablaUxNablaV) {
            //Optimized.Start();
            //this.EdgeForm_Perf(inp, Koeff_UxV, Koeff_NablaUxV, Koeff_UxNablaV, Koeff_NablaUxNablaV);
            //base.EdgeForm(inp, Koeff_UxV, Koeff_NablaUxV, Koeff_UxNablaV, Koeff_NablaUxNablaV);
            //Optimized.Stop();
            
            /*
            {
                var chk_Koeff_UxV = Koeff_UxV.CloneAs();
                var chk_Koeff_NablaUxV = Koeff_NablaUxV.CloneAs();
                var chk_Koeff_UxNablaV = Koeff_UxNablaV.CloneAs();
                var chk_Koeff_NablaUxNablaV = Koeff_NablaUxNablaV.CloneAs();
                chk_Koeff_UxV.Clear();
                chk_Koeff_UxNablaV.Clear();
                chk_Koeff_NablaUxV.Clear();
                chk_Koeff_NablaUxNablaV.Clear();

                Classic.Start();
                base.EdgeForm(inp, chk_Koeff_UxV, chk_Koeff_NablaUxV, chk_Koeff_UxNablaV, chk_Koeff_NablaUxNablaV);
                Classic.Stop();

                chk_Koeff_NablaUxNablaV.Acc(-1.0, Koeff_NablaUxNablaV);
                chk_Koeff_UxNablaV.Acc(-1.0, Koeff_UxNablaV);
                chk_Koeff_NablaUxV.Acc(-1.0, Koeff_NablaUxV);
                chk_Koeff_UxV.Acc(-1.0, Koeff_UxV);

                if (chk_Koeff_NablaUxNablaV.L2Norm() > 1.0e-10)
                    throw new ApplicationException();
                if (chk_Koeff_UxNablaV.L2Norm() > 1.0e-10)
                    throw new ApplicationException();
                if (chk_Koeff_NablaUxV.L2Norm() > 1.0e-10)
                    throw new ApplicationException();
                if (chk_Koeff_UxV.L2Norm() > 1.0e-10)
                    throw new ApplicationException();
            } // */
        //}



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
            get { return VariableNames.VelocityVector(this.m_D); }
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
            get { return null;  }
        }

        /*
        /// <summary>
        /// Performance implementation
        /// </summary>
        public void EdgeForm_Perf(LevSetIntParams inp, MultidimensionalArray Koeff_UxV, MultidimensionalArray Koeff_NablaUxV, MultidimensionalArray Koeff_UxNablaV, MultidimensionalArray Koeff_NablaUxNablaV) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_UxV.CheckLengths(Len, N, NoOfVars, 2, 2);
            Koeff_NablaUxV.CheckLengths(Len, N, NoOfVars, 2, 2, D);
            Koeff_UxNablaV.CheckLengths(Len, N, NoOfVars, 2, 2, D);
            Koeff_NablaUxNablaV.CheckLengths(Len, N, NoOfVars, 2, 2, D, D);

            var h_max = this.m_LsTrk.GridDat.Cells.h_max;
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;

            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            
            if (D == 2) {
                switch (this.component) {
                    case 0: EdgeForm_Perf_21(inp, Koeff_UxV, Koeff_NablaUxV, Koeff_UxNablaV, Len, N, lsTrk, h_min, ref posSpc, ref negSpc); break;
                    case 1: EdgeForm_Perf_22(inp, Koeff_UxV, Koeff_NablaUxV, Koeff_UxNablaV, Len, N, lsTrk, h_min, ref posSpc, ref negSpc); break;
                }
            } else if (D == 3) {
                throw new NotImplementedException("todo");
            } else {
                throw new NotSupportedException("unknown spatial dimension");
            }
            Koeff_NablaUxNablaV.Clear();

            throw new ApplicationException("untested code");
        }

        private void EdgeForm_Perf_21(LevSetIntParams inp, MultidimensionalArray Koeff_UxV, MultidimensionalArray Koeff_NablaUxV, MultidimensionalArray Koeff_UxNablaV, int Len, int N, LevelSetTracker lsTrk, double[] h_min, ref SpeciesId posSpc, ref SpeciesId negSpc) {
            double N_mu = this.muA;
            double P_mu = this.muB;
            for (int j = 0; j < Len; j++) { // loop over items...

                ReducedRegionCode rrc;
                int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                Debug.Assert(NoOf == 2);
                int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                int jCell = j + inp.i0;
                double eta = (penalty/h_min[jCell])*(Math.Abs(muA) > Math.Abs(muB) ? muA : muB);

                for (int n = 0; n < N; n++) { // loop over nodes...
                    double N1 = inp.Normal[j, n, 0];
                    double N2 = inp.Normal[j, n, 1];

                    Koeff_UxV[j, n, 0, iSpcNeg, iSpcNeg] = -1.0*eta;
                    Koeff_UxV[j, n, 0, iSpcPos, iSpcNeg] = eta;
                    Koeff_UxV[j, n, 0, iSpcNeg, iSpcPos] = eta;
                    Koeff_UxV[j, n, 0, iSpcPos, iSpcPos] = -1.0*eta;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcNeg, 0] = -1.50*N_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcNeg, 0] = 1.50*N_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcPos, 0] = .500*P_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcPos, 0] = -.500*P_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcNeg, 1] = -.500*N_mu*N2;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcNeg, 1] = .500*N_mu*N2;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcPos, 1] = -.500*P_mu*N2;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcPos, 1] = .500*P_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcNeg, 0] = -1.50*N_mu*N1;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcNeg, 0] = .500*P_mu*N1;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcPos, 0] = 1.50*N_mu*N1;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcPos, 0] = -.500*P_mu*N1;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcNeg, 1] = -.500*N_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcNeg, 1] = -.500*P_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcPos, 1] = .500*N_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcPos, 1] = .500*P_mu*N2;
                    Koeff_UxV[j, n, 1, iSpcNeg, iSpcNeg] = 0.0;
                    Koeff_UxV[j, n, 1, iSpcPos, iSpcNeg] = 0.0;
                    Koeff_UxV[j, n, 1, iSpcNeg, iSpcPos] = 0.0;
                    Koeff_UxV[j, n, 1, iSpcPos, iSpcPos] = 0.0;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcNeg, 0] = -1.0*N_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcNeg, 0] = N_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcPos, 0] = P_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcPos, 0] = -1.0*P_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcNeg, 1] = 0.0;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcNeg, 1] = 0.0;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcPos, 1] = 0.0;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcPos, 1] = 0.0;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcNeg, 0] = 0.0;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcNeg, 0] = 0.0;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcPos, 0] = 0.0;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcPos, 0] = 0.0;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcNeg, 1] = -1.0*N_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcNeg, 1] = P_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcPos, 1] = N_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcPos, 1] = -1.0*P_mu*N1;
                }
            }
        }

        private void EdgeForm_Perf_22(LevSetIntParams inp, MultidimensionalArray Koeff_UxV, MultidimensionalArray Koeff_NablaUxV, MultidimensionalArray Koeff_UxNablaV, int Len, int N, LevelSetTracker lsTrk, double[] h_min, ref SpeciesId posSpc, ref SpeciesId negSpc) {
            double N_mu = this.muA;
            double P_mu = this.muB;
            for (int j = 0; j < Len; j++) { // loop over items...

                ReducedRegionCode rrc;
                int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                Debug.Assert(NoOf == 2);
                int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                int jCell = j + inp.i0;
                double eta = (penalty/h_min[jCell])*(Math.Abs(muA) > Math.Abs(muB) ? muA : muB);

                for (int n = 0; n < N; n++) { // loop over nodes...
                    double N1 = inp.Normal[j, n, 0];
                    double N2 = inp.Normal[j, n, 1];


                    Koeff_UxV[j, n, 0, iSpcNeg, iSpcNeg] = 0.0;
                    Koeff_UxV[j, n, 0, iSpcPos, iSpcNeg] = 0.0;
                    Koeff_UxV[j, n, 0, iSpcNeg, iSpcPos] = 0.0;
                    Koeff_UxV[j, n, 0, iSpcPos, iSpcPos] = 0.0;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcNeg, 0] = 0.0;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcNeg, 0] = 0.0;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcPos, 0] = 0.0;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcPos, 0] = 0.0;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcNeg, 1] = -1.0*N_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcNeg, 1] = N_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcNeg, iSpcPos, 1] = P_mu*N1;
                    Koeff_NablaUxV[j, n, 0, iSpcPos, iSpcPos, 1] = -1.0*P_mu*N1;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcNeg, 0] = -1.0*N_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcNeg, 0] = P_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcPos, 0] = N_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcPos, 0] = -1.0*P_mu*N2;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcNeg, 1] = 0.0;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcNeg, 1] = 0.0;
                    Koeff_UxNablaV[j, n, 0, iSpcNeg, iSpcPos, 1] = 0.0;
                    Koeff_UxNablaV[j, n, 0, iSpcPos, iSpcPos, 1] = 0.0;
                    Koeff_UxV[j, n, 1, iSpcNeg, iSpcNeg] = -1.0*eta;
                    Koeff_UxV[j, n, 1, iSpcPos, iSpcNeg] = eta;
                    Koeff_UxV[j, n, 1, iSpcNeg, iSpcPos] = eta;
                    Koeff_UxV[j, n, 1, iSpcPos, iSpcPos] = -1.0*eta;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcNeg, 0] = -.500*N_mu*N1;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcNeg, 0] = .500*N_mu*N1;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcPos, 0] = -.500*P_mu*N1;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcPos, 0] = .500*P_mu*N1;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcNeg, 1] = -1.50*N_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcNeg, 1] = 1.50*N_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcNeg, iSpcPos, 1] = .500*P_mu*N2;
                    Koeff_NablaUxV[j, n, 1, iSpcPos, iSpcPos, 1] = -.500*P_mu*N2;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcNeg, 0] = -.500*N_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcNeg, 0] = -.500*P_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcPos, 0] = .500*N_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcPos, 0] = .500*P_mu*N1;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcNeg, 1] = -1.50*N_mu*N2;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcNeg, 1] = .500*P_mu*N2;
                    Koeff_UxNablaV[j, n, 1, iSpcNeg, iSpcPos, 1] = 1.50*N_mu*N2;
                    Koeff_UxNablaV[j, n, 1, iSpcPos, iSpcPos, 1] = -.500*P_mu*N2;
                }
            }
        }

        */
    }


    /// <summary>
    /// 
    /// </summary>
    public class GeneralizedViscosityAtLevelSet_FullySymmetric : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public GeneralizedViscosityAtLevelSet_FullySymmetric(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component, double _rhoA, double _rhoB, double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma, double _pc) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penalty = _penalty;
            this.component = _component;
            this.m_D = lstrk.GridDat.SpatialDimension;

            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            //this.M = _M;
            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.Rint = _Rint;
            //this.TintMin = _TintMin;
            this.Tsat = _Tsat;
            this.sigma = _sigma;
            this.pc = _pc;
        }

        double muA;
        double muB;
        double penalty;
        int component;
        int m_D;

        double rhoA;
        double rhoB;

        double kA;
        double kB;
        double hVapA;   // for the identification of the liquid phase
        double Rint;
        //double TintMin;
        double Tsat;
        double sigma;
        double pc;

        //double M;


        private double ComputeEvaporationMass_Macro(double[] GradT_A, double[] GradT_B, double[] n) {

            double hVap = 0.0;
            double qEvap = 0.0;
            if(hVapA > 0) {
                hVap = hVapA;
                for(int d = 0; d < m_D; d++)
                    qEvap += (kA * GradT_A[d] - kB * GradT_B[d]) * n[d];
            } else {
                hVap = -hVapA;
                for(int d = 0; d < m_D; d++)
                    qEvap += (kB * GradT_B[d] - kA * GradT_A[d]) * n[d];
            }

            return qEvap / hVap;
        }

        private double ComputeEvaporationMass_Micro(double T_A, double T_B, double curv, double p_disp) {

            if(hVapA == 0.0)
                return 0.0;

            double pc0 = (pc < 0.0) ? sigma * curv + p_disp : pc;      // augmented capillary pressure (without nonlinear evaporative masss part)

            double TintMin = 0.0;
            double hVap = 0.0;
            double qEvap = 0.0;
            if(hVapA > 0) {
                hVap = hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rhoA)));
                if(T_A > TintMin)
                    qEvap = -(T_A - TintMin) / Rint;
            } else if(hVapA < 0) {
                hVap = -hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rhoB)));
                if(T_B > TintMin)
                    qEvap = (T_B - TintMin) / Rint;
            }

            return qEvap / hVap;
        }


        private double ComputeEvaporationMass(double[] paramsNeg, double[] paramsPos, double[] N, bool microRegion) {

            double M = 0.0;
            if(microRegion) {
                M = ComputeEvaporationMass_Micro(paramsNeg[m_D], paramsPos[m_D], paramsNeg[m_D + 1], paramsNeg[m_D + 2]);
            } else {
                M = ComputeEvaporationMass_Macro(paramsNeg.GetSubVector(0, m_D), paramsPos.GetSubVector(0, m_D), N);
            }

            return M;

        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            //Debug.Assert(this.ArgumentOrdering.Count == D);
            //Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            //Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            //Debug.Assert(Grad_uA.GetLength(1) == D);
            //Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for(int d = 0; d < D; d++) {
                //Grad_uA_xN += Grad_uA[component, d] * N[d];
                //Grad_uB_xN += Grad_uB[component, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }
            double Ret = 0.0;

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            if(hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            //double M = ComputeEvaporationMass_Macro(inp.ParamsNeg.GetSubVector(0, m_D), inp.ParamsPos.GetSubVector(0, m_D), N);
            //double M = ComputeEvaporationMass_Micro(inp.ParamsNeg[m_D], inp.ParamsPos[m_D], inp.ParamsNeg[m_D + 1], inp.ParamsNeg[m_D + 2]);
            double M = -0.1; // ComputeEvaporationMass(inp.ParamsNeg, inp.ParamsPos, N, evapMicroRegion[inp.jCell]);
            if(M == 0.0)
                return 0.0;

            Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
            Debug.Assert(uB.Length == this.ArgumentOrdering.Count);
            //switch (m_ViscosityImplementation) {
            //    // old Form (H-Implementation)
            //    case ViscosityImplementation.H: {
            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret += 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * M * ((1/rhoA) - (1/rhoB)) * N[component];     // symmetry term
            Ret -= (penalty / hCutCellMin) * M * ((1 / rhoA) - (1 / rhoB)) * N[component] * (vA - vB) * muMax; // penalty term
                                                                                                               // Transpose Term
            for(int i = 0; i < D; i++) {
                //Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
                Ret += 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * N[component] * M * ((1 / rhoA) - (1 / rhoB)) * N[i];
            }

            return -Ret;
        }


        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        BitArray evapMicroRegion;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;

            if(csA.UserDefinedValues.Keys.Contains("EvapMicroRegion"))
                evapMicroRegion = (BitArray)csA.UserDefinedValues["EvapMicroRegion"];

        }

        //private static bool rem = true;

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.GradV | TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(new string[] { "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, m_D), VariableNames.Temperature, "Curvature", "DisjoiningPressure"); //;
            }
        }


    }



    //public class ViscosityAtLevelSet_Slip : BoSSS.Foundation.XDG.ILevelSetForm, ILevelSetEquationComponentCoefficient {

    //    LevelSetTracker m_LsTrk;

    //    public ViscosityAtLevelSet_Slip(LevelSetTracker lstrk, double _muA, double _muB, double _penalty, int _component) {
    //        this.m_LsTrk = lstrk;
    //        this.muA = _muA;
    //        this.muB = _muB;
    //        this.penalty = _penalty;
    //        this.component = _component;
    //        this.m_D = lstrk.GridDat.SpatialDimension;
    //    }

    //    double muA;
    //    double muB;
    //    double penalty;
    //    int component;
    //    int m_D;


    //    /// <summary>
    //    /// 
    //    /// </summary>
    //    public double LevelSetForm(ref CommonParamsLs inp,
    //        double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
    //        double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
    //        double[] N = inp.n;
    //        double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

    //        int D = N.Length;
    //        Debug.Assert(this.ArgumentOrdering.Count == D);
    //        Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
    //        Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
    //        Debug.Assert(Grad_uA.GetLength(1) == D);
    //        Debug.Assert(Grad_uB.GetLength(1) == D);


    //        double Ret = 0.0;

    //        double PosCellLengthScale = PosLengthScaleS[inp.jCell];
    //        double NegCellLengthScale = NegLengthScaleS[inp.jCell];

    //        double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
    //        if(hCutCellMin <= 1.0e-10 * hCellMin)
    //            // very small cell -- clippling
    //            hCutCellMin = hCellMin;

    //        Debug.Assert(uA.Length == this.ArgumentOrdering.Count);
    //        Debug.Assert(uB.Length == this.ArgumentOrdering.Count);

    //        double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;

    //        if(NegSlipLengths[inp.jCell] == 0.0 && NegSlipLengths[inp.jCell] == 0.0) {

    //            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
    //            for(int d = 0; d < D; d++) {
    //                Grad_uA_xN += Grad_uA[component, d] * N[d];
    //                Grad_uB_xN += Grad_uB[component, d] * N[d];
    //                Grad_vA_xN += Grad_vA[d] * N[d];
    //                Grad_vB_xN += Grad_vB[d] * N[d];
    //            }

    //            Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
    //            Ret -= 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * (uA[component] - uB[component]);     // symmetry term
    //            Ret += (penalty / hCutCellMin) * (uA[component] - uB[component]) * (vA - vB) * muMax; // penalty term
    //                                                                                                  // Transpose Term
    //            for(int i = 0; i < D; i++) {
    //                Ret -= 0.5 * (muA * Grad_uA[i, component] + muB * Grad_uB[i, component]) * (vA - vB) * N[i];  // consistency term
    //                Ret -= 0.5 * (muA * Grad_vA[i] + muB * Grad_vB[i]) * (uA[i] - uB[i]) * N[component];  // symmetry term
    //            }

    //        } else {

    //            for(int dN = 0; dN < m_D; dN++) {
    //                for(int dD = 0; dD < m_D; dD++) {
    //                    Ret -= 0.5 * N[dN] * (muA * Grad_uA[dD, dN] + muB * Grad_uB[dN, dD]) * N[dD] * (vA - vB) * N[component];     // consistency term
    //                    Ret -= 0.5 * N[dN] * (muA * Grad_vA[dN] + muB * Grad_vB[dD]) * N[component] * (uA[dD] - uB[dN]) * N[dD];    // symmetry term
    //                                                                                                                                // transposed term
    //                    Ret -= 0.5 * N[dN] * (muA * Grad_uA[dN, dD] + muB * Grad_uB[dN, dD]) * N[dD] * (vA - vB) * N[component];    // consistency term
    //                    Ret -= 0.5 * N[component] * (muA * Grad_vA[dD] + muB * Grad_vB[dD]) * N[dD] * (uA[dN] - uB[dN]) * N[dN];    // symmetry term
    //                }
    //                Ret += (penalty / hCutCellMin) * (uA[dN] - uB[dN]) * (vA - vB) * N[dN] * muMax;     // penalty term
    //            }
    //        }

    //        return Ret;
    //    }


    //    MultidimensionalArray PosLengthScaleS;
    //    MultidimensionalArray NegLengthScaleS;

    //    MultidimensionalArray PosSlipLengths;
    //    MultidimensionalArray NegSlipLengths;

    //    public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
    //        NegLengthScaleS = csA.CellLengthScales;
    //        PosLengthScaleS = csB.CellLengthScales;

    //        NegSlipLengths = (MultidimensionalArray)csA.UserDefinedValues["SlipLengths"];
    //        PosSlipLengths = (MultidimensionalArray)csB.UserDefinedValues["SlipLengths"];
    //    }


    //    public int LevelSetIndex {
    //        get { return 0; }
    //    }

    //    public IList<string> ArgumentOrdering {
    //        get { return VariableNames.VelocityVector(this.m_D); }
    //    }

    //    public SpeciesId PositiveSpecies {
    //        get { return m_LsTrk.GetSpeciesId("B"); }
    //    }

    //    public SpeciesId NegativeSpecies {
    //        get { return m_LsTrk.GetSpeciesId("A"); }
    //    }

    //    public TermActivationFlags LevelSetTerms {
    //        get {
    //            return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
    //        }
    //    }

    //    public IList<string> ParameterOrdering {
    //        get { return null; }
    //    }


    //}

}
