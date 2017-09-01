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
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    

    class LinearLevelSetComponentVectorizer :  
        ILinearLevelSetComponent_V, 
        ILinearLevelSetComponent_GradV, 
        ILinearLevelSetComponent_GradUxGradV, 
        ILinearLevelSetComponent_GradUxV, 
        ILinearLevelSetComponent_UxGradV, 
        ILinearLevelSetComponent_UxV //
    {

        /// <summary>
        /// ctor.
        /// </summary>
        public LinearLevelSetComponentVectorizer(LevelSetTracker lsTrk, ILevelSetComponent _OrgComponent) {
            this.m_LsTrk = lsTrk;
            this.ArgumentOrdering = _OrgComponent.ArgumentOrdering.ToArray();
            this.ParameterOrdering = _OrgComponent.ParameterOrdering != null ? _OrgComponent.ParameterOrdering.ToArray() : null;
            this.LevelSetIndex = _OrgComponent.LevelSetIndex;
            this.PositiveSpecies = _OrgComponent.PositiveSpecies;
            this.NegativeSpecies = _OrgComponent.NegativeSpecies;
            this.LevelSetTerms = _OrgComponent.LevelSetTerms;
            this.OrgComponent = _OrgComponent;
        }

        /// <summary>
        /// The original component that is beeing vectorized.
        /// </summary>
        ILevelSetComponent OrgComponent;


        LevelSetTracker m_LsTrk;

        ///// <summary>
        ///// override this method to implement the bilinear form on the edge.
        ///// </summary>
        //public abstract double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp, 
        //    double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
        //    double vA, double vB, double[] Grad_vA, double[] Grad_vB);


        private double GetCoeff(ref double d1, ref double d2, ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, ref double vA, ref double vB, double[] Grad_vA, double[] Grad_vB) {
            //if(this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) != 0.0) {
            //    double funky = this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
            //    Console.WriteLine(funky);
            //}
            Debug.Assert(this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) == 0.0);

            d2 = 1; // set test function
            double a1 = this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the source/affine part
            Debug.Assert((this.LevelSetTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);
            d1 = 1; // set trial function
            double a = this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the bilinear+affine part
            d1 = 0; // reset
            d2 = 0; // reset
            return a - a1;
        }

        private double GetSourceCoeff(ref double d2, ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, ref double vA, ref double vB, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) == 0.0);

            d2 = 1; // set test function
            double a1 = this.OrgComponent.LevelSetForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the source part
            Debug.Assert((this.LevelSetTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);
            d2 = 0; // reset
            return a1;
        }
        
        public int LevelSetIndex {
            get; 
            private set;
        }

        public SpeciesId PositiveSpecies {
            get;
            private set;
        }

        public SpeciesId NegativeSpecies {
            get; 
            private set;
        }



        public TermActivationFlags LevelSetTerms {
            get;
            private set;
        }

        public IList<string> ArgumentOrdering {
            get;
            private set;
        }

        public IList<string> ParameterOrdering {
            get;
            private set;
        }

        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            throw new NotSupportedException("Should not be called.");
        }

        void ILinearLevelSetComponent_UxV.LevelSetForm_UxV(LevSetIntParams inp, MultidimensionalArray Koeff_UxV) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_UxV.CheckLengths(Len, N, NoOfVars, 2, 2);


            // create temp mem:
            double[] node = new double[D];
            double[] normal = new double[D];
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParamsLs cp = default(CommonParamsLs);
            cp.n = normal;
            cp.x = node;
            cp.ParamsNeg = ParamsNeg;
            cp.ParamsPos = ParamsPos;
            cp.time = inp.time;


            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;


            //LevelSetSignCode pos;
            //LevelSetSignCode neg;
            //GetSignCode(out neg, out pos);
            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...

                    ReducedRegionCode rrc;
                    int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                    Debug.Assert(NoOf == 2);
                    //int iSpcPos = lsTrk.GetSpeciesIndex(rrc, pos);
                    //int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, neg);
                    int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                    int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                    cp.jCell = j + inp.i0;
                    if (inp.PosCellLengthScale != null)
                        cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                    else
                        cp.PosCellLengthScale = double.NaN;
                    if (inp.NegCellLengthScale != null)
                        cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                    else
                        cp.NegCellLengthScale = double.NaN;
                    
                    for (int n = 0; n < N; n++) { // loop over nodes...
                        inp.Normal.CopyTo(normal, j, n, -1);
                        inp.X.CopyTo(node, j, n, -1);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParamsPos[i][j, n];
                            ParamsNeg[i] = inp.ParamsNeg[i][j, n];
                        }

                        Koeff_UxV[j, n, c, iSpcNeg, iSpcNeg] = GetCoeff(ref uA[c], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_UxV[j, n, c, iSpcPos, iSpcNeg] = GetCoeff(ref uA[c], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_UxV[j, n, c, iSpcNeg, iSpcPos] = GetCoeff(ref uB[c], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_UxV[j, n, c, iSpcPos, iSpcPos] = GetCoeff(ref uB[c], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                    }
                }
            }

        }

        void ILinearLevelSetComponent_UxGradV.LevelSetForm_UxGradV(LevSetIntParams inp, MultidimensionalArray Koeff_UxNablaV) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_UxNablaV.CheckLengths(Len, N, NoOfVars, 2, 2, D);
            

            // create temp mem:
            double[] node = new double[D];
            double[] normal = new double[D];
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParamsLs cp = default(CommonParamsLs);
            cp.n = normal;
            cp.x = node;
            cp.ParamsNeg = ParamsNeg;
            cp.ParamsPos = ParamsPos;
            cp.time = inp.time;


            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;


            //LevelSetSignCode pos;
            //LevelSetSignCode neg;
            //GetSignCode(out neg, out pos);
            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...

                    ReducedRegionCode rrc;
                    int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                    Debug.Assert(NoOf == 2);
                    //int iSpcPos = lsTrk.GetSpeciesIndex(rrc, pos);
                    //int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, neg);
                    int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                    int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                    cp.jCell = j + inp.i0;
                    if (inp.PosCellLengthScale != null)
                        cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                    else
                        cp.PosCellLengthScale = double.NaN;
                    if (inp.NegCellLengthScale != null)
                        cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                    else
                        cp.NegCellLengthScale = double.NaN;

                    for (int n = 0; n < N; n++) { // loop over nodes...
                        inp.Normal.CopyTo(normal, j, n, -1);
                        inp.X.CopyTo(node, j, n, -1);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParamsPos[i][j, n];
                            ParamsNeg[i] = inp.ParamsNeg[i][j, n];
                        }

                        

                        for (int d = 0; d < D; d++) {
                            Koeff_UxNablaV[j, n, c, iSpcNeg, iSpcNeg, d] = GetCoeff(ref uA[c], ref Grad_vA[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_UxNablaV[j, n, c, iSpcPos, iSpcNeg, d] = GetCoeff(ref uA[c], ref Grad_vB[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_UxNablaV[j, n, c, iSpcNeg, iSpcPos, d] = GetCoeff(ref uB[c], ref Grad_vA[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_UxNablaV[j, n, c, iSpcPos, iSpcPos, d] = GetCoeff(ref uB[c], ref Grad_vB[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        }
                    }
                }
            }
        }

        void ILinearLevelSetComponent_GradUxV.LevelSetForm_GradUxV(LevSetIntParams inp, MultidimensionalArray Koeff_NablaUxV) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_NablaUxV.CheckLengths(Len, N, NoOfVars, 2, 2, D);


            // create temp mem:
            double[] node = new double[D];
            double[] normal = new double[D];
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParamsLs cp = default(CommonParamsLs);
            cp.n = normal;
            cp.x = node;
            cp.ParamsNeg = ParamsNeg;
            cp.ParamsPos = ParamsPos;
            cp.time = inp.time;


            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;


            //LevelSetSignCode pos;
            //LevelSetSignCode neg;
            //GetSignCode(out neg, out pos);
            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...

                    ReducedRegionCode rrc;
                    int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                    Debug.Assert(NoOf == 2);
                    //int iSpcPos = lsTrk.GetSpeciesIndex(rrc, pos);
                    //int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, neg);
                    int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                    int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                    cp.jCell = j + inp.i0;
                    if (inp.PosCellLengthScale != null)
                        cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                    else
                        cp.PosCellLengthScale = double.NaN;
                    if (inp.NegCellLengthScale != null)
                        cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                    else
                        cp.NegCellLengthScale = double.NaN;

                    for (int n = 0; n < N; n++) { // loop over nodes...
                        inp.Normal.CopyTo(normal, j, n, -1);
                        inp.X.CopyTo(node, j, n, -1);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParamsPos[i][j, n];
                            ParamsNeg[i] = inp.ParamsNeg[i][j, n];
                        }

                        
                        for (int d = 0; d < D; d++) {
                            Koeff_NablaUxV[j, n, c, iSpcNeg, iSpcNeg, d] = GetCoeff(ref Grad_uA[c, d], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_NablaUxV[j, n, c, iSpcPos, iSpcNeg, d] = GetCoeff(ref Grad_uA[c, d], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_NablaUxV[j, n, c, iSpcNeg, iSpcPos, d] = GetCoeff(ref Grad_uB[c, d], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_NablaUxV[j, n, c, iSpcPos, iSpcPos, d] = GetCoeff(ref Grad_uB[c, d], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        }
                    }
                }
            }
        }

        void ILinearLevelSetComponent_GradUxGradV.LevelSetForm_GradUxGradV(LevSetIntParams inp, MultidimensionalArray Koeff_NablaUxNablaV) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_NablaUxNablaV.CheckLengths(Len, N, NoOfVars, 2, 2, D, D);


            // create temp mem:
            double[] node = new double[D];
            double[] normal = new double[D];
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParamsLs cp = default(CommonParamsLs);
            cp.n = normal;
            cp.x = node;
            cp.ParamsNeg = ParamsNeg;
            cp.ParamsPos = ParamsPos;
            cp.time = inp.time;


            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;


            //LevelSetSignCode pos;
            //LevelSetSignCode neg;
            //GetSignCode(out neg, out pos);
            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...

                    ReducedRegionCode rrc;
                    int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                    Debug.Assert(NoOf == 2);
                    //int iSpcPos = lsTrk.GetSpeciesIndex(rrc, pos);
                    //int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, neg);
                    int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                    int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                    cp.jCell = j + inp.i0;
                    if (inp.PosCellLengthScale != null)
                        cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                    else
                        cp.PosCellLengthScale = double.NaN;
                    if (inp.NegCellLengthScale != null)
                        cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                    else
                        cp.NegCellLengthScale = double.NaN;

                    for (int n = 0; n < N; n++) { // loop over nodes...
                        inp.Normal.CopyTo(normal, j, n, -1);
                        inp.X.CopyTo(node, j, n, -1);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParamsPos[i][j, n];
                            ParamsNeg[i] = inp.ParamsNeg[i][j, n];
                        }
                        
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                Koeff_NablaUxNablaV[j, n, c, iSpcNeg, iSpcNeg, d1, d2] = GetCoeff(ref Grad_uA[c, d1], ref Grad_vA[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                                Koeff_NablaUxNablaV[j, n, c, iSpcPos, iSpcNeg, d1, d2] = GetCoeff(ref Grad_uA[c, d1], ref Grad_vB[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                                Koeff_NablaUxNablaV[j, n, c, iSpcNeg, iSpcPos, d1, d2] = GetCoeff(ref Grad_uB[c, d1], ref Grad_vA[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                                Koeff_NablaUxNablaV[j, n, c, iSpcPos, iSpcPos, d1, d2] = GetCoeff(ref Grad_uB[c, d1], ref Grad_vB[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            }
                        }
                    }
                }
            }
        }

        void ILinearLevelSetComponent_GradV.LevelSetForm_GradV(LevSetIntParams inp, MultidimensionalArray Koeff_NablaV) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_NablaV.CheckLengths(Len, N, 2, D);


            // create temp mem:
            double[] node = new double[D];
            double[] normal = new double[D];
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParamsLs cp = default(CommonParamsLs);
            cp.n = normal;
            cp.x = node;
            cp.ParamsNeg = ParamsNeg;
            cp.ParamsPos = ParamsPos;
            cp.time = inp.time;


            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;
                        
            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            for (int j = 0; j < Len; j++) { // loop over items...

                ReducedRegionCode rrc;
                int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                Debug.Assert(NoOf == 2);
                //int iSpcPos = lsTrk.GetSpeciesIndex(rrc, pos);
                //int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, neg);
                int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                cp.jCell = j + inp.i0;
                if (inp.PosCellLengthScale != null)
                    cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                else
                    cp.PosCellLengthScale = double.NaN;
                if (inp.NegCellLengthScale != null)
                    cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                else
                    cp.NegCellLengthScale = double.NaN;

                for (int n = 0; n < N; n++) { // loop over nodes...
                    inp.Normal.CopyTo(normal, j, n, -1);
                    inp.X.CopyTo(node, j, n, -1);
                    for (int i = 0; i < NP; i++) {
                        ParamsPos[i] = inp.ParamsPos[i][j, n];
                        ParamsNeg[i] = inp.ParamsNeg[i][j, n];
                    }
                    
                    for (int d = 0; d < D; d++) {
                        Koeff_NablaV[j, n, iSpcNeg, d] = GetSourceCoeff(ref Grad_vA[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_NablaV[j, n, iSpcPos, d] = GetSourceCoeff(ref Grad_vB[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                    }
                }
            }
        }

        void ILinearLevelSetComponent_V.LevelSetForm_V(LevSetIntParams inp, MultidimensionalArray Koeff_V) {
            int j0 = inp.i0;
            int Len = inp.Len;
            int N = inp.X.GetLength(1); // nodes per cell
            int D = inp.X.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            LevelSetTracker lsTrk = m_LsTrk;

            // check dimension of input array
            Koeff_V.CheckLengths(Len, N, 2);


            // create temp mem:
            double[] node = new double[D];
            double[] normal = new double[D];
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParamsLs cp = default(CommonParamsLs);
            cp.n = normal;
            cp.x = node;
            cp.ParamsNeg = ParamsNeg;
            cp.ParamsPos = ParamsPos;
            cp.time = inp.time;


            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            var h_min = this.m_LsTrk.GridDat.Cells.h_min;


            //LevelSetSignCode pos;
            //LevelSetSignCode neg;
            //GetSignCode(out neg, out pos);
            SpeciesId posSpc = this.PositiveSpecies;
            SpeciesId negSpc = this.NegativeSpecies;

            for (int j = 0; j < Len; j++) { // loop over items...
                ReducedRegionCode rrc;
                int NoOf = lsTrk.GetNoOfSpecies(j + inp.i0, out rrc);
                Debug.Assert(NoOf == 2);
                //int iSpcPos = lsTrk.GetSpeciesIndex(rrc, pos);
                //int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, neg);
                int iSpcPos = lsTrk.GetSpeciesIndex(rrc, posSpc);
                int iSpcNeg = lsTrk.GetSpeciesIndex(rrc, negSpc);
                cp.jCell = j + inp.i0;
                if (inp.PosCellLengthScale != null)
                    cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                else
                    cp.PosCellLengthScale = double.NaN;
                if (inp.NegCellLengthScale != null)
                    cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                else
                    cp.NegCellLengthScale = double.NaN;

                for (int n = 0; n < N; n++) { // loop over nodes...
                    inp.Normal.CopyTo(normal, j, n, -1);
                    inp.X.CopyTo(node, j, n, -1);
                    for (int i = 0; i < NP; i++) {
                        ParamsPos[i] = inp.ParamsPos[i][j, n];
                        ParamsNeg[i] = inp.ParamsNeg[i][j, n];
                    }

                    Koeff_V[j, n, iSpcNeg] = GetSourceCoeff(ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                    Koeff_V[j, n, iSpcPos] = GetSourceCoeff(ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                }
            }
        }
    }
}
