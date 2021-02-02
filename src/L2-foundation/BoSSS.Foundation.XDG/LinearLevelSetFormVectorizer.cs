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
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG {

    

    class LinearLevelSetFormVectorizer :  
        ILevelSetForm_V, 
        ILevelSetForm_GradV, 
        ILevelSetForm_GradUxGradV, 
        ILevelSetForm_GradUxV, 
        ILevelSetForm_UxGradV, 
        ILevelSetForm_UxV,
        ILevelSetFormSetup //
    {

        //LevelSetTracker lsTrk;

        public void Setup(LevelSetTracker _lsTrk) {
            //this.lsTrk = _lsTrk;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        public LinearLevelSetFormVectorizer(ILevelSetForm _OrgComponent, LevelSetTracker _lsTrk) {
            this.ArgumentOrdering = _OrgComponent.ArgumentOrdering.ToArray();
            this.ParameterOrdering = _OrgComponent.ParameterOrdering != null ? _OrgComponent.ParameterOrdering.ToArray() : null;
            this.LevelSetIndex = _OrgComponent.LevelSetIndex;
            this.PositiveSpecies = _OrgComponent.PositiveSpecies;
            this.NegativeSpecies = _OrgComponent.NegativeSpecies;
            this.LevelSetTerms = _OrgComponent.LevelSetTerms;
            this.OrgComponent = _OrgComponent;
            //this.lsTrk = _lsTrk;
        }

        /// <summary>
        /// The original component that is beeing vectorized.
        /// </summary>
        ILevelSetForm OrgComponent;

        private double GetCoeff(ref double d1, ref double d2, ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, ref double vA, ref double vB, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(this.OrgComponent.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) == 0.0);

            d2 = 1; // set test function
            double a1 = this.OrgComponent.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the source/affine part
            Debug.Assert((this.LevelSetTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);
            d1 = 1; // set trial function
            double a = this.OrgComponent.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the bilinear+affine part
            d1 = 0; // reset
            d2 = 0; // reset
            return a - a1;
        }

        public override string ToString() {
            return this.OrgComponent.GetType().FullName;
        }

        private double GetSourceCoeff(ref double d2, ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, ref double vA, ref double vB, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(this.OrgComponent.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) == 0.0);

            d2 = 1; // set test function
            double a1 = this.OrgComponent.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the source part
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

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.None; }
        }
        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public IList<string> ArgumentOrdering {
            get;
            private set;
        }

        public IList<string> ParameterOrdering {
            get;
            private set;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            throw new NotSupportedException("Should not be called.");
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException("Should not be called.");
        }

        void IInnerEdgeform_UxV.InternalEdge_UxV(ref EdgeFormParams inp, MultidimensionalArray Koeff_UxV) {
            //(LevSetIntParams inp, MultidimensionalArray Koeff_UxV) {
            int j0 = inp.e0;
            int Len = inp.Len;
            int N = inp.Nodes.GetLength(1); // nodes per cell
            int D = inp.Nodes.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            
            // check dimension of input array
            Koeff_UxV.CheckLengths(Len, N, 2, 2, NoOfVars);


            // create temp mem:
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParams cp;
            cp.Normal = new Vector(D);
            cp.X = new Vector(D);
            cp.Parameters_IN = ParamsNeg;
            cp.Parameters_OUT = ParamsPos;
            cp.time = inp.time;
            cp.iEdge = -123456;
            cp.GridDat = inp.GridDat;

            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];

             //var Reg = lsTrk.Regions;
            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...

                    int iSpcNeg = 0;
                    int iSpcPos = 1;

                    cp.jCellIn = j + inp.e0;
                    cp.jCellOut = cp.jCellIn;
                    
                    for (int n = 0; n < N; n++) { // loop over nodes...
                        cp.Normal.SetFrom(inp.Normals, j, n);
                        cp.X.SetFrom(inp.Nodes, j, n);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParameterVars_OUT[i][j, n];
                            ParamsNeg[i] = inp.ParameterVars_IN[i][j, n];
                        }

                        Koeff_UxV[j, n, iSpcNeg, iSpcNeg, c] = GetCoeff(ref uA[c], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_UxV[j, n, iSpcPos, iSpcNeg, c] = GetCoeff(ref uA[c], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_UxV[j, n, iSpcNeg, iSpcPos, c] = GetCoeff(ref uB[c], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_UxV[j, n, iSpcPos, iSpcPos, c] = GetCoeff(ref uB[c], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                    }
                }
            }

        }

        void IInnerEdgeform_UxGradV.InternalEdge_UxGradV(ref EdgeFormParams inp, MultidimensionalArray Koeff_UxNablaV) { 
            int j0 = inp.e0;
            int Len = inp.Len;
            int N = inp.Nodes.GetLength(1); // nodes per cell
            int D = inp.Nodes.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            
            // check dimension of input array
            Koeff_UxNablaV.CheckLengths(Len, N, 2, 2, NoOfVars, D);
            

            // create temp mem:
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParams cp;
            cp.Normal = new Vector(D);
            cp.X = new Vector(D);
            cp.Parameters_IN = ParamsNeg;
            cp.Parameters_OUT = ParamsPos;
            cp.time = inp.time;
            cp.iEdge = -123456;
            cp.GridDat = inp.GridDat;

            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];

            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...
                    int iSpcNeg = 0;
                    int iSpcPos = 1;

                    cp.jCellIn = j + inp.e0;
                    cp.jCellOut = cp.jCellIn;
                    //if (inp.PosCellLengthScale != null)
                    //    cp.PosCellLengthScale = inp.PosCellLengthScale[cp.jCell];
                    //else
                    //    cp.PosCellLengthScale = double.NaN;
                    //if (inp.NegCellLengthScale != null)
                    //    cp.NegCellLengthScale = inp.NegCellLengthScale[cp.jCell];
                    //else
                    //    cp.NegCellLengthScale = double.NaN;

                    for (int n = 0; n < N; n++) { // loop over nodes...
                        cp.Normal.SetFrom(inp.Normals, j, n);
                        cp.X.SetFrom(inp.Nodes, j, n);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParameterVars_OUT[i][j, n];
                            ParamsNeg[i] = inp.ParameterVars_OUT[i][j, n];
                        }

                        

                        for (int d = 0; d < D; d++) {
                            Koeff_UxNablaV[j, n, iSpcNeg, iSpcNeg, c, d] = GetCoeff(ref uA[c], ref Grad_vA[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_UxNablaV[j, n, iSpcPos, iSpcNeg, c, d] = GetCoeff(ref uA[c], ref Grad_vB[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_UxNablaV[j, n, iSpcNeg, iSpcPos, c, d] = GetCoeff(ref uB[c], ref Grad_vA[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_UxNablaV[j, n, iSpcPos, iSpcPos, c, d] = GetCoeff(ref uB[c], ref Grad_vB[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        }
                    }
                }
            }
        }

        void IInnerEdgeform_GradUxV.InternalEdge_GradUxV(ref EdgeFormParams inp, MultidimensionalArray Koeff_NablaUxV) { 
            //.LevelSetForm_GradUxV(LevSetIntParams inp, MultidimensionalArray Koeff_NablaUxV) {
            int j0 = inp.e0;
            int Len = inp.Len;
            int N = inp.Nodes.GetLength(1); // nodes per cell
            int D = inp.Nodes.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            
            // check dimension of input array
            Koeff_NablaUxV.CheckLengths(Len, N, 2, 2, NoOfVars, D);


            // create temp mem:
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParams cp;
            cp.Normal = new Vector(D);
            cp.X = new Vector(D);
            cp.Parameters_IN = ParamsNeg;
            cp.Parameters_OUT = ParamsPos;
            cp.time = inp.time;
            cp.iEdge = -123456;
            cp.GridDat = inp.GridDat;

            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];

            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...
                    int iSpcNeg = 0;
                    int iSpcPos = 1;

                    cp.jCellIn = j + inp.e0;
                    cp.jCellOut = cp.jCellIn;

                    for (int n = 0; n < N; n++) { // loop over nodes...
                        cp.Normal.SetFrom(inp.Normals, j, n);
                        cp.X.SetFrom(inp.Nodes, j, n);
                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParameterVars_OUT[i][j, n];
                            ParamsNeg[i] = inp.ParameterVars_IN[i][j, n];
                        }

                        
                        for (int d = 0; d < D; d++) {
                            Koeff_NablaUxV[j, n, iSpcNeg, iSpcNeg, c, d] = GetCoeff(ref Grad_uA[c, d], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_NablaUxV[j, n, iSpcPos, iSpcNeg, c, d] = GetCoeff(ref Grad_uA[c, d], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_NablaUxV[j, n, iSpcNeg, iSpcPos, c, d] = GetCoeff(ref Grad_uB[c, d], ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            Koeff_NablaUxV[j, n, iSpcPos, iSpcPos, c, d] = GetCoeff(ref Grad_uB[c, d], ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        }
                    }
                }
            }
        }

        void IInnerEdgeform_GradUxGradV.InternalEdge_GradUxGradV(ref EdgeFormParams inp, MultidimensionalArray Koeff_NablaUxNablaV) { 
        //void ILevelSetForm_GradUxGradV.LevelSetForm_GradUxGradV(LevSetIntParams inp, MultidimensionalArray Koeff_NablaUxNablaV) {
            int j0 = inp.e0;
            int Len = inp.Len;
            int N = inp.Nodes.GetLength(1); // nodes per cell
            int D = inp.Nodes.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            
            // check dimension of input array
            Koeff_NablaUxNablaV.CheckLengths(Len, N, 2, 2, NoOfVars, D, D);


            // create temp mem:
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParams cp;
            cp.Normal = new Vector(D);
            cp.X = new Vector(D);
            cp.Parameters_IN = ParamsNeg;
            cp.Parameters_OUT = ParamsPos;
            cp.time = inp.time;
            cp.iEdge = -123456;
            cp.GridDat = inp.GridDat;

            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];

           
            for (int c = 0; c < NoOfVars; c++) { // loop over variables...
                for (int j = 0; j < Len; j++) { // loop over items...

                    int iSpcNeg = 0;
                    int iSpcPos = 1;

                    cp.jCellIn = j + inp.e0;
                    cp.jCellOut = cp.jCellIn;

                    for (int n = 0; n < N; n++) { // loop over nodes...
                        cp.Normal.SetFrom(inp.Normals, j, n);
                        cp.X.SetFrom(inp.Nodes, j, n);

                        for (int i = 0; i < NP; i++) {
                            ParamsPos[i] = inp.ParameterVars_OUT[i][j, n];
                            ParamsNeg[i] = inp.ParameterVars_IN[i][j, n];
                        }
                        
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                Koeff_NablaUxNablaV[j, n, iSpcNeg, iSpcNeg, c, d1, d2] = GetCoeff(ref Grad_uA[c, d1], ref Grad_vA[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                                Koeff_NablaUxNablaV[j, n, iSpcPos, iSpcNeg, c, d1, d2] = GetCoeff(ref Grad_uA[c, d1], ref Grad_vB[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                                Koeff_NablaUxNablaV[j, n, iSpcNeg, iSpcPos, c, d1, d2] = GetCoeff(ref Grad_uB[c, d1], ref Grad_vA[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                                Koeff_NablaUxNablaV[j, n, iSpcPos, iSpcPos, c, d1, d2] = GetCoeff(ref Grad_uB[c, d1], ref Grad_vB[d2], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                            }
                        }
                    }
                }
            }
        }

        void IInnerEdgeSource_GradV.InternalEdge_GradV(ref EdgeFormParams inp, MultidimensionalArray Koeff_NablaV) { 
            //void ILevelSetForm_GradV.LevelSetForm_GradV(LevSetIntParams inp, MultidimensionalArray Koeff_NablaV) {
            int j0 = inp.e0;
            int Len = inp.Len;
            int N = inp.Nodes.GetLength(1); // nodes per cell
            int D = inp.Nodes.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
            
            // check dimension of input array
            Koeff_NablaV.CheckLengths(Len, N, 2, D);


            // create temp mem:
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParams cp;
            cp.Normal = new Vector(D);
            cp.X = new Vector(D);
            cp.Parameters_IN = ParamsNeg;
            cp.Parameters_OUT = ParamsPos;
            cp.time = inp.time;
            cp.iEdge = -123456;
            cp.GridDat = inp.GridDat;

            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];
            //var h_min = this.m_LsTrk.GridDat.Cells.h_min;
                        
            for (int j = 0; j < Len; j++) { // loop over items...

                int iSpcNeg = 0;
                int iSpcPos = 1;

                cp.jCellIn = j + inp.e0;
                cp.jCellOut = cp.jCellIn;


                for (int n = 0; n < N; n++) { // loop over nodes...
                    cp.Normal.SetFrom(inp.Normals, j, n);
                    cp.X.SetFrom(inp.Nodes, j, n);
                    for (int i = 0; i < NP; i++) {
                        ParamsPos[i] = inp.ParameterVars_OUT[i][j, n];
                        ParamsNeg[i] = inp.ParameterVars_IN[i][j, n];
                    }
                    
                    for (int d = 0; d < D; d++) {
                        Koeff_NablaV[j, n, iSpcNeg, d] = GetSourceCoeff(ref Grad_vA[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                        Koeff_NablaV[j, n, iSpcPos, d] = GetSourceCoeff(ref Grad_vB[d], ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                    }
                }
            }
        }

        void IInnerEdgeSource_V.InternalEdge_V(ref EdgeFormParams inp, MultidimensionalArray Koeff_V) {
            int j0 = inp.e0;
            int Len = inp.Len;
            int N = inp.Nodes.GetLength(1); // nodes per cell
            int D = inp.Nodes.GetLength(2); // spatial dim.
            int NoOfVars = this.ArgumentOrdering.Count;
           
            // check dimension of input array
            Koeff_V.CheckLengths(Len, N, 2);


            // create temp mem:
            int NP = (this.ParameterOrdering != null) ? this.ParameterOrdering.Count : 0;

            double[] ParamsPos = new double[NP];
            double[] ParamsNeg = new double[NP];
            CommonParams cp;
            cp.Normal = new Vector(D);
            cp.X = new Vector(D);
            cp.Parameters_IN = ParamsNeg;
            cp.Parameters_OUT = ParamsPos;
            cp.time = inp.time;
            cp.iEdge = -123456;
            cp.GridDat = inp.GridDat;

            // temp mem.
            double[] uA = new double[NoOfVars];
            double[] uB = new double[NoOfVars];
            double[,] Grad_uA = new double[NoOfVars, D];
            double[,] Grad_uB = new double[NoOfVars, D];
            double vA = 0;
            double vB = 0;
            double[] Grad_vA = new double[D];
            double[] Grad_vB = new double[D];

            for (int j = 0; j < Len; j++) { // loop over items...
                int iSpcNeg = 0;
                int iSpcPos = 1;
                cp.jCellIn = j + inp.e0;
                cp.jCellOut = cp.jCellIn;
  
                for (int n = 0; n < N; n++) { // loop over nodes...
                    cp.Normal.SetFrom(inp.Normals, j, n);
                    cp.X.SetFrom(inp.Nodes, j, n);
                    for (int i = 0; i < NP; i++) {
                        ParamsPos[i] = inp.ParameterVars_OUT[i][j, n];
                        ParamsNeg[i] = inp.ParameterVars_IN[i][j, n];
                    }

                    Koeff_V[j, n, iSpcNeg] = GetSourceCoeff(ref vA, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                    Koeff_V[j, n, iSpcPos] = GetSourceCoeff(ref vB, ref cp, uA, uB, Grad_uA, Grad_uB, ref vA, ref vB, Grad_vA, Grad_vB);
                }
            }
        }
    }
}
