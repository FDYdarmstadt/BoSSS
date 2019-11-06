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
using BoSSS.Foundation;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using ilPSP;

namespace BoSSS.Foundation.Quadrature.Linear {


    class LinearVolumeFormVectorizer : IVolumeForm_GradUxGradV, IVolumeForm_UxV, IVolumeForm_UxGradV, IVolumeForm_GradUxV, IVolumeSource_V, IVolumeSource_GradV {

        public TermActivationFlags VolTerms{ 
            get;
            private set;
        }
        
        public LinearVolumeFormVectorizer(IVolumeForm _volForm) {
            this.VolTerms = _volForm.VolTerms;
            this.volForm = _volForm;
        }


        IVolumeForm volForm;

        /// <summary>
        /// Returns name of nested form.
        /// </summary>
        public override string ToString() {
            return volForm.GetType().FullName;
        }


        /// <summary>
        /// tests for components which are bilinear in u and v
        /// </summary>
        private double GetCoeff(ref double TrialVar, ref double TestVar, ref CommonParamsVol inp) {
            int D = inp.GridDat.SpatialDimension;
            int NoArgs = this.ArgumentOrdering.Count;

            Debug.Assert(u.Length == NoArgs);
            Debug.Assert(Grad_u.GetLength(0) == NoArgs);
            Debug.Assert(Grad_u.GetLength(1) == D);
            Debug.Assert(Grad_v.Length == D);

            Debug.Assert(this.volForm.VolumeForm(ref inp, this.u, this.Grad_u, this.v, Grad_v) == 0.0);

            TestVar = 1;
            double a1 = this.volForm.VolumeForm(ref inp, this.u, this.Grad_u, this.v, Grad_v);
            Debug.Assert((this.VolTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);
            TrialVar = 1;
            double a = this.volForm.VolumeForm(ref inp, this.u, this.Grad_u, this.v, Grad_v);
            TrialVar = 0;
            TestVar = 0;
            return a - a1;
        }

        /// <summary>
        /// tests for components which are linear in v and independent of u
        /// </summary>
        private double GetCoeff(ref double TestVar, ref CommonParamsVol inp) {
            int D = inp.GridDat.SpatialDimension;
            int NoArgs = this.ArgumentOrdering.Count;

            Debug.Assert(u.Length == NoArgs);
            Debug.Assert(Grad_u.GetLength(0) == NoArgs);
            Debug.Assert(Grad_u.GetLength(1) == D);
            Debug.Assert(Grad_v.Length == D);

            Debug.Assert(this.volForm.VolumeForm(ref inp, this.u, this.Grad_u, this.v, Grad_v) == 0.0);

            TestVar = 1;
            double a1 = this.volForm.VolumeForm(ref inp, this.u, this.Grad_u, this.v, Grad_v);
            Debug.Assert((this.VolTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);
            TestVar = 0;
            return a1;
        }
        


        /// <summary>
        /// ordering of argument variables
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return this.volForm.ArgumentOrdering;
            }
        }

        /// <summary>
        /// ordering of parameter variables
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return this.volForm.ParameterOrdering;
            }
        }
                
        // "global" variables
        int D;
        int K;
        int L;
        int NoArgs;
        int NoParams;
        double[] u;
        double[,] Grad_u;
        double[] Grad_v;
        double v = 0.0;


        private void InitGlobals(VolumFormParams efp, int __K) {
            D = efp.GridDat.SpatialDimension;
            K = __K;
            L = efp.Len;
            NoArgs = this.ArgumentOrdering.Count;
            NoParams = this.ParameterOrdering != null ? this.ParameterOrdering.Count : 0;

            u = new double[NoArgs];
            Grad_u = new double[NoArgs, D];
            Grad_v = new double[D];
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this.volForm.VolumeForm(ref cpv, U, GradU, V, GradV);
        }

        void IVolumeForm_UxV.Form(ref VolumFormParams prm, MultidimensionalArray UxV) {
            Debug.Assert(prm.Len == UxV.GetLength(0));
            int L = prm.Len;
            int __K = UxV.GetLength(1); // no of nodes
            Debug.Assert(UxV.GetLength(2) == this.ArgumentOrdering.Count);
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            int _NOargs = this.ArgumentOrdering.Count;
            this.InitGlobals(prm, __K);
            int D = prm.GridDat.SpatialDimension;


            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new double[D];
            cpv.time = prm.time;

            for(int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for(int k = 0; k < __K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }
                    for(int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }


                    for(int c = 0; c < _NOargs; c++) { // loop over trial (codomain) variable components...
                        UxV[l, k, c] = GetCoeff(ref this.u[c], ref this.v, ref cpv);
                    }
                }
            }
        }

        void IVolumeForm_UxGradV.Form(ref VolumFormParams prm, MultidimensionalArray UxGradV) {
            Debug.Assert(prm.Len == UxGradV.GetLength(0));
            int L = prm.Len;
            int __K = UxGradV.GetLength(1); // no of nodes
            Debug.Assert(UxGradV.GetLength(2) == this.ArgumentOrdering.Count);
            int D = UxGradV.GetLength(3);
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            int _NOargs = this.ArgumentOrdering.Count;
            this.InitGlobals(prm, __K);

            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new double[D];
            cpv.time = prm.time;

            for(int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for(int k = 0; k < __K; k++) { // loop over nodes...
                    for(int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }
                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }

                    for(int c = 0; c < _NOargs; c++) { // loop over trial (codomain) variable components...

                        for(int d = 0; d < D; d++) { // test variable spatial direction
                            UxGradV[l, k, c, d] = GetCoeff(ref this.u[c], ref this.Grad_v[d], ref cpv);
                        }
                    }
                }
            }
        }

        void IVolumeForm_GradUxV.Form(ref VolumFormParams prm, MultidimensionalArray GradUxV) {
            Debug.Assert(prm.Len == GradUxV.GetLength(0));
            int L = prm.Len;
            int __K = GradUxV.GetLength(1); // no of nodes
            Debug.Assert(GradUxV.GetLength(2) == this.ArgumentOrdering.Count);
            int D = GradUxV.GetLength(3);
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            int _NOargs = this.ArgumentOrdering.Count;
            this.InitGlobals(prm, __K);
            
            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new double[D];
            cpv.time = prm.time;

            for(int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for(int k = 0; k < __K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }
                    for(int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }

                    for(int c = 0; c < _NOargs; c++) { // loop over trial (codomain) variable components...

                        for(int e = 0; e < D; e++) { // trial variable spatial direction
                            GradUxV[l, k, c, e] = GetCoeff(ref this.Grad_u[c, e], ref this.v, ref cpv);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IVolumeForm_GradUxGradV.Form"/>
        /// </summary>
        void IVolumeForm_GradUxGradV.Form(ref VolumFormParams prm, MultidimensionalArray GradUxGradV) {
            Debug.Assert(prm.Len == GradUxGradV.GetLength(0));
            int L = prm.Len;
            int __K = GradUxGradV.GetLength(1); // no of nodes
            Debug.Assert(GradUxGradV.GetLength(2) == this.ArgumentOrdering.Count);
            int D = GradUxGradV.GetLength(3);
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            int _NOargs = this.ArgumentOrdering.Count;
            this.InitGlobals(prm, __K);

            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new double[D];
            cpv.time = prm.time;

            for(int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for(int k = 0; k < __K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }
                    for(int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }

                    for(int c = 0; c < _NOargs; c++) { // loop over trial (codomain) variable components...

                        for(int e = 0; e < D; e++) { // trial variable spatial direction
                            for(int d = 0; d < D; d++) { // test variable spatial direction
                                GradUxGradV[l, k, c, d, e] = GetCoeff(ref this.Grad_u[c, e], ref this.Grad_v[d], ref cpv);
                            }
                        }
                    }
                }
            }
        }

        void IVolumeSource_V.Form(ref VolumFormParams prm, MultidimensionalArray V) {
            Debug.Assert(prm.Len == V.GetLength(0));
            int L = prm.Len;
            int __K = V.GetLength(1); // no of nodes
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            this.InitGlobals(prm, __K);
            int D = prm.GridDat.SpatialDimension;


            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new double[D];
            cpv.time = prm.time;

            for(int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for(int k = 0; k < __K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }
                    for(int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }


                    V[l, k] += GetCoeff(ref this.v, ref cpv);
                    
                }
            }
        }

        void IVolumeSource_GradV.Form(ref VolumFormParams prm, MultidimensionalArray GradV) {
            Debug.Assert(prm.Len == GradV.GetLength(0));
            int L = prm.Len;
            int __K = GradV.GetLength(1); // no of nodes
            int D = prm.GridDat.SpatialDimension;
            Debug.Assert(D == GradV.GetLength(2));
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            this.InitGlobals(prm, __K);

            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new double[D];
            cpv.time = prm.time;

            for(int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for(int k = 0; k < __K; k++) { // loop over nodes...
                    for(int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }
                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }

                    for(int d = 0; d < D; d++) { // test variable spatial direction
                        GradV[l, k, d] += GetCoeff(ref this.Grad_v[d], ref cpv);
                    }
                }
            }
        }
    }


    class LinearEdgeFormVectorizer : IEdgeform_GradUxV, IEdgeform_UxGradV, IEdgeSource_V, IEdgeSource_GradV, IEdgeform_UxV, IEdgeform_GradUxGradV {

        
        public TermActivationFlags BoundaryEdgeTerms { 
            get;
            private set;
        }

        public TermActivationFlags InnerEdgeTerms{ 
            get;
            private set;
        }


        
        public LinearEdgeFormVectorizer(IEdgeForm _edgeForm) {
            this.InnerEdgeTerms = _edgeForm.InnerEdgeTerms;
            this.BoundaryEdgeTerms = _edgeForm.BoundaryEdgeTerms;
            this.edgeForm = _edgeForm;
        }


        
        IEdgeForm edgeForm;

        /// <summary>
        /// Returns name of nested form.
        /// </summary>
        public override string ToString() {
            return edgeForm.GetType().FullName;
        }


        /// <summary>
        /// ordering of argument variables
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return this.edgeForm.ArgumentOrdering;
            }
        }

        /// <summary>
        /// ordering of parameter variables
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return this.edgeForm.ParameterOrdering;
            }
        }


        
        // "global" variables
        double[] uA;
        double[] uB;
        double[,] Grad_uA;
        double[,] Grad_uB;
        double vA = 0.0;
        double vB = 0.0;
        double[] Grad_vA;
        double[] Grad_vB;
        int D;
        int K;
        int L;
        int NoArgs;
        int NoParams;


        private double GetCoeff(ref double d1, ref double d2, ref CommonParams inp) {
            int D = inp.Normale.Length;
            Debug.Assert(this.NoArgs == this.ArgumentOrdering.Count);

            Debug.Assert(uA.Length == NoArgs);
            Debug.Assert(uB.Length == NoArgs);
            Debug.Assert(Grad_uA.GetLength(0) == NoArgs);
            Debug.Assert(Grad_uB.GetLength(0) == NoArgs);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);
            Debug.Assert(Grad_vA.Length == D);
            Debug.Assert(Grad_vB .Length == D);

            Debug.Assert(this.edgeForm.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) == 0.0);

            d2 = 1; // set test function
            double a1 = this.edgeForm.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the source/affine part
            Debug.Assert((this.InnerEdgeTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);
            d1 = 1; // set trial function
            double a = this.edgeForm.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB); // probe for the bilinear+affine part
            d1 = 0; // reset
            d2 = 0; // reset
            return a - a1;
        }

        private double GetSourceCoeff(ref double d2, ref CommonParams inp) {
            Debug.Assert(this.edgeForm.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB) == 0.0);

            d2 = 1;
            double a = this.edgeForm.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
            Debug.Assert((this.InnerEdgeTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a == 0);

            d2 = 0;
            return a;
        }

        private double GetCoeffBnd(ref double d1, ref double d2, ref CommonParamsBnd inp) {
            Debug.Assert(this.edgeForm.BoundaryEdgeForm(ref inp, uA, Grad_uA, vA, Grad_vA) == 0.0);

            d2 = 1;
            double a1 = this.edgeForm.BoundaryEdgeForm(ref inp, uA, Grad_uA, vA, Grad_vA);
            Debug.Assert((this.BoundaryEdgeTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a1 == 0);

            d1 = 1;
            double a = this.edgeForm.BoundaryEdgeForm(ref inp, uA, Grad_uA, vA, Grad_vA);
            d1 = 0;
            d2 = 0;
            return a - a1;
        }

        private double GetSourceCoeffBnd(ref double d2, ref CommonParamsBnd inp) {
            Debug.Assert(this.edgeForm.BoundaryEdgeForm(ref inp, uA, Grad_uA, vA, Grad_vA) == 0.0);

            d2 = 1;
            double a = this.edgeForm.BoundaryEdgeForm(ref inp, uA, Grad_uA, vA, Grad_vA);
            Debug.Assert((this.BoundaryEdgeTerms & (TermActivationFlags.GradV | TermActivationFlags.V)) != 0 || a == 0);

            d2 = 0;
            return a;
        }

        /// <summary>
        /// see <see cref="IEdgeform_GradUxV.InternalEdge"/>
        /// </summary>
        void IEdgeform_GradUxV.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray GradUxV) {
            InitGlobals(efp);

            Debug.Assert(efp.ParameterVars_IN.Length == NoParams);
            Debug.Assert(efp.ParameterVars_OUT.Length == NoParams);


            CommonParams cp = default(CommonParams);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.Parameters_OUT = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            for (int l = 0; l < L; l++) { // loop over edges 
                cp.iEdge = efp.e0 + l;
                
                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cp.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        for (int c = 0; c < NoArgs; c++) {
                            GradUxV[l, k, 0, 0, c, d] = GetCoeff(ref Grad_uA[c, d], ref vA, ref cp);
                            GradUxV[l, k, 0, 1, c, d] = GetCoeff(ref Grad_uB[c, d], ref vA, ref cp);
                            GradUxV[l, k, 1, 0, c, d] = GetCoeff(ref Grad_uA[c, d], ref vB, ref cp);
                            GradUxV[l, k, 1, 1, c, d] = GetCoeff(ref Grad_uB[c, d], ref vB, ref cp);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeform_GradUxV.BoundaryEdge"/>
        /// </summary>
        void IEdgeform_GradUxV.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray GradUxV) {
            InitGlobals(efp);

            CommonParamsBnd cp = default(CommonParamsBnd);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;
                cp.EdgeTag = _EdgeTags[cp.iEdge];

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        for (int c = 0; c < NoArgs; c++) {
                            GradUxV[l, k, c, d] = GetCoeffBnd(ref Grad_uA[c, d], ref vA, ref cp);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeform_UxGradV.InternalEdge"/>
        /// </summary>
        void IEdgeform_UxGradV.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray UxGradV) {
            InitGlobals(efp);

            Debug.Assert(L == efp.Len);
            Debug.Assert(efp.ParameterVars_IN.Length == NoParams);
            Debug.Assert(efp.ParameterVars_OUT.Length == NoParams);
            Debug.Assert(efp.Len == UxGradV.GetLength(0));
            Debug.Assert(2 == UxGradV.GetLength(2));
            Debug.Assert(2 == UxGradV.GetLength(3));
            Debug.Assert(NoArgs == UxGradV.GetLength(4));
            Debug.Assert(D == UxGradV.GetLength(5));

            CommonParams cp = default(CommonParams);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.Parameters_OUT = new double[NoParams];
            cp.GridDat = efp.GridDat;

            for (int l = 0; l < L; l++) { // loop over edges 
                cp.iEdge = efp.e0 + l;
                
                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cp.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        for (int c = 0; c < NoArgs; c++) {
                            UxGradV[l, k, 0, 0, c, d] = GetCoeff(ref uA[c], ref Grad_vA[d], ref cp);
                            UxGradV[l, k, 0, 1, c, d] = GetCoeff(ref uB[c], ref Grad_vA[d], ref cp);
                            UxGradV[l, k, 1, 0, c, d] = GetCoeff(ref uA[c], ref Grad_vB[d], ref cp);
                            UxGradV[l, k, 1, 1, c, d] = GetCoeff(ref uB[c], ref Grad_vB[d], ref cp);
                        }
                    }
                }
            }
        }

        void IEdgeform_GradUxGradV.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray GradUxGradV) {
            InitGlobals(efp);

            Debug.Assert(L == efp.Len);
            Debug.Assert(efp.ParameterVars_IN.Length == NoParams);
            Debug.Assert(efp.ParameterVars_OUT.Length == NoParams);
            Debug.Assert(efp.Len == GradUxGradV.GetLength(0));
            Debug.Assert(2 == GradUxGradV.GetLength(2));
            Debug.Assert(2 == GradUxGradV.GetLength(3));
            Debug.Assert(NoArgs == GradUxGradV.GetLength(4));
            Debug.Assert(D == GradUxGradV.GetLength(6));
            Debug.Assert(D == GradUxGradV.GetLength(6));

            CommonParams cp = default(CommonParams);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.Parameters_OUT = new double[NoParams];
            cp.GridDat = efp.GridDat;

            for(int l = 0; l < L; l++) { // loop over edges 
                cp.iEdge = efp.e0 + l;

                for(int k = 0; k < K; k++) { // loop over quadrature nodes

                    for(int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for(int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cp.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            for(int c = 0; c < NoArgs; c++) {
                                GradUxGradV[l, k, 0, 0, c, d1, d2] = GetCoeff(ref Grad_uA[c, d1], ref Grad_vA[d2], ref cp);
                                GradUxGradV[l, k, 0, 1, c, d1, d2] = GetCoeff(ref Grad_uB[c, d1], ref Grad_vA[d2], ref cp);
                                GradUxGradV[l, k, 1, 0, c, d1, d2] = GetCoeff(ref Grad_uA[c, d1], ref Grad_vB[d2], ref cp);
                                GradUxGradV[l, k, 1, 1, c, d1, d2] = GetCoeff(ref Grad_uB[c, d1], ref Grad_vB[d2], ref cp);
                            }
                        }
                    }
                }
            }
        }

        private void InitGlobals(EdgeFormParams efp) {
            D = efp.GridDat.SpatialDimension;
            K = efp.Normals.GetLength(1);
            L = efp.Len;
            NoArgs = this.ArgumentOrdering.Count;
            NoParams = this.ParameterOrdering != null ? this.ParameterOrdering.Count : 0;

            uA = new double[NoArgs];
            uB = new double[NoArgs];
            Grad_uA = new double[NoArgs, D];
            Grad_uB = new double[NoArgs, D];
            Grad_vA = new double[D];
            Grad_vB = new double[D];
        }

        /// <summary>
        /// see <see cref="IEdgeform_UxGradV.BoundaryEdge"/>
        /// </summary>
        void IEdgeform_UxGradV.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray UxGradV) {
            InitGlobals(efp);

            CommonParamsBnd cp = default(CommonParamsBnd);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;
                cp.EdgeTag = _EdgeTags[cp.iEdge];

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        for (int c = 0; c < NoArgs; c++) {
                            UxGradV[l, k, c, d] = GetCoeffBnd(ref uA[c], ref Grad_vA[d], ref cp);
                        }
                    }
                }
            }
        }

        
        void IEdgeform_GradUxGradV.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray GradUxGradV) {
            InitGlobals(efp);

            CommonParamsBnd cp = default(CommonParamsBnd);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for(int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;
                cp.EdgeTag = _EdgeTags[cp.iEdge];

                for(int k = 0; k < K; k++) { // loop over quadrature nodes

                    for(int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for(int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            for(int c = 0; c < NoArgs; c++) {
                                GradUxGradV[l, k, c, d1, d2] = GetCoeffBnd(ref Grad_uA[c,d1], ref Grad_vA[d2], ref cp);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeSource_V.InternalEdge"/>
        /// </summary>
        void IEdgeSource_V.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray V) {
            InitGlobals(efp);

            CommonParams cp = default(CommonParams);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.Parameters_OUT = new double[NoParams];
            cp.GridDat = efp.GridDat;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cp.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    V[l, k, 0] = GetSourceCoeff(ref vA, ref cp);
                    V[l, k, 1] = GetSourceCoeff(ref vB, ref cp);
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeSource_V.BoundaryEdge"/>
        /// </summary>
        void IEdgeSource_V.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray V) {
            InitGlobals(efp);

            CommonParamsBnd cp = default(CommonParamsBnd);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;
                cp.EdgeTag = _EdgeTags[cp.iEdge];

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    V[l, k] = GetSourceCoeffBnd(ref vA, ref cp);
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeSource_GradV.InternalEdge"/>
        /// </summary>
        void IEdgeSource_GradV.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray GradV) {
            InitGlobals(efp);

            CommonParams cp = default(CommonParams);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.Parameters_OUT = new double[NoParams];
            cp.GridDat = efp.GridDat;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cp.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        GradV[l, k, 0, d] = GetSourceCoeff(ref Grad_vA[d], ref cp);
                        GradV[l, k, 1, d] = GetSourceCoeff(ref Grad_vB[d], ref cp);
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeSource_GradV.BoundaryEdge"/>
        /// </summary>
        void IEdgeSource_GradV.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray GradV) {
            InitGlobals(efp);

            CommonParamsBnd cp = default(CommonParamsBnd);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;
                cp.EdgeTag = _EdgeTags[cp.iEdge];

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        GradV[l, k, d] = GetSourceCoeffBnd(ref Grad_vA[d], ref cp);
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeform_UxV.InternalEdge"/>
        /// </summary>
        void IEdgeform_UxV.InternalEdge(ref EdgeFormParams efp, MultidimensionalArray UxV) {
            InitGlobals(efp);

            Debug.Assert(efp.ParameterVars_IN.Length == NoParams);
            Debug.Assert(efp.ParameterVars_OUT.Length == NoParams);

            CommonParams cp = default(CommonParams);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.Parameters_OUT = new double[NoParams];
            cp.GridDat = efp.GridDat;

            for (int l = 0; l < L; l++) { // loop over edges 
                cp.iEdge = efp.e0 + l;

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cp.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for (int c = 0; c < NoArgs; c++) {
                        UxV[l, k, 0, 0, c] = GetCoeff(ref uA[c], ref vA, ref cp);
                        UxV[l, k, 0, 1, c] = GetCoeff(ref uB[c], ref vA, ref cp);
                        UxV[l, k, 1, 0, c] = GetCoeff(ref uA[c], ref vB, ref cp);
                        UxV[l, k, 1, 1, c] = GetCoeff(ref uB[c], ref vB, ref cp);
                    }
                }
            }
        }

        /// <summary>
        /// see <see cref="IEdgeform_UxV.BoundaryEdge"/>
        /// </summary>
        void IEdgeform_UxV.BoundaryEdge(ref EdgeFormParams efp, MultidimensionalArray UxV) {
            InitGlobals(efp);

            CommonParamsBnd cp = default(CommonParamsBnd);
            cp.Normale = new double[D];
            cp.X = new double[D];
            cp.Parameters_IN = new double[NoParams];
            cp.GridDat = efp.GridDat;
            cp.time = efp.time;

            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges
                cp.iEdge = efp.e0 + l;
                cp.EdgeTag = _EdgeTags[cp.iEdge];

                for (int k = 0; k < K; k++) { // loop over quadrature nodes

                    for (int d = 0; d < D; d++) {
                        cp.Normale[d] = efp.Normals[l, k, d];
                        cp.X[d] = efp.NodesGlobal[l, k, d];
                    }

                    for (int np = 0; np < NoParams; np++) {
                        cp.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }


                    for (int c = 0; c < NoArgs; c++) {
                        UxV[l, k, c] = GetCoeffBnd(ref uA[c], ref vA, ref cp);
                    }
                }
            }
        }


       
        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return this.edgeForm.InnerEdgeForm(ref inp, _uA, _uB, _Grad_uA, _Grad_uB, _vA, _vB, _Grad_vA, _Grad_vB);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return this.edgeForm.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
        }


    }

}
