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
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Quadrature.NonLin {


    class NonlinVolumeFormVectorizer : INonlinVolumeForm_GradV, INonlinVolumeForm_V, IMultitreadSafety {

        public TermActivationFlags VolTerms {
            get;
            private set;
        }

        public NonlinVolumeFormVectorizer(IVolumeForm _volForm) {
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


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return this.volForm.VolumeForm(ref cpv, U, GradU, V, GradV);
        }

        void INonlinVolumeForm_V.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {
            int L = prm.Len;
            Debug.Assert(f.GetLength(0) == L);
            int K = f.GetLength(1); // no of nodes per cell
            int D = prm.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == prm.ParameterVars.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == U.Length);
            Debug.Assert(_NOargs == GradU.Length);

            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new Vector(D);
            cpv.time = prm.time;
            double[] _GradV = new double[D];
            double[,] _GradU = new double[_NOargs, D];
            double[] _U = new double[_NOargs];
            double _V = 1.0;

#if DEBUG
            MultidimensionalArray f_check = null;

            if (volForm is INonlinVolumeForm_V volForm_) {
                f_check = f.CloneAs();
                volForm_.Form(ref prm, U, GradU, f);
                var f_tmp = f;
                f = f_check;
                f_check = f_tmp;
            }
#endif

            for (int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for (int k = 0; k < K; k++) { // loop over nodes...


                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }
                    for (int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        if (U[na] != null) {
                            _U[na] = U[na][l, k];
                        }
                        if (GradU[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU[na, d] = GradU[na][l, k, d];
                            }
                        }
                    }

                    f[l, k] += volForm.VolumeForm(ref cpv, _U, _GradU, _V, _GradV);
                }
            }
#if DEBUG
            if (f_check != null) {
                double f_RelErr = f_check.L2Dist(f) / Math.Max(f.L2Norm(), 1);
                Debug.Assert(f_RelErr < 1e-14);
            }
#endif

        }

        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {
            int L = prm.Len;
            Debug.Assert(f.GetLength(0) == L);
            int K = f.GetLength(1); // no of nodes per cell
            int D = prm.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == prm.ParameterVars.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == U.Length);
            Debug.Assert(_NOargs == GradU.Length);

            CommonParamsVol cpv;
            cpv.GridDat = prm.GridDat;
            cpv.Parameters = new double[_NOParams];
            cpv.Xglobal = new Vector(D);
            cpv.time = prm.time;
            double[] _GradV = new double[D];
            double[,] _GradU = new double[_NOargs, D];
            double[] _U = new double[_NOargs];
            double _V = 0.0;
#if DEBUG
            MultidimensionalArray f_check = null;

            if (volForm is INonlinVolumeForm_GradV volForm_) {
                f_check = f.CloneAs();
                volForm_.Form(ref prm, U, GradU, f);
                var f_tmp = f;
                f = f_check;
                f_check = f_tmp;
            }
#endif
            //unsafe {
            //    fixed(double* pParameters = cpv.Parameters, p_U = _U, p_GradU = _GradU, p_GradV = _GradV) {
            //bool* pUisNull = stackalloc bool[_NOargs];
            //bool* pGradUisNull = stackalloc bool[_NOargs];

            //for(int na = 0; na < _NOargs; na++) {
            //    pUisNull[na] = (U[na] != null);
            //    pGradUisNull[na] = (GradU[na] != null);
            //}

            // unsafe opt bringt hier relativ wenig/nix

            for (int l = 0; l < L; l++) { // loop over cells...
                cpv.jCell = prm.j0 + l;

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int d = 0; d < D; d++) {
                        cpv.Xglobal[d] = prm.Xglobal[l, k, d];
                    }
                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters[np] = prm.ParameterVars[np][l, k];
                    }

                    /*
                    bool* ppUisNull = pUisNull;
                    bool* ppGradUisNull = pGradUisNull;
                    double* pp_U = p_U, pp_GradU = p_GradU;
                    for(int na = 0; na < _NOargs; na++) {
                        if(*ppUisNull) {
                            *pp_U = U[na][l, k];
                        } else {
                            *pp_U = double.NaN;
                        }
                        pp_U++;
                        if(*ppGradUisNull) {
                            for(int d = 0; d < D; d++) {
                                *pp_GradU = GradU[na][l, k, d];
                                pp_GradU++;
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
                                *pp_GradU = double.NaN;
                                pp_GradU++;
                            }
                        }
                        ppUisNull++;
                        ppGradUisNull++;
                    }
                    */
                    for (int na = 0; na < _NOargs; na++) {
                        if (U[na] != null) {
                            _U[na] = U[na][l, k];
                        } else {
                            _U[na] = double.NaN;
                        }
                        if (GradU[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU[na, d] = GradU[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU[na, d] = double.NaN;
                            }
                        }
                    }

                    /*
                    double* pp_GradV = p_GradV;
                    for(int d = 0; d < D; d++) {
                        *pp_GradV = 1.0;
                        f[l, k, d] += volForm.VolumeForm(ref cpv, _U, _GradU, _V, _GradV);
                        *pp_GradV = 0.0;
                        pp_GradV++;
                    }
                     */
                    for (int d = 0; d < D; d++) {
                        _GradV[d] = 1.0;
                        f[l, k, d] += volForm.VolumeForm(ref cpv, _U, _GradU, _V, _GradV);
                        _GradV[d] = 0.0;
                    }
                }
            }
#if DEBUG
            if (f_check != null) {
                double f_RelErr = f_check.L2Dist(f) / Math.Max(f.L2Norm(), 1);
                Debug.Assert(f_RelErr < 1e-14);
            }
#endif
        }

        bool m_IsMultithreadSafe;

        public bool IsMultithreadSafe {
            get {
                SetupMultithread();
                return m_IsMultithreadSafe;
            }
        }

        bool m_SetupMultithread_executed = false;

        void SetupMultithread() {
            if (m_SetupMultithread_executed == true)
                return;
            m_SetupMultithread_executed = true;

            if (this.volForm is IMultitreadSafety ms) {
                if (ms.IsMultithreadSafe) {
                    // nothing to do
                    m_IsMultithreadSafe = true;
                } else {
                    var clone = (IVolumeForm)ms.CloneForThread();
                    if (clone != null) {
                        this.volForm = clone;
                        m_IsMultithreadSafe = true;
                    } else {
                        m_IsMultithreadSafe = false;
                    }
                }

                m_SetupMultithread_executed = true;
            } else {
                // per default we just assume nothing goes wrong.
                m_IsMultithreadSafe = true;
            }



        }


        public IEquationComponent CloneForThread() {
            SetupMultithread();
            if (this.m_IsMultithreadSafe == true)
                return this; // it's not necessary to clone the vectorizer, since these guys are anyway created for each thread
            else
                return null;
        }

        public object GetPadlock() {
            SetupMultithread();
            if (this.volForm is IMultitreadSafety ms)
                return ms.GetPadlock();
            else
                return this.volForm;
        }

    }


    class NonlinEdgeFormVectorizer : INonlinEdgeForm_V, INonlinEdgeForm_GradV, IMultitreadSafety {

        public NonlinEdgeFormVectorizer(IEdgeForm __edgeForm) {
            this.BoundaryEdgeTerms = __edgeForm.BoundaryEdgeTerms;
            this.InnerEdgeTerms = __edgeForm.InnerEdgeTerms;
            this.edgeForm = __edgeForm;
        }


        IEdgeForm edgeForm;

        /// <summary>
        /// Returns name of nested form.
        /// </summary>
        public override string ToString() {
            return edgeForm.GetType().FullName;
        }


        public TermActivationFlags BoundaryEdgeTerms {
            get;
            private set;
        }

        public TermActivationFlags InnerEdgeTerms {
            get;
            private set;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return edgeForm.InnerEdgeForm(ref inp, _uA, _uB, _Grad_uA, _Grad_uB, _vA, _vB, _Grad_vA, _Grad_vB);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return edgeForm.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
        }

        public IList<string> ArgumentOrdering {
            get {
                return edgeForm.ArgumentOrdering;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return edgeForm.ParameterOrdering;
            }
        }

        void INonlinInnerEdgeForm_V.NonlinInternalEdge_V(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            int L = efp.Len;
            Debug.Assert(fin.GetLength(0) == L);
            Debug.Assert(fot.GetLength(0) == L);
            int K = fin.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            var E2C = efp.GridDat.iGeomEdges.CellIndices;
            Debug.Assert(_NOParams == efp.ParameterVars_IN.Length);
            Debug.Assert(_NOParams == efp.ParameterVars_OUT.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == Uin.Length);
            Debug.Assert(_NOargs == GradUin.Length);
            Debug.Assert(_NOargs == Uout.Length);
            Debug.Assert(_NOargs == GradUout.Length);
            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            CommonParams cpv;
            cpv.GridDat = efp.GridDat;
            cpv.Parameters_IN = new double[_NOParams];
            cpv.Parameters_OUT = new double[_NOParams];
            cpv.Normal = new Vector(D);
            cpv.X = new Vector(D);
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            double[] _GradV_ot = new double[D];
            double[,] _GradU_ot = new double[_NOargs, D];
            double[] _U_ot = new double[_NOargs];
            double _V_ot = 0.0;

#if DEBUG
            MultidimensionalArray fin_check = null;
            MultidimensionalArray fot_check = null;
            if (edgeForm is INonlinEdgeForm_V edgeForm_) {
                fin_check = fin.CloneAs();
                fot_check = fot.CloneAs();

                edgeForm_.NonlinInternalEdge_V(ref efp, Uin, Uout, GradUin, GradUout, fin, fot);

                var fin_tmp = fin;
                var fot_tmp = fot;
                fin = fin_check;
                fot = fot_check;
                fin_check = fin_tmp;
                fot_check = fot_tmp;
            }
#endif

            for (int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = _EdgeTags[cpv.iEdge];
                cpv.jCellIn = E2C[cpv.iEdge, 0];
                cpv.jCellOut = E2C[cpv.iEdge, 1];

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cpv.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        Debug.Assert((Uin[na] != null) == (Uout[na] != null));
                        if (Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                            _U_ot[na] = Uout[na][l, k];
                        } else {
                            _U_in[na] = 0;
                            _U_ot[na] = 0;
                        }
                        Debug.Assert((GradUin[na] != null) == (GradUout[na] != null));
                        if (GradUin[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                                _GradU_ot[na, d] = GradUout[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                                _GradU_ot[na, d] = 0;
                            }
                        }
                    }

                    _V_ot = 0;
                    _V_in = 1;
                    fin[l, k] += edgeForm.InnerEdgeForm(ref cpv, _U_in, _U_ot, _GradU_in, _GradU_ot, _V_in, _V_ot, _GradV_in, _GradV_ot);
                    _V_in = 0;
                    _V_ot = 1;
                    fot[l, k] += edgeForm.InnerEdgeForm(ref cpv, _U_in, _U_ot, _GradU_in, _GradU_ot, _V_in, _V_ot, _GradV_in, _GradV_ot);
                }
            }

#if DEBUG
            if (fin_check != null) {
                double fin_RelErr = fin_check.L2Dist(fin) / Math.Max(Math.Max(fin.L2Norm(), fin_check.L2Norm()), 1);
                double fot_RelErr = fot_check.L2Dist(fot) / Math.Max(Math.Max(fot.L2Norm(), fot_check.L2Norm()), 1);
                Debug.Assert(fin_RelErr < 1e-14 && fot_RelErr < 1e-14);

            }
#endif
        }

        void INonlinBoundaryEdgeForm_V.NonlinBoundaryEdge_V(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin, MultidimensionalArray f) {
            int L = efp.Len;
            Debug.Assert(f.GetLength(0) == L);
            int K = f.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == efp.ParameterVars_IN.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == Uin.Length);
            Debug.Assert(_NOargs == GradUin.Length);

            CommonParamsBnd cpv;
            cpv.GridDat = efp.GridDat;
            cpv.Parameters_IN = new double[_NOParams];
            cpv.Normal = new Vector(D);
            cpv.X = new Vector(D);
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;


#if DEBUG
            MultidimensionalArray f_check = null;

            if (edgeForm is INonlinEdgeForm_V edgeForm_) {
                f_check = f.CloneAs();
                edgeForm_.NonlinBoundaryEdge_V(ref efp, Uin, GradUin, f);
                var f_tmp = f;
                f = f_check;
                f_check = f_tmp;
            }
#endif

            for (int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        if (Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if (GradUin[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    _V_in = 1;
                    f[l, k] += edgeForm.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                }
            }

#if DEBUG
            if (f_check != null) { 
                double f_RelErr = f_check.L2Dist(f) / Math.Max(f.L2Norm(), 1); 
                Debug.Assert(f_RelErr < 1e-14);
            }
#endif
        }


        void INonlinInnerEdgeForm_GradV.NonlinInternalEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            int L = efp.Len;
            Debug.Assert(fin.GetLength(0) == L);
            Debug.Assert(fot.GetLength(0) == L);
            int K = fin.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;
            var E2C = efp.GridDat.iGeomEdges.CellIndices;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == efp.ParameterVars_IN.Length);
            Debug.Assert(_NOParams == efp.ParameterVars_OUT.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == Uin.Length);
            Debug.Assert(_NOargs == GradUin.Length);
            Debug.Assert(_NOargs == Uout.Length);
            Debug.Assert(_NOargs == GradUout.Length);
            var _EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            CommonParams cpv;
            cpv.GridDat = efp.GridDat;
            cpv.Parameters_IN = new double[_NOParams];
            cpv.Parameters_OUT = new double[_NOParams];
            cpv.Normal = new Vector(D);
            cpv.X = new Vector(D);
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            double[] _GradV_ot = new double[D];
            double[,] _GradU_ot = new double[_NOargs, D];
            double[] _U_ot = new double[_NOargs];
            double _V_ot = 0.0;

#if DEBUG
            MultidimensionalArray fin_check = null;
            MultidimensionalArray fot_check = null;
            if (edgeForm is INonlinEdgeForm_GradV edgeForm_) {
                fin_check = fin.CloneAs();
                fot_check = fot.CloneAs();

                edgeForm_.NonlinInternalEdge_GradV(ref efp, Uin, Uout, GradUin, GradUout, fin, fot);

                var fin_tmp = fin;
                var fot_tmp = fot;
                fin = fin_check;
                fot = fot_check;
                fin_check = fin_tmp;
                fot_check = fot_tmp;
            }
#endif

            for (int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = _EdgeTags[cpv.iEdge];
                cpv.jCellIn = E2C[cpv.iEdge, 0];
                cpv.jCellOut = E2C[cpv.iEdge, 1];

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        cpv.Parameters_OUT[np] = efp.ParameterVars_OUT[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        Debug.Assert((Uin[na] != null) == (Uout[na] != null));
                        if (Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                            _U_ot[na] = Uout[na][l, k];
                        } else {
                            _U_in[na] = 0;
                            _U_ot[na] = 0;
                        }
                        Debug.Assert((GradUin[na] != null) == (GradUout[na] != null));
                        if (GradUin[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                                _GradU_ot[na, d] = GradUout[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                                _GradU_ot[na, d] = 0;
                            }
                        }
                    }

                    for (int d = 0; d < D; d++) {
                        _GradV_in[d] = 1;
                        fin[l, k, d] += edgeForm.InnerEdgeForm(ref cpv, _U_in, _U_ot, _GradU_in, _GradU_ot, _V_in, _V_ot, _GradV_in, _GradV_ot);
                        _GradV_in[d] = 0;
                        _GradV_ot[d] = 1;
                        fot[l, k, d] += edgeForm.InnerEdgeForm(ref cpv, _U_in, _U_ot, _GradU_in, _GradU_ot, _V_in, _V_ot, _GradV_in, _GradV_ot);
                        _GradV_ot[d] = 0;
                    }
                }
            }

#if DEBUG
            if (fin_check != null) {
                double fin_RelErr = fin_check.L2Dist(fin) / Math.Max(Math.Max(fin.L2Norm(), fin_check.L2Norm()), 1);
                double fot_RelErr = fot_check.L2Dist(fot) / Math.Max(Math.Max(fot.L2Norm(), fot_check.L2Norm()), 1);               
                Debug.Assert(fin_RelErr < 1e-14 && fot_RelErr < 1e-14);

            }
#endif
        }

        void INonlinBoundaryEdgeForm_GradV.NonlinBoundaryEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin, MultidimensionalArray f) {
            int L = efp.Len;
            Debug.Assert(f.GetLength(0) == L);
            int K = f.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == efp.ParameterVars_IN.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == Uin.Length);
            Debug.Assert(_NOargs == GradUin.Length);

            CommonParamsBnd cpv;
            cpv.GridDat = efp.GridDat;
            cpv.Parameters_IN = new double[_NOParams];
            cpv.Normal = new Vector(D);
            cpv.X = new Vector(D);
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

#if DEBUG
            MultidimensionalArray f_check = null;

            if (edgeForm is INonlinEdgeForm_GradV edgeForm_) {
                f_check = f.CloneAs();
                edgeForm_.NonlinBoundaryEdge_GradV(ref efp, Uin, GradUin, f);
                var f_tmp = f;
                f = f_check;
                f_check = f_tmp;
            }
#endif

            for (int l = 0; l < L; l++) { // loop over edges ...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for (int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        if (Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if (GradUin[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    for (int d = 0; d < D; d++) {
                        _GradV_in[d] = 1;
                        f[l, k, d] += edgeForm.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                        _GradV_in[d] = 0;
                    }
                }
            }

#if DEBUG
            if (f_check != null) {
                double f_RelErr = f_check.L2Dist(f) / Math.Max(f.L2Norm(), 1);              
                Debug.Assert(f_RelErr < 1e-14);
            }
#endif
        }

        bool m_IsMultithreadSafe;

        public bool IsMultithreadSafe {
            get {
                SetupMultithread();
                return m_IsMultithreadSafe;
            }
        }

        bool m_SetupMultithread_executed = false;

        void SetupMultithread() {
            if (m_SetupMultithread_executed == true)
                return;
            m_SetupMultithread_executed = true;

            if (this.edgeForm is IMultitreadSafety ms) {
                if (ms.IsMultithreadSafe) {
                    // nothing to do
                    m_IsMultithreadSafe = true;
                } else {
                    var clone = (IEdgeForm)ms.CloneForThread();
                    if (clone != null) {
                        this.edgeForm = clone;
                        m_IsMultithreadSafe = true;
                    } else {
                        m_IsMultithreadSafe = false;
                    }
                }

                m_SetupMultithread_executed = true;
            } else {
                // per default we just assume nothing goes wrong.
                m_IsMultithreadSafe = true;
            }



        }


        public IEquationComponent CloneForThread() {
            SetupMultithread();
            if (this.m_IsMultithreadSafe == true)
                return this; // it's not necessary to clone the vectorizer, since these guys are anyway created for each thread
            else
                return null;
        }

        public object GetPadlock() {
            SetupMultithread();
            if (this.edgeForm is IMultitreadSafety ms)
                return ms.GetPadlock();
            else
                return this.edgeForm;
        }
    }
}
