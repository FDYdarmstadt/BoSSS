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

using System.Diagnostics;
using BoSSS.Platform;
using System;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;

namespace BoSSS.Foundation {

    /// <summary>
    /// This class encapsulates the evaluation of DG fields (value, gradient, mean-value,...)
    /// in a way that it can be commonly used by the single-phase and the XDG framework.
    /// The main purpose of this hack is to eliminate duplicate code and make performance optimizations available to 
    /// all kinds of DG fields.
    /// The methods are not intended for direct use outside of layer L2.
    /// </summary>
    public partial class DGField {

        /// <summary>
        /// evaluation of DG field; may be used in derived classes to implement <see cref="DGField.Evaluate(int,int,NodeSet,MultidimensionalArray,int,double)"/>.
        /// </summary>
        protected static void EvaluateInternal(int j0, int L, NodeSet NS, Basis basis, MultidimensionalArray Coördinates, int coördOffset, MultidimensionalArray ResultAcc, double ResultPreScale) {

            int D, N, NumNodes;
            bool AffineLinear;
            CheckArgs(j0, L, NS, basis, Coördinates, ResultAcc, out D, out N, out NumNodes, out AffineLinear);
            Debug.Assert(ResultAcc.Dimension == 2);
            Debug.Assert(L == ResultAcc.GetLength(0));
            Debug.Assert(NumNodes == ResultAcc.GetLength(1));

            /*
            MultidimensionalArray BasisValues;
            BasisValues = basis.CellEval(NS, j0, L);
            Debug.Assert(BasisValues.GetLength(0) == L, "No. of. cells mismatch");
            Debug.Assert(BasisValues.GetLength(1) == M, "No. of. nodes mismatch");
            Debug.Assert(BasisValues.GetLength(2) == N, "No. of. basis elements mismatch");

            ResultAcc.Multiply(1.0, Coördinates, BasisValues, ResultPreScale, "jm", "jn", "jmn");
            */

            //int[] geom2log = basis.GridDat.iGeomCells.GeomCell2LogicalCell;

            MultidimensionalArray BasisValues = basis.Evaluate(NS);
            
            if(L == 1 && AffineLinear) {
                // Special optimization for single-cell evaluation:
                // this happens very often for edge quadrature, so it is quite relevant.
                double scale0 = basis.GridDat.ChefBasis.Scaling[j0];
                ResultAcc.Multiply(scale0, Coördinates, BasisValues, ResultPreScale, ref mp_jk_jm_km); //"jk", "jm", "km");
            } else {
                int iBuf;
                MultidimensionalArray trfCoördinates = TempBuffer.GetTempMultidimensionalarray(out iBuf, L, N);
                TransformCoördinates(j0, L, basis, Coördinates, coördOffset, N, AffineLinear, trfCoördinates);

                if(ResultAcc.IsContinious && trfCoördinates.IsContinious && BasisValues.IsContinious) {
                    unsafe {
                        fixed(double* _pResultAcc = ResultAcc.Storage, _ptrfCoördinates = trfCoördinates.Storage, _pBasisValues = BasisValues.Storage) {
                            double* pResultAcc = _pResultAcc + ResultAcc.Index(0,0);
                            double* ptrfCoördinates = _ptrfCoördinates + trfCoördinates.Index(0,0);
                            double* pBasisValues = _pBasisValues + BasisValues.Index(0, 0);

                            
//#if DEBUG
//                            MultidimensionalArray check = ResultAcc.CloneAs();
                            
//#endif

                            int _M = ResultAcc.GetLength(1);   // entspricht k   (node    index)
                            int _N = ResultAcc.GetLength(0);   // entspricht j   (cell    index)
                            int _K = BasisValues.GetLength(1); // entspricht m   (DG mode index)

                            // NOTE: dimensions in FORTRAN order:
                            // pBasisValues     :  _K x _M
                            // ptrfCoördinates  :  _K x _N
                            // pResultAcc       :  _M x _N 
                            //
                            // => FORTRAN GEMM
                            // pResultAcc = pBasisValues^T * ptrfCoördinates

                            int TRANSA = 'T';
                            int TRANSB = 'N';

                            BLAS.dgemm(TRANSA, TRANSB, _M, _N, _K, 
                                1.0, 
                                pBasisValues, _K,
                                ptrfCoördinates, _K,
                                ResultPreScale,
                                pResultAcc, _M);
                        
//#if DEBUG
//                            check.Multiply(1.0, trfCoördinates, BasisValues, ResultPreScale, ref mp_jk_jm_km);
//                            check.Acc(-1.0, ResultAcc);
//                            double error = check.L2Norm();
//                            Console.WriteLine("GEMM error: " + error);
//                            Debug.Assert(error < 1.0);
//#endif
                        }

                    }
                } else {
                    ResultAcc.Multiply(1.0, trfCoördinates, BasisValues, ResultPreScale, ref mp_jk_jm_km); //"jk", "jm", "km");
                }
                
                
                TempBuffer.FreeTempBuffer(iBuf);
            }
        }

        static MultidimensionalArray.MultiplyProgram mp_jk_jm_km = MultidimensionalArray.MultiplyProgram.Compile("jk", "jm", "km");

        private static void TransformCoördinates(int j0, int L, 
            Basis basis, 
            MultidimensionalArray Coördinates, int coördOffset, 
            int N, bool AffineLinear, 
            MultidimensionalArray trfCoördinates) {
            int[] geom2log = basis.GridDat.iGeomCells.GeomCell2LogicalCell;

            if (geom2log != null) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // aggregation grid branch -- apply index trafo to coordinates
                // (j0..j0+L) are geometical grid indices
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                // extract trafo
                MultidimensionalArray trafo = basis.GridDat.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(j0, L, basis.Degree);
                Debug.Assert(trafo.GetLength(0) == L);
                if (trafo.GetLength(1) > N)
                    trafo = trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { L - 1, N - 1, N - 1 });

                // apply trafo
                unsafe {
                    fixed (int* p_geom2log = geom2log) {
                        trfCoördinates.Multiply(1.0, trafo, Coördinates, 0.0, ref mp_jm_jmn_Tjn,
                            p_geom2log, p_geom2log,
                            trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0,
                            trfPreOffset_B: coördOffset, trfCycle_B: 1, trfPostOffset_B: 0); // geom cell to logical cell trafo for coördinates

                    }
                }
            } else if (AffineLinear) {
                // +++++++++++++++++++++++++
                // affine-linear grid branch
                // +++++++++++++++++++++++++

                // extract coördinates
                MultidimensionalArray _Coördinates;
                if(coördOffset > 0 || coördOffset + L < Coördinates.GetLength(0)) {
                    _Coördinates = Coördinates.ExtractSubArrayShallow(new[] { j0, 0 }, new[] { j0 + L - 1, N - 1 });
                } else {
                    _Coördinates = Coördinates;
                }
                
                // extract trafo
                MultidimensionalArray scale = basis.GridDat.ChefBasis.Scaling.ExtractSubArrayShallow(new int[] { j0 }, new int[] { j0 + L - 1 });

                // apply trafo
                trfCoördinates.Multiply(1.0, scale, _Coördinates, 0.0, ref mp_jn_j_jn);

            } else {
                // ++++++++++++++++++
                // curved cell branch
                // ++++++++++++++++++

                // extract coördinates
                MultidimensionalArray _Coördinates;
                if(coördOffset > 0 || coördOffset + L < Coördinates.GetLength(0)) {
                    _Coördinates = Coördinates.ExtractSubArrayShallow(new[] { j0, 0 }, new[] { j0 + L - 1, N - 1 });
                } else {
                    _Coördinates = Coördinates;
                }

                // extract trafo
                MultidimensionalArray trafo = basis.GridDat.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(j0, L, basis.Degree);
                Debug.Assert(trafo.GetLength(0) == L);
                if (trafo.GetLength(1) > N)
                    trafo = trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { L - 1, N - 1, N - 1 });

                // apply trafo
                trfCoördinates.Multiply(1.0, trafo, _Coördinates, 0.0, ref mp_jm_jmn_jn);
            }
        }

        static MultidimensionalArray.MultiplyProgram mp_jn_j_jn = MultidimensionalArray.MultiplyProgram.Compile("jn", "j", "jn");
        static MultidimensionalArray.MultiplyProgram mp_jm_jmn_jn = MultidimensionalArray.MultiplyProgram.Compile("jm", "jmn", "jn");
        static MultidimensionalArray.MultiplyProgram mp_jm_jmn_Tjn = MultidimensionalArray.MultiplyProgram.Compile("jm", "jmn", "T(j)n", true);

        private static void CheckArgs(int j0, int L, NodeSet NS, Basis basis, MultidimensionalArray Coördinates, MultidimensionalArray ResultAcc, out int D, out int N, out int M, out bool AffineLinear) {
            int[] g2l = basis.GridDat.iGeomCells.GeomCell2LogicalCell;

            D = basis.GridDat.SpatialDimension; // spatial dimension
            if(g2l == null)
                N = basis.GetLength(j0);      // number of coordinates per cell -- standard grid
            else 
                N = basis.GetLength(g2l[j0]); // number of coordinates per cell -- aggregation grid
            M = NS.NoOfNodes;            // number of nodes
            AffineLinear = basis.GridDat.iGeomCells.IsCellAffineLinear(j0);

            Debug.Assert(basis.GetType() == typeof(Basis));

            Debug.Assert(Coördinates.Dimension == 2);
            Debug.Assert(Coördinates.GetLength(0) >= j0 + L || basis.GridDat.iGeomCells.GeomCell2LogicalCell != null);
            Debug.Assert(Coördinates.GetLength(1) == N);

#if DEBUG
            for(int i = 1; i < L; i++) {
                int jCell = j0 + i;
                Debug.Assert(basis.GridDat.iGeomCells.IsCellAffineLinear(jCell) == AffineLinear);

                if(g2l == null)
                    Debug.Assert(basis.GetLength(jCell) == N);
                else 
                    Debug.Assert(basis.GetLength(g2l[jCell]) == N);
            }
#endif
        }

        /// <summary>
        /// evaluation of DG field gradient; may be used in derived classes to implement <see cref="EvaluateGradient"/>.
        /// </summary>
        protected static void EvaluateGradientInternal(int j0, int L, NodeSet NS, Basis basis, MultidimensionalArray Coördinates, int coördOffset, MultidimensionalArray ResultAcc, double ResultPreScale) {


            int D, N, K;
            bool AffineLinear;
            CheckArgs(j0, L, NS, basis, Coördinates, ResultAcc, out D, out N, out K, out AffineLinear);
            Debug.Assert(ResultAcc.Dimension == 3);
            Debug.Assert(L == ResultAcc.GetLength(0));
            Debug.Assert(K == ResultAcc.GetLength(1));
            Debug.Assert(D == ResultAcc.GetLength(2));

            /*
            MultidimensionalArray BasisGradValues;
            BasisGradValues = basis.CellEvalGradient(NS, j0, L);

            ResultAcc.Multiply(1.0, Coördinates, BasisGradValues, ResultPreScale, "jmd", "jn", "jmnd");
            */

            MultidimensionalArray BasisGradValues = basis.EvaluateGradient(NS);
            int iBuf1, iBuf2;
            if(AffineLinear) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // affine-linear-cell:
                // Inverse Jacobian different for each cell, but constant among nodes
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                if(L == 1) {
                    // Special optimization for single-cell evaluation:
                    // this happens very often for edge quadrature, so it is quite relevant.
                    double scale0 = basis.GridDat.ChefBasis.Scaling[j0];
                    MultidimensionalArray GradientRef = TempBuffer.GetTempMultidimensionalarray(out iBuf2, L, K, D);
                    MultidimensionalArray InvJacobi = basis.GridDat.iGeomCells.InverseTransformation.ExtractSubArrayShallow(j0, -1, -1);

                    GradientRef.Multiply(scale0, Coördinates, BasisGradValues, 0.0, ref mp_jke_jm_kme);  // gradient in reference coördinates
                    ResultAcc.Multiply(1.0, InvJacobi, GradientRef, ResultPreScale, ref mp_jkd_ed_jke);
                    
                    TempBuffer.FreeTempBuffer(iBuf2);
                } else {
                    MultidimensionalArray trfCoördinates = TempBuffer.GetTempMultidimensionalarray(out iBuf1, L, N);
                    MultidimensionalArray GradientRef = TempBuffer.GetTempMultidimensionalarray(out iBuf2, L, K, D);
                    MultidimensionalArray InvJacobi = basis.GridDat.iGeomCells.InverseTransformation.ExtractSubArrayShallow(new int[] { j0, 0, 0 }, new int[] { j0 + L - 1, D - 1, D - 1 });

                    TransformCoördinates(j0, L, basis, Coördinates, coördOffset, N, AffineLinear, trfCoördinates);
                    GradientRef.Multiply(1.0, trfCoördinates, BasisGradValues, 0.0, ref mp_jke_jm_kme);  // gradient in reference coördinates
                    ResultAcc.Multiply(1.0, InvJacobi, GradientRef, ResultPreScale, ref mp_jkd_jed_jke);
                    
                    TempBuffer.FreeTempBuffer(iBuf1);
                    TempBuffer.FreeTempBuffer(iBuf2);
                }
            } else {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // curved-cell:
                // Inverse Jacobian different for each node and each cell
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                MultidimensionalArray trfCoördinates = TempBuffer.GetTempMultidimensionalarray(out iBuf1, L, N);
                MultidimensionalArray GradientRef = TempBuffer.GetTempMultidimensionalarray(out iBuf2, L, K, D);
                MultidimensionalArray InvJacobi = basis.GridDat.InverseJacobian.GetValue_Cell(NS, j0, L);
                
                TransformCoördinates(j0, L, basis, Coördinates, coördOffset, N, AffineLinear, trfCoördinates);
                GradientRef.Multiply(1.0, trfCoördinates, BasisGradValues, 0.0, ref mp_jke_jm_kme);  // gradient in reference coördinates
                ResultAcc.Multiply(1.0, InvJacobi, GradientRef, ResultPreScale, ref mp_jkd_jked_jke);

                TempBuffer.FreeTempBuffer(iBuf1);
                TempBuffer.FreeTempBuffer(iBuf2);
            }
        }

        static MultidimensionalArray.MultiplyProgram mp_jke_jm_kme = MultidimensionalArray.MultiplyProgram.Compile("jke", "jm", "kme");
        static MultidimensionalArray.MultiplyProgram mp_jkd_jked_jke = MultidimensionalArray.MultiplyProgram.Compile("jkd", "jked", "jke");
        static MultidimensionalArray.MultiplyProgram mp_jkd_ed_jke = MultidimensionalArray.MultiplyProgram.Compile("jkd", "ed", "jke");
        static MultidimensionalArray.MultiplyProgram mp_jkd_jed_jke = MultidimensionalArray.MultiplyProgram.Compile("jkd", "jed", "jke");


        protected static void EvaluateEdgeInternal(int e0, int Len, NodeSet NS, Basis _Basis, MultidimensionalArray Coord,
            MultidimensionalArray valIN, MultidimensionalArray valOT,
            MultidimensionalArray meanValIN, MultidimensionalArray meanValOT,
            MultidimensionalArray gradIN, MultidimensionalArray gradOT,
            double ResultPreScale) {

            // checks and init
            // ===============

            var grd = _Basis.GridDat;
            int NoOfNodes = NS.NoOfNodes;
            bool AffineLinear = grd.iGeomEdges.IsEdgeAffineLinear(e0);
            Debug.Assert(NS.GetNodeCoordinateSystem(grd) == NodeCoordinateSystem.EdgeCoord);
            Debug.Assert(valIN == null || valIN.Dimension == 2);
            Debug.Assert(valOT == null || valOT.Dimension == 2);
            Debug.Assert(valIN == null || valIN.GetLength(0) == Len);
            Debug.Assert(valOT == null || valOT.GetLength(0) == Len);
            Debug.Assert(valIN == null || valIN.GetLength(1) == NoOfNodes);
            Debug.Assert(valOT == null || valOT.GetLength(1) == NoOfNodes);

            int[,] trfIdx = grd.iGeomEdges.Edge2CellTrafoIndex;
            int[,] E2Cl = grd.iGeomEdges.LogicalCellIndices;
            int[,] E2Cg = grd.iGeomEdges.CellIndices;

            

            // transform DG coördinates
            // ========================

            int Nin = _Basis.GetLength(E2Cl[e0, 0]);
            int Not = Nin;

#if DEBUG
            for (int e = 0; e < Len; e++) {
                int iEdge = e + e0;
                int jCellIN = E2Cl[iEdge, 0];
                int jCellOT = E2Cl[iEdge, 1];
                Debug.Assert(_Basis.GetLength(jCellIN) == Nin);
                if (jCellOT >= 0)
                    Debug.Assert(_Basis.GetLength(jCellOT) == Not);
            }
#endif
            int iBufIN, iBufOT = 0;
            MultidimensionalArray trfCoördinatesIN = TempBuffer.GetTempMultidimensionalarray(out iBufIN, Len, Nin);
            MultidimensionalArray trfCoördinatesOT = Not > 0 ? TempBuffer.GetTempMultidimensionalarray(out iBufOT, Len, Not) : null;
            TransformCoördinatesEdge(e0, Len, grd, Coord, Nin, Not, _Basis.Degree, AffineLinear, trfCoördinatesIN, trfCoördinatesOT);

            // Evaluate
            // ========

            unsafe {
                fixed(int* pTrfIndex = trfIdx) {

                    MultidimensionalArray BasisValues = null;
                    Debug.Assert((valIN != null) == (valOT != null));
                    if(valIN != null) {
                        // compute the values
                        // -------------------

                        BasisValues = grd.ChefBasis.EdgeEval.GetValues(NS, e0, Len, _Basis.Degree);
                        Debug.Assert(BasisValues.Dimension == 3);

                        MultidimensionalArray BasisValuesIN, BasisValuesOT;
                        if(BasisValues.GetLength(2) > Nin)
                            BasisValuesIN = BasisValues.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { BasisValues.GetLength(0) - 1, NoOfNodes - 1, Nin - 1 });
                        else
                            BasisValuesIN = BasisValues;
                        if(BasisValues.GetLength(2) > Not) {
                            if(Nin == Not)
                                BasisValuesOT = BasisValuesIN;
                            else
                                BasisValuesOT = BasisValues.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { BasisValues.GetLength(0) - 1, NoOfNodes - 1, Nin - 1 });
                        } else {
                            BasisValuesOT = BasisValues;
                        }

                        {
                            valIN.Multiply(1.0, trfCoördinatesIN, BasisValuesIN, ResultPreScale, ref mp_ik_im_Tikm,
                                pTrfIndex, pTrfIndex,
                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0), trfCycle_B: 2, trfPostOffset_B: 0);
                        }
                        if(Not > 0) {
                            valOT.Multiply(1.0, trfCoördinatesOT, BasisValuesOT, ResultPreScale, ref mp_ik_im_Tikm,
                                pTrfIndex, pTrfIndex,
                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);
                        }
                    }

                    Debug.Assert((meanValIN != null) == (meanValOT != null));
                    if(meanValIN != null) {
                        // compute the mean values
                        // -----------------------

                        var _Basis0Values = BasisValues;
                        if(_Basis0Values == null)
                            _Basis0Values = grd.ChefBasis.EdgeEval.GetValues(NS, e0, Len, 0);
                        // assume the 0-th basis polynomial \f$ \phi_0 \f$ is constant!
                        _Basis0Values = _Basis0Values.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { BasisValues.GetLength(0) - 1, -1, -1 });

                        // DG coördinates for the 0-th mode:
                        var _trfCoördinatesIN = trfCoördinatesIN.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Len-1, -1 });
                        var _trfCoördinatesOT = trfCoördinatesOT.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Len-1, -1 });

                        {
                            meanValIN.Multiply(1.0, _trfCoördinatesIN, _Basis0Values, ResultPreScale, ref mp_i_i_Ti,
                                pTrfIndex, pTrfIndex,
                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0), trfCycle_B: 2, trfPostOffset_B: 0);
                        }
                        if(Not > 0) {
                            meanValOT.Multiply(1.0, _trfCoördinatesOT, _Basis0Values, ResultPreScale, ref mp_i_i_Ti,
                                pTrfIndex, pTrfIndex,
                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);
                        }
                    }

                    Debug.Assert((gradIN != null) == (gradOT != null));
                    if(gradIN != null) {
                        // compute gradient values
                        // -----------------------

                        int D = grd.SpatialDimension;
                        int iBuf2;
                        MultidimensionalArray GradientRef = TempBuffer.GetTempMultidimensionalarray(out iBuf2, Len, NoOfNodes, D);

                        MultidimensionalArray BasisGradValues = grd.ChefBasis.EdgeGradientEval.GetValues(NS, e0, Len, _Basis.Degree);
                        Debug.Assert(BasisGradValues.Dimension == 4);

                        MultidimensionalArray BasisGradValuesIN, BasisGradValuesOT;
                        if(BasisGradValues.GetLength(2) > Nin)
                            BasisGradValuesIN = BasisGradValues.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { BasisGradValues.GetLength(0) - 1, NoOfNodes - 1, Nin - 1, D - 1 });
                        else
                            BasisGradValuesIN = BasisGradValues;
                        if(BasisGradValues.GetLength(2) > Not) {
                            if(Nin == Not)
                                BasisGradValuesOT = BasisGradValuesIN;
                            else
                                BasisGradValuesOT = BasisGradValues.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { BasisGradValues.GetLength(0) - 1, NoOfNodes - 1, Nin - 1, D - 1 });
                        } else {
                            BasisGradValuesOT = BasisGradValues;
                        }


                        if(AffineLinear) {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // affine-linear-cell:
                            // Inverse Jacobian different for each cell, but constant among nodes
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            var InvJacobi = grd.iGeomCells.InverseTransformation;

                            Debug.Assert(grd is Grid.Classic.GridData, "implementation only valid for classic grid");
                            Debug.Assert(object.ReferenceEquals(E2Cg, E2Cl));

                            fixed(int* pE2Cl = E2Cl) {

                                {
                                    GradientRef.Multiply(1.0, trfCoördinatesIN, BasisGradValuesIN, 0.0, ref mp_ike_im_Tikme,
                                        pTrfIndex, pTrfIndex,
                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0), trfCycle_B: 2, trfPostOffset_B: 0);  // gradient in reference coördinates

                                    gradIN.Multiply(1.0, InvJacobi, GradientRef, ResultPreScale, ref mp_ikd_Tied_ike,
                                        pE2Cl, pE2Cl, 
                                        trfPreOffset_A: (2 * e0), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                }

                                if(Not > 0) {
                                    GradientRef.Multiply(1.0, trfCoördinatesOT, BasisGradValuesOT, 0.0, ref mp_ike_im_Tikme,
                                        pTrfIndex, pTrfIndex,
                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);  // gradient in reference coördinates

                                    gradOT.Multiply(1.0, InvJacobi, GradientRef, ResultPreScale, ref mp_ikd_Tied_ike,
                                        pE2Cl, pE2Cl,
                                        trfPreOffset_A: (2 * e0 + 1), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                }
                            }


                        } else {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // curved-cell:
                            // Inverse Jacobian different for each node and each cell
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            MultidimensionalArray invJacobiIN, invJacobiOT;
                            var TiJ = grd.InverseJacobian.GetValue_EdgeDV(NS, e0, Len);
                            invJacobiIN = TiJ.Item1;
                            invJacobiOT = TiJ.Item2;


                            {
                                GradientRef.Multiply(1.0, trfCoördinatesIN, BasisGradValuesIN, 0.0, ref mp_ike_im_Tikme,
                                    pTrfIndex, pTrfIndex,
                                    trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0), trfCycle_B: 2, trfPostOffset_B: 0);  // gradient in reference coördinates, i.e. \f$ \nabla_{\vec{xi}} \f$

                                gradIN.Multiply(1.0, invJacobiIN, GradientRef, ResultPreScale, ref mp_ikd_iked_ike);   // gradient in physical coördinates, i.e. \f$ \nabla_{\vec{x}}n \f$
                            }

                            if(Not > 0) {
                                GradientRef.Multiply(1.0, trfCoördinatesOT, BasisGradValuesOT, 0.0, ref mp_ike_im_Tikme,
                                    pTrfIndex, pTrfIndex,
                                    trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);  // gradient in reference coördinates, i.e. \f$ \nabla_{\vec{xi}} \f$

                                gradOT.Multiply(1.0, invJacobiOT, GradientRef, ResultPreScale, ref mp_ikd_iked_ike);   // gradient in physical coördinates, i.e. \f$ \nabla_{\vec{x}} \f$
                            }
                        }

                        TempBuffer.FreeTempBuffer(iBuf2);

                    }
                }
            }


            TempBuffer.FreeTempBuffer(iBufIN);
            if(Not > 0)
                TempBuffer.FreeTempBuffer(iBufOT);

            /*
            {
                var _BasisValues = grd.ChefBasis.EdgeEval.GetValues(NS, e0, Len, _Basis.Degree);
                        

                var resultINAccCheck = valIN.CloneAs();
                var resultOTAccCheck = valOT.CloneAs();

                for(int e = 0; e < Len; e++) {
                    int iEdge = e + e0;
                    int jCellIN = E2C[iEdge, 0];
                    int jCellOT = E2C[iEdge, 1];
                    int iTrfIN = trfIdx[iEdge, 0];
                    int iTrfOT = trfIdx[iEdge, 1];

                    for(int k = 0; k < NoOfNodes; k++) {
                        int BN = _Basis.GetLength(jCellIN);

                        resultINAccCheck[e, k] *= ResultPreScale;
                        resultOTAccCheck[e, k] *= ResultPreScale;

                        for(int n = 0; n < BN; n++) {
                            double Cinerr = trfCoördinatesIN[e, n] - Coord[jCellIN, n];
                            double Coterr = jCellOT >= 0 ? trfCoördinatesOT[e, n] - Coord[jCellOT, n] : 0.0;

                            if(Math.Abs(Cinerr) > 1.0e-5)
                                Console.WriteLine("44fuckIN" + e);
                            if(Math.Abs(Coterr) > 1.0e-5)
                                Console.WriteLine("44fuckOT" + e);


                            resultINAccCheck[e, k] += Coord[jCellIN, n] * _BasisValues[iTrfIN, k, n];
                            if(jCellOT >= 0)
                                resultOTAccCheck[e, k] += Coord[jCellOT, n] * _BasisValues[iTrfOT, k, n];

                        }

                        double Vinerr = resultINAccCheck[e, k] - valIN[e, k];
                        double Voterr = jCellOT >= 0 ? resultOTAccCheck[e, k] - valOT[e, k] : 0.0;

                        

                        if(Math.Abs(Vinerr) > 1.0e-5)
                            Console.WriteLine("44fuckIN" + e);
                        if(Math.Abs(Voterr) > 1.0e-5)
                            Console.WriteLine("44fuckOT" + e);
                    }

                }

                //resultINAcc.Set(resultINAccCheck);
                //resultOTAcc.Set(resultOTAccCheck);
            }  */
        }

        static MultidimensionalArray.MultiplyProgram mp_i_i_Ti = MultidimensionalArray.MultiplyProgram.Compile("i", "i", "T(i)", true);
        static MultidimensionalArray.MultiplyProgram mp_ikd_iked_ike = MultidimensionalArray.MultiplyProgram.Compile("ikd", "iked", "ike", true);
        static MultidimensionalArray.MultiplyProgram mp_ikd_Tied_ike = MultidimensionalArray.MultiplyProgram.Compile("ikd", "T(i)ed", "ike", true);
        static MultidimensionalArray.MultiplyProgram mp_ike_im_Tikme = MultidimensionalArray.MultiplyProgram.Compile("ike", "im", "T(i)kme", true);
        static MultidimensionalArray.MultiplyProgram mp_ik_im_Tikm = MultidimensionalArray.MultiplyProgram.Compile("ik", "im", "T(i)km", true);


        private static void TransformCoördinatesEdge(int e0, int L, IGridData grd,
            MultidimensionalArray Coördinates, int N_IN, int N_OT, int Degree, bool AffineLinear,
            MultidimensionalArray trfCoördinatesIN, MultidimensionalArray trfCoördinatesOT) //
        {

            int[,] Edge2Cell_g = grd.iGeomEdges.CellIndices;
            int[,] Edge2Cell_l = grd.iGeomEdges.LogicalCellIndices;


            int jCellMin_g = int.MaxValue, jCellMax_g = int.MinValue;
            for(int e = 0; e < L; e++) {
                int iEdge = e + e0;
                int jCellIN = Edge2Cell_g[iEdge, 0];
                int jCellOT = Edge2Cell_g[iEdge, 1];
                {
                    //Debug.Assert(BasisValues.GetLength(2) >= _Basis.GetLength(jCellIN));
                    jCellMin_g = Math.Min(jCellMin_g, jCellIN);
                    jCellMax_g = Math.Max(jCellMax_g, jCellIN);
                    
                }
                if(jCellOT >= 0) {
                    //Debug.Assert(BasisValues.GetLength(2) >= _Basis.GetLength(jCellOT));
                    jCellMin_g = Math.Min(jCellMin_g, jCellOT);
                    jCellMax_g = Math.Max(jCellMax_g, jCellOT);
                }
                Debug.Assert(grd.iGeomEdges.IsEdgeAffineLinear(iEdge) == AffineLinear);
            }

            unsafe {
                fixed(int* pEdge2Cell_g = Edge2Cell_g, pEdge2Cell_l = Edge2Cell_l) {

                    if(AffineLinear) {
                        Debug.Assert(grd is BoSSS.Foundation.Grid.Classic.GridData);
                        Debug.Assert(object.ReferenceEquals(Edge2Cell_g, Edge2Cell_l)); 
                        // this branch only for Classic.GridData

                        MultidimensionalArray scale = grd.ChefBasis.Scaling;

                        // trfCoördinatesIN[i,n] = scale[T(i)]*Coördinates[T(i),n], where T(i) = Edge2Cell[i + e0,0]
                        {
                            trfCoördinatesIN.Multiply(1.0, scale, Coördinates, 0.0, ref mp_in_Ti_Tin,
                                pEdge2Cell_g, pEdge2Cell_g,
                                trfPreOffset_A: (2 * e0), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0), trfCycle_B: 2, trfPostOffset_B: 0);
                        }

                        // trfCoördinatesOT[i,n] = scale[T(i)]*Coördinates[T(i),n], where T(i) = Edge2Cell[i + e0,1]
                        if(N_OT > 0) {
                            trfCoördinatesOT.Multiply(1.0, scale, Coördinates, 0.0, ref mp_in_Ti_Tin,
                                pEdge2Cell_g, pEdge2Cell_g,
                                trfPreOffset_A: (2 * e0 + 1), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: (2 * e0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);
                        }


                    } else {
                        MultidimensionalArray trafo = grd.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(jCellMin_g, jCellMax_g - jCellMin_g + 1, Degree);
                        MultidimensionalArray trafoIN, trafoOT;
                        if(trafo.GetLength(1) > N_IN)
                            trafoIN = trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { trafo.GetLength(0) - 1, N_IN - 1, N_IN - 1 });
                        else
                            trafoIN = trafo;
                        if(N_OT <= 0)
                            trafoOT = null;
                        else if(N_IN == N_OT)
                            trafoOT = trafoIN;
                        else
                            trafoOT = trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { trafo.GetLength(0) - 1, N_OT - 1, N_OT - 1 });
                        {
                            Debug.Assert(trfCoördinatesIN.GetLength(1) == N_IN);
                            Debug.Assert(Coördinates.GetLength(1) == N_IN);


                            // trfCoördinatesIN[i,m] = sum_{n} trafo[T(i),m,n]*Coördinates[T(i),m,n], where T(i) = Edge2Cell[i,0]
                            trfCoördinatesIN.Multiply(1.0, trafoIN, Coördinates, 0.0, ref mp_im_Timn_Tin,
                                pEdge2Cell_g, pEdge2Cell_l,
                                trfPreOffset_A: (2 * e0), trfCycle_A: 2, trfPostOffset_A: -jCellMin_g, 
                                trfPreOffset_B: (2 * e0), trfCycle_B: 2, trfPostOffset_B: 0);
                        }
                        if(N_OT > 0) {
                            Debug.Assert(trfCoördinatesOT.GetLength(1) == N_OT);
                            Debug.Assert(Coördinates.GetLength(1) == N_OT);

                            trfCoördinatesOT.Multiply(1.0, trafoOT, Coördinates, 0.0, ref mp_im_Timn_Tin,
                                pEdge2Cell_g, pEdge2Cell_l,
                                trfPreOffset_A: (2 * e0 + 1), trfCycle_A: 2, trfPostOffset_A: -jCellMin_g, 
                                trfPreOffset_B: (2 * e0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);
                        }
                    }
                }
            }
        }
        static MultidimensionalArray.MultiplyProgram mp_in_Ti_Tin = MultidimensionalArray.MultiplyProgram.Compile("in", "T(i)", "T(i)n", true);
        static MultidimensionalArray.MultiplyProgram mp_im_Timn_Tin = MultidimensionalArray.MultiplyProgram.Compile("im", "T(i)mn", "T(i)n", true);
    }
}
