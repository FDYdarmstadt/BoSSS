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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using ilPSP.Utils;
using ilPSP.Tracing;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    public class LevelSetQuadRuleFactory2 {


        class QRF : IQuadRuleFactory<QuadRule> {
            internal LevelSetQuadRuleFactory2 m_Owner;

            /// <summary>
            /// cache<br/>
            /// key: quadrature order; <br/>
            /// value: cached quadrature rule
            /// </summary>
            internal Dictionary<int, ChunkRulePair<QuadRule>[]> Rules;

            /// <summary>
            /// If there are any cached rules, this method returns their order.
            /// </summary>
            public int[] GetCachedRuleOrders() {
                return Rules.Keys.ToArray();
            }


            public RefElement RefElement {
                get { return m_Owner.Kref; }
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");

#if DEBUG
                if (mask.Except(m_Owner.MaxGrid).NoOfItemsLocally > 0)
                    throw new NotSupportedException("'mask' must be a subset of the cut cells, for my reference element.");
#endif
                if (!Rules.ContainsKey(order))
                    m_Owner.GetQuadRuleSet_Internal(order);

                if (mask.NoOfItemsLocally == m_Owner.MaxGrid.NoOfItemsLocally) {
                    // aggressive
                    return Rules[order];
                } else {
                    var Rule = Rules[order];

                    SubGrid S = new SubGrid((CellMask)mask);
                    var jsub2jcell = S.SubgridIndex2LocalCellIndex;
                    var Ret = new ChunkRulePair<QuadRule>[S.LocalNoOfCells_WithExternal];

                    int L = Ret.Length, H = Rule.Length;
                    int h = 0;
                    for (int jsub = 0; jsub < L; jsub++) {
                        int jCell = jsub2jcell[jsub];
                        
                        Debug.Assert(Rule[h].Chunk.Len == 1);

                        while(jCell > Rule[h].Chunk.i0) {
                            h++;
                        }
                        
                        Debug.Assert(jCell == Rule[h].Chunk.i0);
                        Ret[jsub] = Rule[h];
                    }
                    
                    return Ret;
                }

            }
        }

        
     
        public IQuadRuleFactory<QuadRule> GetSurfaceFactory() {
            return new QRF() {
                m_Owner = this,
                Rules = this.m_SurfaceRules
            };
        }

        

        /// <summary>
        /// grid reference element.
        /// </summary>
        RefElement Kref;

        /// <summary>
        /// index of the respective level-set
        /// </summary>
        protected int LevelSetIndex;

        /// <summary>
        /// if true, the nodes for the surface integration are 'projected' onto the zero-level-set
        /// </summary>
        bool SurfaceNodesOnZeroLevset;

        /// <summary>
        /// if true, the accuracy of the quadrature is checked after solution of the system
        /// </summary>
        bool Docheck;

        /// <summary>
        /// if true, Stoke's theorem is used to create the surface rule
        /// </summary>
        bool UseStokes;

        /// <summary>
        /// if true, Gauß theorem is used to create the surface rule
        /// </summary>
        bool UseGauß;

        /// <summary>
        /// the awesome level set tracker
        /// </summary>
        protected LevelSetTracker tracker;
        

        IQuadRuleFactory<CellBoundaryQuadRule> cellFaceFactory;
        IQuadRuleFactory<CellBoundaryQuadRule> LevelSetBoundaryLineFactory;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="tracker"></param>
        /// <param name="iLevSet"></param>
        /// <param name="_SurfaceNodesOnZeroLevset">if true, the nodes for the surface integration are 'projected' onto the zero-level-set</param>
        /// <param name="_DoCheck">
        /// if true, the accuracy of the quadrature is checked after solution of the system
        /// </param>
        /// <param name="_LevelSetBoundaryLineFactory"></param>
        /// <param name="cellFaceFactory"></param>
        /// <param name="_UseGauß"></param>
        /// <param name="_UseStokes"></param>
        public LevelSetQuadRuleFactory2(LevelSetTracker.LevelSetData levelSetData,
            IQuadRuleFactory<CellBoundaryQuadRule> cellFaceFactory,
            IQuadRuleFactory<CellBoundaryQuadRule> _LevelSetBoundaryLineFactory,
            bool _UseStokes, bool _UseGauß,
            bool _SurfaceNodesOnZeroLevset = false,
            bool _DoCheck = false) 
        {
            this.UseStokes = _UseStokes;
            this.UseGauß = _UseGauß;
            this.LevelSetData = levelSetData;

            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.Docheck = _DoCheck;

            //if(_UseGauß && !object.ReferenceEquals(cellFaceFactory.RefElement, _LevelSetBoundaryLineFactory.RefElement))
            //    throw new ArgumentException("boundary factory and boundary lne factory must operate on the same reference element.");
            //if(_UseStokes && _LevelSetBoundaryLineFactory == null)
            //    throw new ArgumentException();
            
            this.Kref = cellFaceFactory.RefElement;
            if (! levelSetData.GridDat.Grid.RefElements.Contains(Kref, ReferenceComparer.Instance)) {
                throw new ArgumentOutOfRangeException(
                    "simplex", "'simplex' must be a volume - reference element");
            }
            
            this.LevelSetIndex = levelSetData.LevelSetIndex;
            this.cellFaceFactory = cellFaceFactory;
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.LevelSetBoundaryLineFactory = _LevelSetBoundaryLineFactory;

            int iKref = this.tracker.GridDat.Grid.RefElements.IndexOf(Kref);
            this.MaxGrid = levelSetData.GridDat.Cells.GetCells4Refelement(iKref).Intersect(
                levelSetData.Region.GetCutCellMask4LevSet(this.LevelSetIndex));
        }

        /// <summary>
        /// Evaluation of the level-set.
        /// </summary>
        LevelSetTracker.LevelSetData LevelSetData;

        /// <summary>
        /// the intersection of the cut cells for Level Set <see cref="LevelSetIndex"/>
        /// and the subgrid for reference element <see cref="Kref"/>
        /// </summary>
        private CellMask MaxGrid;


        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_SurfaceRules = new Dictionary<int,ChunkRulePair<QuadRule>[]>();




        public static Stopwatch stpwGetQuadRuleSet = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_StokesRHS = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_GaussRHS = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_SolveRHS = new Stopwatch();



        void GetQuadRuleSet_Internal(int order) {
            using (new FuncTrace()) {
                int original_order = order;
                stpwGetQuadRuleSet.Start();
                order = Math.Max(order, 2); // Ansatz won't work below 2!
                int D = this.tracker.GridDat.SpatialDimension; ;
                
                // check arguments, init
                // =====================

                CellMask _mask = this.MaxGrid;

                // subgrid on which the volume rule should be constructed
                // ======================================================
                CellBoundaryQuadratureScheme cellBndSchme = new CellBoundaryQuadratureScheme(this.cellFaceFactory, _mask);
                CellBoundaryQuadratureScheme cellBndLineSchme = null;
                if (this.UseStokes)
                    cellBndLineSchme = new CellBoundaryQuadratureScheme(this.LevelSetBoundaryLineFactory, _mask);

                // set up
                // ======

                int NoOfEqTotal, NoOfStokesEq, NoOfGaußEq;
                DivergenceFreeBasis DivFreeBasis = this.UseGauß ? new DivergenceFreeBasis(tracker.GridDat, this.Kref, order + 1) : null;
                Basis ScalarBasis = this.UseStokes ? new Basis(this.tracker.GridDat, order + 1) : null; // we loose a order of 1 for the volume rule due to the divergence operator

                NoOfGaußEq = DivFreeBasis != null ? (DivFreeBasis.Count/D) : 0;
                NoOfStokesEq = (ScalarBasis != null ? ScalarBasis.Length : 0) * D;

                NoOfEqTotal = NoOfGaußEq + NoOfStokesEq;
                
                // define Nodes
                // ============
                NodeSet NodeSet = null;
                {
                    if (this.Kref.GetType() == typeof(Square)) {


                        int K = (int)Math.Ceiling(Math.Sqrt(NoOfEqTotal*1.7)) + 1;
                        
                        var Nodes1D = GenericBlas.Linspace(-1, 1, K);

                        var _NodeSet = MultidimensionalArray.Create(K*K, 2);
                        int n = 0;
                        for (int i = 0; i < K /*&& n <= NoOfEq*1.1*/; i++) {
                            for (int j = 0; j < K /*&& n <= NoOfEq*1.1*/; j++) {
                                _NodeSet[n, 0] = Nodes1D[i];
                                _NodeSet[n, 1] = Nodes1D[j];
                                n++;
                            }
                        }

                        NodeSet = new NodeSet(this.Kref, _NodeSet);
                    } else {
                        for (int o = 1; o < 1000000; o++) {
                            var qr = Kref.GetBruteForceQuadRule(o, 0);
                            if (qr.NoOfNodes >= (NoOfEqTotal*1.1)) {
                                NodeSet = qr.Nodes;
                                break;
                            }
                        }
                    }

                }
                int NoOfNodes = NodeSet.GetLength(0);
                
                // find RHS integrals
                // ==================

                MultidimensionalArray RHS_Gauß = null;
                if(this.UseGauß) {
                    stpwGetQuadRuleSet_GaussRHS.Start();
                    RHS_Gauß = this.GaußAnsatzRHS(DivFreeBasis, cellBndSchme, _mask, order);
                    stpwGetQuadRuleSet_GaussRHS.Stop();
                    Debug.Assert(RHS_Gauß.Dimension == 2);
                    Debug.Assert(RHS_Gauß.GetLength(0) == NoOfGaußEq);
                    Debug.Assert(RHS_Gauß.GetLength(1) == _mask.NoOfItemsLocally);
                }

                MultidimensionalArray RHS_Stokes = null;
                if (this.UseStokes) {
                    stpwGetQuadRuleSet_StokesRHS.Start();
                    RHS_Stokes = this.StokesAnsatzRHS(ScalarBasis, cellBndLineSchme, _mask, order);
                    stpwGetQuadRuleSet_StokesRHS.Stop();
                    Debug.Assert(RHS_Stokes.Dimension == 2);
                    Debug.Assert(RHS_Stokes.GetLength(0) == NoOfStokesEq);
                    Debug.Assert(RHS_Stokes.GetLength(1) == _mask.NoOfItemsLocally);
                }

                // construct da rule!
                // ==================
                ChunkRulePair<QuadRule>[] SurfaceRule;
                {
                    SurfaceRule = new ChunkRulePair<QuadRule>[_mask.NoOfItemsLocally];
                    var grddat = this.tracker.GridDat;
                    
                    // loop over cells in subgrid...
                    int jSub = 0;
                    foreach(int jCell in _mask.ItemEnum) { // loop over cells in the mask


                        // setup System
                        // ============

                        NodeSet surfNodes;
                        if(this.SurfaceNodesOnZeroLevset)
                            surfNodes = ProjectOntoLevset(jCell, NodeSet);
                        else
                            surfNodes = NodeSet;

                        MultidimensionalArray metrics;
                        {
                            int iKref = grddat.Cells.GetRefElementIndex(jCell);
                            metrics = LevelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(surfNodes, jCell, 1);
                        }

                        MultidimensionalArray Mtx_Gauss = null;
                        if(this.UseGauß) {
                            Mtx_Gauss = GaußAnsatzMatrix(DivFreeBasis, surfNodes, jCell);
                            Debug.Assert(Mtx_Gauss.Dimension == 2);
                            Debug.Assert(Mtx_Gauss.GetLength(0) == NoOfGaußEq);
                            Debug.Assert(Mtx_Gauss.GetLength(1) == NoOfNodes);
                        }

                        MultidimensionalArray Mtx_Stokes = null;
                        if(this.UseStokes) {
                            Mtx_Stokes = this.StokesAnsatzMatrix(ScalarBasis, surfNodes, jCell);
                            Debug.Assert(Mtx_Stokes.Dimension == 2);
                            Debug.Assert(Mtx_Stokes.GetLength(0) == NoOfStokesEq);
                            Debug.Assert(Mtx_Stokes.GetLength(1) == NoOfNodes);
                            for(int i = 0; i < NoOfStokesEq; i++) {
                                for(int j = 0; j < NoOfNodes; j++) {
                                    Mtx_Stokes[i, j] /= metrics[0, j];
                                }
                            }
                        }

                        stpwGetQuadRuleSet_SolveRHS.Start();

                        // convert to FORTRAN order
                        MultidimensionalArray _Mtx = MultidimensionalArray.Create(NoOfEqTotal, NoOfNodes);
                        if(this.UseGauß) {
                            _Mtx.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NoOfGaußEq - 1, NoOfNodes - 1 }).Set(Mtx_Gauss);
                        }
                        if(this.UseStokes) {
                            _Mtx.ExtractSubArrayShallow(new int[] { NoOfGaußEq, 0 }, new int[] { NoOfStokesEq + NoOfGaußEq - 1, NoOfNodes - 1 }).Set(Mtx_Stokes);
                        }


                        // convert to FORTRAN order
                        Debug.Assert(NoOfNodes >= NoOfEqTotal);
                        MultidimensionalArray _RHS = MultidimensionalArray.Create(NoOfNodes, 1); // this is also output, so it must be larger!
                        {
                            for(int i = 0; i < NoOfGaußEq; i++) {
                                _RHS[i, 0] = RHS_Gauß[i, jSub];
                            }
                            for(int i = 0; i < NoOfStokesEq; i++) {
                                _RHS[i + NoOfGaußEq, 0] = RHS_Stokes[i, jSub];
                            }
                        }

                        MultidimensionalArray __RHS = null, __Mtx = null;
                        if(this.Docheck) {
                            // values used for testing:
                            __RHS = MultidimensionalArray.Create(NoOfEqTotal);
                            for(int i = 0; i < NoOfEqTotal; i++)
                                __RHS[i] = _RHS[i, 0];
                            __Mtx = MultidimensionalArray.Create(_Mtx.NoOfRows, _Mtx.NoOfCols);
                            __Mtx.SetMatrix(_Mtx);
                        }

                        // solve system
                        // ============

                        //int M = _Mtx.NoOfRows;
                        //int N = _Mtx.NoOfCols;
                        //LAPACK.F77_LAPACK.DGELSY(M, N, _Mtx.Entries, _RHS.Entries, 1, 1.0e-14);
                        _Mtx.LeastSquareSolve(_RHS);


                        if(this.Docheck) {
                            // Probe:
                            MultidimensionalArray X = MultidimensionalArray.Create(NoOfNodes); // weights
                            X.ResizeShallow(NoOfNodes, 1).SetMatrix(_RHS);
                            __RHS.Multiply(-1.0, __Mtx, X, 1.0, "j", "jk", "k");
                            double L2_ERR = __RHS.L2Norm();
                            if(L2_ERR > 1.0e-7)
                                throw new ApplicationException("Quadrature rule in cell " + jCell + " seems to be not very precise: L2_ERR = " + L2_ERR);
                            //Debug.Assert(L2_ERR < 1.0e-8, "Quadrature rule in cell " + jCell + " seems to be not very precise: L2_ERR = " + L2_ERR);
                            //if (L2_ERR > 1.0e-9)
                            //    Console.WriteLine("Warning: Quadrature rule in cell " + jCell + ": L2_ERR = " + L2_ERR);
                        }


                        stpwGetQuadRuleSet_SolveRHS.Stop();

                        // return da rule!
                        // ===============

                        {


                            {
                                // the surface rule
                                // ----------------

                                QuadRule qr_l = new QuadRule() {
                                    OrderOfPrecision = order,
                                    Weights = MultidimensionalArray.Create(NoOfNodes),
                                    Nodes = surfNodes
                                };


                                for(int k = 0; k < NoOfNodes; k++) {
                                    qr_l.Weights[k] = _RHS[k, 0] / metrics[0, k];
                                }

                                SurfaceRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);
                            }
                        }
                        jSub++;
                    }
                }


                stpwGetQuadRuleSet.Stop();

                if(!this.m_SurfaceRules.ContainsKey(original_order)) {
                    this.m_SurfaceRules.Add(original_order, SurfaceRule);
                }
                
            }
        }

        

        protected MultidimensionalArray GaußAnsatzMatrix(DivergenceFreeBasis TestBasis, NodeSet SurfQrNodes, int jCell) {
            int N = TestBasis.Count;
            int D = this.tracker.GridDat.Grid.SpatialDimension;
            Debug.Assert(N % D == 0);
            N /= D;
            int ksurf = SurfQrNodes.GetLength(0);
            Debug.Assert(SurfQrNodes.Dimension == 2);
            Debug.Assert(SurfQrNodes.GetLength(1) == D);
            int iKref = this.tracker.GridDat.Cells.GetRefElementIndex(jCell);

            MultidimensionalArray EdgMatrix = MultidimensionalArray.Create(N, ksurf);

            // evaluate Basis and Gradient of Basis
            var Phi = TestBasis.Values.GetValues(SurfQrNodes); // test function, n
            var LevelSetNormals = this.LevelSetData.GetLevelSetReferenceNormals(SurfQrNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);
            //metrics = tracker.GetLevelSetNormalReferenceToPhysicalMetrics(this.LevelSetIndex, 0, jCell, 1);
            
            // multiply
            for(int k = 0; k < ksurf; k++) { // loop over nodes...
                for(int n = 0; n < N; n++) { // loop over basis polynomials...
                    double acc = 0.0;
                    for(int d = 0; d < D; d++) {// loop over spatial dimension...
                        acc += Phi[k, n*D + d] * LevelSetNormals[k, d];
                    }
                    EdgMatrix[n, k] = acc;
                }
            }
            
            // resize and return:
            return EdgMatrix;
        }


        NodeSet ProjectOntoLevset(int jCell, NodeSet Nodes) {
            int D = Nodes.GetLength(1);
            int NoOfNodes = Nodes.GetLength(0);
            var m_Context = this.tracker.GridDat;
            LevelSet LevSet = (LevelSet)(this.tracker.LevelSets[this.LevelSetIndex]);

            MultidimensionalArray LevSetValues = MultidimensionalArray.Create(1, NoOfNodes);
            MultidimensionalArray LevSetGrad = MultidimensionalArray.Create(1, NoOfNodes, D);

            MultidimensionalArray x0_i_Local = MultidimensionalArray.Create(1, NoOfNodes, D);
            MultidimensionalArray x0_i_Global = MultidimensionalArray.Create(1, NoOfNodes, D); // quadrature nodes in global coordinates

            MultidimensionalArray x0_ip1_Local = MultidimensionalArray.Create(1, NoOfNodes, D);
            MultidimensionalArray x0_ip1_Global = MultidimensionalArray.Create(NoOfNodes, D);

            // set initial value;


            x0_i_Local.SetSubArray(Nodes, 0, -1, -1);

            int NN = NoOfNodes;
            for(int i = 0; i < 10; i++) {

                double radiusError = 0;


                int j = jCell;


                LevSet.Evaluate(j, 1, Nodes, LevSetValues, 0, 0.0);
                LevSet.EvaluateGradient(j, 1, Nodes, LevSetGrad);


                m_Context.TransformLocal2Global(new NodeSet(this.Kref, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)), j, 1, x0_i_Global, 0);

                for(int nn = 0; nn < NN; nn++) {

                    double sc = 0;
                    for(int d = 0; d < D; d++) {
                        sc += LevSetGrad[0, nn, d].Pow2();
                    }


                    for(int d = 0; d < D; d++) {
                        double xd = x0_i_Global[0, nn, d] - LevSetGrad[0, nn, d] * LevSetValues[0, nn] / sc;
                        x0_ip1_Global[nn, d] = xd;
                    }

                    radiusError += Math.Abs(LevSetValues[0, nn]);

                }

                m_Context.TransformGlobal2Local(x0_ip1_Global, x0_ip1_Local, j, 1, 0);


                // next iter: x0_i <- x0_{i+1}
                x0_i_Local.Set(x0_ip1_Local);
                Nodes = (new NodeSet(this.Kref, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)));
            }

            return Nodes;
        }

        protected MultidimensionalArray GaußAnsatzRHS(DivergenceFreeBasis TestBasis, CellBoundaryQuadratureScheme cellBndScheme, CellMask _mask, int order) {
            var _Context = this.tracker.GridDat;
            int D = this.tracker.GridDat.Grid.SpatialDimension;
            int N = TestBasis.Count;
            var coordSys = CoordinateSystem.Reference;
            var LsTrk = this.tracker;
            int Nrhs = _mask.NoOfItemsLocally;

            Debug.Assert(N % D == 0);
            N /= D;
            MultidimensionalArray RHS = MultidimensionalArray.Create(N, Nrhs);

            var splx = this.Kref;
            int NoOfFaces = splx.NoOfFaces;
            //var normals = _Context.GridDat.Normals;
            CellBoundaryQuadrature<CellBoundaryQuadRule> qBnd = null;

            int jSgrd = 0;
            qBnd = CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { N },
                _Context, cellBndScheme.Compile(_Context, order),
                delegate(int i0, int Length, CellBoundaryQuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet Nodes = QR.Nodes;
                    MultidimensionalArray BasisValues;
                    if (coordSys == CoordinateSystem.Physical) {
                        //BasisValues = TestBasis.CellEval(Nodes, i0, Length);
                        throw new NotImplementedException("todo");
                    } else if (coordSys == CoordinateSystem.Reference) {
                        BasisValues = TestBasis.Values.GetValues(Nodes);
                    } else
                        throw new NotImplementedException();

                    for (int i = 0; i < Length; i++) { // loop over cells

                        CellBoundaryQuadRule cR = qBnd.CurrentRule;
                        int[] NodesPerEdge = cR.NumbersOfNodesPerFace;
                        Debug.Assert(object.ReferenceEquals(splx, cR.RefElement));

                        int iNode = 0;

                        Debug.Assert(NoOfFaces == NodesPerEdge.Length);
                        for (int e = 0; e < NoOfFaces; e++) { // loop over the faces of the cell
                            for (int _n = 0; _n < NodesPerEdge[e]; _n++) { // loop over nodes in one edge
                                for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
                                    double acc = 0;
                                    for (int d = 0; d < D; d++) { // loop over spatial directions

                                        if (coordSys == CoordinateSystem.Physical) {
                                            throw new NotImplementedException("todo");
                                            //int q = _Context.GridDat.LocalCellIndexToEdges[i+i0, e];
                                            //int iEdge = Math.Abs(q) - 1;
                                            //double Nsign = Math.Sign(q);
                                            //double Nd = normals[iEdge, d];
                                            //EvalResult[i, iNode, n, d] = BasisValues[i, iNode, n]*Nd*Nsign;
                                        } else {
                                            Debug.Assert(coordSys == CoordinateSystem.Reference);
                                            double Nd = splx.FaceNormals[e, d];
                                            //Debug.Assert(Nd == normals[iEdge, d]*Nsign);
                                            acc += BasisValues[iNode, n*D + d]*Nd;
                                        }
                                    }
                                    EvalResult[i, iNode, n] = acc;
                                }

                                iNode++;
                            }
                        }
                        Debug.Assert(iNode == EvalResult.GetLength(1));
                    }
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, jSgrd }, new int[] { N-1, jSgrd-1 });

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0 }, new int[] { i-1, e-1, N-1 });
                            ResPart.Acc(1.0, ip);
                        }
                        jSgrd++;
                    }
                },
                cs: coordSys);
            qBnd.Execute();

            
            return RHS;
        }

        MultidimensionalArray StokesAnsatzMatrix(Basis TestBasis, NodeSet surfaceNodes, int jCell) {
            int N = TestBasis.Length;
            int NoOfNodes = surfaceNodes.GetLength(0);
            int D = surfaceNodes.GetLength(1);
            var GridDat = this.tracker.GridDat;
            Debug.Assert(D == GridDat.SpatialDimension);
            int iKref = GridDat.Cells.GetRefElementIndex(jCell);
            var scalings = GridDat.Cells.JacobiDet;
            int iLevSet = this.LevelSetIndex;

            if (!GridDat.Cells.IsCellAffineLinear(jCell))
                throw new NotSupportedException();

            
            //var Phi = TestBasis.Evaluate(0);              // reference
            //var GradPhi = TestBasis.EvaluateGradient(0);  // reference
            var Phi = TestBasis.CellEval(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);                  // physical
            var GradPhi = TestBasis.CellEvalGradient(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1, -1);  // physical

            //var LevsetNormal = this.LsTrk.GetLevelSetReferenceNormals(iLevSet, 0, jCell, 1);  // reference
            //var Curvature = this.LsTrk.GetLevelSetReferenceCurvature(iLevSet, 0, jCell, 1);   // reference
            var LevsetNormal = this.LevelSetData.GetLevelSetNormals(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1); // physical
            var Curvature = MultidimensionalArray.Create(1, NoOfNodes);                                  // physical
            ((LevelSet)(this.tracker.LevelSets[iLevSet])).EvaluateTotalCurvature(jCell, 1, surfaceNodes, Curvature);  // physical

            var Coeffs = MultidimensionalArray.Create(D, N, NoOfNodes);

            if (D == 2) {
                for (int k = 0; k < NoOfNodes; k++) { // loop over nodes
                    double Nx = LevsetNormal[k, 0];
                    double Ny = LevsetNormal[k, 1];
                    double kappa = Curvature[0, k];
                    
                    
                    double Prj_11 = 1.0 - Nx*Nx, Prj_12 = -Nx*Ny,
                        Prj_21 = -Ny*Nx, Prj_22 = 1.0 - Ny*Ny;

                    for (int n = 0; n < N; n++) {
                        double Phi_kn = Phi[k, n];
                        double dPhi_dx_kn = GradPhi[k, n, 0];
                        double dPhi_dy_kn = GradPhi[k, n, 1];

                        Coeffs[0, n, k] = -Phi_kn*kappa*Nx + Prj_11*dPhi_dx_kn + Prj_12*dPhi_dy_kn;
                        Coeffs[1, n, k] = -Phi_kn*kappa*Ny + Prj_21*dPhi_dx_kn + Prj_22*dPhi_dy_kn;
                    }
                }
            } else if (D==3) {
                throw new NotImplementedException("to do.");
            } else {
                throw new NotSupportedException("Unknown spatial dimension.");
            }
            
            Coeffs.Scale(scalings[jCell]);  // physical

            return Coeffs.ResizeShallow(N*D, NoOfNodes);
        }

        static void tangente(double[] SurfN, double[] EdgeN, double[] tan) {
            Debug.Assert(SurfN.Length == EdgeN.Length);
            if (SurfN.Length != 2)
                throw new NotSupportedException();

            tan[0] = -SurfN[1];
            tan[1] = SurfN[0];

            if (GenericBlas.InnerProd(tan, EdgeN) < 0.0)
                tan.ScaleV(-1.0);
        }

        MultidimensionalArray StokesAnsatzRHS(Basis TestBasis, CellBoundaryQuadratureScheme cellBndSchme, CellMask _mask, int order) {
            var GridDat = this.tracker.GridDat;
            CellBoundaryQuadrature<CellBoundaryQuadRule> qBnd = null;

            int N = TestBasis.Length;
            int D = GridDat.SpatialDimension;
            MultidimensionalArray RHS = MultidimensionalArray.Create(D, N, _mask.NoOfItemsLocally);

            double[] CellN = new double[D];       // cell normal
            double[] SurfN = new double[D];       // level-set normal
            double[] OutwardTang = new double[D]; // level-set tangent, outward of cell

            if (D != 2)
                throw new NotSupportedException("Currently only supported for spatial dimension of 2.");

            //MultidimensionalArray Nudes = null;
            int jSgrd = 0;
            qBnd = CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { D, N },
                GridDat, cellBndSchme.Compile(GridDat, order),
                delegate(int i0, int Length, CellBoundaryQuadRule NS, MultidimensionalArray EvalResult) { // Evaluate

                    //MultidimensionalArray BasisValues = TestBasis.Evaluate(0);                 // reference
                    //var LSNormals = LsTrk.GetLevelSetReferenceNormals(iLevSet, 0, i0, Length); // reference
                    MultidimensionalArray BasisValues = TestBasis.CellEval(NS.Nodes, i0, Length);         // physical
                    MultidimensionalArray LSNormals = this.LevelSetData.GetLevelSetNormals(NS.Nodes, i0, Length);  // physical

                    for (int i = 0; i < Length; i++) { // loop over cells
                        //if(i0 + i == 1) {
                        //    EvalResult.ExtractSubArrayShallow(i, -1, -1, -1).Clear();
                        //    continue;
                        //}

                        CellBoundaryQuadRule cR = qBnd.CurrentRule;
                        int[] NodesPerEdge = cR.NumbersOfNodesPerFace;
                        var Kref = cR.RefElement;
                        int NoOfFaces = Kref.NoOfFaces;
                        int iNode = 0;

                        Debug.Assert(NoOfFaces == NodesPerEdge.Length);
                        for (int e = 0; e < NoOfFaces; e++) { // loop over the faces of the cell

                            if (NodesPerEdge[e] <= 0)
                                continue;

                            // reference: 
                            //for (int d = 0; d < D; d++) {
                            //    CellN[d] = Kref.FaceNormals[e, d];
                            //}
                            // ~~~~

                            // physical:
                            var FaceNodes = new NodeSet(Kref, cR.Nodes.ExtractSubArrayShallow(new int[] { iNode, 0 }, new int[] { iNode + NodesPerEdge[e] - 1, D - 1 }));
                            var FaceNormals = MultidimensionalArray.Create(NodesPerEdge[e], D);
                            GridDat.Edges.GetNormalsForCell(FaceNodes, i0, e, FaceNormals);
                            // ~~~~

                            for (int _n = 0; _n < NodesPerEdge[e]; _n++) { // loop over nodes in one edge
                                for (int d = 0; d < D; d++) {
                                    SurfN[d] = LSNormals[i, iNode, d];
                                    CellN[d] = FaceNormals[_n, d]; // physical
                                }
                                tangente(SurfN, CellN, OutwardTang);

                                for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
                                    for (int d = 0; d < D; d++) { // loop over spatial direction
                                        EvalResult[i, iNode, d, n] = BasisValues[i, iNode, n]*OutwardTang[d]; // physical
                                        //EvalResult[i, iNode, d, n] = BasisValues[iNode, n]*OutwardTang[d]; // reference
                                    }
                                }

                                iNode++;
                            }
                        }
                        Debug.Assert(iNode == EvalResult.GetLength(1));
                    }
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, 0, jSgrd }, new int[] { D - 1, N - 1, jSgrd - 1 });
                        int NoOfFaces = ResultsOfIntegration.GetLength(1);

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0, 0 }, new int[] { i-1, e-1, D-1, N-1 });                            
                            ResPart.Acc(1.0, ip);
                        }
                        jSgrd++;
                    }
                },
                cs: CoordinateSystem.Physical);
            qBnd.Execute();

            var ret = RHS.ResizeShallow(N*D, _mask.NoOfItemsLocally);
            return ret;
        }


    }
}

