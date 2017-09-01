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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// HMF-surface integrals in 2D, based on Stokes integral theorem.
    /// </summary>
    public class SurfaceStokes_2D {



        class QRF : IQuadRuleFactory<QuadRule> {
            internal SurfaceStokes_2D m_Owner;

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
                get {
                    return m_Owner.RefElement;
                }
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");

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

                        while (jCell > Rule[h].Chunk.i0) {
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
        RefElement RefElement;

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
        /// the awesome level set tracker
        /// </summary>
        protected LevelSetTracker tracker;

        IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory;

        IQuadRuleFactory<CellBoundaryQuadRule> LevelSetBoundaryLineFactory;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tracker"></param>
        /// <param name="iLevSet"></param>
        /// <param name="_SurfaceNodesOnZeroLevset">if true, the nodes for the surface integration are 'projected' onto the zero-level-set</param>
        /// <param name="_DoCheck">
        /// if true, the accuracy of the quadrature is checked after solution of the system
        /// </param>
        /// <param name="_LevelSetBoundaryLineFactory"></param>
        /// <param name="cellBoundaryFactory"></param>
        public SurfaceStokes_2D(LevelSetTracker tracker, int iLevSet,
            IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory,
            IQuadRuleFactory<CellBoundaryQuadRule> _LevelSetBoundaryLineFactory,
            bool _SurfaceNodesOnZeroLevset = false,
            bool _DoCheck = false) {
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.Docheck = _DoCheck;

            this.RefElement = cellBoundaryFactory.RefElement;
            if (!tracker.GridDat.Grid.RefElements.Contains(RefElement, ReferenceComparer.Instance)) {
                throw new ArgumentOutOfRangeException(
                    "simplex", "'simplex' must be a volume - reference element");
            }

            this.tracker = tracker;
            this.LevelSetIndex = iLevSet;
            this.cellBoundaryFactory = cellBoundaryFactory;
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.LevelSetBoundaryLineFactory = _LevelSetBoundaryLineFactory;

            int iKref = this.tracker.GridDat.Grid.RefElements.IndexOf(RefElement);
            this.MaxGrid = this.tracker.GridDat.Cells.GetCells4Refelement(iKref).Intersect(
                tracker._Regions.GetCutCellMask4LevSet(this.LevelSetIndex));
        }

        /// <summary>
        /// the intersection of the cut cells for Level Set <see cref="LevelSetIndex"/>
        /// and the subgrid for reference element <see cref="RefElement"/>
        /// </summary>
        private CellMask MaxGrid;


        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_SurfaceRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();


        public static Stopwatch stpwGetQuadRuleSet = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_StokesRHS = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_SolveRHS = new Stopwatch();


        private int OrderToInternalOrder(int o) {
            /*
            if (o < 2)
                o = 1;
            o += 1;
            return o;
            */
            return o;
        }


        void GetQuadRuleSet_Internal(int ReqOrder) {
            using (new FuncTrace()) {

                int IntOrder = OrderToInternalOrder(ReqOrder);

                stpwGetQuadRuleSet.Start();

                Debug.Assert(IntOrder >= 2);

                if (this.m_SurfaceRules.ContainsKey(ReqOrder))
                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky
                if (this.m_SurfaceRules.ContainsKey(IntOrder))
                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky

                // check arguments, init
                // =====================

                CellMask _mask = this.MaxGrid;

                // subgrid on which the volume rule should be constructed
                // ======================================================
                CellBoundaryQuadratureScheme cellBndSchme = new CellBoundaryQuadratureScheme(this.cellBoundaryFactory, _mask);
                CellBoundaryQuadratureScheme cellBndLineSchme = new CellBoundaryQuadratureScheme(this.LevelSetBoundaryLineFactory, _mask);

                // set up
                // ======

                int NoOfEqTotal;
                var TestBasis = new Basis(this.tracker.GridDat, IntOrder); // we loose a order of 1 for the volume rule due to the divergence operator
                NoOfEqTotal = TestBasis.Length;

                // define Nodes
                // ============
                NodeSet NodeSet = null;
                {
                    if (this.RefElement.GetType() == typeof(Square)) {

                        int K = (int)Math.Ceiling(Math.Sqrt(NoOfEqTotal * 1.75)) + 1;

                        var Nodes1D = GenericBlas.Linspace(-1, 1, K);

                        var _NodeSet = new NodeSet(this.RefElement, K * K, 2);
                        int n = 0;
                        for (int i = 0; i < K /*&& n <= NoOfEq*1.1*/; i++) {
                            for (int j = 0; j < K /*&& n <= NoOfEq*1.1*/; j++) {
                                _NodeSet[n, 0] = Nodes1D[i];
                                _NodeSet[n, 1] = Nodes1D[j];
                                n++;
                            }
                        }


                        //NodeSet = MultidimensionalArray.Create(n, 2);
                        //NodeSet.Set(_NodeSet.ExtractSubArrayShallow(new int[] { 0, 0}, new int[] { n - 1, 1 }));
                        NodeSet = _NodeSet;
                    } else {
                        for (int o = 1; o < 1000000; o++) {
                            var qr = RefElement.GetBruteForceQuadRule(o, 0);
                            if (qr.NoOfNodes >= (NoOfEqTotal * 1.1)) {
                                NodeSet = qr.Nodes;
                                break;
                            }
                        }
                    }

                    //double ratio = ((double)NodeSet.GetLength(0)) / ((double)NoOfEqTotal);
                    //Console.WriteLine("nodes over eq " + ratio);

                    NodeSet.LockForever();
                }
                int NoOfNodes = NodeSet.GetLength(0);
                Debug.Assert(NoOfNodes >= NoOfEqTotal);


                // find RHS integrals
                // ==================


                MultidimensionalArray RHS_Stokes = null;
                stpwGetQuadRuleSet_StokesRHS.Start();
                RHS_Stokes = this.StokesAnsatzRHS(TestBasis, cellBndLineSchme, _mask, IntOrder);
                stpwGetQuadRuleSet_StokesRHS.Stop();
                Debug.Assert(RHS_Stokes.Dimension == 2);
                Debug.Assert(RHS_Stokes.GetLength(0) == NoOfEqTotal);
                Debug.Assert(RHS_Stokes.GetLength(1) == _mask.NoOfItemsLocally);

                // construct da rule!
                // ==================
                ChunkRulePair<QuadRule>[] SurfaceRule;
                {
                    SurfaceRule = new ChunkRulePair<QuadRule>[_mask.NoOfItemsLocally];
                    var grddat = this.tracker.GridDat;
                    int D = grddat.SpatialDimension;

                    // loop over cells in subgrid...
                    int jSub = 0;
                    foreach (int jCell in _mask.ItemEnum) { // loop over cells in the mask


                        // setup System
                        // ============

                        NodeSet VolNodes = NodeSet;
                        NodeSet surfNodes;
                        if (this.SurfaceNodesOnZeroLevset)
                            surfNodes = ProjectOntoLevset(jCell, VolNodes);
                        else
                            surfNodes = VolNodes;


                        MultidimensionalArray Mtx_Stokes = null;
                        {
                            Mtx_Stokes = this.StokesAnsatzMatrix(TestBasis, surfNodes, jCell);
                            Debug.Assert(Mtx_Stokes.Dimension == 2);
                            Debug.Assert(Mtx_Stokes.GetLength(0) == NoOfEqTotal);
                            Debug.Assert(Mtx_Stokes.GetLength(1) == NoOfNodes);

                        }

                        stpwGetQuadRuleSet_SolveRHS.Start();

                        Debug.Assert(NoOfNodes >= NoOfEqTotal);
                        MultidimensionalArray _RHS = MultidimensionalArray.Create(NoOfNodes, 1); // this is also output, so it must be larger!
                        {
                            for (int i = 0; i < NoOfEqTotal; i++) {
                                _RHS[i, 0] = RHS_Stokes[i, jSub];
                            }

                        }

                        MultidimensionalArray __RHS = null, __Mtx = null;
                        if (this.Docheck) {
                            // values used for testing:
                            __RHS = MultidimensionalArray.Create(NoOfEqTotal);
                            for (int i = 0; i < NoOfEqTotal; i++)
                                __RHS[i] = _RHS[i, 0];
                            __Mtx = MultidimensionalArray.Create(Mtx_Stokes.NoOfRows, Mtx_Stokes.NoOfCols);
                            __Mtx.SetMatrix(Mtx_Stokes);
                        }

                        // solve system
                        // ============

                        Mtx_Stokes.LeastSquareSolve(_RHS);


                        if (this.Docheck) {
                            // Probe:
                            MultidimensionalArray X = MultidimensionalArray.Create(NoOfNodes); // weights
                            X.ResizeShallow(NoOfNodes, 1).SetMatrix(_RHS);

                            double MaxWeight = X.Max(wi => wi.Pow2());
                            double MinWeight = X.Min(wi => wi.Pow2());
                            double wr = MaxWeight / MinWeight;
                            //Console.WriteLine("MaxWight: " + MaxWeight  + ",  Weight Ratio: " + wr);

                            __RHS.Multiply(-1.0, __Mtx, X, 1.0, "j", "jk", "k");
                            double L2_ERR = __RHS.L2Norm();
                            if (L2_ERR > 1.0e-7) {
                                Console.WriteLine("Un-precise quadrature order " + IntOrder + " rule in cell " + jCell + ": L2_ERR = " + L2_ERR);
                                //throw new ApplicationException("Quadrature rule in cell " + jCell + " seems to be not very precise: L2_ERR = " + L2_ERR);
                            }
                        }


                        stpwGetQuadRuleSet_SolveRHS.Stop();

                        // return da rule!
                        // ===============
                        var metrics = tracker.GetLevelSetNormalReferenceToPhysicalMetrics(this.LevelSetIndex, surfNodes, jCell, 1);

                        {
                            {
                                // the surface rule
                                // ----------------

                                QuadRule qr_l = new QuadRule() {
                                    OrderOfPrecision = IntOrder,
                                    Weights = MultidimensionalArray.Create(NoOfNodes),
                                    Nodes = surfNodes
                                };

                                for (int k = 0; k < NoOfNodes; k++) {
                                    qr_l.Weights[k] = _RHS[k, 0] / metrics[0, k];
                                }

                                SurfaceRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);
                            }
                        }
                        jSub++;
                    }
                }


                stpwGetQuadRuleSet.Stop();

                this.m_SurfaceRules.Add(IntOrder, SurfaceRule);
                if (ReqOrder != IntOrder) {
                    this.m_SurfaceRules.Add(ReqOrder, SurfaceRule);
                }
            }
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
            for (int i = 0; i < 10; i++) {

                double radiusError = 0;


                int j = jCell;


                LevSet.Evaluate(j, 1, Nodes, LevSetValues, 0, 0.0);
                LevSet.EvaluateGradient(j, 1, Nodes, LevSetGrad);


                m_Context.TransformLocal2Global(new NodeSet(this.RefElement, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)), j, 1, x0_i_Global, 0);

                for (int nn = 0; nn < NN; nn++) {

                    double sc = 0;
                    for (int d = 0; d < D; d++) {
                        sc += LevSetGrad[0, nn, d].Pow2();
                    }


                    for (int d = 0; d < D; d++) {
                        double xd = x0_i_Global[0, nn, d] - LevSetGrad[0, nn, d] * LevSetValues[0, nn] / sc;
                        x0_ip1_Global[nn, d] = xd;
                    }

                    radiusError += Math.Abs(LevSetValues[0, nn]);

                }

                m_Context.TransformGlobal2Local(x0_ip1_Global, x0_ip1_Local, j, 1, 0);


                // next iter: x0_i <- x0_{i+1}
                x0_i_Local.Set(x0_ip1_Local);
                Nodes = (new NodeSet(this.RefElement, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)));
            }

            return Nodes;
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


            var Phi = TestBasis.Evaluate(surfaceNodes);              // reference
            var GradPhi = TestBasis.EvaluateGradient(surfaceNodes);  // reference
            var LevsetNormal = this.tracker.GetLevelSetReferenceNormals(iLevSet, surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);   // reference

            var Coeffs = MultidimensionalArray.Create(N, NoOfNodes);

            if (D == 2) {
                for (int k = 0; k < NoOfNodes; k++) { // loop over nodes
                    double N1 = LevsetNormal[k, 0];
                    double N2 = LevsetNormal[k, 1];



                    for (int n = 0; n < N; n++) {
                        double Phi_kn = Phi[k, n];
                        double dPhi_dx_kn = GradPhi[k, n, 0];
                        double dPhi_dy_kn = GradPhi[k, n, 1];

                        Coeffs[n, k] = dPhi_dy_kn * N1 - dPhi_dx_kn * N2;
                    }
                }
            } else if (D == 3) {
                throw new NotImplementedException("to do.");
            } else {
                throw new NotSupportedException("Unknown spatial dimension.");
            }

            //Coeffs.Scale(scalings[jCell]);  // physical

            return Coeffs;
        }

        

        MultidimensionalArray StokesAnsatzRHS(Basis TestBasis, CellBoundaryQuadratureScheme cellBndSchme, CellMask _mask, int order) {
            var GridDat = this.tracker.GridDat;
            CellBoundaryQuadrature<CellBoundaryQuadRule> qBnd = null;
            int iLevSet = this.LevelSetIndex;
            int N = TestBasis.Length;
            int D = GridDat.SpatialDimension;
            MultidimensionalArray RHS = MultidimensionalArray.Create(N, _mask.NoOfItemsLocally);

            double CellN0, CellN1; // cell normal
            double Tang0, Tang1;  // level-set tangent, outward of cell

            if (D != 2)
                throw new NotSupportedException("Currently only supported for spatial dimension of 2.");

            //MultidimensionalArray Nudes = null;
            int jSgrd = 0;
            qBnd = CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { N },
                GridDat, cellBndSchme.Compile(GridDat, order),
                delegate (int i0, int Length, CellBoundaryQuadRule NS, MultidimensionalArray EvalResult) { // Evaluate
                    int NoOfNodes = NS.NoOfNodes;
                    MultidimensionalArray BasisValues = TestBasis.Evaluate(NS.Nodes);                 // reference
                    var LSNormals = this.tracker.GetLevelSetReferenceNormals(iLevSet, NS.Nodes, i0, Length); // reference

                    for (int i = 0; i < Length; i++) { // loop over cells


                        CellBoundaryQuadRule cR = qBnd.CurrentRule;
                        int[] NodesPerEdge = cR.NumbersOfNodesPerFace;
                        var Kref = cR.RefElement;
                        int NoOfFaces = Kref.NoOfFaces;
                        int iNode = 0;

                        Debug.Assert(NoOfFaces == NodesPerEdge.Length);
                        for (int e = 0; e < NoOfFaces; e++) { // loop over the faces of the cell

                            if (NodesPerEdge[e] <= 0)
                                continue;

                            CellN0 = Kref.FaceNormals[e, 0];
                            CellN1 = Kref.FaceNormals[e, 1];

                            for (int _n = 0; _n < NodesPerEdge[e]; _n++) { // loop over nodes in one edge
                                Tang0 = -LSNormals[i, iNode, 1];
                                Tang1 = LSNormals[i, iNode, 0];
                                double Sign = Math.Sign(CellN0 * Tang0 + CellN1 * Tang1);
                                
                                for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
                                    EvalResult[i, iNode, n] = BasisValues[iNode, n] * Sign;
                                }

                                iNode++;
                            }
                        }
                        Debug.Assert(iNode == EvalResult.GetLength(1));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, jSgrd }, new int[] { N - 1, jSgrd - 1 });
                        int NoOfFaces = ResultsOfIntegration.GetLength(1);

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0 }, new int[] { i - 1, e - 1, N - 1 });
                            ResPart.Acc(1.0, ip);
                        }
                        jSgrd++;
                    }
                },
                cs: CoordinateSystem.Physical);
            qBnd.Execute();

            var ret = RHS.ResizeShallow(N, _mask.NoOfItemsLocally);
            return ret;
        }


    }


    /// <summary>
    /// HMF-surface integrals in 2D, based on Stokes integral theorem.
    /// </summary>
    public class SurfaceStokes_2D_Curvature {

        class QRF : IQuadRuleFactory<QuadRule> {
            internal SurfaceStokes_2D_Curvature m_Owner;

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
                get {
                    return m_Owner.RefElement;
                }
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");

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

                        while (jCell > Rule[h].Chunk.i0) {
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
        RefElement RefElement;

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
        /// the awesome level set tracker
        /// </summary>
        protected LevelSetTracker tracker;

        IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory;

        IQuadRuleFactory<CellBoundaryQuadRule> LevelSetBoundaryLineFactory;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tracker"></param>
        /// <param name="iLevSet"></param>
        /// <param name="_SurfaceNodesOnZeroLevset">if true, the nodes for the surface integration are 'projected' onto the zero-level-set</param>
        /// <param name="_DoCheck">
        /// if true, the accuracy of the quadrature is checked after solution of the system
        /// </param>
        /// <param name="_LevelSetBoundaryLineFactory"></param>
        /// <param name="cellBoundaryFactory"></param>
        public SurfaceStokes_2D_Curvature(LevelSetTracker tracker, int iLevSet,
            IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory,
            IQuadRuleFactory<CellBoundaryQuadRule> _LevelSetBoundaryLineFactory,
            bool _SurfaceNodesOnZeroLevset = false,
            bool _DoCheck = false) {
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.Docheck = _DoCheck;

            this.RefElement = cellBoundaryFactory.RefElement;
            if (!tracker.GridDat.Grid.RefElements.Contains(RefElement, ReferenceComparer.Instance)) {
                throw new ArgumentOutOfRangeException(
                    "simplex", "'simplex' must be a volume - reference element");
            }

            this.tracker = tracker;
            this.LevelSetIndex = iLevSet;
            this.cellBoundaryFactory = cellBoundaryFactory;
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.LevelSetBoundaryLineFactory = _LevelSetBoundaryLineFactory;

            int iKref = this.tracker.GridDat.Grid.RefElements.IndexOf(RefElement);
            this.MaxGrid = this.tracker.GridDat.Cells.GetCells4Refelement(iKref).Intersect(
                tracker._Regions.GetCutCellMask4LevSet(this.LevelSetIndex));
        }

        /// <summary>
        /// the intersection of the cut cells for Level Set <see cref="LevelSetIndex"/>
        /// and the subgrid for reference element <see cref="RefElement"/>
        /// </summary>
        private CellMask MaxGrid;


        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_SurfaceRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();


        public static Stopwatch stpwGetQuadRuleSet = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_StokesRHS = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_SolveRHS = new Stopwatch();

               

        private int OrderToInternalOrder(int _o) {
            int o = _o;
            if (o < 2)
                o = 1;
            o += 1; // for the volume rules, we loose one order of precision 

            //Console.WriteLine("Rem: even-order elim deaktiviert: Stokes-Surface-internal {0} -> {1} ", _o, o);
            if (o % 2 == 0)
                o += 1; // for some reason, the even orders fail in the one-step construction
            return o;
        }


        void GetQuadRuleSet_Internal(int ReqOrder) {
            using (new FuncTrace()) {

                int IntOrder = OrderToInternalOrder(ReqOrder);

                stpwGetQuadRuleSet.Start();

                Debug.Assert(IntOrder >= 2);

                if (this.m_SurfaceRules.ContainsKey(ReqOrder))
                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky
                if (this.m_SurfaceRules.ContainsKey(IntOrder))
                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky

                // check arguments, init
                // =====================

                CellMask _mask = this.MaxGrid;

                // subgrid on which the volume rule should be constructed
                // ======================================================
                CellBoundaryQuadratureScheme cellBndSchme = new CellBoundaryQuadratureScheme(this.cellBoundaryFactory, _mask);
                CellBoundaryQuadratureScheme cellBndLineSchme = new CellBoundaryQuadratureScheme(this.LevelSetBoundaryLineFactory, _mask);

                // set up
                // ======

                int NoOfEqTotal;
                var TestBasis = new Basis(this.tracker.GridDat, IntOrder); // we loose a order of 1 for the volume rule due to the divergence operator
                NoOfEqTotal = TestBasis.Length * this.tracker.GridDat.Grid.SpatialDimension;

                // define Nodes
                // ============
                NodeSet NodeSet = null;
                {
                    if (this.RefElement.GetType() == typeof(Square)) {

                        int K = (int)Math.Ceiling(Math.Sqrt(NoOfEqTotal * 1.75)) + 1;

                        var Nodes1D = GenericBlas.Linspace(-1, 1, K);

                        var _NodeSet = new NodeSet(this.RefElement, K * K, 2);
                        int n = 0;
                        for (int i = 0; i < K /*&& n <= NoOfEq*1.1*/; i++) {
                            for (int j = 0; j < K /*&& n <= NoOfEq*1.1*/; j++) {
                                _NodeSet[n, 0] = Nodes1D[i];
                                _NodeSet[n, 1] = Nodes1D[j];
                                n++;
                            }
                        }
                                                
                        NodeSet = _NodeSet;
                    } else {
                        for (int o = 1; o < 1000000; o++) {
                            var qr = RefElement.GetBruteForceQuadRule(o, 0);
                            if (qr.NoOfNodes >= (NoOfEqTotal * 1.1)) {
                                NodeSet = qr.Nodes;
                                break;
                            }
                        }
                    }

                    NodeSet.LockForever();
                }
                int NoOfNodes = NodeSet.GetLength(0);
                Debug.Assert(NoOfNodes >= NoOfEqTotal);


                // find RHS integrals
                // ==================


                stpwGetQuadRuleSet_StokesRHS.Start();
                MultidimensionalArray RHS_Stokes = this.StokesAnsatzRHS(TestBasis, cellBndLineSchme, _mask, IntOrder);
                stpwGetQuadRuleSet_StokesRHS.Stop();
                Debug.Assert(RHS_Stokes.Dimension == 2);
                Debug.Assert(RHS_Stokes.GetLength(0) == NoOfEqTotal);
                Debug.Assert(RHS_Stokes.GetLength(1) == _mask.NoOfItemsLocally);

      
                // construct da rule!
                // ==================
                ChunkRulePair<QuadRule>[] SurfaceRule;
                {
                    SurfaceRule = new ChunkRulePair<QuadRule>[_mask.NoOfItemsLocally];
                    var grddat = this.tracker.GridDat;
                    int D = grddat.SpatialDimension;

                    // loop over cells in subgrid...
                    int jSub = 0;
                    foreach (int jCell in _mask.ItemEnum) { // loop over cells in the mask


                        // setup System
                        // ============

                        NodeSet VolNodes = NodeSet;
                        NodeSet surfNodes;
                        if (this.SurfaceNodesOnZeroLevset)
                            surfNodes = ProjectOntoLevset(jCell, VolNodes);
                        else
                            surfNodes = VolNodes;


                        MultidimensionalArray Mtx_Stokes = null;
                        {
                            Mtx_Stokes = this.StokesAnsatzMatrix(TestBasis, surfNodes, jCell);
                            Debug.Assert(Mtx_Stokes.Dimension == 2);
                            Debug.Assert(Mtx_Stokes.GetLength(0) == NoOfEqTotal);
                            Debug.Assert(Mtx_Stokes.GetLength(1) == NoOfNodes);
                        }

                        stpwGetQuadRuleSet_SolveRHS.Start();

                        Debug.Assert(NoOfNodes >= NoOfEqTotal);
                        MultidimensionalArray _RHS = MultidimensionalArray.Create(NoOfNodes, 1); // this is also output, so it must be larger!
                        {
                            for (int i = 0; i < NoOfEqTotal; i++) {
                                _RHS[i, 0] = RHS_Stokes[i, jSub];
                            }

                        }

                        MultidimensionalArray __RHS = null, __Mtx = null;
                        if (this.Docheck) {
                            // values used for testing:
                            __RHS = MultidimensionalArray.Create(NoOfEqTotal);
                            for (int i = 0; i < NoOfEqTotal; i++)
                                __RHS[i] = _RHS[i, 0];
                            __Mtx = MultidimensionalArray.Create(Mtx_Stokes.NoOfRows, Mtx_Stokes.NoOfCols);
                            __Mtx.SetMatrix(Mtx_Stokes);
                        }

                        // solve system
                        // ============

                        //int M = _Mtx.NoOfRows;
                        //int N = _Mtx.NoOfCols;
                        //LAPACK.F77_LAPACK.DGELSY(M, N, _Mtx.Entries, _RHS.Entries, 1, 1.0e-14);
                        Mtx_Stokes.LeastSquareSolve(_RHS);


                        if (this.Docheck) {
                            // Probe:
                            MultidimensionalArray X = MultidimensionalArray.Create(NoOfNodes); // weights
                            X.ResizeShallow(NoOfNodes, 1).SetMatrix(_RHS);

                            double MaxWeight = X.Max(wi => wi.Pow2());
                            double MinWeight = X.Min(wi => wi.Pow2());
                            double wr = MaxWeight / MinWeight;
                            //Console.WriteLine("MaxWight: " + MaxWeight  + ",  Weight Ratio: " + wr);

                            __RHS.Multiply(-1.0, __Mtx, X, 1.0, "j", "jk", "k");
                            double L2_ERR = __RHS.L2Norm();
                            if (L2_ERR > 1.0e-7) {
                                Console.WriteLine("Un-precise quadrature order " + IntOrder + " rule in cell " + jCell + ": L2_ERR = " + L2_ERR);
                                //throw new ApplicationException("Quadrature rule in cell " + jCell + " seems to be not very precise: L2_ERR = " + L2_ERR);
                            }
                        }


                        stpwGetQuadRuleSet_SolveRHS.Stop();

                        // return da rule!
                        // ===============
                        var metrics = tracker.GetLevelSetNormalReferenceToPhysicalMetrics(this.LevelSetIndex, surfNodes, jCell, 1);

                        {
                            QuadRule qr_l = new QuadRule() {
                                OrderOfPrecision = IntOrder,
                                Weights = MultidimensionalArray.Create(NoOfNodes),
                                Nodes = surfNodes
                            };

                            for (int k = 0; k < NoOfNodes; k++) {
                                qr_l.Weights[k] = _RHS[k, 0] / metrics[0, k];
                            }

                            SurfaceRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);

                        }
                        jSub++;
                    }
                }


                stpwGetQuadRuleSet.Stop();

                this.m_SurfaceRules.Add(IntOrder, SurfaceRule);
                if (ReqOrder != IntOrder) {
                    this.m_SurfaceRules.Add(ReqOrder, SurfaceRule);
                }
            }
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
            for (int i = 0; i < 10; i++) {

                double radiusError = 0;


                int j = jCell;


                LevSet.Evaluate(j, 1, Nodes, LevSetValues, 0, 0.0);
                LevSet.EvaluateGradient(j, 1, Nodes, LevSetGrad);


                m_Context.TransformLocal2Global(new NodeSet(this.RefElement, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)), j, 1, x0_i_Global, 0);

                for (int nn = 0; nn < NN; nn++) {

                    double sc = 0;
                    for (int d = 0; d < D; d++) {
                        sc += LevSetGrad[0, nn, d].Pow2();
                    }


                    for (int d = 0; d < D; d++) {
                        double xd = x0_i_Global[0, nn, d] - LevSetGrad[0, nn, d] * LevSetValues[0, nn] / sc;
                        x0_ip1_Global[nn, d] = xd;
                    }

                    radiusError += Math.Abs(LevSetValues[0, nn]);

                }

                m_Context.TransformGlobal2Local(x0_ip1_Global, x0_ip1_Local, j, 1, 0);


                // next iter: x0_i <- x0_{i+1}
                x0_i_Local.Set(x0_ip1_Local);
                Nodes = (new NodeSet(this.RefElement, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)));
            }

            return Nodes;
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

            //uint lh = GridDat.NSC.LockNodeSetFamily(GridDat.NSC.CreateContainer(surfaceNodes, iKref, -1.0));

            var Phi = TestBasis.Evaluate(surfaceNodes);              // reference
            var GradPhi = TestBasis.EvaluateGradient(surfaceNodes);  // reference
            //var Phi = TestBasis.CellEval(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);                  // physical
            //var GradPhi = TestBasis.CellEvalGradient(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1, -1);  // physical

            var LevsetNormal = this.tracker.GetLevelSetReferenceNormals(iLevSet, surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);   // reference
            var Curvature = this.tracker.GetLevelSetReferenceCurvature(iLevSet, surfaceNodes, jCell, 1);   // reference
            //var LevsetNormal = this.tracker.GetLevelSetNormals(iLevSet, surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1); // physical
            //var Curvature = MultidimensionalArray.Create(1, NoOfNodes);                                  // physical
            //this.tracker.LevelSets[iLevSet].EvaluateTotalCurvature(jCell, 1, surfaceNodes, Curvature);  // physical

            var Coeffs = MultidimensionalArray.Create(D, N, NoOfNodes);

            if (D == 2) {
                for (int k = 0; k < NoOfNodes; k++) { // loop over nodes
                    double Nx = LevsetNormal[k, 0];
                    double Ny = LevsetNormal[k, 1];
                    double kappa = Curvature[0, k];


                    double Prj_11 = 1.0 - Nx * Nx, Prj_12 = -Nx * Ny,
                        Prj_21 = -Ny * Nx, Prj_22 = 1.0 - Ny * Ny;

                    for (int n = 0; n < N; n++) {
                        double Phi_kn = Phi[k, n];
                        double dPhi_dx_kn = GradPhi[k, n, 0];
                        double dPhi_dy_kn = GradPhi[k, n, 1];

                        Coeffs[0, n, k] = -Phi_kn * kappa * Nx + Prj_11 * dPhi_dx_kn + Prj_12 * dPhi_dy_kn;
                        Coeffs[1, n, k] = -Phi_kn * kappa * Ny + Prj_21 * dPhi_dx_kn + Prj_22 * dPhi_dy_kn;
                    }
                }
            } else if (D == 3) {
                throw new NotImplementedException("to do.");
            } else {
                throw new NotSupportedException("Unknown spatial dimension.");
            }

            //Coeffs.Scale(scalings[jCell]);  // physical

            return Coeffs.ResizeShallow(N * D, NoOfNodes);
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
            int iLevSet = this.LevelSetIndex;
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
                delegate (int i0, int Length, CellBoundaryQuadRule NS, MultidimensionalArray EvalResult) { // Evaluate
                    int NoOfNodes = NS.NoOfNodes;
                    MultidimensionalArray BasisValues = TestBasis.Evaluate(NS.Nodes);                 // reference
                    var LSNormals = this.tracker.GetLevelSetReferenceNormals(iLevSet, NS.Nodes, i0, Length); // reference
                    //MultidimensionalArray BasisValues = TestBasis.CellEval(NS.Nodes, i0, Length);         // physical
                    //MultidimensionalArray LSNormals = this.tracker.GetLevelSetNormals(0, NS.Nodes, i0, Length);  // physical

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
                            for (int d = 0; d < D; d++) {
                                CellN[d] = Kref.FaceNormals[e, d];
                            }
                            // ~~~~

                            //// physical:
                            //NodeSet FaceNodes = new NodeSet(this.RefElement, cR.Nodes.ExtractSubArrayShallow(new int[] { iNode, 0 }, new int[] { iNode + NodesPerEdge[e] - 1, D - 1 }));
                            //var FaceNormals = MultidimensionalArray.Create(NodesPerEdge[e], D);
                            //GridDat.Edges.GetNormalsForCell(FaceNodes, i0, e, FaceNormals);
                            //// ~~~~

                            for (int _n = 0; _n < NodesPerEdge[e]; _n++) { // loop over nodes in one edge
                                for (int d = 0; d < D; d++) {
                                    SurfN[d] = LSNormals[i, iNode, d];
                                    //CellN[d] = FaceNormals[_n, d]; // physical
                                }
                                tangente(SurfN, CellN, OutwardTang);

                                for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
                                    for (int d = 0; d < D; d++) { // loop over spatial direction
                                        //EvalResult[i, iNode, d, n] = BasisValues[i, iNode, n] * OutwardTang[d]; // physical
                                        EvalResult[i, iNode, d, n] = BasisValues[iNode, n] * OutwardTang[d]; // reference
                                    }
                                }

                                iNode++;
                            }
                        }
                        Debug.Assert(iNode == EvalResult.GetLength(1));
                    }
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, 0, jSgrd }, new int[] { D - 1, N - 1, jSgrd - 1 });
                        int NoOfFaces = ResultsOfIntegration.GetLength(1);

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0, 0 }, new int[] { i - 1, e - 1, D - 1, N - 1 });
                            ResPart.Acc(1.0, ip);
                        }
                        jSgrd++;
                    }
                },
                cs: CoordinateSystem.Physical);
            qBnd.Execute();

            var ret = RHS.ResizeShallow(N * D, _mask.NoOfItemsLocally);
            return ret;
        }


    }
}
