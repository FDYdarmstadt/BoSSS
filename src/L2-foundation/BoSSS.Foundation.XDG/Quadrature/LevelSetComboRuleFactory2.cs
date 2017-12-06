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
    /// One-step rules (creating 
    /// surface rules for integrals \f$ \oint_{\frakI \cap K_j } \ldots \dS \f$ and 
    /// volume  rules for integrals \f$ \int_{ \frakB \cap K_j } \ldots \dV \f$ 
    /// in one step), using Gauss and optionally, the Stokes theorem.
    /// Supports only 2D.
    /// </summary>
    public class LevelSetComboRuleFactory2 {

       
        class QRF : IQuadRuleFactory<QuadRule> {
            internal LevelSetComboRuleFactory2 m_Owner;

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

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");


                int InternalSurfaceOrder = m_Owner.OrderToInternalOrder(RequestedOrder);
                int InternalVolumeOrder = InternalSurfaceOrder - 1;

                int FiledOrder;
                Debug.Assert(object.ReferenceEquals(this.Rules, m_Owner.m_VolumeRules) || object.ReferenceEquals(this.Rules, m_Owner.m_SurfaceRules));
                if (object.ReferenceEquals(this.Rules, m_Owner.m_VolumeRules))
                    // I'm a volume rule factory
                    FiledOrder = InternalVolumeOrder;
                else
                    // I'm a surface rule factory
                    FiledOrder = InternalSurfaceOrder;

#if DEBUG
                if (mask.Except(m_Owner.MaxGrid).NoOfItemsLocally > 0)
                    throw new NotSupportedException("'mask' must be a subset of the cut cells, for my reference element.");
#endif
                if (!Rules.ContainsKey(FiledOrder))
                    m_Owner.GetQuadRuleSet_Internal(InternalSurfaceOrder);

                if (mask.NoOfItemsLocally == m_Owner.MaxGrid.NoOfItemsLocally) {
                    // aggressive
                    return Rules[FiledOrder];
                } else {
                    var Rule = Rules[FiledOrder];

                    int L = mask.NoOfItemsLocally, H = Rule.Length;
                    var Ret = new ChunkRulePair<QuadRule>[L];
                    int h = 0;
                    //for (int jsub = 0; jsub < L; jsub++) {
                    //    int jCell = jsub2jcell[jsub];
                    int jsub = 0;
                    foreach(int jCell in mask.ItemEnum) { 

                        Debug.Assert(Rule[h].Chunk.Len == 1);

                        while (jCell > Rule[h].Chunk.i0) {
                            h++;
                        }

                        Debug.Assert(jCell == Rule[h].Chunk.i0);
                        Ret[jsub] = Rule[h];
#if DEBUG
                        Ret[jsub].Rule.Weights.CheckForNanOrInf();
                        Ret[jsub].Rule.Nodes.CheckForNanOrInf();
#endif
                        jsub++;
                    }
                    Debug.Assert(jsub == L);

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

        public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            return new QRF() {
                m_Owner = this,
                Rules = this.m_VolumeRules
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
        /// if true, also Stoke's theorem is used to create the surface rule
        /// </summary>
        bool UseAlsoStokes;

        ///// <summary>
        ///// the awesome level set tracker
        ///// </summary>
        //protected LevelSetTracker tracker;


        IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory;
        IQuadRuleFactory<CellBoundaryQuadRule> LevelSetBoundaryLineFactory;

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="_SurfaceNodesOnZeroLevset">if true, the nodes for the surface integration are 'projected' onto the zero-level-set</param>
        /// <param name="_DoCheck">
        /// if true, the accuracy of the quadrature is checked after solution of the system
        /// </param>
        /// <param name="_LevelSetBoundaryLineFactory"></param>
        /// <param name="cellBoundaryFactory"></param>
        /// <param name="_UseAlsoStokes">if true, also Stoke's theorem is used to create the surface rule</param>
        /// <param name="lsData"></param>
        public LevelSetComboRuleFactory2(LevelSetTracker.LevelSetData lsData,
            IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory,
            IQuadRuleFactory<CellBoundaryQuadRule> _LevelSetBoundaryLineFactory,
            bool _SurfaceNodesOnZeroLevset = false,
            bool _UseAlsoStokes = true,
            bool _DoCheck = false) {
            this.UseAlsoStokes = _UseAlsoStokes;
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.Docheck = _DoCheck;
            this.LevelSetData = lsData;

            if (_UseAlsoStokes && !object.ReferenceEquals(cellBoundaryFactory.RefElement, _LevelSetBoundaryLineFactory.RefElement))
                throw new ArgumentException("boundary factory and boundary lne factory must operate on the same reference element.");
            if (!_UseAlsoStokes && _LevelSetBoundaryLineFactory != null)
                throw new ArgumentException();

            this.RefElement = cellBoundaryFactory.RefElement;
            if (!lsData.GridDat.Grid.RefElements.Contains(RefElement, ReferenceComparer.Instance)) {
                throw new ArgumentOutOfRangeException(
                    "simplex", "'simplex' must be a volume - reference element");
            }

            this.LevelSetIndex = lsData.LevelSetIndex;
            this.cellBoundaryFactory = cellBoundaryFactory;
            this.SurfaceNodesOnZeroLevset = _SurfaceNodesOnZeroLevset;
            this.LevelSetBoundaryLineFactory = _LevelSetBoundaryLineFactory;

            int iKref = lsData.GridDat.Grid.RefElements.IndexOf(RefElement);
            this.MaxGrid = lsData.GridDat.Cells.GetCells4Refelement(iKref).Intersect(
                this.LevelSetData.Region.GetCutCellMask4LevSet(this.LevelSetIndex));
        }

        LevelSetTracker.LevelSetData LevelSetData;

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


        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_VolumeRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();



        public static Stopwatch stpwGetQuadRuleSet = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_StokesRHS = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_GaussRHS = new Stopwatch();
        public static Stopwatch stpwGetQuadRuleSet_SolveRHS = new Stopwatch();


        private int OrderToInternalOrder(int _o) {
            int o = _o;
            if (o < 2)
                o = 2;
            o += 1; // for the volume rules, we loose one order of precision 

            if (o % 2 == 0)
                o += 1; // for some reason, the even orders fail in the one-step construction
            return o;
        }
        
        void GetQuadRuleSet_Internal(int IntOrder) {
            using (new FuncTrace()) {

                //int IntOrder = OrderToInternalOrder(ReqOrder);
                
                stpwGetQuadRuleSet.Start();

                if ((IntOrder <= 2) || (IntOrder % 2 == 0))
                    throw new ArgumentOutOfRangeException();
                

//                if (this.m_VolumeRules.ContainsKey(ReqOrder))
//                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky
                if (this.m_VolumeRules.ContainsKey(IntOrder - 1))
                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky
//                if (this.m_SurfaceRules.ContainsKey(ReqOrder))
//                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky
                if (this.m_SurfaceRules.ContainsKey(IntOrder))
                    throw new ApplicationException("illegal call"); // rule should be in cache, request for re-computation indicates that something is funky



                // check arguments, init
                // =====================

                CellMask _mask = this.MaxGrid;

                // subgrid on which the volume rule should be constructed
                // ======================================================
                CellBoundaryQuadratureScheme cellBndSchme = new CellBoundaryQuadratureScheme(this.cellBoundaryFactory, _mask);
                CellBoundaryQuadratureScheme cellBndLineSchme = null;
                if (this.UseAlsoStokes)
                    cellBndLineSchme = new CellBoundaryQuadratureScheme(this.LevelSetBoundaryLineFactory, _mask);

                // set up
                // ======

                int NoOfEqPerSy, NoOfEqTotal;
                var TestBasis = new Basis(this.LevelSetData.GridDat, IntOrder);
                NoOfEqPerSy = TestBasis.Length * this.LevelSetData.GridDat.Grid.SpatialDimension;
                NoOfEqTotal = this.UseAlsoStokes ? NoOfEqPerSy * 2 : NoOfEqPerSy;
                
                // define Nodes
                // ============
                NodeSet NodeSet = null;
                {
                    if (this.RefElement.GetType() == typeof(Square)) {


                        int K = (int)Math.Ceiling(Math.Sqrt(NoOfEqTotal * 1.75)) + 1;

                        var Nodes1D = GenericBlas.Linspace(-1, 1, K);

                        var _NodeSet = new NodeSet(this.RefElement, K*K, 2);
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
                Debug.Assert(NoOfNodes * 2 >= NoOfEqTotal);


                // find RHS integrals
                // ==================

                stpwGetQuadRuleSet_GaussRHS.Start();
                var RHS_Gauss = this.GaußAnsatzRHS(TestBasis, cellBndSchme, _mask, IntOrder);
                stpwGetQuadRuleSet_GaussRHS.Stop();
                Debug.Assert(RHS_Gauss.Dimension == 2);
                Debug.Assert(RHS_Gauss.GetLength(0) == NoOfEqPerSy);
                Debug.Assert(RHS_Gauss.GetLength(1) == _mask.NoOfItemsLocally);
#if DEBUG
                RHS_Gauss.CheckForNanOrInf();
#endif

                MultidimensionalArray RHS_Stokes = null;
                if (this.UseAlsoStokes) {
                    stpwGetQuadRuleSet_StokesRHS.Start();
                    RHS_Stokes = this.StokesAnsatzRHS_RefBasis(TestBasis, cellBndLineSchme, _mask, IntOrder);
                    //RHS_Stokes = this.StokesAnsatzRHS_PhysBasis(TestBasis, cellBndLineSchme, _mask, IntOrder);
                    stpwGetQuadRuleSet_StokesRHS.Stop();
                    Debug.Assert(RHS_Stokes.Dimension == 2);
                    Debug.Assert(RHS_Stokes.GetLength(0) == NoOfEqPerSy);
                    Debug.Assert(RHS_Stokes.GetLength(1) == _mask.NoOfItemsLocally);
#if DEBUG
                    RHS_Stokes.CheckForNanOrInf();
#endif
                }

                // construct da rule!
                // ==================
                ChunkRulePair<QuadRule>[] VolumeRule;
                ChunkRulePair<QuadRule>[] SurfaceRule;
                {
                    VolumeRule = new ChunkRulePair<QuadRule>[_mask.NoOfItemsLocally];
                    SurfaceRule = new ChunkRulePair<QuadRule>[_mask.NoOfItemsLocally];
                    var grddat = this.LevelSetData.GridDat;
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

                        MultidimensionalArray metrics;
                        var Mtx_Gauss = GaußAnsatzMatrix(TestBasis, VolNodes, surfNodes, jCell, out metrics);
                        Debug.Assert(Mtx_Gauss.Dimension == 2);
                        Debug.Assert(Mtx_Gauss.GetLength(0) == NoOfEqPerSy);
                        Debug.Assert(Mtx_Gauss.GetLength(1) == NoOfNodes * 2);
#if DEBUG
                        Mtx_Gauss.CheckForNanOrInf();
#endif


                        MultidimensionalArray Mtx_Stokes = null;
                        if (this.UseAlsoStokes) {
                            Mtx_Stokes = this.StokesAnsatzMatrix_RefBasis(TestBasis, surfNodes, jCell);
                            //Mtx_Stokes = this.StokesAnsatzMatrix_PhysBasis(TestBasis, surfNodes, jCell);
                            Debug.Assert(Mtx_Stokes.Dimension == 2);
                            Debug.Assert(Mtx_Stokes.GetLength(0) == NoOfEqPerSy);
                            Debug.Assert(Mtx_Stokes.GetLength(1) == NoOfNodes);

                            /*
                            // Nur für PhysBasis
                            for (int i = 0; i < NoOfEqPerSy; i++) {
                                for (int j = 0; j < NoOfNodes; j++) {
                                    Mtx_Stokes[i, j] /= metrics[0, j];
                                }
                            }
                            */
#if DEBUG
                            Mtx_Stokes.CheckForNanOrInf();
#endif
                        }

                        stpwGetQuadRuleSet_SolveRHS.Start();

                        // convert to FORTRAN order
                        MultidimensionalArray _Mtx = MultidimensionalArray.Create(NoOfEqTotal, NoOfNodes * 2);
                        _Mtx.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NoOfEqPerSy - 1, NoOfNodes * 2 - 1 }).Set(Mtx_Gauss);
                        if (this.UseAlsoStokes) {
                            _Mtx.ExtractSubArrayShallow(new int[] { NoOfEqPerSy, NoOfNodes }, new int[] { 2 * NoOfEqPerSy - 1, 2 * NoOfNodes - 1 }).Set(Mtx_Stokes);

                        }
                        
                        // convert to FORTRAN order
                        Debug.Assert(NoOfNodes * 2 >= NoOfEqTotal);
                        MultidimensionalArray RHSandSolution = MultidimensionalArray.Create(NoOfNodes * 2, 1); // this is also output, so it must be larger!
                        {
                            int I = NoOfEqPerSy, J = _mask.NoOfItemsLocally;
                            for (int i = 0; i < I; i++) {
                                RHSandSolution[i, 0] = RHS_Gauss[i, jSub];
                            }
                            if (this.UseAlsoStokes) {
                                for (int i = 0; i < I; i++) {
                                    RHSandSolution[i + NoOfEqPerSy, 0] = RHS_Stokes[i, jSub];
                                }
                            }
                        }

                        MultidimensionalArray __RHS = null, __Mtx = null;
                        if (this.Docheck) {
                            // values used for testing:
                            __RHS = MultidimensionalArray.Create(NoOfEqTotal);
                            for (int i = 0; i < NoOfEqTotal; i++)
                                __RHS[i] = RHSandSolution[i, 0];
                            __Mtx = MultidimensionalArray.Create(_Mtx.NoOfRows, _Mtx.NoOfCols);
                            __Mtx.SetMatrix(_Mtx);
                        }

                        // solve system
                        // ============

                        _Mtx.LeastSquareSolve(RHSandSolution);
#if DEBUG
                        RHSandSolution.CheckForNanOrInf();
#endif

                        

                        if (this.Docheck) {
                            // Probe:
                            MultidimensionalArray X = MultidimensionalArray.Create(NoOfNodes * 2); // weights
                            X.ResizeShallow(NoOfNodes * 2, 1).SetMatrix(RHSandSolution);

                            double MaxWeight = X.Max(wi => wi.Pow2());
                            double MinWeight = X.Min(wi => wi.Pow2());
                            double wr = MaxWeight / MinWeight;
                            
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

                        {
                            {
                                // the volume rule!
                                // ----------------

                                QuadRule qr_l = new QuadRule() {
                                    OrderOfPrecision = IntOrder - 1, 
                                    Weights = MultidimensionalArray.Create(NoOfNodes),
                                    Nodes = VolNodes
                                };

                                int Kend = VolNodes.GetLength(0);

                                for (int k = 0; k < Kend; k++) {
                                    qr_l.Weights[k] = RHSandSolution[k, 0];
                                }

                                VolumeRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);
                                
                            }

                            {
                                // the surface rule
                                // ----------------

                                QuadRule qr_l = new QuadRule() {
                                    OrderOfPrecision = IntOrder,
                                    Weights = MultidimensionalArray.Create(NoOfNodes),
                                    Nodes = surfNodes
                                };

                                int Kend = VolNodes.GetLength(0) + surfNodes.GetLength(0);
                                int Kstr = VolNodes.GetLength(0);

                                for (int k = Kstr; k < Kend; k++) {
                                    qr_l.Weights[k - Kstr] = RHSandSolution[k, 0] / metrics[0, k - Kstr];
                                }

                                SurfaceRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), qr_l);
                            }
                        }
                        jSub++;
                    }
                }


                stpwGetQuadRuleSet.Stop();

                // cache the rules
                // ================
                this.m_VolumeRules.Add(IntOrder - 1, VolumeRule);
                this.m_SurfaceRules.Add(IntOrder, SurfaceRule);
            }
        }



        protected MultidimensionalArray GaußAnsatzMatrix(Basis TestBasis, NodeSet VolQrNodes, NodeSet SurfQrNodes, int jCell, out MultidimensionalArray metrics) {
            int N = TestBasis.Length;
            int D = this.LevelSetData.GridDat.Grid.SpatialDimension;
            int Kvol = VolQrNodes.GetLength(0);
            int ksurf = SurfQrNodes.GetLength(0);
            Debug.Assert(VolQrNodes.Dimension == 2);
            Debug.Assert(VolQrNodes.GetLength(1) == D);
            Debug.Assert(SurfQrNodes.Dimension == 2);
            Debug.Assert(SurfQrNodes.GetLength(1) == D);
            int iKref = this.LevelSetData.GridDat.Cells.GetRefElementIndex(jCell);

            MultidimensionalArray Matrix = MultidimensionalArray.Create(N, D, Kvol + ksurf);

            var VolMatrix = Matrix.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { N - 1, D - 1, Kvol - 1 });
            var EdgMatrix = Matrix.ExtractSubArrayShallow(new int[] { 0, 0, Kvol }, new int[] { N - 1, D - 1, ksurf + Kvol - 1 });


            // evaluate Basis and Gradient of Basis
            //uint lh = NSC.LockNodeSetFamily(NSC.CreateContainer(VolQrNodes, 0, 0.0), NSC.CreateContainer(SurfQrNodes, iKref, -1));
            var GradPhi = TestBasis.EvaluateGradient(VolQrNodes); // test function, n
            var Phi = TestBasis.Evaluate(SurfQrNodes); // test function, n
            var LevelSetNormals = this.LevelSetData.GetLevelSetReferenceNormals(SurfQrNodes, jCell, 1).ExtractSubArrayShallow(0,-1,-1);

            metrics = this.LevelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(SurfQrNodes, jCell, 1);
            
            // multiply
            for (int d = 0; d < D; d++) // loop over spatial dimension...
                for (int k = 0; k < Kvol; k++) // loop over nodes...
                    for (int n = 0; n < N; n++) // loop over basis polynomials...
                        VolMatrix[n, d, k] = GradPhi[k, n, d];

            // multiply
            //for (int d = 0; d < D; d++) // loop over spatial dimension...
            //    for (int k = 0; k < ksurf; k++) // loop over nodes...
            //        for (int n = 0; n < N; n++) // loop over basis polynomials...
            //            EdgMatrix[n, d, k] = Phi[k, n]*LevelSetNormals[0, k, d];
            EdgMatrix.Multiply(1.0, Phi, LevelSetNormals, 0.0, "ndk", "kn", "kd");


            // resize and return:
            return Matrix.ResizeShallow(N * D, Kvol + ksurf);
        }

        
        NodeSet ProjectOntoLevset(int jCell, NodeSet Nodes) {

            int D = Nodes.GetLength(1);
            int NoOfNodes = Nodes.GetLength(0);
            var m_Context = this.LevelSetData.GridDat;
            //LevelSet LevSet = (LevelSet)(this.tracker.LevelSets[this.LevelSetIndex]);

            MultidimensionalArray LevSetValues;// = MultidimensionalArray.Create(1, NoOfNodes);
            MultidimensionalArray LevSetGrad;// = MultidimensionalArray.Create(1, NoOfNodes, D);

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

                
                LevSetValues = this.LevelSetData.GetLevSetValues(Nodes, j, 1);
                LevSetGrad = this.LevelSetData.GetLevelSetGradients(Nodes, j, 1);
                
                
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

        protected MultidimensionalArray GaußAnsatzRHS(Basis TestBasis, CellBoundaryQuadratureScheme cellBndScheme, CellMask _mask, int order) {
            var _Context = this.LevelSetData.GridDat;
            int N = TestBasis.Length;
            int D = this.LevelSetData.GridDat.Grid.SpatialDimension;
            var coordSys = CoordinateSystem.Reference;
            var b = TestBasis;
            int Nrhs = _mask.NoOfItemsLocally;

            MultidimensionalArray RHS = MultidimensionalArray.Create(N, D, Nrhs);

            var splx = this.RefElement;
            int NoOfFaces = splx.NoOfFaces;
            //var normals = _Context.GridDat.Normals;
            ICompositeQuadRule<CellBoundaryQuadRule> bndRule = cellBndScheme.Compile(_Context, order);
            CellBoundaryQuadrature<CellBoundaryQuadRule> qBnd = null;
            int jSgrd = 0;
            qBnd = CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { N, D },
                _Context, bndRule,
                delegate(int i0, int Length, CellBoundaryQuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    int NoOfNodes = QR.NoOfNodes;
                    MultidimensionalArray BasisValues;
                    if (coordSys == CoordinateSystem.Physical) {
                        BasisValues = b.CellEval(QR.Nodes, i0, Length);
                    } else if (coordSys == CoordinateSystem.Reference) {
                        BasisValues = b.Evaluate(QR.Nodes);
                    } else
                        throw new NotImplementedException();

                    for (int i = 0; i < Length; i++) { // loop over cell

                        CellBoundaryQuadRule cR = qBnd.CurrentRule;
                        int[] NodesPerEdge = cR.NumbersOfNodesPerFace;
                        Debug.Assert(object.ReferenceEquals(splx, cR.RefElement));

                        int iNode = 0;

                        Debug.Assert(NoOfFaces == NodesPerEdge.Length);
                        for (int e = 0; e < NoOfFaces; e++) { // loop over the faces of the cell
                            for (int _n = 0; _n < NodesPerEdge[e]; _n++) { // loop over nodes in one edge
                                for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
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
                                            EvalResult[i, iNode, n, d] = BasisValues[iNode, n] * Nd;
                                        }
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
                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, 0, jSgrd }, new int[] { N - 1, D - 1, jSgrd - 1 });

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0, 0 }, new int[] { i - 1, e - 1, N - 1, D - 1 });
                            ResPart.Acc(1.0, ip);
                        }
                        jSgrd++;
                    }
                },
                cs: coordSys);
            qBnd.Execute();

            var ret = RHS.ResizeShallow(N * D, Nrhs);
            return ret;
        }

        /// <summary>
        /// Matrix (LHS) for the Stokes/curvature Ansatz, in the _reference_ coordinate system.
        /// </summary>
        MultidimensionalArray StokesAnsatzMatrix_RefBasis(Basis TestBasis, NodeSet surfaceNodes, int jCell) {
            int N = TestBasis.Length;
            int NoOfNodes = surfaceNodes.GetLength(0);
            int D = surfaceNodes.GetLength(1);
            var GridDat = this.LevelSetData.GridDat;
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

            var LevsetNormal = this.LevelSetData.GetLevelSetReferenceNormals(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1); ;  // reference
            var Curvature = this.LevelSetData.GetLevelSetReferenceCurvature(surfaceNodes, jCell, 1);   // reference
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
                        Debug.Assert(!(double.IsNaN(Coeffs[0, n, k]) || double.IsInfinity(Coeffs[0, n, k])));
                        Debug.Assert(!(double.IsNaN(Coeffs[1, n, k]) || double.IsInfinity(Coeffs[1, n, k])));
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

        /// <summary>
        /// Matrix (LHS) for the Stokes/curvature Ansatz, in the _physical_ coordinate system.
        /// </summary>
        MultidimensionalArray StokesAnsatzMatrix_PhysBasis(Basis TestBasis, NodeSet surfaceNodes, int jCell) {
            int N = TestBasis.Length;
            int NoOfNodes = surfaceNodes.GetLength(0);
            int D = surfaceNodes.GetLength(1);
            var GridDat = this.LevelSetData.GridDat;
            Debug.Assert(D == GridDat.SpatialDimension);
            int iKref = GridDat.Cells.GetRefElementIndex(jCell);
            var scalings = GridDat.Cells.JacobiDet;
            int iLevSet = this.LevelSetIndex;

            if (!GridDat.Cells.IsCellAffineLinear(jCell))
                throw new NotSupportedException();

            //uint lh = GridDat.NSC.LockNodeSetFamily(GridDat.NSC.CreateContainer(surfaceNodes, iKref, -1.0));

            //var Phi = TestBasis.Evaluate(0);              // reference
            //var GradPhi = TestBasis.EvaluateGradient(0);  // reference
            var Phi = TestBasis.CellEval(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);                  // physical
            var GradPhi = TestBasis.CellEvalGradient(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1, -1);  // physical

            //var LevsetNormal = this.LsTrk.GetLevelSetReferenceNormals(iLevSet, 0, jCell, 1);  // reference
            //var Curvature = this.LsTrk.GetLevelSetReferenceCurvature(iLevSet, 0, jCell, 1);   // reference
            var LevsetNormal = this.LevelSetData.GetLevelSetNormals(surfaceNodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1); // physical
            var Curvature = MultidimensionalArray.Create(1, NoOfNodes);                                  // physical
            this.LevelSetData.LevelSet.EvaluateTotalCurvature(jCell, 1, surfaceNodes, Curvature);  // physical

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

            Coeffs.Scale(scalings[jCell]);  // physical

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

        /// <summary>
        /// RHS for the Stokes/curvature Ansatz, in the _reference_ coordinate system.
        /// </summary>
        MultidimensionalArray StokesAnsatzRHS_RefBasis(Basis TestBasis, CellBoundaryQuadratureScheme cellBndSchme, CellMask _mask, int order) {
            var GridDat = this.LevelSetData.GridDat;
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
                delegate(int i0, int Length, CellBoundaryQuadRule NS, MultidimensionalArray EvalResult) { // Evaluate
                    int NoOfNodes = NS.NoOfNodes;
                    MultidimensionalArray BasisValues = TestBasis.Evaluate(NS.Nodes);                 // reference
                    var LSNormals = this.LevelSetData.GetLevelSetReferenceNormals(NS.Nodes, i0, Length); // reference
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
                                        EvalResult[i, iNode, d, n] = BasisValues[iNode, n]*OutwardTang[d]; // reference
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


        /// <summary>
        /// RHS for the Stokes/curvature Ansatz, in the _physical_ coordinate system.
        /// </summary>
        MultidimensionalArray StokesAnsatzRHS_PhysBasis(Basis TestBasis, CellBoundaryQuadratureScheme cellBndSchme, CellMask _mask, int order) {
            var GridDat = this.LevelSetData.GridDat;
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
                delegate (int i0, int Length, CellBoundaryQuadRule NS, MultidimensionalArray EvalResult) { // Evaluate
                    int NoOfNodes = NS.NoOfNodes;
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
                            NodeSet FaceNodes = new NodeSet(this.RefElement, cR.Nodes.ExtractSubArrayShallow(new int[] { iNode, 0 }, new int[] { iNode + NodesPerEdge[e] - 1, D - 1 }));
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
                                        EvalResult[i, iNode, d, n] = BasisValues[i, iNode, n] * OutwardTang[d]; // physical
                                        //EvalResult[i, iNode, d, n] = BasisValues[iNode, n]*OutwardTang[d]; // reference
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

