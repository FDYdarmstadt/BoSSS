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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// This factory produces quadrature rules which are,
    /// for each cell \f$ K\f$  in a volume mask,
    /// capable of computing (an approximation of)
    /// \f[ 
    ///    \int\limits_{\{ \vec{x}; \varphi(\vec{x}) {\leq \atop \geq} 0 \} \cap K}  f \ d \vec{x},
    /// \f]
    /// where \f$ \varphi\f$  denotes the level set function.
    /// </summary>
    abstract public class LevelSetVolumeQuadRuleFactory2 : IQuadRuleFactory<QuadRule> {


        /// <summary>
        /// grid reference element.
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        private JumpTypes jumpType;


        private IQuadRuleFactory<QuadRule> edgeRuleFactory;

        IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory;

        private IQuadRuleFactory<QuadRule> surfaceRuleFactory;


        /// <summary>
        /// constructor.
        /// </summary>
        /// <param name="edgeRuleFactory">
        /// Some factory that provides quadrature rulz for the integration 
        /// over 
        /// \f[ 
        ///  \partial K \cap \{ \vec{x}; \varphi(\vec{x}) {\leq \atop \geq} 0 \}.
        /// \f]
        /// (Here, \f$ \partial K\f$   the boundary of some cell 
        /// \f$ K\f$  and 
        /// \f$ \varphi\f$  denotes the level set function.)
        /// </param>
        /// <param name="surfaceRuleFactory">
        /// Some factory that provides quadrature rulz for the integration 
        /// over the zero level set, i.e.
        /// \f[ 
        ///   \{ \vec{x}; \varphi(\vec{x}) = 0 \} \cap K.
        ///\f]
        /// (Here, \f$ \partial K\f$   the boundary of some cell 
        /// \f$ K\f$  and 
        /// \f$ \varphi\f$  denotes the level set function.)
        /// </param>
        /// <param name="simplex"></param>
        /// <param name="jumpType"></param>
        /// <param name="LevSetData"></param>
        /// <param name="cellBoundaryFactory"></param>
        internal LevelSetVolumeQuadRuleFactory2(RefElement simplex, LevelSetTracker.LevelSetData LevSetData, IQuadRuleFactory<QuadRule> edgeRuleFactory, IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory, IQuadRuleFactory<QuadRule> surfaceRuleFactory, JumpTypes jumpType) {
            if (jumpType != JumpTypes.Heaviside) {
                throw new NotSupportedException();
            }
            if ((cellBoundaryFactory != null) == (edgeRuleFactory != null))
                throw new ArgumentException("only one of them can be used at once.");

            if (edgeRuleFactory != null && (LevSetData.GridDat.Edges.EdgeRefElements.Contains(edgeRuleFactory.RefElement, (a, b) => object.ReferenceEquals(a, b))))
                throw new ArgumentException("Illegal edge rule.");
            if (cellBoundaryFactory != null && (!object.ReferenceEquals(simplex, cellBoundaryFactory.RefElement)))
                throw new ArgumentException("Illegal cell boundary rule.");
            if (!object.ReferenceEquals(surfaceRuleFactory.RefElement, simplex))
                throw new ArgumentException("Illegal level-set surface rule.");


            this.RefElement = simplex;
            this.edgeRuleFactory = edgeRuleFactory;
            this.surfaceRuleFactory = surfaceRuleFactory;
            this.jumpType = jumpType;
            this.cellBoundaryFactory = cellBoundaryFactory;
            this.LevelSetData = LevSetData;
        }

        /// <summary>
        /// Evaluation of level-set fields
        /// </summary>
        protected LevelSetTracker.LevelSetData LevelSetData {
            get;
            private set;
        }


        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            order = Math.Max(order, 2); // Ansatz won't work below 2!



            // check arguments, init
            // =====================

            if (!(mask is CellMask)) 
                throw new ArgumentException("Must be a cell mask.", "mask");
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");

            CellMask _mask = mask as CellMask;

#if DEBUG
            if (mask.Except(LevelSetData.Region.GetCutCellSubgrid4LevSet(LevelSetData.LevelSetIndex).VolumeMask).NoOfItemsLocally > 0)
                throw new NotSupportedException("'mask' must be a subset of the cut cells.");
#endif

            // subgrid on which the volume rule should be constructed
            // ======================================================
            SubGrid sgrd = new SubGrid(_mask);
            CellQuadratureScheme surfaceScheme = new CellQuadratureScheme(surfaceRuleFactory, sgrd.VolumeMask);
            EdgeQuadratureScheme edgeScheme = null;
            CellBoundaryQuadratureScheme cellBndSchme = null;
            if (edgeRuleFactory != null)
                edgeScheme = new EdgeQuadratureScheme(edgeRuleFactory, sgrd.AllEdgesMask);
            if (cellBoundaryFactory != null)
                cellBndSchme = new CellBoundaryQuadratureScheme(cellBoundaryFactory, sgrd.VolumeMask);
            if ((edgeScheme == null) == (cellBndSchme == null))
                throw new Exception("internal error");

            // set up
            // ======

            int NoOfEq;
            this.SetUp(order, out NoOfEq);

            // find RHS integrals
            // ==================

            var RHS = GetRHS(edgeScheme, cellBndSchme, surfaceScheme, sgrd, order);
            if (RHS.Dimension != 2)
                throw new ApplicationException();
            if (RHS.GetLength(0) != NoOfEq)
                throw new ApplicationException();
            if (RHS.GetLength(1) != sgrd.LocalNoOfCells)
                throw new ApplicationException();

            // construct da rule!
            // ==================

            var rule = ConstructRule_AllTheSameNodes(order, sgrd, NoOfEq, RHS);
            //var rule = ConstructRule_Experimental(order, sgrd, NoOfEq, RHS);

            // return
            // ======
            return rule;
        }

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        /// <summary>
        /// uses the same quadrature nodes in all cells
        /// </summary>
        private ChunkRulePair<QuadRule>[] ConstructRule_AllTheSameNodes(int order, SubGrid sgrd, int NoOfEq, MultidimensionalArray RHS) {
            // define Nodes
            // ============
            NodeSet NodeSet = null;
            for (int o = 0; o < 1000000; o++) {
                var qr = RefElement.GetQuadratureRule(o);
                if (qr.NoOfNodes >= (NoOfEq)) {
                    NodeSet = qr.Nodes;
                    break;
                }
            }
            int NoOfNodes = NodeSet.GetLength(0);

            // setup System
            // ============

            var Mtx = GetMatrix(NodeSet, int.MinValue);
            if (Mtx.Dimension != 2)
                throw new ApplicationException();
            if (Mtx.GetLength(0) != NoOfEq)
                throw new ApplicationException();
            if (Mtx.GetLength(1) != NoOfNodes)
                throw new ApplicationException();

            //// convert to FORTRAN order
            //FullMatrix _Mtx = new FullMatrix(Mtx);
            ////_Mtx.ToTxtFile("C:\\tmptest\\Msd.txt");




            ////(new FullMatrix(RHS)).ToTxtFile("C:\\tmptest\\RHSsd.txt");

            //// convert to FORTRAN order
            var _RHS = MultidimensionalArray.Create(NoOfNodes, sgrd.LocalNoOfCells); // this is also output, so it must be larger!
            {
                int I = NoOfEq, J = sgrd.LocalNoOfCells;
                //for(int i = 0; i < I; i++)
                //    for(int j = 0; j < J; j++)
                //        _RHS[i, j] = RHS[i, j];
                _RHS.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { I - 1, J - 1 }).Set(RHS);

            }


            // solve system
            // ============

            int M = Mtx.NoOfRows;
            int N = Mtx.NoOfCols;
            //LAPACK.F77_LAPACK.DGELSY(M, N, _Mtx.Entries, _RHS.Entries, sgrd.LocalNoOfCells, 1.0e-14);
            Mtx.LeastSquareSolve(_RHS); // solve for multiple rhs!

            #if DEBUG
            {
                // Probe:
                MultidimensionalArray X = MultidimensionalArray.Create(NoOfNodes, sgrd.LocalNoOfCells);
                X.AccMatrix(1.0, _RHS);

                MultidimensionalArray Residual = RHS.CloneAs();
                Residual.Multiply(-1.0, Mtx, X, 1.0, "ik", "ij", "jk");
                //double L2_ERR = 0;
                //RHS.ApplyAll(delegate(double x) {
                //    L2_ERR += x * x;
                //    return x;
                //});
                //if (L2_ERR >= 1.0e-8)
                //    throw new ArithmeticException("Volume rule seems to be not very precise: L2_ERR = " + L2_ERR);
                System.IO.StringWriter errlog = null;
                System.Collections.BitArray errBitMask = null;
                for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                    double RHSnorm = 0;
                    double ResiNorm = 0;
                    for (int i = 0; i < NoOfEq; i++) {
                        ResiNorm += Residual[i, j].Pow2();
                        RHSnorm += RHS[i, j].Pow2();
                    }

                    if (ResiNorm / RHSnorm > 1.0e-8) {
                        int jCell = sgrd.SubgridIndex2LocalCellIndex[j];
                        if (errlog == null) {
                            errlog = new System.IO.StringWriter();
                            errBitMask = new System.Collections.BitArray(sgrd.GridData.iLogicalCells.NoOfLocalUpdatedCells);
                        }
                        errBitMask[jCell] = true;
                        errlog.WriteLine("High residual in cell {0} {1}, Residual norm = {2}, RHS norm = {3}", jCell, j, RHSnorm, ResiNorm);
                    }
                }
                
                if (errlog != null) {
                    errlog.Flush();
                    string msg = errlog.ToString();
                    errlog.Dispose();

                    Console.Error.WriteLine(errlog);

                    CellMask errMask = new CellMask(sgrd.GridData, errBitMask);
                    errMask.SaveToTextFile("error.csv", WriteHeader:false);
                    throw new ArithmeticException(msg);
                }
            }
            #endif
            
            // return da rule!
            // ===============


            var rule = new ChunkRulePair<QuadRule>[sgrd.LocalNoOfCells];
            int L = rule.Length;
            int[] Sgrd2Local = sgrd.SubgridIndex2LocalCellIndex;
            for (int l = 0; l < L; l++) {
                QuadRule qr_l = new QuadRule();
                qr_l.Nodes = NodeSet;
                qr_l.OrderOfPrecision = order;
                qr_l.Weights = MultidimensionalArray.Create(NoOfNodes);
                for (int k = 0; k < NoOfNodes; k++)
                    qr_l.Weights[k] = _RHS[k, l];

                Chunk c;
                c.i0 = Sgrd2Local[l];
                c.Len = 1;

                rule[l] = new ChunkRulePair<QuadRule>(c, qr_l);
            }
            return rule;
        }

        /*

        /// <summary>
        /// use only quadrature nodes within the species domain
        /// </summary>
        private ChunkRulePair<QuadRule>[] ConstructRule_Experimental(int order, SubGrid sgrd, int NoOfEq, MultidimensionalArray RHS) {

            var rule = new ChunkRulePair<QuadRule>[sgrd.LocalNoOfCells];
            int[] Sgrd2Local = sgrd.SubgridIndex2LocalCellIndex;
            var grddat = this.tracker.Ctx.GridDat;
            int D = grddat.SpatialDimension;
            
            // loop over cells in subgrid...
            int jSub = 0;
            foreach (var chk in sgrd.VolumeMask) {
                for (int jCell = chk.i0; jCell < chk.JE; jCell++) {

                    // define Nodes
                    // ============
                    MultidimensionalArray NodeSet = null;
                    MultidimensionalArray GuessWeight = null;
                    for (int o = 1; o < 1000000; o++) {
                        var qr = Simplex.GetBruteForceQuadRule(o,0);

                        List<int> goodNodes = new List<int>();

                        var GlobalNodes = MultidimensionalArray.Create(1, qr.Nodes.GetLength(0), qr.Nodes.GetLength(1));
                        grddat.TransformLocal2Global(qr.Nodes, GlobalNodes, jCell, 1, 0);

                        for (int k = 0; k < qr.NoOfNodes; k++) {
                            double Radius = 0;
                            for (int d = 0; d < D; d++)
                                Radius += GlobalNodes[0, k, d].Pow2();

                            Radius = Math.Sqrt(Radius);

                            if (Radius < 0.85)
                                // node [k] is in the "B"-region of the specific example that i am testing.
                                goodNodes.Add(k);
                        }

                       

                        if (goodNodes.Count >= (NoOfEq*10)) {
                            NodeSet = MultidimensionalArray.Create(goodNodes.Count, D);
                            GuessWeight = MultidimensionalArray.Create(goodNodes.Count);
                            int kk = 0;
                            foreach( int k in goodNodes) {
                                for (int d = 0; d < D; d++) {
                                    NodeSet[kk, d] = qr.Nodes[k, d];
                                }
                                GuessWeight[kk] = qr.Weights[k];
                                kk++;
                            }
                            break;
                        }
                    }
                    int NoOfNodes = NodeSet.GetLength(0);

                    // setup System
                    // ============

                    var Mtx = GetMatrix(NodeSet);
                    if (Mtx.Dimension != 2)
                        throw new ApplicationException();
                    if (Mtx.GetLength(0) != NoOfEq)
                        throw new ApplicationException();
                    if (Mtx.GetLength(1) != NoOfNodes)
                        throw new ApplicationException();

                    // convert to FORTRAN order
                    FullMatrix _Mtx = new FullMatrix(Mtx);
                    //_Mtx.ToTxtFile("C:\\tmptest\\Msd.txt");




                    //(new FullMatrix(RHS)).ToTxtFile("C:\\tmptest\\RHSsd.txt");

                    // convert to FORTRAN order
                    //GuessWeight.Clear();
                    FullMatrix _RHS = new FullMatrix(NoOfNodes, 1); // this is also output, so it must be larger!
                    {
                        int I = NoOfEq, J = sgrd.LocalNoOfCells;
                        for (int i = 0; i < I; i++) {
                            double acc = 0;
                            for (int k = 0; k < NoOfNodes; k++)
                                acc += _Mtx[i, k]*GuessWeight[k];
                            
                            _RHS[i, 0] = -acc + RHS[i, jSub];
                        }

                    }


                    // solve system
                    // ============

                    int M = _Mtx.NoOfRows;
                    int N = _Mtx.NoOfCols;
                    LAPACK.F77_LAPACK.DGELSY(M, N, _Mtx.Entries, _RHS.Entries, 1, 1.0e-14);

                    


//#if DEBUG
//                    // Probe:
//                    MultidimensionalArray X = MultidimensionalArray.Create(NoOfNodes, 1);
//                    X.AccM(1.0, _RHS);

//                    RHS.Multiply(-1.0, Mtx, X, 1.0, "ik", "ij", "jk");
//                    double L2_ERR = 0;
//                    RHS.DoAll(delegate(double x) { L2_ERR += x*x; return x; });
//                    Debug.Assert(L2_ERR < 1.0e-8, "Volume rule seems to be not very precise: L2_ERR = " + L2_ERR);

//                    //Console.WriteLine("Volume rule precision: L2_ERR = " + L2_ERR.ToString("0.###E-00") + ", order = " + order);

//#endif
                    // return da rule!
                    // ===============
                    
                    {
                        QuadRule qr_l = new QuadRule();
                        qr_l.Nodes = NodeSet;
                        qr_l.OrderOfPrecision = order;
                        qr_l.Weights = MultidimensionalArray.Create(NoOfNodes);

                        double mini = double.MaxValue;
                        double maxi = double.MinValue;

                        for (int k = 0; k < NoOfNodes; k++) {
                            qr_l.Weights[k] = GuessWeight[k] + _RHS[k, 0];

                            mini = Math.Min(mini, qr_l.Weights[k]);
                            maxi = Math.Max(maxi, qr_l.Weights[k]);
                        }

                        Chunk c;
                        c.i0 = Sgrd2Local[jSub];
                        c.Len = 1;

                        rule[jSub] = new ChunkRulePair<QuadRule>(c, qr_l);
                    }
                    jSub++;
                }
            }
            return rule;
        }

        */

        abstract protected void SetUp(int order, out int NoOfEqs);


        abstract protected MultidimensionalArray GetMatrix(NodeSet Nodes, int dummy);


        abstract protected MultidimensionalArray GetRHS(EdgeQuadratureScheme edgeScheme, CellBoundaryQuadratureScheme cellBndScheme, CellQuadratureScheme surfScheme, SubGrid sgrd, int order);

    }


    public class LevelSetVolumeQuadRuleFactory2a : LevelSetVolumeQuadRuleFactory2 {


        public LevelSetVolumeQuadRuleFactory2a(RefElement simplex, LevelSetTracker.LevelSetData lsData, IQuadRuleFactory<QuadRule> edgeRuleFactory, IQuadRuleFactory<QuadRule> surfaceRuleFactory, JumpTypes jumpType) :
            base(simplex, lsData, edgeRuleFactory, null, surfaceRuleFactory, jumpType) {
        }

        public LevelSetVolumeQuadRuleFactory2a(RefElement simplex, LevelSetTracker.LevelSetData lsData, IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory, IQuadRuleFactory<QuadRule> surfaceRuleFactory, JumpTypes jumpType) :
            base(simplex, lsData, null, cellBoundaryFactory, surfaceRuleFactory, jumpType) {
        }


        Basis TestBasis;

        protected override void SetUp(int order, out int NoOfEqs) {
            var Gdat = base.LevelSetData.GridDat;
            TestBasis = new Basis(base.LevelSetData.GridDat, order / 2);
            NoOfEqs = TestBasis.Length * TestBasis.Length * Gdat.SpatialDimension;
        }

        protected override MultidimensionalArray GetMatrix(NodeSet Nodes, int dummy) {
            int N = TestBasis.Length;
            int D = base.LevelSetData.GridDat.Grid.SpatialDimension;
            int K = Nodes.GetLength(0);
            Debug.Assert(Nodes.Dimension == 2);
            Debug.Assert(Nodes.GetLength(1) == D);
            
            MultidimensionalArray Matrix = MultidimensionalArray.Create(N, N, D, K);

            // evaluate Basis and Gradient of Basis
            var GradPhi = TestBasis.EvaluateGradient(Nodes); // test function, n
            var Phi = TestBasis.Evaluate(Nodes);         // basis function, m
            
            // multiply
            Matrix.Multiply(1.0, GradPhi, Phi, 0.0, "nmdk", "kmd", "kn");
            Matrix.Multiply(1.0, Phi, GradPhi, 1.0, "nmdk", "km", "knd");


            // resize and return:
            return Matrix.ResizeShallow(N * N * D, K);
        }

        protected override MultidimensionalArray GetRHS(EdgeQuadratureScheme edgeScheme, CellBoundaryQuadratureScheme cellBndScheme, CellQuadratureScheme surfScheme, SubGrid sgrd, int order) {
            var _Context = base.LevelSetData.GridDat;
            int N = TestBasis.Length;
            int D = _Context.SpatialDimension;
            var coordSys = CoordinateSystem.Reference;
            var b = this.TestBasis;
            var LsData = base.LevelSetData;
            var AllowedCells = sgrd.VolumeMask.GetBitMask();
            int Nrhs = sgrd.LocalNoOfCells;
            int[] j2Sgrd = sgrd.LocalCellIndex2SubgridIndex;

            MultidimensionalArray RHS = MultidimensionalArray.Create(N, N, D, Nrhs);

            Debug.Assert((edgeScheme != null) ^ (cellBndScheme != null));
            if (edgeScheme != null)
                RhsCellEdgeA(edgeScheme, order, _Context, N, D, coordSys, b, AllowedCells, j2Sgrd, RHS);
            if (cellBndScheme != null)
                RhsCellEdgeB(cellBndScheme, order, _Context, N, D, coordSys, b, AllowedCells, j2Sgrd, RHS);


            // level set
            {
                double NormalSign = -1;
                //var scaling = _Context.Cells.AbsJacobiDet;

                CellQuadrature.GetQuadrature(new int[] { N, N, D },
                    _Context, surfScheme.Compile(_Context, order),
                    delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                        int NoOfNodes = QR.NoOfNodes;
                        if (coordSys == CoordinateSystem.Physical) {
                            var BasisValues = b.CellEval(QR.Nodes, i0, Length);
                            var Normals = LsData.GetLevelSetNormals(QR.Nodes, i0, Length);

                            // EvalResult[i,k,n,m,d] = BasisValues[i,k,n]*BasisValues[i,k,m]*Normals[i,k,d]
                            EvalResult.Multiply(NormalSign, BasisValues, BasisValues, Normals, 0.0, "iknmd", "ikn", "ikm", "ikd");

                        } else if (coordSys == CoordinateSystem.Reference) {
                            //var BasisValues = b.CellEval(0, i0, Length);
                            var BasisValues = b.Evaluate(QR.Nodes);
                            var Normals = LsData.GetLevelSetReferenceNormals(QR.Nodes, i0, Length).ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, D - 1 });


                            // EvalResult[i,k,n,m,d] = BasisValues[i,k,n]*BasisValues[i,k,m]*Normals[i,k,d]
                            EvalResult.Multiply(NormalSign, BasisValues, BasisValues, Normals, 0.0, "iknmd", "kn", "km", "ikd");

                            var metrics = LsData.GetLevelSetNormalReferenceToPhysicalMetrics(QR.Nodes, i0, Length);
                            for (int i = 0; i < Length; i++) {
                                for (int k = 0; k < NoOfNodes; k++) {

                                    double magic_number = metrics[i, k];


                                    var EV_ik = EvalResult.ExtractSubArrayShallow(i, k, -1, -1, -1);
                                    EV_ik.Scale(magic_number);
                                }

                            }
                        } else
                            throw new NotImplementedException();

                    },
                    delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                        for (int i = 0; i < Length; i++) {
                            int jCell = i + i0;
                            int jSgrd = j2Sgrd[jCell];

                            var A = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, 0, 0, 0 }, new int[] { i - 1, N - 1, N - 1, D - 1 });
                            var B = RHS.ExtractSubArrayShallow(new int[] { 0, 0, 0, jSgrd }, new int[] { N - 1, N - 1, D - 1, jSgrd - 1 });
                            B.Acc(1.0, A);
                        }
                    },
                    cs: coordSys).Execute();
            }


            return RHS.ResizeShallow(N * N * D, Nrhs);

        }

        /// <summary>
        /// compute RHS/ cell boundary terms/ from CELL BOUNDARY rule (<paramref name="cellBndSchme"/>)
        /// </summary>
        private void RhsCellEdgeB(CellBoundaryQuadratureScheme cellBndSchme, int order, GridData _Context, int N, int D, CoordinateSystem coordSys, Basis b, System.Collections.BitArray AllowedCells, int[] j2Sgrd, MultidimensionalArray RHS) {
            RefElement splx = this.RefElement;
            int _NoOfFaces = splx.NoOfFaces;
            var FaceNormals = splx.FaceNormals;
            var EdgeNormals = base.LevelSetData.GridDat.Edges.NormalsForAffine;
            CellBoundaryQuadrature<CellBoundaryQuadRule> qBnd = null;


            qBnd = CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { N, N, D },
                _Context, cellBndSchme.Compile(_Context, order),
                delegate(int i0, int Length, CellBoundaryQuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    var NS = QR.Nodes;
                    Debug.Assert(NS.GetNodeCoordinateSystem(_Context) == NodeCoordinateSystem.CellCoord);
                    
                    MultidimensionalArray BasisValues;
                    if (coordSys == CoordinateSystem.Physical) {
                        BasisValues = b.CellEval(NS, i0, Length);
                    } else if (coordSys == CoordinateSystem.Reference) {
                        BasisValues = b.Evaluate(NS);
                    } else
                        throw new NotImplementedException();

                    var metrix = splx.FaceTrafoGramianSqrt;
                    //var Cell2Edges = _Context.Cells.Cells2Edges;
                    //var FaceIndices = _Context.Edges.FaceIndices;

                    for (int i = 0; i < Length; i++) { // loop over cell

                        CellBoundaryQuadRule cR = qBnd.CurrentRule;
                        int[] NodesPerFace = cR.NumbersOfNodesPerFace;

                        int iNode = 0;
                        int NoOfFaces = _NoOfFaces;
                        Debug.Assert(NoOfFaces == NodesPerFace.Length);
                        //var Cell2Edges_j = Cell2Edges[i+i0];

                        Debug.Assert(NoOfFaces == NodesPerFace.Length);
                        for (int iface = 0; iface < NoOfFaces; iface++) { // loop over edges of the cell...

                            double trafo = metrix[iface];

                            //double Nsign = Math.Sign(q);
                            //if (AllowedEdges[iEdge] == false)
                            //    Debug.Assert(NodesPerEdge[e] == 0);
                            for (int _n = 0; _n < NodesPerFace[iface]; _n++) { // loop over nodes in one edge

                                for (int m = 0; m < N; m++) { // loop over Basis polynomials
                                    for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
                                        for (int d = 0; d < D; d++) { // loop over spatial directions

                                            if (coordSys == CoordinateSystem.Physical) {
                                                //if (!_Context.Edges.IsEdgeAffineLinear(iEdge))
                                                //    throw new NotImplementedException("curved-element-edge: todo.");
                                                //double Nd = EdgeNormals[iEdge, d];
                                                //EvalResult[i, iNode, n, m, d] = BasisValues[i, iNode, n]*BasisValues[i, iNode, m]*Nd*Nsign;
                                                throw new NotImplementedException("todo");
                                            } else {
                                                Debug.Assert(coordSys == CoordinateSystem.Reference);
                                                double Nd = FaceNormals[iface, d];
                                                //Debug.Assert(Nd == normals[iEdge, d]*Nsign);
                                                EvalResult[i, iNode, n, m, d] = BasisValues[iNode, n] * BasisValues[iNode, m] * Nd;
                                            }
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
                    //var A = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, u, 0, 0, 0 }, new int[] { i-1, u-1, N-1, N-1, D-1 });
                    //var B = RHS.ExtractSubArrayShallow(new int[] { 0, 0, 0, jSgrd }, new int[] { N-1, N-1, D-1, jSgrd-1 });
                    ////var B = SurfResultA1.ExtractSubArrayShallow(new int[] { jCell, 0, 0, 0 }, new int[] { jCell-1, N-1, N-1, D-1 });
                    //B.Acc(1.0, A);


                    int NoOfFaces = _NoOfFaces;
                    for (int i = 0; i < Length; i++) {
                        int jCell = i0 + i;
                        int jSgrd = j2Sgrd[jCell];

                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, 0, 0, jSgrd }, new int[] { N - 1, N - 1, D - 1, jSgrd - 1 });

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0, 0, 0 }, new int[] { i - 1, e - 1, N - 1, N - 1, D - 1 });
                            ResPart.Acc(1.0, ip);
                        }

                    }
                },
                cs: coordSys);
            qBnd.Execute();
        }

        /// <summary>
        /// compute RHS/ cell boundary terms/ from EDGE rule (<paramref name="edgeScheme"/>)
        /// </summary>
        private void RhsCellEdgeA(EdgeQuadratureScheme edgeScheme, int order, GridData _Context, int N, int D, CoordinateSystem coordSys, Basis b, System.Collections.BitArray AllowedCells, int[] j2Sgrd, MultidimensionalArray RHS) {

            var normals = _Context.Edges.NormalsForAffine;
            GridData grd = _Context;
            RefElement volSimplx = this.RefElement;
            int NoOfFaces = volSimplx.NoOfFaces;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;
            MultidimensionalArray Normal = MultidimensionalArray.Create(D);

            var edgeRule = edgeScheme.Compile(_Context, order);

            EdgeQuadrature.GetQuadrature(new int[] { 2, N, N, D },
                _Context, edgeRule,
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet edgeNodeSet = QR.Nodes;
                    Debug.Assert(edgeNodeSet.GetNodeCoordinateSystem(grd) == NodeCoordinateSystem.EdgeCoord);

                    for (int i = 0; i < Length; i++) { // loop over edges
                        int iEdge = i + i0;

                        for (int u = 0; u < 2; u++) { // loop over IN,OUT cells
                            int jCell, iTrafo, iFace;
                            iTrafo = grd.Edges.Edge2CellTrafoIndex[iEdge, u];
                            jCell = grd.Edges.CellIndices[iEdge, u];
                            iFace = grd.Edges.FaceIndices[iEdge, u];
                            if (jCell < 0)
                                break;

                            double Nsign;
                            if (u == 0) {
                                Nsign = 1.0;
                            } else {
                                Nsign = -1.0;
                            }

                            // normal
                            if (coordSys == CoordinateSystem.Physical) {
                                Normal.Set(normals.ExtractSubArrayShallow(iEdge, -1));
                            } else {
                                for (int d = 0; d < D; d++)
                                    Normal[d] = volSimplx.FaceNormals[iFace, d];
                                Nsign = 1.0;
                            }

                            // evaluate Basis:
                            MultidimensionalArray BasisValues;
                            NodeSet volNodeSet = edgeNodeSet.GetVolumeNodeSet(grd, iTrafo);
                            if (coordSys == CoordinateSystem.Physical) {
                                BasisValues = b.CellEval(volNodeSet, jCell, Length).ExtractSubArrayShallow(0, -1, -1);
                            } else {
                                BasisValues = b.Evaluate(volNodeSet);
                            }

                            // where to store the result of the current edge?
                            var ER_iu = EvalResult.ExtractSubArrayShallow(i, -1, u, -1, -1, -1);


                            // Compute
                            ER_iu.Multiply(Nsign, BasisValues, BasisValues, Normal, 0.0, "knmd", "kn", "km", "d");

                        }

                    }
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) { // loop over edges
                        int iEdge = i + i0;

                        for (int u = 0; u < 2; u++) { // loop over IN,OUT cells
                            int jCell, iFace;
                            iFace = grd.Edges.FaceIndices[iEdge, u];
                            jCell = grd.Edges.CellIndices[iEdge, u];
                            if (u == 0) {
                            } else {
                                if (jCell < 0)
                                    break;
                            }

                            if (AllowedCells[jCell] != true)
                                continue;

                            int jSgrd = j2Sgrd[jCell];


                            // put the result to where it should be
                            var A = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, u, 0, 0, 0 }, new int[] { i - 1, u - 1, N - 1, N - 1, D - 1 });
                            var B = RHS.ExtractSubArrayShallow(new int[] { 0, 0, 0, jSgrd }, new int[] { N - 1, N - 1, D - 1, jSgrd - 1 });
                            //var B = SurfResultA1.ExtractSubArrayShallow(new int[] { jCell, 0, 0, 0 }, new int[] { jCell-1, N-1, N-1, D-1 });
                            B.Acc(1.0, A);

                        }
                    }

                },
                cs: coordSys).Execute();
        }
    }

    public class LevelSetVolumeQuadRuleFactory2b : LevelSetVolumeQuadRuleFactory2 {

        public LevelSetVolumeQuadRuleFactory2b(RefElement simplex, LevelSetTracker.LevelSetData lsData, IQuadRuleFactory<QuadRule> edgeRuleFactory, IQuadRuleFactory<QuadRule> surfaceRuleFactory, JumpTypes jumpType) :
            base(simplex, lsData, edgeRuleFactory, null, surfaceRuleFactory, jumpType) {
        }

        public LevelSetVolumeQuadRuleFactory2b(RefElement simplex, LevelSetTracker.LevelSetData lsData, IQuadRuleFactory<CellBoundaryQuadRule> cellBoundaryFactory, IQuadRuleFactory<QuadRule> surfaceRuleFactory, JumpTypes jumpType) :
            base(simplex, lsData, null, cellBoundaryFactory, surfaceRuleFactory, jumpType) {
        }


        Basis TestBasis;

        protected override void SetUp(int order, out int NoOfEqs) {
            TestBasis = new Basis(base.LevelSetData.GridDat, order);
            NoOfEqs = TestBasis.Length * base.LevelSetData.GridDat.SpatialDimension;
        }

        protected override MultidimensionalArray GetMatrix(NodeSet Nodes, int jCell) {
            int N = TestBasis.Length;
            int D = base.LevelSetData.GridDat.SpatialDimension;
            int K = Nodes.GetLength(0);
            Debug.Assert(Nodes.Dimension == 2);
            Debug.Assert(Nodes.GetLength(1) == D);
            
            MultidimensionalArray Matrix = MultidimensionalArray.Create(N, D, K);

            // evaluate Basis and Gradient of Basis
            var GradPhi = TestBasis.EvaluateGradient(Nodes); // test function, n
            
            // phi                        dphi_dx            dphi_dy
            // 1                          0                  0
            // x y                        1 y                x 1
            // xx xy yy                   x y 0              0 y x
            // xxx xxy xyy yyy            xx xy yy 0         0 xx xy yy
            // xxxx xxxy xxyy xyyy yyyy   xxx xxy xyy yyy 0  0 xxx xxy xyy yyy

            // multiply
            for (int d = 0; d < D; d++) // loop over spatial dimension...
                for (int k = 0; k < K; k++) // loop over nodes...
                    for (int n = 0; n < N; n++) // loop over basis polynomials...
                        Matrix[n, d, k] = GradPhi[k, n, d];



            // resize and return:
            return Matrix.ResizeShallow(N * D, K);
        }

 
        protected override MultidimensionalArray GetRHS(EdgeQuadratureScheme edgeScheme, CellBoundaryQuadratureScheme cellBndScheme, CellQuadratureScheme surfScheme, SubGrid sgrd, int order) {
            var LsData = base.LevelSetData;
            var _Context = LsData.GridDat;
            int N = TestBasis.Length;
            int D = LsData.GridDat.SpatialDimension;
            var coordSys = CoordinateSystem.Reference;
            var b = this.TestBasis;
            var AllowedCells = sgrd.VolumeMask.GetBitMask();
            int Nrhs = sgrd.LocalNoOfCells;
            int[] j2Sgrd = sgrd.LocalCellIndex2SubgridIndex;

            MultidimensionalArray RHS = MultidimensionalArray.Create(N, D, Nrhs);

            Debug.Assert((edgeScheme != null) ^ (cellBndScheme != null));
            if (edgeScheme != null)
                RhsCellEdgeA(edgeScheme, order, _Context, N, D, coordSys, b, AllowedCells, j2Sgrd, RHS);
            if (cellBndScheme != null)
                RhsCellEdgeB(cellBndScheme, order, _Context, N, D, coordSys, b, AllowedCells, j2Sgrd, RHS);


            // level set
            {
                double NormalSign = -1;
                //var scaling = _Context.Cells.AbsJacobiDet;

                //var A1 = surfScheme.Compile(_Context.GridDat, order);
                //var A2 = this.FakeSurfaceScheme(_Context, surfScheme.Domain, order);


                CellQuadrature.GetQuadrature(new int[] { N, D },
                    _Context,
                    surfScheme.Compile(_Context, order + 1),
                    delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                        NodeSet NS = QR.Nodes;
                        int NoOfNodes = NS.NoOfNodes;
                        if (coordSys == CoordinateSystem.Physical) {
                            var BasisValues = b.CellEval(NS, i0, Length);
                            var Normals = LsData.GetLevelSetNormals(NS, i0, Length);

                            // EvalResult[i,k,n,d] = BasisValues[i,k,n]*Normals[i,k,d]
                            EvalResult.Multiply(NormalSign, BasisValues, Normals, 0.0, "iknd", "ikn", "ikd");

                        } else if (coordSys == CoordinateSystem.Reference) {
                            //var BasisValues = b.CellEval(0, i0, Length);
                            var BasisValues = b.Evaluate(NS);
                            var Normals = LsData.GetLevelSetReferenceNormals(NS, i0, Length).ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, D - 1 });


                            // EvalResult[i,k,n,m,d] = BasisValues[i,k,n]*Normals[i,k,d]
                            EvalResult.Multiply(NormalSign, BasisValues, Normals, 0.0, "iknd", "kn", "ikd");




                            var metrics = LsData.GetLevelSetNormalReferenceToPhysicalMetrics(NS, i0, Length);
                            for (int i = 0; i < Length; i++) {
                                for (int k = 0; k < NoOfNodes; k++) {

                                    double magic_number = metrics[i, k];


                                    var EV_ik = EvalResult.ExtractSubArrayShallow(i, k, -1, -1);
                                    EV_ik.Scale(magic_number);
                                }

                            }

                        } else
                            throw new NotImplementedException();

                    },
                    delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                        //SurfResultB.ExtractSubArrayShallow(new int[] { i0, 0, 0, 0 }, new int[] { i0+Length-1, N-1, N-1, D-1 })
                        // .Acc(1.0, ResultsOfIntegration);

                        for (int i = 0; i < Length; i++) {
                            int jCell = i + i0;
                            int jSgrd = j2Sgrd[jCell];

                            var A = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, 0, 0 }, new int[] { i - 1, N - 1, D - 1 });
                            var B = RHS.ExtractSubArrayShallow(new int[] { 0, 0, jSgrd }, new int[] { N - 1, D - 1, jSgrd - 1 });
                            B.Acc(1.0, A);
                        }
                    },
                    cs: coordSys).Execute();
            }

            return RHS.ResizeShallow(N * D, Nrhs);
        }


        /// <summary>
        /// compute RHS/ cell boundary terms/ from EDGE rule (<paramref name="edgeScheme"/>)
        /// </summary>
        private void RhsCellEdgeA(EdgeQuadratureScheme edgeScheme, int order, GridData _Context, int N, int D, CoordinateSystem coordSys, Basis b, System.Collections.BitArray AllowedCells, int[] j2Sgrd, MultidimensionalArray RHS) {

            var normals = _Context.Edges.NormalsForAffine;
            GridData grd = _Context;
            RefElement volSimplx = base.RefElement;
            int NoOfFaces = volSimplx.NoOfFaces;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;
            MultidimensionalArray Normal = MultidimensionalArray.Create(D);

            var edgeRule = edgeScheme.Compile(_Context, order + 1);

            EdgeQuadrature.GetQuadrature(new int[] { 2, N, D },
                _Context, edgeRule,
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    int NoOfNodes = NS.NoOfNodes;
                    Debug.Assert(NS.GetNodeCoordinateSystem(grd) == NodeCoordinateSystem.EdgeCoord);

                    for (int i = 0; i < Length; i++) { // loop over edges
                        int iEdge = i + i0;

                        for (int u = 0; u < 2; u++) { // loop over IN,OUT cells
                            int jCell, iTrafo, iNodeSet, iFace;
                            iTrafo = grd.Edges.Edge2CellTrafoIndex[iEdge, u];
                            jCell = grd.Edges.CellIndices[iEdge, u];
                            iFace = grd.Edges.FaceIndices[iEdge, u];
                            iNodeSet = iTrafo;
                            if (jCell < 0)
                                break;

                            double Nsign;
                            if (u == 0) {
                                Nsign = 1.0;
                            } else {
                                Nsign = -1.0;
                            }

                            // normal
                            if (coordSys == CoordinateSystem.Physical) {
                                if (!grd.Edges.IsEdgeAffineLinear(iEdge))
                                    throw new NotImplementedException("Nonlinear elements: todo.");
                                Normal.Set(normals.ExtractSubArrayShallow(iEdge, -1));
                            } else {
                                for (int d = 0; d < D; d++)
                                    Normal[d] = volSimplx.FaceNormals[iFace, d];
                                Nsign = 1.0;
                            }

                            // evaluate Basis:
                            MultidimensionalArray BasisValues;
                            NodeSet volNodeSet = NS.GetVolumeNodeSet(grd,iNodeSet);
                            if (coordSys == CoordinateSystem.Physical) {
                                BasisValues = b.CellEval(volNodeSet, jCell, Length).ExtractSubArrayShallow(0, -1, -1);
                            } else {
                                BasisValues = b.Evaluate(volNodeSet);
                            }

                            // where to store the result of the current edge?
                            var ER_iu = EvalResult.ExtractSubArrayShallow(i, -1, u, -1, -1);


                            // Compute
                            ER_iu.Multiply(Nsign, BasisValues, Normal, 0.0, "knd", "kn", "d");

                        }

                    }
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) { // loop over edges
                        int iEdge = i + i0;

                        for (int u = 0; u < 2; u++) { // loop over IN,OUT cells
                            int jCell;
                            jCell = grd.Edges.CellIndices[iEdge, u];
                            if (jCell < 0)
                                break;

                            if (AllowedCells[jCell] != true)
                                continue;

                            int jSgrd = j2Sgrd[jCell];


                            // put the result to where it should be
                            var A = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, u, 0, 0 }, new int[] { i - 1, u - 1, N - 1, D - 1 });
                            var B = RHS.ExtractSubArrayShallow(new int[] { 0, 0, jSgrd }, new int[] { N - 1, D - 1, jSgrd - 1 });
                            //var B = SurfResultA1.ExtractSubArrayShallow(new int[] { jCell, 0, 0, 0 }, new int[] { jCell-1, N-1, N-1, D-1 });
                            B.Acc(1.0, A);

                        }
                    }

                },
                cs: coordSys).Execute();
        }

        /// <summary>
        /// compute RHS/ cell boundary terms/ from CELL BOUNDARY rule (<paramref name="cellBndSchme"/>)
        /// </summary>
        private void RhsCellEdgeB(CellBoundaryQuadratureScheme cellBndSchme, int order, GridData _Context, int N, int D, CoordinateSystem coordSys, Basis b, System.Collections.BitArray AllowedCells, int[] j2Sgrd, MultidimensionalArray RHS) {
            RefElement splx = base.RefElement;
            int _NoOfFaces = splx.NoOfFaces;
            var normals = _Context.Edges.NormalsForAffine;
            CellBoundaryQuadrature<CellBoundaryQuadRule> qBnd = null;
            var metrix = splx.FaceTrafoGramianSqrt;

            qBnd = CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(new int[] { N, D },
                _Context, cellBndSchme.Compile(_Context, order + 1),
                delegate(int i0, int Length, CellBoundaryQuadRule NS, MultidimensionalArray EvalResult) { // Evaluate
                    int NoOfNodes = NS.NoOfNodes;
                    MultidimensionalArray BasisValues;
                    if (coordSys == CoordinateSystem.Physical) {
                        BasisValues = b.CellEval(NS.Nodes, i0, Length);
                    } else if (coordSys == CoordinateSystem.Reference) {
                        BasisValues = b.Evaluate(NS.Nodes);
                    } else
                        throw new NotImplementedException();

                    //var Cell2Edges = _Context.Cells.Cells2Edges;
                    //var FaceIndices = _Context.Edges.FaceIndices;
                    int NoOfFaces = _NoOfFaces;

                    for (int i = 0; i < Length; i++) { // loop over cell

                        CellBoundaryQuadRule cR = qBnd.CurrentRule;
                        int[] NodesPerFace = cR.NumbersOfNodesPerFace;
                        Debug.Assert(NoOfFaces == NodesPerFace.Length);
                        int iNode = 0;



                        Debug.Assert(NoOfFaces == NodesPerFace.Length);
                        for (int iFace = 0; iFace < NoOfFaces; iFace++) { // loop over faces of the cell...
                            //int q = Cell2Edges_j[_e];
                            //int iEdge = Math.Abs(q) - 1;
                            //int e = FaceIndices[iEdge, q >= 0 ? 0 : 1];
                            //e = _e;

                            double trafo = metrix[iFace];

                            //double Nsign = Math.Sign(q);
                            //if (AllowedEdges[iEdge] == false)
                            //    Debug.Assert(NodesPerEdge[e] == 0);
                            for (int _n = 0; _n < NodesPerFace[iFace]; _n++) { // loop over nodes in one edge


                                for (int n = 0; n < N; n++) { // loop over Test polynomials (the same as the basis polynomials)
                                    for (int d = 0; d < D; d++) { // loop over spatial directions

                                        if (coordSys == CoordinateSystem.Physical) {
                                            //if (!_Context.Edges.IsEdgeAffineLinear(iEdge))
                                            //    throw new NotImplementedException("curved elm: todo");
                                            //double Nd = normals[iEdge, d];
                                            //EvalResult[i, iNode, n, d] = BasisValues[i, iNode, n]*Nd*Nsign;
                                            throw new NotImplementedException("todo");


                                        } else {
                                            Debug.Assert(coordSys == CoordinateSystem.Reference);
                                            double Nd = splx.FaceNormals[iFace, d];
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
                    //var A = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, u, 0, 0, 0 }, new int[] { i-1, u-1, N-1, N-1, D-1 });
                    //var B = RHS.ExtractSubArrayShallow(new int[] { 0, 0, 0, jSgrd }, new int[] { N-1, N-1, D-1, jSgrd-1 });
                    ////var B = SurfResultA1.ExtractSubArrayShallow(new int[] { jCell, 0, 0, 0 }, new int[] { jCell-1, N-1, N-1, D-1 });
                    //B.Acc(1.0, A);
                    int NoOfFaces = _NoOfFaces;


                    for (int i = 0; i < Length; i++) {
                        int jCell = i0 + i;
                        int jSgrd = j2Sgrd[jCell];

                        var ResPart = RHS.ExtractSubArrayShallow(new int[] { 0, 0, jSgrd }, new int[] { N - 1, D - 1, jSgrd - 1 });

                        for (int e = 0; e < NoOfFaces; e++) {
                            var ip = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { i, e, 0, 0 }, new int[] { i - 1, e - 1, N - 1, D - 1 });
                            ResPart.Acc(1.0, ip);
                        }

                    }
                },
                cs: coordSys);
            qBnd.Execute();
        }
    }

}
