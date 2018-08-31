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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using System.IO;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    public class LineAndPointQuadratureFactory {

        public LineAndPointQuadratureFactory(RefElement Kref, 
            LevelSetTracker.LevelSetData levelSetData, 
            bool SupportPointrule, LineSegment.IRootFindingAlgorithm rootFindingAlgorithm = null) {
            this.m_RefElement = Kref;
            this.Tolerance = 1e-13;
            this.LevelSetData = levelSetData;
            this.SupportPointrule = SupportPointrule;
            this.RootFindingAlgorithm = rootFindingAlgorithm ?? new LineSegment.SafeGuardedNewtonMethod(this.Tolerance*0.1);
            this.referenceLineSegments = GetReferenceLineSegments(out this.segmentSorting, this.m_RefElement, this.RootFindingAlgorithm, levelSetData, this.LevelSetIndex);
           

            if (this.LevelSetData.GridDat.Grid.SpatialDimension != 2)
                throw new NotSupportedException();

            this.iKref = this.LevelSetData.GridDat.Grid.RefElements.IndexOf(this.m_RefElement, (A, B) => object.ReferenceEquals(A, B));
            if (this.iKref < 0)
                throw new ArgumentException("Reference element cannot be found in the provided grid.");

            

            this.MaxGrid = this.LevelSetData.GridDat.Cells.GetCells4Refelement(this.iKref).Intersect(
                LevelSetData.Region.GetCutCellMask4LevSet(this.LevelSetIndex).ToGeometicalMask());
        }

        int LevelSetIndex {
            get {
                return LevelSetData.LevelSetIndex;
            }
        }

        /// <summary>
        /// the intersection of the cut cells for Level Set <see cref="LevelSetIndex"/>
        /// and the subgrid for reference element <see cref="m_RefElement"/>
        /// </summary>
        private CellMask MaxGrid;

        public double Tolerance {
            get;
            private set;
        }

        private static Line lineSimplex = Line.Instance;

        private LineSegment[] referenceLineSegments;

        private int[] segmentSorting;

        public LineSegment.IRootFindingAlgorithm RootFindingAlgorithm {
            get;
            private set;
        }

        int iKref;

        /// <summary>
        /// grid reference element.
        /// </summary>
        RefElement m_RefElement;

        /// <summary>
        ///Node-wise evaluation of the level-set field.
        /// </summary>
        protected LevelSetTracker.LevelSetData LevelSetData;

        bool SupportPointrule;


        ///// <summary>
        ///// the awesome level set tracker
        ///// </summary>
        //protected LevelSetTracker tracker;

        abstract internal class QRF : IQuadRuleFactory<CellBoundaryQuadRule> {
            internal LineAndPointQuadratureFactory m_Owner;

            /// <summary>
            /// cache<br/>
            /// key: quadrature order; <br/>
            /// value: cached quadrature rule
            /// </summary>
            internal Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> Rules;

            /// <summary>
            /// If there are any cached rules, this method returns their order.
            /// </summary>
            public int[] GetCachedRuleOrders() {
                return Rules.Keys.ToArray();
            }


            public RefElement RefElement {
                get {
                    return m_Owner.m_RefElement;
                }
            }


            private ChunkRulePair<CellBoundaryQuadRule> GetUncutRule(int jCell, int order) {
                int iLs = this.m_Owner.LevelSetData.LevelSetIndex;
                int iDist = LevelSetTracker.DecodeLevelSetDist(this.m_Owner.LevelSetData.Region.m_LevSetRegions[jCell], iLs);
                if ( this.m_Owner.LevelSetData.GridDat.Cells.GetRefElementIndex(jCell) != m_Owner.iKref)
                    throw new ArgumentException("illegal cell mask.");
                if (iDist > 0)
                    return new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), GetPosRule(order));
                else if (iDist < 0)
                    return new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), GetNegRule(order));
                else
                    throw new ApplicationException();
            }

            /// <summary>
            /// Rule for an _uncut_ cell in the positive level-set region.
            /// </summary>
            protected abstract CellBoundaryQuadRule GetPosRule(int order);

            /// <summary>
            /// Rule for an _uncut_ cell in the negative level-set region.
            /// </summary>
            protected abstract CellBoundaryQuadRule GetNegRule(int order);


            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");

                if (!Rules.ContainsKey(order))
                    m_Owner.GetQuadRuleSet_Internal(order);

                var Ret = new ChunkRulePair<CellBoundaryQuadRule>[mask.NoOfItemsLocally];
                var Rule = Rules[order];

                int iR = 0;
                int iA = 0;
                var enuB = mask.GetItemEnumerator();

                while (enuB.MoveNext()) {
                    int jCell = enuB.Current;

                    if (iA >= Rule.Length) {
                        Ret[iR] = GetUncutRule(jCell, order); // outside
                        Debug.Assert(Ret[iR].Rule.Nodes.IsLocked);
                        iR++;
                    } else {
                        if (Rule[iA].Chunk.i0 > jCell) {
                            Ret[iR] = GetUncutRule(jCell, order); // outside
                            Debug.Assert(Ret[iR].Rule.Nodes.IsLocked);
                            iR++;
                        } else {
                            while (iA < Rule.Length && Rule[iA].Chunk.i0 < jCell) {
                                Debug.Assert(Rule[iA].Chunk.Len == 1);
                                iA++;
                            }

                            if (iA >= Rule.Length) {
                                Ret[iR] = GetUncutRule(jCell, order); // outside
                                Debug.Assert(Ret[iR].Rule.Nodes.IsLocked);
                                iR++;
                            } else {
                                if (Rule[iA].Chunk.i0 == jCell) {
                                    Debug.Assert(Rule[iA].Chunk.Len == 1);
                                    Ret[iR] = Rule[iA]; // inside
                                    Debug.Assert(Ret[iR].Rule.Nodes.IsLocked);
                                    iR++;
                                } else {
                                    Ret[iR] = GetUncutRule(jCell, order); // outside
                                    Debug.Assert(Ret[iR].Rule.Nodes.IsLocked);
                                    iR++;
                                }
                            }
                        }
                    }

                    

                }


                return Ret;


                /*
#if DEBUG
                if (mask.Except(m_Owner.MaxGrid).NoOfItemsLocally > 0)
                    throw new NotSupportedException("'mask' must be a subset of the cut cells, for my reference element.");
#endif
                if (!Rules.ContainsKey(order))
                    m_Owner.GetQuadRuleSet_Internal(order);

                if (mask.NoOfItemsLocally == m_Owner.MaxGrid.NoOfItemsLocally) {
                    // aggressive
                    return this.Rules[order];
                } else {
                    var Rule = Rules[order];

                    var Ret = new ChunkRulePair<CellBoundaryQuadRule>[mask.NoOfItemsLocally];

                    int L = Ret.Length, H = Rule.Length;
                    int h = 0;
                    int jsub = -1;
                    foreach(int jCell in mask.ItemEnum) {
                        jsub++;

                        Debug.Assert(Rule[h].Chunk.Len == 1);
                        while (jCell > Rule[h].Chunk.i0) {
                            h++;
                        }

                        Debug.Assert(jCell == Rule[h].Chunk.i0);
                        Ret[jsub] = Rule[h];
                    }

                    return Ret;
                }
                 * 
                 */
            }
        }

        internal class LineQRF : QRF {
            public LineQRF(LineAndPointQuadratureFactory o)
                : base() {
                base.m_Owner = o;
                base.Rules = o.m_LineMeasure;
                int D = o.LevelSetData.GridDat.SpatialDimension;
                this.empty = CellBoundaryQuadRule.CreateEmpty(o.m_RefElement, 1, D, o.referenceLineSegments.Length);
                this.empty.Nodes.LockForever();
            }

            CellBoundaryQuadRule empty;

            int fullOrder = int.MinValue;
            CellBoundaryQuadRule full;

            protected override CellBoundaryQuadRule GetPosRule(int order) {

                if (full == null || fullOrder != order) {
                    full = this.RefElement.GetBoundaryQuadRule(order);
                }
                Debug.Assert(full.Nodes.IsLocked);
                return full;
            }

            protected override CellBoundaryQuadRule GetNegRule(int order) {
                Debug.Assert(empty.Nodes.IsLocked);
                return empty;
            }
        }

        public IQuadRuleFactory<CellBoundaryQuadRule> GetLineFactory() {
            int D = this.LevelSetData.GridDat.SpatialDimension;
            return new LineQRF(this);
        }


        class PointQRF : QRF {
            public PointQRF(LineAndPointQuadratureFactory o)
                : base() {
                base.m_Owner = o;
                base.Rules = o.m_PointMeasure;
                int D = o.LevelSetData.GridDat.SpatialDimension;
                if (D != 2)
                    throw new NotSupportedException("the point measure is only supported in 2D");
                this.empty = CellBoundaryQuadRule.CreateEmpty(o.m_RefElement, 1, D, o.referenceLineSegments.Length);
                this.empty.Nodes.LockForever();
            }

            CellBoundaryQuadRule empty;

            protected override CellBoundaryQuadRule GetPosRule(int order) {
                Debug.Assert(empty.Nodes.IsLocked);
                return empty;
            }

            protected override CellBoundaryQuadRule GetNegRule(int order) {
                Debug.Assert(empty.Nodes.IsLocked);
                return empty;
            }
        }


        public IQuadRuleFactory<CellBoundaryQuadRule> GetPointFactory() {
            if(!this.SupportPointrule)
                throw new NotSupportedException("point measure creation is turned off");
            return new PointQRF(this);
        }

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> m_PointMeasure = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]>();


        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> m_LineMeasure = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]>();




        void GetQuadRuleSet_Internal(int order) {
            var mask = this.MaxGrid;
            var grdDat = this.LevelSetData.GridDat;
            int D = grdDat.SpatialDimension;
            var _Cells = grdDat.Cells;
            var scalings = grdDat.Edges.SqrtGramian;
            var cell2Edge = grdDat.Cells.Cells2Edges;
            //var edge2Cell = tracker.GridDat.Edges.CellIndices;
            //var FaceIdx = tracker.GridDat.Edges.FaceIndices;
            var EdgeData = grdDat.Edges;
            int levSetIndex = this.LevelSetIndex;
            var Vertices = this.m_RefElement.Vertices;
            var iKref = this.iKref;
            double[] EdgeToVolumeTransformationDeterminants = this.m_RefElement.FaceTrafoGramianSqrt;
            QuadRule baseRule = lineSimplex.GetQuadratureRule(order);
            double baseRuleWeightSum = baseRule.Weights.Sum();

            int[][] Vtx2Face;
            switch (this.m_RefElement.SpatialDimension) {
                case 2: Vtx2Face = this.m_RefElement.VertexIndicesToFaces; break;
                case 3: Vtx2Face = this.m_RefElement.VertexIndicesToCoFaces; throw new NotSupportedException("total ungetestet! => untersuchen!");
                default: throw new NotSupportedException();
            }

            double isVertexTol = this.Tolerance*10;
            double isVertexSaveTol = isVertexTol*4;
            double rootFilterTol = isVertexSaveTol*4;

#if DEBUG
            for (int e = 0; e < referenceLineSegments.Length; e++) {
                LineSegment referenceSegment = referenceLineSegments[e];

                double[] P0 = referenceSegment.GetPointOnSegment(-1.0);
                double[] P1 = referenceSegment.GetPointOnSegment(+1.0);

                int iP0 = this.m_RefElement.Vertices.MindistRow(P0);
                int iP1 = this.m_RefElement.Vertices.MindistRow(P1);

                int _iP0 = this.m_RefElement.Vertices.MindistRow(P0);
                int _iP1 = this.m_RefElement.Vertices.MindistRow(P1);

                Debug.Assert(GenericBlas.L2Dist(P0, this.m_RefElement.Vertices.GetRow(referenceSegment.iVertexStart)) <= 1.0e-8);
                Debug.Assert(GenericBlas.L2Dist(P1, this.m_RefElement.Vertices.GetRow(referenceSegment.iVertexEnd)) <= 1.0e-8);
                Debug.Assert(_iP0 == iP0);
                Debug.Assert(_iP1 == iP1);
            }
#endif

            // required for the point rule
            // ===========================



            var PointMeasure_result = new List<ChunkRulePair<CellBoundaryQuadRule>>(mask.NoOfItemsLocally);
            var LineMeasure_result = new List<ChunkRulePair<CellBoundaryQuadRule>>(mask.NoOfItemsLocally);
            foreach (Chunk chunk in mask) {            // loop over cells
                for (int i = 0; i < chunk.Len; i++) {  // loop over cells...
                    int jCell = i + chunk.i0;
                    Debug.Assert(iKref == _Cells.GetRefElementIndex(jCell));

                    var cell2Edge_j = cell2Edge[jCell];



                    //if (cache.ContainsKey(cell)) {
                    //    result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                    //        Chunk.GetSingleElementChunk(cell), cache[cell]));
                    //    continue;
                    //}

                    List<double[]> LnMeas_nodes = new List<double[]>();
                    List<double> LnMeas_weights = new List<double>();
                    int[] LnMeas_noOfNodesPerEdge = new int[referenceLineSegments.Length];
                    int[] PtMeas_noOfNodesPerEdge = new int[referenceLineSegments.Length];

                    // ==================
                    // find roots
                    // ==================

                    double[][] _roots = new double[referenceLineSegments.Length][];
                    int TotalNoOfRoots = 0;
                    for (int e = 0; e < referenceLineSegments.Length; e++) {
                        LineSegment referenceSegment = referenceLineSegments[e];
                        double[] roots = referenceSegment.GetRoots(LevelSetData.LevelSet, jCell, iKref);
                        _roots[e] = FilterRoots(roots, rootFilterTol);
                        TotalNoOfRoots += _roots[e].Length;
                    }




                    // ============================
                    // create line measure
                    // ============================
                    bool LineMeasureEmptyOrFull = false;
                    List<LineSegment>[] ActiveSegments = new List<LineSegment>[_roots.Length];
                    double nodesSum = 0;
                    {
                        double Fullsum = 0;
                        for (int e = 0; e < referenceLineSegments.Length; e++) {
                            ActiveSegments[e] = new List<LineSegment>();
                            LineSegment referenceSegment = referenceLineSegments[e];
                            double edgeDet = EdgeToVolumeTransformationDeterminants[e];
                            double[] roots = _roots[e];
                            Fullsum += referenceSegment.Length*edgeDet;
                            LineSegment[] subSegments = referenceSegment.Split(roots);
                            

                            for (int k = 0; k < subSegments.Length; k++) {
                                // Evaluate sub segment at center to determine sign
                                NodeSet _point = new NodeSet(this.m_RefElement, subSegments[k].GetPointOnSegment(0.0));
                                
                                double weightFactor = subSegments[k].Length / referenceSegment.Length;

                                if (weightFactor < this.Tolerance)
                                    // segment has a length of approximately 0.0 => no need to care about it.
                                    continue;

                                weightFactor *= edgeDet;


                                //bool center = false;
                                MultidimensionalArray levelSetValue = LevelSetData.GetLevSetValues(_point, jCell, 1);
                                if (levelSetValue[0, 0] <= -Tolerance) {
                                    continue;
                                }

                                

                                ActiveSegments[e].Add(subSegments[k]);

                                //bool start, end;
                                //double[] StPoint = subSegments[k].GetPointOnSegment(-1.0);
                                //using (tracker.GridDat.NSC.CreateLock(MultidimensionalArray.CreateWrapper(point, 1, D), this.iKref, -1.0)) {
                                //    MultidimensionalArray levelSetValue = tracker.GetLevSetValues(levSetIndex, 0, jCell, 1);
                                //    start = levelSetValue[0, 0] > -Tolerance;
                                //}
                                //double[] EnPoint = subSegments[k].GetPointOnSegment(+1.0);
                                //using (tracker.GridDat.NSC.CreateLock(MultidimensionalArray.CreateWrapper(point, 1, D), this.iKref, -1.0)) {
                                //    MultidimensionalArray levelSetValue = tracker.GetLevSetValues(levSetIndex, 0, jCell, 1);
                                //    end = levelSetValue[0, 0] > -Tolerance;
                                //}

                                for (int m = 0; m < baseRule.NoOfNodes; m++) {
                                    // Base rule _always_ is a line rule, thus Nodes[*, _0_]
                                    double[] point = subSegments[k].GetPointOnSegment(baseRule.Nodes[m, 0]);

                                    LnMeas_weights.Add(weightFactor * baseRule.Weights[m]);
                                    LnMeas_nodes.Add(point);

                                    LnMeas_noOfNodesPerEdge[e]++;
                                }
                            }
                        }

                        if (LnMeas_weights.Count == 0) {
                            var emptyrule = CellBoundaryQuadRule.CreateEmpty(this.m_RefElement, 1, D, referenceLineSegments.Length);
                            // create a rule with just one node and weight zero;
                            // this should avoid some special-case handling for empty rules
                            emptyrule.NumbersOfNodesPerFace[0] = 1;
                            emptyrule.Nodes.LockForever();

                            LineMeasure_result.Add(new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), emptyrule));

                            LineMeasureEmptyOrFull = true;
                        } else {
                            NodeSet localNodes = new NodeSet(this.m_RefElement, LnMeas_nodes.Count, D);
                            for (int j = 0; j < LnMeas_nodes.Count; j++) {
                                for (int d = 0; d < D; d++) {
                                    localNodes[j, d] = LnMeas_nodes[j][d];
                                }
                            }
                            localNodes.LockForever();

                            CellBoundaryQuadRule subdividedRule = new CellBoundaryQuadRule() {
                                OrderOfPrecision = order,
                                Weights = MultidimensionalArray.Create(LnMeas_weights.Count),
                                Nodes = localNodes,
                                NumbersOfNodesPerFace = LnMeas_noOfNodesPerEdge
                            };
                            subdividedRule.Weights.SetSubVector(LnMeas_weights, -1);

                            LineMeasure_result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                                Chunk.GetSingleElementChunk(jCell), subdividedRule));

                            nodesSum = subdividedRule.Weights.Sum();
                            LineMeasureEmptyOrFull = (Math.Abs(nodesSum) <= this.Tolerance*10) || (Math.Abs(nodesSum - Fullsum) <= this.Tolerance*10);
                        }
                    }


                    // ============================
                    // create point measure
                    // ============================
                    if (D == 2 && this.SupportPointrule) { // only in 2D

                        double[][] _newRoots = new double[_roots.Length][];
                        List<double> newRoot = new List<double>(5);
                        int[][] _rootIsAtVertex = new int[_roots.GetLength(0)][]; // indices correspond with '_roots':
                        //                                                           1st index: reference line segment; 
                        //                                                           2nd index: root enumeration
                        //                                              content: negative, if the root is not located at a vertex (of the ref element),
                        //                                                       otherwise the vertex index.
                        bool bRootAtVertex = false; 
                        List<int> rootIsAtVertex = new List<int>(5);
                        
                        int newRootCounter = 0;
                        for (int _e = 0; _e < _roots.Length; _e++) {
                            newRoot.Clear();
                            rootIsAtVertex.Clear();

                            int e_curr = this.segmentSorting[_e];
                            int e_next = this.segmentSorting[(_e + 1)%_roots.Length];
                            int e_prev = this.segmentSorting[_e > 0 ? _e - 1 : _roots.Length - 1];

                            var segs_curr = ActiveSegments[Math.Abs(e_curr)];
                            if (segs_curr.Count == 0) {
                                _newRoots[Math.Abs(e_curr)] = new double[0];
                                _rootIsAtVertex[Math.Abs(e_curr)] = new int[0];
                                continue;
                            }
                            var segs_prev = ActiveSegments[Math.Abs(e_prev)];
                            var segs_next = ActiveSegments[Math.Abs(e_next)];

                            int iVtxStart, iVtxEnd;
                            if(e_curr >= 0) {
                                iVtxStart = this.referenceLineSegments[Math.Abs(e_curr)].iVertexStart;
                                iVtxEnd = this.referenceLineSegments[Math.Abs(e_curr)].iVertexEnd;
                            } else {
                                iVtxStart = this.referenceLineSegments[Math.Abs(e_curr)].iVertexEnd;
                                iVtxEnd = this.referenceLineSegments[Math.Abs(e_curr)].iVertexStart;
                            }

                            bool prevOnEdge;
                            if (segs_prev.Count > 0) {
                                double endCoord;
                                if (e_prev >= 0) {
                                    endCoord = segs_prev.Last().EndCoord;
                                } else {
                                    endCoord = -segs_prev.First().StartCoord;
                                }

                                prevOnEdge = (endCoord == 1.0);
                            } else {
                                prevOnEdge = false;
                            }

                            bool nextOnEdge;
                            if (segs_next.Count > 0) {
                                double endCoord;
                                if (e_next >= 0) {
                                    endCoord = segs_next.First().StartCoord;
                                } else {
                                    endCoord = -segs_next.Last().EndCoord;
                                }

                                nextOnEdge = (endCoord == -1.0);
                            } else {
                                nextOnEdge = false;
                            }

                            if (e_curr >= 0) {
                                int I = segs_curr.Count;
                                for (int j = 0; j < I; j++) {
                                    var seg = segs_curr[j];

                                    if (j <= 0) {
                                        bool startAtVertex = seg.StartCoord == -1.0;
                                        if (!(startAtVertex && prevOnEdge)) {
                                            newRoot.Add(seg.StartCoord);
                                            if (startAtVertex) {
                                                rootIsAtVertex.Add(iVtxStart);
                                                bRootAtVertex = true;
                                            } else {
                                                rootIsAtVertex.Add(int.MinValue);
                                            }
                                        }
                                    } else {
                                        if (!(seg.StartCoord == segs_curr[j - 1].EndCoord)) {
                                            newRoot.Add(seg.StartCoord);
                                            rootIsAtVertex.Add(int.MinValue);
                                        }
                                    }

                                    if (j >= (I - 1)) {
                                        bool endAtVertex = seg.EndCoord == +1.0;
                                        if (!(endAtVertex && nextOnEdge)) {
                                            newRoot.Add(seg.EndCoord);
                                            if (endAtVertex) {
                                                rootIsAtVertex.Add(iVtxEnd);
                                                bRootAtVertex = true;
                                            } else {
                                                rootIsAtVertex.Add(int.MinValue);
                                            }
                                        }
                                    } else {
                                        if (!(seg.EndCoord == segs_curr[j + 1].StartCoord)) {
                                            newRoot.Add(seg.EndCoord);
                                            rootIsAtVertex.Add(int.MinValue);
                                        }
                                    }
                                }
                            } else {
                                int I = segs_curr.Count;
                                for (int j = I - 1; j >= 0; j--) {
                                    var seg = segs_curr[j];

                                    if (j >= (I - 1)) {
                                        bool startAtVertex = (-seg.EndCoord == -1.0);
                                        if (!(startAtVertex && prevOnEdge)) {
                                            newRoot.Add(seg.EndCoord);
                                            if (startAtVertex) {
                                                rootIsAtVertex.Add(iVtxStart);
                                                bRootAtVertex = true;
                                            } else {
                                                rootIsAtVertex.Add(int.MinValue);
                                            }
                                        }
                                    } else {
                                        if (!(seg.EndCoord == segs_curr[j + 1].StartCoord)) {
                                            newRoot.Add(seg.EndCoord);
                                            rootIsAtVertex.Add(int.MinValue);
                                        }
                                    }

                                    if (j <= 0) {
                                        bool endAtVertex = (-seg.StartCoord == +1.0);
                                        if (!(endAtVertex && nextOnEdge)) {
                                            newRoot.Add(seg.StartCoord);
                                            if (endAtVertex) {
                                                rootIsAtVertex.Add(iVtxEnd);
                                                bRootAtVertex = true;
                                            } else {
                                                rootIsAtVertex.Add(int.MinValue);
                                            }
                                        }
                                    } else {
                                        if (!(seg.StartCoord == segs_curr[j - 1].EndCoord)) {
                                            newRoot.Add(seg.StartCoord);
                                            rootIsAtVertex.Add(int.MinValue);
                                        }
                                    }
                                }
                            }

                            Debug.Assert(newRoot.Count == rootIsAtVertex.Count);
                            _newRoots[Math.Abs(e_curr)] = newRoot.ToArray();
                            _rootIsAtVertex[Math.Abs(e_curr)] = rootIsAtVertex.ToArray();
                            newRootCounter += newRoot.Count();
                        }

                        if (newRootCounter%2 != 0)
                            throw new ArgumentException("error in alg");

                        _roots = _newRoots;


                        if (LineMeasureEmptyOrFull) {
                            // ++++++++++++++
                            // the empty case
                            // ++++++++++++++

                            var emptyrule = CellBoundaryQuadRule.CreateEmpty(this.m_RefElement, 1, D, referenceLineSegments.Length);
                            // create a rule with just one node and weight zero;
                            // this should avoid some special-case handling for empty rules
                            emptyrule.NumbersOfNodesPerFace[0] = 1;
                            emptyrule.Nodes.LockForever();

                            PointMeasure_result.Add(new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), emptyrule));
                        }
                        else  
                        {
                            // identify if roots are at vertices of the RefElement
                            // ===================================================

                            List<double[]> PtMeas_nodes = new List<double[]>();
                            for (int e = 0; e < referenceLineSegments.Length; e++) {
                                LineSegment referenceSegment = referenceLineSegments[e];
                                double[] roots = _roots[e];
                                
                                /*
                                 _rootIsAtVertex[e] = new int[roots.Length];
                                //_rootIsAtVertex[e].SetAll(int.MinValue);
                                 */

                                // check whether the first and/or last root correspond with a vertex
                                for (int h = 0; h < roots.Length; h += Math.Max(roots.Length - 1, 1)) {
                                    PtMeas_nodes.Add(referenceSegment.GetPointOnSegment(roots[h]));

                                    /*
                                    if (Math.Abs(roots[h] - (-1.0)) < isVertexTol) {
                                        // found a root at a RefElement-Vertex
                                        // this happens seldom and needs special treatment (see below)
                                        
                                        _rootIsAtVertex[e][h] = referenceSegment.iVertexStart;
                                        bRootAtVertex = true;
                                    }

                                    if (Math.Abs(roots[h] - (+1.0)) < isVertexTol) {
                                        // found a root at a RefElement-Vertex
                                        // this happens seldom and needs special treatment (see below)
                                        
                                        _rootIsAtVertex[e][h] = referenceSegment.iVertexEnd;
                                        bRootAtVertex = true;
                                    }
                                     */
                                }
                            }

                            if (TotalNoOfRoots == 2 && !bRootAtVertex) {
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                                // this will hopefully cover most of the cases;
                                // we try to omit expensive handling of special cases
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++

                                Debug.Assert(PtMeas_nodes.Count == 2);

                                CellBoundaryQuadRule subdividedRule = new CellBoundaryQuadRule() {
                                    OrderOfPrecision = order,
                                    Weights = MultidimensionalArray.Create(PtMeas_nodes.Count),
                                    Nodes = new NodeSet(this.m_RefElement, PtMeas_nodes.Count, D),
                                    NumbersOfNodesPerFace = PtMeas_noOfNodesPerEdge
                                };

                                var PtMeas_weights = NewMethod(scalings, EdgeData, jCell, cell2Edge_j, _roots);

                                subdividedRule.Weights.SetVector(PtMeas_weights);

                                for (int j = 0; j < PtMeas_nodes.Count; j++) {
                                    for (int d = 0; d < D; d++) {
                                        subdividedRule.Nodes[j, d] = PtMeas_nodes[j][d];
                                    }
                                }

                                for (int e = 0; e < _roots.Length; e++) {
                                    subdividedRule.NumbersOfNodesPerFace[e] = _roots[e].Length;
                                }

                                subdividedRule.Nodes.LockForever();

                                PointMeasure_result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                                    Chunk.GetSingleElementChunk(jCell), subdividedRule));
                            } else {
                                // +++++++++++++++++++++++++++++++++++
                                // the general case -- a bit tricky
                                // +++++++++++++++++++++++++++++++++++

                                // ----
                                // this is one of the most fucked-up pieces of code that I have ever written.
                                // ----

                                // ensure that roots at Refelement-vertices are present on all faces
                                // =================================================================
                                List<Tuple<int, int>>[] NodesAtVertices = null;
                                if (bRootAtVertex) {
                                    NodesAtVertices = new List<Tuple<int, int>>[this.m_RefElement.NoOfVertices];

                                    for (int e = 0; e < referenceLineSegments.Length; e++) {
                                        LineSegment referenceSegment = referenceLineSegments[e];
                                        double[] roots = _roots[e];

                                        for (int h = roots.Length - 1; h >= 0; h--) {
                                            int iVertex = _rootIsAtVertex[e][h];

                                            double[] pt = null;

                                            if (iVertex >= 0) {
                                                // ensure that this root is present for all other faces
                                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++


                                                if (NodesAtVertices[iVertex] == null)
                                                    NodesAtVertices[iVertex] = new List<Tuple<int, int>>();

                                                var nodeFace = new Tuple<int, int>(e, h);
                                                if (!NodesAtVertices[iVertex].Contains(nodeFace))
                                                    NodesAtVertices[iVertex].Add(nodeFace);

                                                foreach (int iF in Vtx2Face[iVertex]) {
                                                    if (iF == e)
                                                        continue;

                                                    if (!_rootIsAtVertex[iF].Contains(iVertex)) {
                                                        if (pt == null)
                                                            pt = referenceSegment.GetPointOnSegment(roots[h]);

                                                        var alpha = referenceLineSegments[iF].GetSegmentCoordinateForPoint(pt);
                                                        Debug.Assert(alpha >= -1.00000001);
                                                        Debug.Assert(alpha <= +1.00000001);
                                                        Debug.Assert(GenericBlas.L2DistPow2(pt, referenceLineSegments[iF].GetPointOnSegment(alpha)) < this.Tolerance*1000);


                                                        alpha.AddToArray(ref _roots[iF]);
                                                        iVertex.AddToArray(ref _rootIsAtVertex[iF]);


                                                        var nodeFace2 = new Tuple<int, int>(iF, _roots[iF].Length - 1);
                                                        if (!NodesAtVertices[iVertex].Contains(nodeFace2))
                                                            NodesAtVertices[iVertex].Add(nodeFace2);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                TotalNoOfRoots = _roots.Sum(r => r.Length);

                                // compute inner products between cell normals and level-set normals
                                // ==================================================================


                                NodeSet[] Nodes = new NodeSet[referenceLineSegments.Length];
                                MultidimensionalArray[] LevSetNormals = new MultidimensionalArray[referenceLineSegments.Length];
                                var simplexNormals = this.m_RefElement.FaceNormals;
                                double[][] Innerproducts = new double[referenceLineSegments.Length][];
                                for (int e = 0; e < referenceLineSegments.Length; e++) {
                                    double[] roots = _roots[e];
                                    int K = roots.Length;
                                    Nodes[e] = new NodeSet(this.m_RefElement, K, D);
                                    var Nodes_e = Nodes[e];

                                    LineSegment referenceSegment = referenceLineSegments[e];
                                    for (int k = 0; k < K; k++) {
                                        Nodes_e.SetRow(k, referenceSegment.GetPointOnSegment(roots[k]));
                                    }
                                    Nodes_e.LockForever();
                                    Innerproducts[e] = new double[K];

                                    if (K > 0) {
                                        LevSetNormals[e] = this.LevelSetData.GetLevelSetReferenceNormals(Nodes_e, jCell, 1);
                                        var LevSetNormals_e = LevSetNormals[e];


                                        var Innerproducts_e =  Innerproducts[e];
                                        for (int k = 0; k < K; k++) {
                                            double acc = 0;
                                            acc += simplexNormals[e, 0]*LevSetNormals_e[0, k, 0];
                                            acc += simplexNormals[e, 1]*LevSetNormals_e[0, k, 1];
                                            Innerproducts_e[k] = acc;
                                        }
                                    }
                                }


                                // treatment of vertex-nodes
                                // =========================
                                bool[][] KilledRoots = _roots.Select(rr => new bool[rr.Length]).ToArray();
                                if (bRootAtVertex) {
                                    for (int iVtx = 0; iVtx < this.m_RefElement.NoOfVertices; iVtx++) {
                                        var Nodes_iVtx = NodesAtVertices[iVtx];
                                        if (Nodes_iVtx == null || Nodes_iVtx.Count == 0) {
                                            continue;

                                        } else if (Nodes_iVtx.Count == 2) {

                                            var t0 = Nodes_iVtx[0];
                                            var t1 = Nodes_iVtx[1];

                                            double ip0 = Innerproducts[t0.Item1][t0.Item2];
                                            double ip1 = Innerproducts[t1.Item1][t1.Item2];
                                            
                                            if (Math.Abs(ip0) < Math.Abs(ip1)) {
                                                // choose node 0
                                                KilledRoots[t1.Item1][t1.Item2] = true;
                                            } else {
                                                // choose node 1
                                                KilledRoots[t0.Item1][t0.Item2] = true;
                                            }

                                        } else {
                                            throw new ApplicationException("error in alg");
                                        }
                                    }
                                }
                                
                                // wackelkandidaten finden
                                // =======================
                                List<Tuple<int, int>> SaveNodes = new List<Tuple<int, int>>();
                                Tuple<int, int> Least_UncertainNode = null;
                                double Least_UncertainNode_flatness = 0;

                                for (int e = 0; e < _roots.Length; e++) {
                                    var Innerproducts_e =  Innerproducts[e];
                                    int K = Innerproducts_e.Length;
                                    for (int k = 0; k < K; k++) {
                                        if (KilledRoots[e][k])
                                            continue;


                                        double flatNess = Math.Abs(Math.Abs(Innerproducts_e[k]) - 1.0);

                                        if (flatNess <= this.Tolerance*10) {
                                            // level-set at node 'k' is almost tangential to face/co-face 'e'

                                            if (Least_UncertainNode == null) {
                                                Least_UncertainNode = new Tuple<int, int>(e, k);
                                                Least_UncertainNode_flatness = flatNess;
                                            } else {
                                                if (flatNess > Least_UncertainNode_flatness) {
                                                    Least_UncertainNode = new Tuple<int, int>(e, k);
                                                    Least_UncertainNode_flatness = flatNess;
                                                }
                                            }
                                        } else {
                                            SaveNodes.Add(new Tuple<int, int>(e, k));
                                        }
                                    }
                                }
                                
                                // create point measure
                                // ====================
                                if (SaveNodes.Count == 0) {
                                    var emptyrule = CellBoundaryQuadRule.CreateEmpty(this.m_RefElement, 1, D, referenceLineSegments.Length);
                                    // create a rule with just one node and weight zero;
                                    // this should avoid some special-case handling for empty rules
                                    emptyrule.NumbersOfNodesPerFace[0] = 1;
                                    emptyrule.Nodes.LockForever();

                                    PointMeasure_result.Add(new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(jCell), emptyrule));
                                    
                                } else {

                                    if ((SaveNodes.Count % 2) != 0) {
                                        if (Least_UncertainNode == null)
                                            throw new ApplicationException();
                                        
                                        if(Least_UncertainNode.Item1 < SaveNodes[0].Item1) {
                                            SaveNodes.Insert(0,Least_UncertainNode);
                                        } else {
                                            bool notatend = false;
                                            for(int jj = 0; jj < SaveNodes.Count; jj++) {
                                                var t = SaveNodes[jj];
                                                if(t.Item1 == Least_UncertainNode.Item1) {
                                                    SaveNodes.Insert(jj, Least_UncertainNode);
                                                    notatend = true;
                                                    break;
                                                }
                                            }
                                            if(!notatend)
                                                SaveNodes.Add(Least_UncertainNode);
                                        }
                                    }

                                    CellBoundaryQuadRule subdividedRule = new CellBoundaryQuadRule() {
                                        OrderOfPrecision = order,
                                        Weights = MultidimensionalArray.Create(SaveNodes.Count),
                                        Nodes = new NodeSet(this.m_RefElement, SaveNodes.Count, D),
                                        NumbersOfNodesPerFace = new int[_roots.Length]
                                    };

                                    var PtMeas_weights = NewMethod(scalings, grdDat.Edges, jCell, cell2Edge_j, _roots);
                                    double[][] _PtMeas_weights = new double[_roots.Length][];
                                    int cnt = 0;
                                    for (int e = 0; e < _roots.Length; e++) {
                                        _PtMeas_weights[e] = PtMeas_weights.GetSubVector(cnt, _roots[e].Length);
                                        cnt += _roots[e].Length;
                                    }

                                    int k = 0;
                                    foreach (var t in SaveNodes) {
                                        int e = t.Item1;
                                        int kk = t.Item2;
                                        subdividedRule.NumbersOfNodesPerFace[e]++;
                                        var pt = referenceLineSegments[e].GetPointOnSegment(_roots[e][kk]);
                                        for (int d = 0; d < D; d++) {
                                            subdividedRule.Nodes[k, d] = pt[d];
                                        }
                                        subdividedRule.Weights[k] = _PtMeas_weights[e][kk];
                                        k++;
                                    }

                                    subdividedRule.Nodes.LockForever();

                                    PointMeasure_result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                                        Chunk.GetSingleElementChunk(jCell), subdividedRule));
                                }

                            }
                        }
                    }
                }
            }
            // cache result
            Debug.Assert(PointMeasure_result.Any(crp => crp.Rule.Nodes.IsLocked == false) == false);
            Debug.Assert(LineMeasure_result.Any(crp => crp.Rule.Nodes.IsLocked == false) == false);
            this.m_PointMeasure.Add(order, PointMeasure_result.ToArray());
            this.m_LineMeasure.Add(order, LineMeasure_result.ToArray());
        }

        private static List<double> NewMethod(MultidimensionalArray scalings, GridData.EdgeData EdgeData, int jCell, int[] cell2Edge_j, double[][] _roots) {
            var edge2Cell = EdgeData.CellIndices;
            var FaceIdx = EdgeData.FaceIndices;
            
            var PtMeas_weights = new List<double>();
            for (int e = 0; e < _roots.Length; e++) { // loop over faces...

                double[] roots = _roots[e];


                int iEdge = -1;
                int _inOut = -1;
                {
                    // finding the edge: this code is provisory

                    //int _inOut = -1;
                    foreach (int em in cell2Edge_j) {
                        int _iedge = Math.Abs(em) - 1;
                        _inOut = em > 0 ? 0 : 1;

                        if (edge2Cell[_iedge, _inOut] == jCell && FaceIdx[_iedge, _inOut] == e) {
                            iEdge = _iedge;
                            break;
                        }
                    }
                    if (iEdge < 0)
                        throw new ApplicationException();
                    //if (!EdgeData.IsEdgeConformal(iEdge, _inOut))
                    //    throw new NotSupportedException("hanging nodes not supported");
                    if (!EdgeData.IsEdgeAffineLinear(iEdge))
                        throw new NotSupportedException("no curved element support yet.");
                }

                //for (int l = 0; l < roots.Length; l++) {
                //    PtMeas_weights.Add(1.0 / scalings[iEdge]);
                //}

                if (!EdgeData.IsEdgeConformal(iEdge, _inOut)) {
                    // compute new scaling for non-conforming edges

                    double scaling = EdgeData.GetSqrtGramianForNonConformEdge(iEdge, _inOut);

                    for (int l = 0; l < roots.Length; l++) {
                        PtMeas_weights.Add(1.0 / scaling);
                    }

                } else {
                    // use standard scaling of the edges
                    for (int l = 0; l < roots.Length; l++) {
                        PtMeas_weights.Add(1.0 / scalings[iEdge]);
                    }
                }
            }
            return PtMeas_weights;
        }


        static private double[] FilterRoots(double[] roots, double tol) {
            if (roots.Length <= 1) {
                return roots;
            } else {
                Array.Sort(roots);
                List<double> newRoots = new List<double>(roots.Length);
                newRoots.Add(roots[0]);
                int L = roots.Length;
                for (int i = 1; i < L; i++) {
                    double dist = roots[i] - roots[i-1];
                    Debug.Assert(dist >= 0);

                    if (dist < tol*10) {
                        newRoots[newRoots.Count - 1] = roots[i];
                        continue;
                    } else {
                        newRoots.Add(roots[i]);
                    }
                }
                return newRoots.ToArray();
            }
        }


        static private LineSegment[] GetReferenceLineSegments(out int[] segmentSort, RefElement Simplex, LineSegment.IRootFindingAlgorithm RootFindingAlgorithm, LevelSetTracker.LevelSetData levelSetData, int LevelSetIndex) {
            Stack<RefElement> simplexHierarchy = new Stack<RefElement>();
            RefElement currentSimplex = Simplex;
            int spatialDimension = Simplex.SpatialDimension;

            int n = 0;
            while (currentSimplex.GetType() != lineSimplex.GetType()) {
                // If n > 2, the edge of the edge of a simplex is not a line.
                // This is hardly likely in 2d/3d.
                if (n > 2) {
                    throw new ApplicationException("Something went terribly wrong.");
                }

                simplexHierarchy.Push(currentSimplex);

                currentSimplex = currentSimplex.FaceRefElement;
                n++;
            }

            MultidimensionalArray vertexCoordinates = MultidimensionalArray.Create(2, 1);
            vertexCoordinates[0, 0] = -1.0;
            vertexCoordinates[1, 0] = 1.0;

            while (simplexHierarchy.Count > 0) {
                currentSimplex = simplexHierarchy.Pop();

                int noOfVertices = vertexCoordinates.GetLength(0);
                int D = currentSimplex.SpatialDimension;
                MultidimensionalArray volumeCoordinates = MultidimensionalArray.Create(
                    noOfVertices * currentSimplex.NoOfFaces, currentSimplex.SpatialDimension);

                for (int e = 0; e < currentSimplex.NoOfFaces; e++) {
                    MultidimensionalArray coordinates = MultidimensionalArray.Create(noOfVertices, D);
                    currentSimplex.TransformFaceCoordinates(e, vertexCoordinates, coordinates);

                    for (int i = 0; i < noOfVertices; i++) {
                        for (int d = 0; d < D; d++) {
                            volumeCoordinates[e * noOfVertices + i, d] = coordinates[i, d];
                        }
                    }
                }

                vertexCoordinates = volumeCoordinates;
            }

            Debug.Assert(
                vertexCoordinates.GetLength(0) % 2 == 0,
                "Even number of vertices expected");
            int initialNumberOfLineSegments = vertexCoordinates.GetLength(0) / 2;

            List<LineSegment> lineSegments = new List<LineSegment>(initialNumberOfLineSegments);
            for (int i = 0; i < initialNumberOfLineSegments; i++) {
                var p0 = vertexCoordinates.GetRow(2 * i + 0);
                int iP0 = Simplex.Vertices.FindRow(p0, 1.0e-8);
                var p1 = vertexCoordinates.GetRow(2 * i + 1);
                int iP1 = Simplex.Vertices.FindRow(p1, 1.0e-8);

                LineSegment newSegment = new LineSegment(spatialDimension, Simplex,
                    p0, p1, iP0, iP1,
                    RootFindingAlgorithm);
                //for (int d = 0; d < spatialDimension; d++) {
                //    newSegment.Start[d] = vertexCoordinates[, d];
                //    newSegment.End[d] = vertexCoordinates[2 * i + 1, d];
                //}

                if (!lineSegments.Contains(newSegment)) {
                    lineSegments.Add(newSegment);
                    //tracker.Subscribe(newSegment);
                }
            }

            foreach (LineSegment segment in lineSegments) {
                LevelSet levelSetField = levelSetData.LevelSet as LevelSet;
                if (levelSetField != null) {
                    segment.ProjectBasisPolynomials(levelSetField.Basis);
                }
            }


            // in 2D, permute the line segments so that they form a continuous line-stroke
            // (required to obtain a 'stable' point measure)

            if (spatialDimension == 2) {
                var sotierung = new List<int>();
                sotierung.Add(0);

                for (int i = 1; i < lineSegments.Count; i++) {
                    var sortSeg = lineSegments[Math.Abs(sotierung[i-1])];
                    int Soll_iVtxStart = sotierung[i-1] >= 0 ? sortSeg.iVertexEnd : sortSeg.iVertexStart;
                    int Soll_iVtxEnd = (i == lineSegments.Count - 1) ? (lineSegments[sotierung[0]].iVertexStart) : int.MinValue;

                    int NoOfhits = 0;
                    int j = -1;
                    foreach (var segment in lineSegments) {
                        j++;
                        if (sotierung
                            .Contains(segment, delegate(int a, LineSegment b) {
                            LineSegment _a = lineSegments[Math.Abs(a)];
                            if (a >= 0) {
                                return (_a.iVertexStart == b.iVertexStart && _a.iVertexEnd == b.iVertexEnd);
                            } else {
                                return (_a.iVertexEnd == b.iVertexStart && _a.iVertexStart == b.iVertexEnd);
                            }
                        }))
                            continue;

                        if (Soll_iVtxEnd < 0) {

                            if (segment.iVertexStart == Soll_iVtxStart) {
                                sotierung.Add(j);
                                NoOfhits++;
                            }

                            if (segment.iVertexEnd == Soll_iVtxStart) {
                                var si = segment.Inverse;
                                //tracker.Subscribe(si);
                                sotierung.Add(-j);
                                NoOfhits++;
                            }

                        } else {
                            if (segment.iVertexStart == Soll_iVtxStart && segment.iVertexEnd == Soll_iVtxEnd) {
                                sotierung.Add(j);
                                NoOfhits++;
                            }

                            if (segment.iVertexEnd == Soll_iVtxStart && segment.iVertexStart == Soll_iVtxEnd) {
                                sotierung.Add(-j);
                                NoOfhits++;
                            }
                        }
                    }
                    if (NoOfhits != 1)
                        throw new ApplicationException("error in algorithm");

                }

#if DEBUG
                for (int i = 0; i < lineSegments.Count; i++) {
                    var seg_i = lineSegments[Math.Abs(sotierung[i])];
                    var seg_in = lineSegments[Math.Abs(sotierung[(i+1)%lineSegments.Count])];

                    int iVtx1 = (sotierung[i] >= 0) ? seg_i.iVertexEnd : seg_i.iVertexStart;
                    int iVtx2 = (sotierung[(i+1)%lineSegments.Count] >= 0) ? seg_in.iVertexStart : seg_in.iVertexEnd;
                    Debug.Assert(iVtx1 == iVtx2);
                }
#endif

                segmentSort = sotierung.ToArray();
            } else {
                segmentSort = null;
            }
            
            return lineSegments.ToArray();
        }

    }
}
