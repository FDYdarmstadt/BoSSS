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

using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid;
using System.Diagnostics;

namespace BoSSS.Foundation.Caching {

    /// <summary>
    /// Caches monomials ( e.g. \f$ x \f$, \f$ x y \f$, \f$ x y^2 \f$) for node sets. 
    /// A global singleton, see <see cref="Instance"/>.
    /// </summary>
    public class MonomialCache : CacheLogic_NsP {
        private MonomialCache() {

        }

        private static MonomialCache instance = null;
        private static readonly object padlock = new object();
        
        /// <summary>
        /// Access to the single, global instance.
        /// </summary>
        public static MonomialCache Instance {
            get {
                lock(padlock) {
                    if(instance == null) {
                        instance = new MonomialCache();
                    }
                    return instance;
                }
            }
        }

        /// <summary>
        /// Un-cached evaluation of monomials at nodes <paramref name="Nodes"/>
        /// up to degree <paramref name="Degree"/>.
        /// </summary>
        protected override MultidimensionalArray Evaluate(NodeSet Nodes, int Degree) {
            return GetMonomialsImpl(Nodes, Degree);
        }

        /// <summary>
        /// Evaluates the monomials at <paramref name="NS"/> up to the
        /// <paramref name="degree"/> of this polynomial
        /// </summary>
        /// <param name="NS">
        /// Points at which the monomials should be evaluated;
        /// 1st index: Point index;
        /// 2nd index: Spatial coordinate index;
        /// </param>
        /// <param name="degree">
        /// Required degree.
        /// </param>
        /// <returns>
        /// Evaluated monomials at <paramref name="NS"/>;
        /// 1st index: Point index;
        /// 2nd index: Spatial coordinate index;
        /// 3rd index: Exponent
        /// If monomials are already cached for degree higher than <paramref name="degree"/>, these are returned, so the is 
        /// range of the third index may be higher than expected.
        /// </returns>
        public MultidimensionalArray GetMonomials(NodeSet NS, int degree) {
            return this.Evaluate(NS, degree);

        }

                
        static MultidimensionalArray GetMonomialsImpl(MultidimensionalArray LocalPoints, int degree) {
            if(LocalPoints.Dimension != 2)
                throw new ArgumentException();
            int N = LocalPoints.GetLength(0);
            int spatialDimension = LocalPoints.GetLength(1);

            MultidimensionalArray monomials = MultidimensionalArray.Create(N, spatialDimension, degree + 1);
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < spatialDimension; j++) {
                    monomials[i, j, 0] = 1.0;
                    double xk = LocalPoints[i, j];
                    for(int k = 0; k < degree; k++) {
                        monomials[i, j, k + 1] = monomials[i, j, k] * xk;
                    }
                }
            }

            return monomials;
        }

    }
    
    /// <summary>
    /// Cache-Logic for values which depend on <see cref="NodeSet"/> and polynomial degree.
    /// </summary>
    public abstract class CacheLogic_NsP {


        /// <summary>
        /// Evaluation of Values that should subsequently be cached.
        /// </summary>
        /// <param name="Nodes">
        /// Points/Nodes to evaluate in the reference coordinate system;
        /// <list type="bullet">
        ///   <item>1st index: node index</item>
        ///   <item>2nd index: spatial dimension; (0 for 1D; 0,1 for 2D; 0,1,2 for 3D)</item>
        /// </list>
        /// </param>
        /// <param name="Degree">polynomial degree</param>
        /// <returns></returns>
        protected abstract MultidimensionalArray Evaluate(NodeSet Nodes, int Degree);

        /// <summary>
        /// Returns the (hopefully) cached value.
        /// </summary>
        public virtual MultidimensionalArray GetValues(NodeSet NS, int Degree) {
            
            // lookup table garbage collection, if necessary
            // ---------------------------------------------
            if(CacheLookup.Count > Cache.NoOfUsedBanks) {
                // definitely to many items in the lockup-table

                List<int> ToRemove = new List<int>();
                foreach(var kv in this.CacheLookup) {
                    if(!Cache.IsAlive(kv.Value.CacheRef))
                        ToRemove.Add(kv.Key);
                }

                foreach(int nsr in ToRemove) {
                    CacheLookup.Remove(nsr);
                    if(nsr == lastNodesetRef)
                        lastNodesetRef = -1;
                }
            }

            // transform node-set reference into cache reference
            // -------------------------------------------------
            Stuple CacheRef;
            bool present = false;
            if(this.lastNodesetRef == NS.Reference) {
                // lucky punch: we don't need to ask the dict
                CacheRef = this.lastCacheRef;
                present = true;
                Debug.Assert(CacheLookup.ContainsKey(lastNodesetRef));
            } else {
                present = this.CacheLookup.TryGetValue(NS.Reference, out CacheRef);
                if(!present) {
                    CacheRef.CacheRef = 0;
                    CacheRef.Degree = -1;
                }
            }

            // try to get value from cache
            // ---------------------------
            MultidimensionalArray vals;
            if(CacheRef.CacheRef <= 0) {
                // value not present in cache -> recompute
                vals = null;
            } else {
                vals = (MultidimensionalArray) Cache.GetItem(CacheRef.CacheRef);
            }

            // check the degree
            // ----------------
            bool reCache = false;
            if(CacheRef.Degree < Degree) {
                // unfortunately, insufficient degree
                reCache = true;
                vals = null;
            }
            
            // return
            // ------

            if(vals != null) {
                //return vals;
            } else {
                Debug.Assert(vals == null || Degree >= CacheRef.Degree);
                vals = this.Evaluate(NS, Degree);
                if(reCache) {
                    CacheRef.CacheRef = Cache.ReCacheItem(vals, vals.Length*sizeof(double), CacheRef.CacheRef);
                } else {
                    CacheRef.CacheRef = Cache.CacheItem(vals, vals.Length * sizeof(double));
                }
                CacheRef.Degree = Degree;
                
                if(present)
                    CacheLookup[NS.Reference] = CacheRef;
                else
                    CacheLookup.Add(NS.Reference, CacheRef);
            }

            this.lastNodesetRef = NS.Reference;
            this.lastCacheRef = CacheRef;
            Debug.Assert(object.ReferenceEquals(vals, Cache.GetItem(CacheRef.CacheRef)));

            return vals;
        }

        int lastNodesetRef = -1;
        Stuple lastCacheRef;
        

        struct Stuple {
            public ulong CacheRef;
            public int Degree;
        }

        Dictionary<int,Stuple> CacheLookup = new Dictionary<int, Stuple>();
       
    }

    /// <summary>
    ///  Delegate-supporting version of <see cref="CacheLogic_Ns"/>.
    /// </summary>
    public class CacheLogicImpl_Ns : CacheLogic_Ns {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="evalfunc">
        ///  The evaluation function, see <see cref="CacheLogic_Ns.Evaluate"/>
        /// </param>
        public CacheLogicImpl_Ns(Func<NodeSet, MultidimensionalArray> evalfunc) {
            this.m_evalfunc = evalfunc;
        }

        Func<NodeSet, MultidimensionalArray> m_evalfunc;

        /// <summary>
        /// Evaluation fo Values if they are not in the cache.
        /// </summary>
        /// <param name="Nodes">node set</param>
        override protected MultidimensionalArray Evaluate(NodeSet Nodes) {
            return m_evalfunc(Nodes);
        }
     
    }

    /// <summary>
    /// Cache-Logic for values which depend on <see cref="NodeSet"/>.
    /// </summary>
    public abstract class CacheLogic_Ns {

        /// <summary>
        /// ctor.
        /// </summary>
        public CacheLogic_Ns() {
            m_internal = new CacheLogicImpl_NsP(Evaluate2);
        }

        CacheLogicImpl_NsP m_internal;

        MultidimensionalArray Evaluate2(NodeSet Nodes, int P) {
            Debug.Assert(P == 0);
            return this.Evaluate(Nodes);
        }


        /// <summary>
        /// Evaluation of Values that should subsequently be cached.
        /// </summary>
        /// <param name="Nodes">
        /// Points/Nodes to evaluate in the reference coordinate system;
        /// </param>
        /// <returns></returns>
        protected abstract MultidimensionalArray Evaluate(NodeSet Nodes);

        /// <summary>
        /// Returns the (hopefully) cached value.
        /// </summary>
        public virtual MultidimensionalArray GetValues(NodeSet N) {
            return m_internal.GetValues(N, 0);
        }
    }

    /// <summary>
    ///  Delegate-supporting version of <see cref="CacheLogicImpl_NsP"/>.
    /// </summary>
    public class CacheLogicImpl_NsP : CacheLogic_NsP {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="evalfunc">
        ///  The evaluation function, see <see cref="CacheLogicImpl_NsP.Evaluate(NodeSet, int)"/>
        /// </param>
        public CacheLogicImpl_NsP(Func<NodeSet, int, MultidimensionalArray> evalfunc) {
            this.m_evalfunc = evalfunc;
        }

        Func<NodeSet, int, MultidimensionalArray> m_evalfunc;

        /// <summary>
        /// Evaluation fo Values if they are not in the cache.
        /// </summary>
        /// <param name="Degree">minimum required degree.</param>
        /// <param name="Nodes">node set</param>
        override protected MultidimensionalArray Evaluate(NodeSet Nodes, int Degree) {
            return m_evalfunc(Nodes, Degree);
        }

    }


    /// <summary>
    /// Delegate-supporting version of <see cref="CacheLogic_CNs"/>.
    /// </summary>
    public class CacheLogicImplBy_CNs : CacheLogic_CNs {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="cv">Compute-value function, see <see cref="CacheLogic_CNs.ComputeValues"/>. </param>
        /// <param name="ca">Memory allocation function, see <see cref="CacheLogic_CNs.Allocate"/>. </param>
        public CacheLogicImplBy_CNs(Grid.Classic.GridData g, Action<NodeSet, int, int, MultidimensionalArray> cv, Func<int,int,NodeSet,MultidimensionalArray> ca)
            : base(g) //
        {
            m_ComputeValues = cv;
            m_alloc = ca;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="cv">Compute-value function, see <see cref="CacheLogic_CNs.ComputeValues"/>. </param>
        /// <param name="bl">Buffer length function.</param>
        public CacheLogicImplBy_CNs(Grid.Classic.GridData g, Action<NodeSet, int, int, MultidimensionalArray> cv, Func<int, int, int[]> bl)
            : base(g) //
        {
            m_ComputeValues = cv;
            m_bufferlengths = bl;
            m_alloc = this.FFReAllocate;
        }

        Action<NodeSet, int, int, MultidimensionalArray> m_ComputeValues;

        Func<int,int,int[]> m_bufferlengths;
        Func<int,int,NodeSet,MultidimensionalArray> m_alloc;


        MultidimensionalArray FFReAllocate(int j0, int Len, NodeSet N) {
            int[] Lengths = this.m_bufferlengths(Len, N.NoOfNodes);
            var buffer = new MultidimensionalArray(Lengths.Length);
            buffer.Allocate(Lengths);
            return buffer;
        }


        /// <summary>
        /// Uses the delegate handed to constructor.
        /// </summary>
        protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
            this.m_ComputeValues(N, j0, Len, output);
        }

        /// <summary>
        /// Uses the delegate handed to constructor.
        /// </summary>
        protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
            return this.m_alloc(i0, Len, NS);
        }
    }

    /// <summary>
    /// Cache-logic for data which depends on cell, <see cref="NodeSet"/> and polynomial degree
    /// </summary>
    public class CacheLogicImpl_CNsP  {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="cv">Compute-value function. </param>
        /// <param name="ca">Memory allocation function. </param>
        public CacheLogicImpl_CNsP(IGridData g, Action<NodeSet, int, int, int, MultidimensionalArray> cv, Func<int, int, int, NodeSet, MultidimensionalArray> ca) {
            m_ComputeValues = cv;
            m_alloc = ca;
            this.GridData = g;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="cv">Compute-value function, see <see cref="CacheLogic_CNs.ComputeValues"/>. </param>
        /// <param name="bl">Buffer length function.</param>
        public CacheLogicImpl_CNsP(IGridData g, Action<NodeSet, int, int, int, MultidimensionalArray> cv, Func<int, int, int, int[]> bl) {
            m_ComputeValues = cv;
            m_bufferlengths = bl;
            m_alloc = this.FFReAllocate;
            this.GridData = g;
        }

        Action<NodeSet, int, int, int, MultidimensionalArray> m_ComputeValues;

        Func<int,int,int,int[]> m_bufferlengths;
        Func<int,int, int,NodeSet,MultidimensionalArray> m_alloc;
        
        MultidimensionalArray FFReAllocate(int j0, int Len, int Degree, NodeSet N) {
            int[] Lengths = this.m_bufferlengths(Len, N.NoOfNodes, Degree);
            var buffer = new MultidimensionalArray(Lengths.Length);
            buffer.Allocate(Lengths);
            return buffer;
        }


        // lucky-punch cache
        int lastCell_j0 = -1;
        int lastCell_Len = -1;
        int lastCell_Degree = -1;
        int lastCell_NSref = -1;
        ulong lastCell_CashRef = 0;


        /// <summary>
        /// returns the values for node set <paramref name="NS"/>,
        /// in the cells (with local index) <paramref name="j0"/> (including)
        /// to <paramref name="j0"/>+<paramref name="Len"/> (excluding).
        /// </summary>
        /// <param name="NS">the node set.</param>
        /// <param name="j0">local cell index of the first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <param name="Degree">Polynomial degree</param>
        /// <returns></returns>
        /// <remarks>
        /// This method returns, if available, cached values  or it re-computes them by
        /// calling <see cref="m_ComputeValues"/>.
        /// </remarks>
        public MultidimensionalArray GetValue_Cell(NodeSet NS, int j0, int Len, int Degree) {
            if(j0 < 0) throw new ArgumentOutOfRangeException("j0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((j0 + Len) > this.GridData.iGeomCells.NoOfCells)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of local cells");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.CellCoord)
                throw new ArgumentException("Expecting a node set in local cell coordinate system");
#if DEBUG
            int iKref = NS.GetVolumeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not a volume node set");
            for(int j = 0; j < Len; j++) {
                if(GridData.iGeomCells.GetRefElementIndex(j + j0) != iKref)
                    throw new ArgumentException("node set ref. element/cell mismatch");
            }
            Debug.Assert(NS.SpatialDimension == this.GridData.SpatialDimension,
                "Mismatch between number of spatial directions in node set and reference element.");
#endif

            MultidimensionalArray R = null;
            if(lastCell_j0 == j0 && lastCell_Len == Len && lastCell_Degree == Degree && lastCell_NSref == NS.Reference) {
                R = (MultidimensionalArray)Cache.GetItem(lastCell_CashRef);
                Debug.Assert((R == null) || (R.GetLength(0) == Len));
                Debug.Assert((R == null) || (R.GetLength(1) == NS.NoOfNodes));
            }

            lastCell_j0 = j0;
            lastCell_Len = Len;
            lastCell_Degree = Degree;
            lastCell_NSref = NS.Reference;

            if(R == null) {
                R = this.m_alloc(j0, Len, Degree, NS);
                this.m_ComputeValues(NS, j0, Len, Degree, R);
                Debug.Assert(R.GetLength(0) == Len);
                Debug.Assert(R.GetLength(1) == NS.NoOfNodes);
                lastCell_CashRef = Cache.CacheItem(R, R.Length * sizeof(double));
            }
            return R;
        }

        public MultidimensionalArray GetValue_EdgeSV(NodeSet NS, int e0, int Len, int Degree) {
            if(e0 < 0) throw new ArgumentOutOfRangeException("e0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((e0 + Len) > this.GridData.iGeomEdges.Count)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of edges");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentException("Expecting a node set in edge coordinate system");
#if DEBUG
            int iKref = NS.GetEdgeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not an edge node set");
            for(int j = 0; j < Len; j++) {
                if(this.GridData.iGeomEdges.GetRefElementIndex(j + e0) != iKref)
                    throw new ArgumentException("node set ref. element/edge missmatch");
            }
#endif

            MultidimensionalArray Rin = null;
            if(lastEdge_e0 == e0 && lastEdge_Len == Len && lastEdge_Degree == Degree && lastEdge_NSref == NS.Reference) {
                Rin = (MultidimensionalArray)Cache.GetItem(lastEdge_CashRefIn);
                
                Debug.Assert((Rin == null) || (Rin.GetLength(0) == Len));
                Debug.Assert((Rin == null) || (Rin.GetLength(1) == NS.NoOfNodes));
            }
            
            bool RinReq = Rin == null;

            lastEdge_e0 = e0;
            lastEdge_Len = Len;
            lastEdge_Degree = Degree;
            lastEdge_NSref = NS.Reference;
            
            if(RinReq) {
                Rin = m_alloc(e0, Len, Degree, NS);
                Debug.Assert(Rin.GetLength(0) == Len);
                Debug.Assert(Rin.GetLength(1) == NS.NoOfNodes);

            
                int[] i0 = new int[Rin.Dimension];
                int[] iE = new int[i0.Length];
                for(int k = 2; k < i0.Length; k++) {
                    iE[k] = Rin.GetLength(k) - 1;
                }
                iE[1] = NS.NoOfNodes - 1;

                int[,] E2C = this.GridData.iLogicalEdges.CellIndices;
                int[,] TrIdx = this.GridData.iGeomEdges.Edge2CellTrafoIndex;

                for(int e = 0; e < Len; e++) {
                    int iEdge = e + e0;
                    int jCell1, edge1;
                    jCell1 = E2C[iEdge, 0];
                    edge1 = TrIdx[iEdge, 0];

                    i0[0] = e;
                    iE[0] = e;

                    this.m_ComputeValues(NS.GetVolumeNodeSet(this.GridData, edge1), jCell1, 1, Degree, Rin.ExtractSubArrayShallow(i0, iE));
                }

                lastEdge_CashRefIn = Cache.CacheItem(Rin, Rin.Length * sizeof(double));
            }

            return Rin;
        }


        int lastEdge_e0 = -1;
        int lastEdge_Len = -1;
        int lastEdge_Degree = -1;
        int lastEdge_NSref = -1;
        ulong lastEdge_CashRefIn = 0;
        ulong lastEdge_CashRefOt = 0;


        /// <summary>
        /// Double-value (i.e. values for IN- and OUT-cell) evaluation at edges,
        /// used for properties which ar discontinuous at edges.
        /// </summary>
        public Tuple<MultidimensionalArray, MultidimensionalArray> GetValue_EdgeDV(NodeSet NS, int e0, int Len, int Degree) {
            if(e0 < 0) throw new ArgumentOutOfRangeException("e0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((e0 + Len) > this.GridData.iGeomEdges.Count)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of edges");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentException("Expecting a node set in edge coordinate system");
#if DEBUG
            int iKref = NS.GetEdgeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not an edge node set");
            for(int j = 0; j < Len; j++) {
                if(this.GridData.iGeomEdges.GetRefElementIndex(j + e0) != iKref)
                    throw new ArgumentException("node set ref. element/edge missmatch");
            }
#endif

            MultidimensionalArray Rin = null, Rot = null;

            if(lastEdge_e0 == e0 && lastEdge_Len == Len && lastEdge_Degree == Degree && lastEdge_NSref == NS.Reference) {
                Rin = (MultidimensionalArray)Cache.GetItem(lastEdge_CashRefIn);
                Rot = (MultidimensionalArray)Cache.GetItem(lastEdge_CashRefOt);

                Debug.Assert((Rin == null) || (Rin.GetLength(0) == Len));
                Debug.Assert((Rot == null) || (Rot.GetLength(0) == Len));
                Debug.Assert((Rin == null) || (Rin.GetLength(1) == NS.NoOfNodes));
                Debug.Assert((Rot == null) || (Rot.GetLength(1) == NS.NoOfNodes));
            }

            bool RinReq = Rin == null;
            bool RotReq = Rin == null;

            lastEdge_e0 = e0;
            lastEdge_Len = Len;
            lastEdge_Degree = Degree;
            lastEdge_NSref = NS.Reference;
            
            
            if(RinReq)
                Rin = this.m_alloc(e0, Len, Degree, NS);
            if(RotReq)
                Rot = this.m_alloc(e0, Len, Degree, NS);
            Debug.Assert(Rin.GetLength(0) == Len);
            Debug.Assert(Rot.GetLength(0) == Len);
            Debug.Assert(Rin.GetLength(1) == NS.NoOfNodes);
            Debug.Assert(Rot.GetLength(1) == NS.NoOfNodes);

            if(RinReq || RotReq) {
                int[] i0 = new int[Rin.Dimension];
                int[] iE = new int[i0.Length];
                for(int k = 2; k < i0.Length; k++) {
                    iE[k] = Rin.GetLength(k) - 1;
                    Debug.Assert(Rin.GetLength(k) == Rot.GetLength(k));
                }
                iE[1] = NS.NoOfNodes - 1;

                int[,] E2C = this.GridData.iLogicalEdges.CellIndices;
                int[,] TrIdx = this.GridData.iGeomEdges.Edge2CellTrafoIndex;

                for(int e = 0; e < Len; e++) {
                    int iEdge = e + e0;
                    int jCell1, jCell2, edge1, edge2;
                    jCell1 = E2C[iEdge, 0];
                    jCell2 = E2C[iEdge, 1];
                    edge1 = TrIdx[iEdge, 0];
                    edge2 = TrIdx[iEdge, 1];

                    i0[0] = e;
                    iE[0] = e;

                    if(RinReq) {
                        this.m_ComputeValues(NS.GetVolumeNodeSet(this.GridData, edge1), jCell1, 1, Degree, Rin.ExtractSubArrayShallow(i0, iE));
                    }
                    if(RotReq) {
                        if(jCell2 >= 0) {
                            this.m_ComputeValues(NS.GetVolumeNodeSet(this.GridData, edge2), jCell2, 1, Degree, Rot.ExtractSubArrayShallow(i0, iE));
                        } else {
                            Rot.ExtractSubArrayShallow(i0, iE).Clear();
                        }
                    }
                }

                if(RinReq)
                    lastEdge_CashRefIn = Cache.CacheItem(Rin, Rin.Length * sizeof(double));
                if(RotReq)
                    lastEdge_CashRefOt = Cache.CacheItem(Rot, Rot.Length * sizeof(double));
            }

            return new Tuple<MultidimensionalArray, MultidimensionalArray>(Rin, Rot);
        }

        /// <summary>
        /// Grid reference.
        /// </summary>
        protected IGridData GridData {
            get;
            private set;
        }

        /// <summary>
        /// Disposes all chached values.
        /// </summary>
        public void Clear() {
            lastEdge_e0 = -1;
            lastEdge_Len = -1;
            lastEdge_Degree = -1;
            lastEdge_NSref = -1;
            lastEdge_CashRefIn = 0;
            lastEdge_CashRefOt = 0;

            lastCell_j0 = -1;
            lastCell_Len = -1;
            lastCell_Degree = -1;
            lastCell_NSref = -1;
            lastCell_CashRef = 0;
        }
    }
    

    /// <summary>
    /// Cache-logic for data which depends on cell and <see cref="NodeSet"/>.
    /// </summary>
    public abstract class CacheLogic_CNs {

        /// <summary>
        /// constructor
        /// </summary>
        protected CacheLogic_CNs(IGridData g) {
            this.GridData = g;
            m_impl = new CacheLogicImpl_CNsP(g, this.__ComputeValues, this.__Allocate);
        }

        CacheLogicImpl_CNsP m_impl;

        void __ComputeValues(NodeSet N, int j0, int Len, int Degree, MultidimensionalArray output) {
            Debug.Assert(Degree == 0);
            this.ComputeValues(N, j0, Len, output);
        }

        MultidimensionalArray __Allocate(int i0, int Len, int Degree, NodeSet NS) {
            Debug.Assert(Degree == 0);
            return this.Allocate(i0, Len, NS);
        }


        /// <summary>
        /// Implementers should override this method for the computation of values which are not in the cache.
        /// </summary>
        /// <param name="N">Nodes at which evaluation should be done.</param>
        /// <param name="j0">local cell index of the first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <param name="output">
        /// reference to storage buffer for result
        /// </param>
        protected abstract void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output);

        /// <summary>
        /// Allocates memory to store results from <see cref="ComputeValues"/>.
        /// </summary>
        /// <param name="i0">Index of first item (edge or cell).</param>
        /// <param name="Len">Number of items.</param>
        /// <param name="NS">Node set.</param>
        /// <returns>
        /// Some buffer, 
        /// 1st length is expected to be equal to <paramref name="Len"/>,
        /// 2nd length is expected to be equal to number of nodes in <paramref name="NS"/>.
        /// </returns>
        protected abstract MultidimensionalArray Allocate(int i0, int Len, NodeSet NS);
                
        /// <summary>
        /// returns the values for nodeset <paramref name="NS"/>,
        /// in the cells (with local index) <paramref name="j0"/> (including)
        /// to <paramref name="j0"/>+<paramref name="Len"/> (excluding).
        /// </summary>
        /// <param name="NS">the node set.</param>
        /// <param name="j0">local cell index of the first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns></returns>
        /// <remarks>
        /// This method returns, if available, cached values  or it re-computes them by
        /// calling <see cref="ComputeValues"/>.
        /// </remarks>
        public MultidimensionalArray GetValue_Cell(NodeSet NS, int j0, int Len) {
            if(j0 < 0) throw new ArgumentOutOfRangeException("j0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((j0 + Len) > this.GridData.iLogicalCells.NoOfCells)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of local cells");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.CellCoord)
                throw new ArgumentException("Expecting a node set in local cell coordinate system");
#if DEBUG
            int iKref = NS.GetVolumeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not a volume node set");
            for(int j = 0; j < Len; j++) {
                if(GridData.iGeomCells.GetRefElementIndex(j + j0) != iKref)
                    throw new ArgumentException("node set ref. element/cell missmatch");
            }
            Debug.Assert(NS.SpatialDimension == this.GridData.SpatialDimension,
                "Mismatch between number of spatial directions in node set and reference element.");
#endif

            return this.m_impl.GetValue_Cell(NS, j0, Len, 0);
        }


        /// <summary>
        /// Single-value evaluation at edges,
        /// used for properties which are unique at edges (e.g. edge normals or global coordinates).
        /// </summary>
        public MultidimensionalArray GetValue_EdgeSV(NodeSet NS, int e0, int Len) {
            if(e0 < 0) throw new ArgumentOutOfRangeException("e0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((e0 + Len) > this.GridData.iLogicalEdges.Count)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of edges");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentException("Expecting a node set in edge coordinate system");
#if DEBUG
            int iKref = NS.GetEdgeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not an edge node set");
            for(int j = 0; j < Len; j++) {
                if(this.GridData.iGeomEdges.GetRefElementIndex(j + e0) != iKref)
                    throw new ArgumentException("node set ref. element/edge missmatch");
            }
#endif

            return this.m_impl.GetValue_EdgeSV(NS, e0, Len, 0);
        }

        /// <summary>
        /// Double-value (i.e. values for IN- and OUT-cell) evaluation at edges,
        /// used for properties which are discontinuous at edges (e.g. DG-fields).
        /// </summary>
        public Tuple<MultidimensionalArray, MultidimensionalArray> GetValue_EdgeDV(NodeSet NS, int e0, int Len) {
            if(e0 < 0) throw new ArgumentOutOfRangeException("e0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((e0 + Len) > this.GridData.iLogicalEdges.Count)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of edges");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentException("Expecting a node set in edge coordinate system");
#if DEBUG
            int iKref = NS.GetEdgeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not an edge node set");
            for(int j = 0; j < Len; j++) {
                if(this.GridData.iGeomEdges.GetRefElementIndex(j + e0) != iKref)
                    throw new ArgumentException("node set ref. element/edge missmatch");
            }
#endif

            return m_impl.GetValue_EdgeDV(NS, e0, Len, 0);
        }

        /// <summary>
        /// The grid; required since one key into the cache logic is the cell index.
        /// </summary>
        protected IGridData GridData {
            get;
            private set;
        }

        /// <summary>
        /// Disposes all chached values.
        /// </summary>
        public void Clear() {
            m_impl.Clear();
        }
    }

    class CorEP_logic {
        public CorEP_logic(int _JE, Func<int, int, int, MultidimensionalArray> _ComputeValues) {
            m_JE = _JE;
            this.m_ComputeValues = _ComputeValues;
        }


        int NoOfCachedCells = 0;
        int NoOfChunks = 0;

        CChunk[] cells2cchunk = null;


        public MultidimensionalArray GetValue(int j0, int Len, int degree) {
            int JE = this.m_JE;

            if(cells2cchunk == null) {

                MultidimensionalArray R = m_ComputeValues(j0, Len, degree);
                ulong cacheRef = Cache.CacheItem(R, R.Length * sizeof(double));

                cells2cchunk = new CChunk[JE];
                CChunk cc = new CChunk() {
                    i0 = j0,
                    iE = j0 + Len - 1,
                    Degree = degree,
                    CacheRef = cacheRef
                };

                for(int j = cc.i0; j <= cc.iE; j++) {
                    cells2cchunk[j] = cc;
                }

                Debug.Assert(NoOfChunks == 0);
                Debug.Assert(NoOfCachedCells == 0);
                NoOfCachedCells = Len;
                NoOfChunks = 1;

                Verify();
                return R;
            }

            Verify();

            CChunk fstCC = cells2cchunk[j0];
            if(fstCC != null && fstCC.i0 <= j0 && fstCC.iE >= (j0 + Len - 1)) {
                // desired result is maybe completly contained in first chunk

                MultidimensionalArray R;
                if(fstCC.Degree >= degree) {
                    R = (MultidimensionalArray)Cache.GetItem(fstCC.CacheRef);
                    Debug.Assert((R == null) || (R.GetLength(0) == fstCC.Len));
                } else {
                    R = m_ComputeValues(fstCC.i0, fstCC.Len, degree);
                    Debug.Assert(R.GetLength(0) == fstCC.Len);
                    fstCC.CacheRef = Cache.ReCacheItem(R, R.Length * sizeof(double), fstCC.CacheRef);
                    fstCC.Degree = degree;
                }

                if(R == null) {
                    KillChunk(fstCC);

                    R = m_ComputeValues(j0, Len, degree);
                    CChunk cc = new CChunk() {
                        i0 = j0,
                        iE = j0 + Len - 1,
                        Degree = degree
                    };
                    cc.CacheRef = Cache.CacheItem(R, R.Length * sizeof(double));

                    NoOfChunks += 1;
                    NoOfCachedCells += cc.Len;

                    for(int j = cc.i0; j <= cc.iE; j++) {
                        cells2cchunk[j] = cc;
                    }
                    Verify();
                    return R;

                } else {
                    if(fstCC.i0 == j0 && fstCC.Len == Len) {
                        Verify();
                        return R;
                    } else {
                        int[] _i0 = new int[R.Dimension];
                        _i0[0] = j0 - fstCC.i0;
                        int[] _iE = R.Lengths;
                        for(int i = _iE.Length - 1; i > 0; i--) {
                            _iE[i]--;
                        }
                        _iE[0] = Len - 1 + _i0[0];

                        Verify();
                        return R.ExtractSubArrayShallow(_i0, _iE);
                    }

                }
            }

            int iE = j0 + Len - 1;
            List<int> foundchunks = null;
            int maxDegree = degree;
            for(int j = j0; j < JE; j++) {
                CChunk cc = cells2cchunk[j];
                if(cc != null) {
                    if(Cache.IsAlive(cc.CacheRef)) {
                        if(foundchunks == null)
                            foundchunks = new List<int>();
                        foundchunks.Add(j);
                        j = cc.iE;
                        maxDegree = Math.Max(maxDegree, cc.Degree);
                    } else {
                        Verify();
                        KillChunk(cc);
                    }
                }
            }

            if(foundchunks == null) {
                MultidimensionalArray R = m_ComputeValues(j0, Len, degree);
                CChunk cc = new CChunk() {
                    i0 = j0,
                    iE = j0 + Len - 1,
                    Degree = degree
                };
                cc.CacheRef = Cache.CacheItem(R, R.Length * sizeof(double));

                for(int j = cc.i0; j <= cc.iE; j++) {
                    cells2cchunk[j] = cc;
                }

                NoOfChunks += 1;
                NoOfCachedCells += cc.Len;
                Verify();
                return R;
            } else {
                CChunk cc = new CChunk() {
                    i0 = Math.Min(j0, cells2cchunk[foundchunks[0]].i0),
                    iE = Math.Max(j0 + Len - 1, cells2cchunk[foundchunks[foundchunks.Count - 1]].iE),
                    Degree = maxDegree
                };
                MultidimensionalArray R = m_ComputeValues(cc.i0, cc.Len, maxDegree);

                foreach(int j in foundchunks) {
                    CChunk kcc = cells2cchunk[j];
                    KillChunk(kcc);
                }

                cc.CacheRef = Cache.CacheItem(R, R.Length * sizeof(double));

                for(int j = cc.i0; j <= cc.iE; j++) {
                    cells2cchunk[j] = cc;
                }

                NoOfChunks += 1;
                NoOfCachedCells += cc.Len;

                if(cc.i0 == j0 && cc.Len == Len) {
                    Verify();
                    return R;
                } else {
                    int[] _i0 = new int[R.Dimension];
                    _i0[0] = j0 - cc.i0;
                    int[] _iE = R.Lengths;
                    for(int i = _iE.Length - 1; i > 0; i--) {
                        _iE[i]--;
                    }
                    _iE[0] = Len - 1 + _i0[0];

                    Verify();
                    return R.ExtractSubArrayShallow(_i0, _iE);
                }
            }
        }

        [Conditional("DEBUG")]
        void Verify() {
#if DEBUG
            int JE = cells2cchunk.Length;
            for(int j = 0; j < JE; j++) {
                CChunk cc = cells2cchunk[j];
                if(cc != null) {
                    for(; j < cc.iE; j++) {
                        Debug.Assert(object.ReferenceEquals(cc, cells2cchunk[j]));
                    }
                    j = cc.iE;


                    MultidimensionalArray stuff = (MultidimensionalArray)Cache.GetItem(cc.CacheRef);
                    Debug.Assert((stuff == null) || (stuff.GetLength(0) == cc.Len));
                }


            }
#endif
        }



        /// <summary>
        /// Removes cache refences to objects that were disposed from the <see cref="Cache"/>.
        /// </summary>
        /// <returns>
        /// True if nothing is in the cache anymore;
        /// </returns>
        public bool GarbageCollect() {

            int JE = cells2cchunk.Length;

            for(int j = 0; j < JE; j++) {
                CChunk cc = cells2cchunk[j];
                if(cc != null) {
                    if(!Cache.IsAlive(cc.CacheRef)) {
                        // chunk is not present in cache anymore

                        KillChunk(cc);
                    }
                    j = cc.iE;
                }
            }

            if(NoOfChunks == 0)
                cells2cchunk = null;

            return (NoOfChunks == 0);
        }

        private void KillChunk(CChunk cc) {
            int iE = cc.iE;
            for(int i = cc.i0; i <= iE; i++) {
                Debug.Assert(object.ReferenceEquals(cells2cchunk[i], cc));
                cells2cchunk[i] = null;
            }

            NoOfChunks--;
            NoOfCachedCells -= (cc.iE - cc.i0 + 1);

            Debug.Assert(NoOfChunks >= 0);
            Debug.Assert(NoOfCachedCells >= 0);

            //if(cc.prev != null)
            //    cc.prev.next = cc.next;
            //if(cc.next != null)
            //    cc.next.prev = cc.prev;
        }

        class CChunk {
            public int i0;
            public int iE;
            public int Len {
                get {
                    Debug.Assert(iE >= i0);
                    return (iE - i0 + 1);
                }
            }
            public int Degree;
            public ulong CacheRef;

            //public CChunk next;
            //public CChunk prev;
        }
  

        int m_JE;

        protected Func<int,int,int,MultidimensionalArray> m_ComputeValues;
    }

    /// <summary>
    /// Cache logic for values which depend on cell and polynomial degree.
    /// </summary>
    public class CacheLogicImpl_CP {

        /// <summary>
        /// constructor
        /// </summary>
        public CacheLogicImpl_CP(IGridData g, Func<int, int, int, MultidimensionalArray> _ComputeValues) {
            this.GridData = g;
            this.m_impl = new CorEP_logic(g.iGeomCells.NoOfCells, _ComputeValues);
        }


        CorEP_logic m_impl;

        
        /// <summary>
        /// returns the values in the cells (with local index)
        /// <paramref name="j0"/> (including) to
        /// <paramref name="j0"/>+<paramref name="Len"/> (excluding).
        /// </summary>
        /// <param name="j0">local cell index of the first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <param name="degree">Polynomial degree</param>
        /// <returns></returns>
        /// <remarks>
        /// This method returns, if available, cached values  or it re-computes them by
        /// calling <see cref="CorEP_logic.GetValue(int, int, int)"/>.
        /// </remarks>
        public MultidimensionalArray GetValue_Cell(int j0, int Len, int degree) {
            int JE = this.GridData.iGeomCells.NoOfCells;
            if(j0 < 0) throw new ArgumentOutOfRangeException("j0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((j0 + Len) > JE)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of local cells");

#if DEBUG
            int iKref = this.GridData.iGeomCells.GetRefElementIndex(j0);
            if(iKref < 0)
                throw new ArgumentException("not a volume node set");
            for(int j = 0; j < Len; j++) {
                if(GridData.iGeomCells.GetRefElementIndex(j + j0) != iKref)
                    throw new ArgumentException("Node set ref. element/cell mismatch.");
            }
#endif

            return m_impl.GetValue(j0, Len, degree);
        }

      
        protected IGridData GridData {
            get;
            private set;
        }

    }

    /// <summary>
    /// Cache logic for values which depend on cell.
    /// </summary>
    public class CacheLogicImpl_C {

        /// <summary>
        /// constructor
        /// </summary>
        public CacheLogicImpl_C(Grid.Classic.GridData g, Del_ComputeValues _ComputeValues) {
            this.GridData = g;
            this.m_ComputeValues = _ComputeValues;
            this.m_impl = new CacheLogicImpl_CP(g, this.RedictComputeValues);
        }

        MultidimensionalArray RedictComputeValues(int j0, int Len, int Degree) {
            Debug.Assert(Degree == 0);
            return m_ComputeValues(j0, Len);
        }
        
        CacheLogicImpl_CP m_impl;

        
        /// <summary>
        /// Implementers should override this method for the computation of
        /// values which are not in the cache.
        /// </summary>
        /// <param name="j0">local cell index of the first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns></returns>
        public delegate MultidimensionalArray Del_ComputeValues(int j0, int Len);


        Del_ComputeValues m_ComputeValues;

        /// <summary>
        /// returns the values in the cells (with local index)
        /// <paramref name="j0"/> (including) to
        /// <paramref name="j0"/>+<paramref name="Len"/> (excluding).
        /// </summary>
        /// <param name="j0">local cell index of the first cell to evaluate</param>
        /// <param name="Len">number of cells to evaluate</param>
        /// <returns></returns>
        /// <remarks>
        /// This method returns, if available, cached values  or it re-computes them by
        /// calling <see cref="CacheLogicImpl_CP.GetValue_Cell(int, int, int)"/>.
        /// </remarks>
        public MultidimensionalArray GetValue_Cell(int j0, int Len) {
            if(j0 < 0) throw new ArgumentOutOfRangeException("j0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((j0 + Len) > this.GridData.Cells.NoOfCells)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of local cells");

#if DEBUG
            int iKref = this.GridData.Cells.GetRefElementIndex(j0);
            if(iKref < 0)
                throw new ArgumentException("not a volume node set");
            for(int j = 0; j < Len; j++) {
                if(GridData.Cells.GetRefElementIndex(j + j0) != iKref)
                    throw new ArgumentException("node set ref. element/cell missmatch");
            }
#endif

            return m_impl.GetValue_Cell(j0, Len, 0);
        }

        /// <summary>
        /// The current grid
        /// </summary>
        protected Grid.Classic.GridData GridData {
            get;
            private set;
        }

    }

    /// <summary>
    /// Cache-logic for data which depends on cell, <see cref="NodeSet"/> and cell-face.
    /// This one is only used for edge normals and edge integration metrics.
    /// </summary>
    public class EdgeNormalsCacheLogic_CNsFace {

        /// <summary>
        /// constructor
        /// </summary>
        public EdgeNormalsCacheLogic_CNsFace(Grid.Classic.GridData g) {
            this.GridData = g;
        }


        int last_e0 = -1;
        int last_Len = -1;
        int last_NSref = -1;
        MultidimensionalArray last_Normals;
        MultidimensionalArray last_metrics;


        /// <summary>
        /// Evaluation of a single-valued property on edges.
        /// </summary>
        public MultidimensionalArray GetNormals_Edge(NodeSet NS, int e0, int Len) {
            if(e0 == last_e0 && Len == last_Len && NS.Reference == last_NSref) {
                // lucky punch successful
                Debug.Assert(last_Normals.GetLength(0) == Len);
                Debug.Assert(last_Normals.GetLength(1) == NS.NoOfNodes);
            } else {
                EdgeEval(NS, e0, Len, out last_Normals, out last_metrics);
                last_e0 = e0;
                last_Len = Len;
                last_NSref = NS.Reference;
            }
            return last_Normals;
        }

        private void EdgeEval(NodeSet NS, int e0, int Len, out MultidimensionalArray Normals, out MultidimensionalArray Metrics) {
            if(e0 < 0) throw new ArgumentOutOfRangeException("e0", "must be greater or equal than zero.");
            if(Len < 1) throw new ArgumentOutOfRangeException("Len", "must be greater or equal than one.");
            if((e0 + Len) > this.GridData.Edges.Count)
                throw new ArgumentOutOfRangeException("j0 + Len exceeds the number of edges");
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentException("Expecting a node set in edge coordinate system");
            bool affineEdge = this.GridData.Edges.IsEdgeAffineLinear(e0);
#if DEBUG
            int iKref = NS.GetEdgeRefElementIndex(this.GridData);
            if(iKref < 0)
                throw new ArgumentException("not an edge node set");
            for(int j = 0; j < Len; j++) {
                if(this.GridData.Edges.GetRefElementIndex(j + e0) != iKref)
                    throw new ArgumentException("node set ref. element/edge missmatch");
                Debug.Assert(affineEdge == this.GridData.Edges.IsEdgeAffineLinear(j + e0));
            }
#endif
            int D = this.GridData.SpatialDimension;
            int NoOfNodes = NS.NoOfNodes;


            Normals = MultidimensionalArray.Create(Len, NoOfNodes, D);
            Metrics = MultidimensionalArray.Create(Len, NoOfNodes);

            var affNrm = this.GridData.Edges.NormalsForAffine;

            if(affineEdge && affNrm != null) {
                for(int i = 0; i < Len; i++) {
                    int iEdge = i + e0;
                    for(int k = 0; k < NoOfNodes; k++) {
                        for(int d = 0; d < D; d++) {
                            Normals[i, k, d] = affNrm[iEdge, d];
                        }
                    }
                }


                if(Metrics != null) {
                    var sqrtGrm = this.GridData.Edges.SqrtGramian;

                    for(int i = 0; i < Len; i++) {
                        int iEdge = i + e0;
                        for(int k = 0; k < NoOfNodes; k++) {
                            Metrics[i,k] = sqrtGrm[iEdge];
                        }
                    }
                }

            } else {
                int[] i0 = new int[3];
                int[] iE = new int[3];
                iE[2] = D - 1;
                iE[1] = NoOfNodes - 1;

                int[,] E2C = this.GridData.Edges.CellIndices;
                byte[,] FcIdx = this.GridData.Edges.FaceIndices;
                int[,] TrIdx = this.GridData.Edges.Edge2CellTrafoIndex;

                for(int e = 0; e < Len; e++) {
                    int iEdge = e + e0;
                    int jCell1, edge1;
                    jCell1 = E2C[iEdge, 0];
                    edge1 = TrIdx[iEdge, 0];
                    int iFace1 = FcIdx[iEdge, 0];

                    i0[0] = e;
                    iE[0] = e - 1;

                    this.GridData.Edges.GetNormalsForCell(NS.GetVolumeNodeSet(this.GridData, edge1), jCell1, iFace1, Normals, Metrics, e);
                }
            }
        }


        /// <summary>
        /// Evaluation of a single-valued property on edges.
        /// </summary>
        public MultidimensionalArray GetIntegrationMetric(NodeSet NS, int e0, int Len) {
            if(e0 == last_e0 && Len == last_Len && NS.Reference == last_NSref) {
                // lucky punch successful
                Debug.Assert(last_metrics.GetLength(0) == Len);
                Debug.Assert(last_metrics.GetLength(1) == NS.NoOfNodes);
            } else {
                EdgeEval(NS, e0, Len, out last_Normals, out last_metrics);
                last_e0 = e0;
                last_Len = Len;
                last_NSref = NS.Reference;
            }
            return last_metrics;
        }

 

        protected Grid.Classic.GridData GridData {
            get;
            private set;
        }
    }


    public class BasisEdgeValuesCacheLogic {

        internal BasisEdgeValuesCacheLogic(IGridData g, bool gradient) {
            this.GridData = g;
            this.m_gradient = gradient;
        }

        bool m_gradient;

        class CacheEntry {
            public ulong CacheRef_Value;
            public int Degree;
            public bool[] ValueAvail;
        }

        Dictionary<int,CacheEntry> m_LinkToCache = new Dictionary<int, CacheEntry>();

        public MultidimensionalArray GetValues(NodeSet NS, int e0, int Len, int Degree) {
            // argcheck
            // ========
            
            if(NS.GetNodeCoordinateSystem(this.GridData) != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentException();
            if(Len < 0)
                throw new ArgumentOutOfRangeException();
            if(e0 < 0 || e0 + Len > this.GridData.iGeomEdges.Count)
                throw new ArgumentOutOfRangeException();

            var E2Ctrafo = this.GridData.iGeomEdges.Edge2CellTrafos;
            int NoOfTrafos = E2Ctrafo.Count;
            int[,] trafoIdx = this.GridData.iGeomEdges.Edge2CellTrafoIndex;
            IList<int> kredIdx = this.GridData.iGeomEdges.Edge2CellTrafosRefElementIndices;
            int[,] cellIdx = this.GridData.iGeomEdges.CellIndices;

            // try to access cache
            // ===================

            // try to find CacheEntry ....
            CacheEntry CE;
            if(!this.m_LinkToCache.TryGetValue(NS.Reference, out CE)) {
                CE = new CacheEntry() {
                    ValueAvail = new bool[E2Ctrafo.Count]
                };
                this.m_LinkToCache.Add(NS.Reference, CE);
            } else {
                if(CE.Degree < Degree) {
                    // We might have some cached values, but with insufficient degree;
                    // throw away everything in this case.

                    CE.Degree = Degree;
                    Array.Clear(CE.ValueAvail, 0, CE.ValueAvail.Length);

                    if(CE.CacheRef_Value > 0)
                        Cache.ForgetItem(CE.CacheRef_Value);
                    CE.CacheRef_Value = 0;
                }
            }

            // try to obtain values from cache...
            MultidimensionalArray Value = null;
            if(CE.CacheRef_Value > 0) {
                Value = (MultidimensionalArray)Cache.GetItem(CE.CacheRef_Value);
                if(Value == null)
                    Array.Clear(CE.ValueAvail, 0, CE.ValueAvail.Length);
            }
            
            // allocate memmory (if necessary) and store it in cache...
            int Nmax = -1;
            Degree = Math.Max(CE.Degree, Degree);
            CE.Degree = Degree;
            PolynomialList[] Polynomials = this.GridData.ChefBasis.GetOrthonormalPolynomials(Degree);
            if(Value == null) {
                //Debug.Assert(N_kref.Length == N_kref.Length);

                for(int iKref = 0; iKref < Polynomials.Length; iKref++) {
                    int N_kref = Polynomials[iKref].Count;
                    Nmax = Math.Max(Nmax, N_kref);
                }

                if(!this.m_gradient) {
                    Value = MultidimensionalArray.Create(NoOfTrafos, NS.NoOfNodes, Nmax);
                } else {
                    Value = MultidimensionalArray.Create(NoOfTrafos, NS.NoOfNodes, Nmax, this.GridData.SpatialDimension);
                }
                CE.CacheRef_Value = Cache.CacheItem(Value, Value.Length * sizeof(double));
            }
            

            // recompute values if necessary
            // =============================
                        
            for(int e = 0; e < Len; e++) { // for all edges in the chung from 'e0' to 'e0 + Len - 1'...
                int iEdge = e + e0;

                // compute the value for the respective edge-to-cell trafo, if necessary:

                int iTrafoIN = trafoIdx[iEdge, 0];
                int iTrafoOT = trafoIdx[iEdge, 1];

                if(CE.ValueAvail[iTrafoIN] == false) {
                    // missing value for (edge-to-cell trafo for) IN-cell
                    int N_kref = Polynomials[kredIdx[iTrafoIN]].Count;
                    if(this.m_gradient)
                        NewMethod2(NS, Degree, Value, iTrafoIN, N_kref);
                    else 
                        NewMethod1(NS, Degree, Value, iTrafoIN, N_kref);
                    CE.ValueAvail[iTrafoIN] = true;
                }

                if(cellIdx[iEdge, 1] >= 0 && CE.ValueAvail[iTrafoOT] == false) {
                    // missing value for (edge-to-cell trafo for) OUT-cell
                    int N_kref = Polynomials[kredIdx[iTrafoOT]].Count;
                    if(this.m_gradient)
                        NewMethod2(NS, Degree, Value, iTrafoOT, N_kref);
                    else
                        NewMethod1(NS, Degree, Value, iTrafoOT, N_kref);
                    CE.ValueAvail[iTrafoOT] = true;
                }
            }
            
            // return
            // ======

            return Value;
        }

        private void NewMethod1(NodeSet NS, int Degree, MultidimensionalArray Value, int iTrafo, int N) {

            // values for node set 'NS' under edge-to-cell trafo #iKrefIN
            MultidimensionalArray value_iTrafo = this.GridData.ChefBasis.BasisValues.GetValues(NS.GetVolumeNodeSet(this.GridData, iTrafo), Degree);
            Debug.Assert(value_iTrafo.GetLength(1) >= N);
            if(value_iTrafo.GetLength(1) > N) {
                value_iTrafo = value_iTrafo.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NS.NoOfNodes - 1, N - 1 });
            }

            Value.ExtractSubArrayShallow(new int[] { iTrafo, 0, 0 }, new int[] { iTrafo - 1, NS.NoOfNodes - 1, N - 1 }).Set(value_iTrafo);
        }

        private void NewMethod2(NodeSet NS, int Degree, MultidimensionalArray Value, int iTrafo, int N) {

            // values for node set 'NS' under edge-to-cell trafo #iKrefIN
            MultidimensionalArray value_iTrafo = this.GridData.ChefBasis.BasisGradientValues.GetValues(NS.GetVolumeNodeSet(this.GridData, iTrafo), Degree);
            Debug.Assert(value_iTrafo.GetLength(1) >= N);
            int D = value_iTrafo.GetLength(2);
            Debug.Assert(D == this.GridData.SpatialDimension);
            if(value_iTrafo.GetLength(1) > N) {
                value_iTrafo = value_iTrafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { NS.NoOfNodes - 1, N - 1, D - 1 });
            }

            Value.ExtractSubArrayShallow(new int[] { iTrafo, 0, 0, 0 }, new int[] { iTrafo - 1, NS.NoOfNodes - 1, N - 1, D - 1 }).Set(value_iTrafo);
        }

        public IGridData GridData {
            get;
            private set;
        }

    }

}
