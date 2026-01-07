
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.XDG.Quadrature {

    /// <summary>
    /// Performs two purposes, mainly for quadrature rules which are expensive to generate:
    /// - rules, once created, should be re-used as much as possible
    /// - multi-threaded creation of these rules
    /// </summary>
    class CachingQuadRuleFactory<TQuadRule> : IQuadRuleFactory<TQuadRule> where TQuadRule : QuadRule {


        public CachingQuadRuleFactory(Func<int, IQuadRuleFactory<TQuadRule>> OriginalRuleFactoryFactory) {
            int NumThreads = ilPSP.Environment.NumThreads;
            OriginalRuleFactory = NumThreads.ForLoop(iThread => OriginalRuleFactoryFactory(iThread));
        }

        /// <summary>
        /// upper limit for number of integration items this rule can actually support
        /// </summary>
        int IMax {
            get {
                var grd = m_MaxDomain.GridData;
                int iMax;
                if(m_MaxDomain is EdgeMask) {
                    iMax = grd.iGeomEdges.Count;
                } else if(m_MaxDomain is CellMask) {
                    iMax = grd.iGeomCells.NoOfLocalUpdatedCells;
                } else {
                    throw new NotImplementedException("ähmmmm");
                }
                return iMax;
            }
        }

        ExecutionMask m_MaxDomain;

        IChunkRulePair<TQuadRule>[] m_QuadRule;

        readonly IQuadRuleFactory<TQuadRule>[] OriginalRuleFactory;

        /// <summary>
        /// mapping from cell/edge index into <see cref="m_QuadRule"/>
        /// </summary>
        int[] m_ItemsToChunk;

        int Order = -1;

        public RefElement RefElement {
            get {
                return OriginalRuleFactory[0].RefElement;
            }
        }

        public int[] GetCachedRuleOrders() {
            if(Order < 0)
                return [];
            else 
                return [Order];
        }

        void Initialize(ExecutionMask __Domain, int __Order) {
            if(Order < 0) {
                // 1st-time-call
                Order = __Order;
                m_MaxDomain = __Domain;
            } else {
                // later call, cells/edges are added
                m_MaxDomain = m_MaxDomain.Union(__Domain);

                if(__Order != Order)
                    throw new ArgumentException("cannot mix orders");
            }

            //// serial version, to be multi-threaded
            //var Rule = OriginalRuleFactory[0];
            //var newQuadRule = Rule.GetQuadRuleSet(__Domain, __Order).ToArray();
            IChunkRulePair<TQuadRule>[] newQuadRuleCombined;
            {
                int NumThreads = ilPSP.Environment.NumThreads;
                IEnumerable<IChunkRulePair<TQuadRule>>[] newQuadRules = new IEnumerable<IChunkRulePair<TQuadRule>>[NumThreads];
                var DomainParts = __Domain.SplitUp(NumThreads);
                bool restore = ilPSP.Environment.StdOut.surpressStream0;
                ilPSP.Environment.ParallelFor(0, NumThreads, delegate (int iPart) {
                    var Domain_i = DomainParts[iPart];
                    if(Domain_i.NoOfItemsLocally > 0)
                        newQuadRules[iPart] = OriginalRuleFactory[iPart].GetQuadRuleSet(Domain_i, __Order);
                    else 
                        newQuadRules[iPart] = [ ];
                }, enablePar: true);
                ilPSP.Environment.StdOut.surpressStream0 = restore;
                newQuadRuleCombined = newQuadRules[0].ToArray();
                for(int iPart = 1; iPart < NumThreads; iPart++) {
                    newQuadRuleCombined = newQuadRuleCombined.Cat(newQuadRules[iPart]);
                }
            }

            if(m_QuadRule == null) {
                m_QuadRule = newQuadRuleCombined;
            } else {
                m_QuadRule = m_QuadRule.MergeDisjointRules(newQuadRuleCombined);
            }

            //
            m_ItemsToChunk = new int[IMax];
            m_ItemsToChunk.SetAll(int.MinValue);
            int NoOfChunks = m_QuadRule.Length;
            for(int iChunk = 0; iChunk < NoOfChunks; iChunk++) {
                var crp = m_QuadRule[iChunk];
                for(int i = 0; i < crp.Chunk.Len; i++) {
                    m_ItemsToChunk[crp.Chunk.i0 + i] = iChunk;
                }
            }

        }

        /*
        static IChunkRulePair<TQuadRule>[] MergeRules(IChunkRulePair<TQuadRule>[] A, IChunkRulePair<TQuadRule>[] B) {
            int iA = 0;
            int iB = 0;
            int AL = A.Length;
            int BL = B.Length;
            int ML = AL + BL;
            IChunkRulePair<TQuadRule>[] M = new IChunkRulePair<TQuadRule>[ML];
            
            for(int iM = 0; iM < ML; iM++) {
                if(iA < AL && (iB >= B.Length || A[iA].Chunk.i0 < B[iB].Chunk.i0)) {
                    M[iM] = A[iA];
                    iA++;
                } else if(iB < BL && (iA >= A.Length || B[iB].Chunk.i0 < A[iA].Chunk.i0)) {
                    M[iM] = B[iB];
                    iB++;
                } else {
                    throw new ArgumentException("cannot merge overlapping quadrature rules");
                }

                if(iM > 0) {
                    if(M[iM].Chunk.i0 < M[iM - 1].Chunk.JE)
                        throw new ArgumentException("rules to merge seem to have overlap - not supported");
                }
            }
            return M;
        }
        */


        public IEnumerable<IChunkRulePair<TQuadRule>> GetQuadRuleSet(ExecutionMask domain, int order) {
            var missingPart = m_MaxDomain != null ? domain.Except(m_MaxDomain) : domain;
            if(this.Order < 0 || missingPart.NoOfItemsLocally > 0)
                Initialize(missingPart, order);

            // init & checks
            // =============
            if(order > this.Order)
                throw new ArgumentException($"requested quadrature order to high (requested {order}, maximum allowed {this.Order})");
            if(domain.GridData != m_MaxDomain.GridData)
                throw new ArgumentException("grid mismatch");
            if(domain.GetType() != m_MaxDomain.GetType())
                throw new ArgumentException($"Expecting an {m_MaxDomain.GetType()}.");
            if(domain.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");
            //if(!domain.IsSubMaskOf(m_MaxDomain)) {
            //    throw new ArgumentException("Cannot return quadrature rule for mask larger than initially specified");
            //}
            var grd = domain.GridData;
#if DEBUG

            var checkBitMask = new BitArray(IMax);
#endif

            // Extract rule for `domain` form the `m_QuadRule`
            // ===============================================
            if(domain.NoOfItemsLocally == m_MaxDomain.NoOfItemsLocally)
                // * `domain` is a subdomain of `m_MaxDomain`, and 
                // * number of elements are equal
                // => domain == m_MaxDomain
                // => no extraction required, return quadrature rule directly
                return m_QuadRule;

            var Ret = new List<IChunkRulePair<TQuadRule>>();

            foreach(var chunk in domain) {
                for(int i = chunk.i0; i < chunk.JE; i++) {
                    int iQuadRule = m_ItemsToChunk[i];
                    var pair = m_QuadRule[iQuadRule];
                    if(pair.Chunk.i0 == i && pair.Chunk.Len == (chunk.JE - i)) {
                        // the pair is an exact fit for the rest of the chunk
                        Ret.Add(pair);
                        i += pair.Chunk.Len - 1;
                    } else {
                        var newChunk = new Chunk() {
                            i0 = i,
                            Len = Math.Min(pair.Chunk.JE, chunk.JE) - i
                        };
                        Ret.Add(new ChunkRulePair<TQuadRule>(newChunk, pair.Rule));
                        i += newChunk.Len - 1;
                    }
#if DEBUG
                    for(int i2 = Ret.Last().Chunk.i0; i2 < Ret.Last().Chunk.JE; i2++) {
                        Debug.Assert(checkBitMask[i2] == false, "rule for respective cell/edge is already set.");
                        checkBitMask[i2] = true;
                    }
#endif
                }

            }

            // return
            // ======

#if DEBUG
            var maskBitMask = domain.GetBitMask();
            for(int i = 0; i < IMax; i++) {
                Debug.Assert(checkBitMask[i] == maskBitMask[i], "not all cells/edges are covered by the rule.");
            }
#endif

            return Ret;
        }

        public override string ToString() {
            return "Caching-" + (this.OriginalRuleFactory?.First()?.GetType()?.Name ?? "null");
        }
    }


}