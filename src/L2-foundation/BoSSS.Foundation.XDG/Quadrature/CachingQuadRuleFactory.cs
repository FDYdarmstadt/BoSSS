
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.XDG.Quadrature {

    /// <summary>
    /// Performs two purposes, mainly for quadrature rules which are expensive to generate:
    /// - rules, once created, should be re-used as much as possible
    /// - multi-threaded creation of these rules
    /// </summary>
    class CachingQuadRuleFactory : IQuadRuleFactory<QuadRule> {


        public CachingQuadRuleFactory(Func<int, IQuadRuleFactory<QuadRule>> __OriginalRuleFactory, ExecutionMask __Domain, int __Order) {
            //OriginalRuleFactory = __OriginalRuleFactory;
            Order = __Order;
            m_MaxDomain = __Domain;

            // serial version, to be multi-threaded
            var Rule = __OriginalRuleFactory(0);
            RefElement = Rule.RefElement;
            m_QuadRule = Rule.GetQuadRuleSet(__Domain, __Order).ToArray();


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

        readonly ExecutionMask m_MaxDomain;



        IChunkRulePair<QuadRule>[] m_QuadRule;

        //readonly Func<int, IQuadRuleFactory<QuadRule>> OriginalRuleFactor

        /// <summary>
        /// mapping from cell/edge index into <see cref="m_QuadRule"/>
        /// </summary>
        int[] m_ItemsToChunk;

        readonly int Order;

        public RefElement RefElement {
            get;
            private set;
        }

        public int[] GetCachedRuleOrders() {
            return [Order];
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask domain, int order) {
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
            if(!domain.IsSubMaskOf(m_MaxDomain))
                throw new ArgumentException("Cannot return quadrature rule for mask larger than initially specified");
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

            var Ret = new List<IChunkRulePair<QuadRule>>();

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
                        Ret.Add(new ChunkRulePair<QuadRule>(newChunk, pair.Rule));
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
            for(int i = 0; i < grd.iLogicalEdges.Count; i++) {
                Debug.Assert(checkBitMask[i] == maskBitMask[i], "not all cells/edges are covered by the rule.");
            }
#endif

            return Ret;
        }
    }


}