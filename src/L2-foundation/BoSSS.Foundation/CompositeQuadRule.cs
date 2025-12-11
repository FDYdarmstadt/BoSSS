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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;

namespace BoSSS.Foundation.Quadrature {

   
	/// <summary>
	/// Instances of this class describe which quadrature rule should be used
	/// at which quadrature item (cell or edge).
	/// </summary>
	/// <remarks>
	/// The main motivation behind this class is the re-use of quadrature nodes
	/// and/or weights.
	/// </remarks>
	public sealed class CompositeQuadRule<TQuadRule> : ICompositeQuadRule<TQuadRule>
        where TQuadRule : QuadRule {

       
        /// <summary>
        /// Pairs of quadrature rules and the domain on which they should be
        /// applied
        /// </summary>
        public List<IChunkRulePair<TQuadRule>> chunkRulePairs = new List<IChunkRulePair<TQuadRule>>();

        public CompositeQuadRule(IIntegrationMetric integrationMetric, IGridData gridData) {
            IntegrationMetric = integrationMetric;
            GridData = gridData;
        }

        #region ICompositeQuadRule<TQuadRule> Members

        /// <summary>
        /// The number of quadrature items (cells or edges)
        /// </summary>
        public int NumberOfItems {
            get {
                return chunkRulePairs.Sum(p => p.Chunk.Len);
            }
        }

        /// <summary>
        /// <see cref="ICompositeQuadRule{TQuadRule}.IntegrationMetric"/>
        /// </summary>
        public IIntegrationMetric IntegrationMetric {
            get;
            set;
        }

        public IGridData GridData {
            get;
            private set;
        }


        #endregion

        #region IEnumerable<IChunkRulePair<TQuadRule>> Members

        /// <summary>
        /// Defines the actual mapping of the integration sub-domains to the
        /// respective quadrature rules.
        /// </summary>
        /// <returns>
        /// An enumerator over all pairs of chunks of the integration domain
        /// and the associated quadrature rules.
        /// </returns>
        public IEnumerator<IChunkRulePair<TQuadRule>> GetEnumerator() {
            return chunkRulePairs.GetEnumerator();
        }

        #endregion

        #region IEnumerable Members

        /// <summary>
        /// <see cref="GetEnumerator"/>
        /// </summary>
        /// <returns></returns>
        IEnumerator IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }

        #endregion

        /// <summary>
        /// Creates a new composite quadrature rules of minimal order
        /// <paramref name="order"/> from a single quadrature rule factory and
        /// an associated domain of integration. 
        /// </summary>
        /// <typeparam name="TDomain">
        /// The type of the domain to be integrated over.
        /// </typeparam>
        /// <param name="ruleFactory">
        /// The quadrature rule factory to be used to create the quadrature
        /// rules for the given domain.
        /// </param>
        /// <param name="order">
        /// The minimal order of the quadrature rules to be used. See
        /// <see cref="IQuadRuleFactory{T}.GetQuadRuleSet"/>.
        /// </param>
        /// <param name="domain">
        /// The domain to be integrated over.
        /// </param>
        /// <param name="scaling">
        /// </param>
        /// <returns>
        /// A composite rule containing the quadrature rules associated with
        /// all elements of the given domain.
        /// </returns>
        /// <remarks>
        /// The main specialty about the resulting composite rule is that the
        /// implementation guarantees that the sets of nodes and weights
        /// associated to the respective quadrature rules share the <b>same</b>
        /// objects of their contents are equal. This is allows for the
        /// optimization of the execution of the quadrature itself.
        /// </remarks>
        public static CompositeQuadRule<TQuadRule> Create<TDomain>(
            IQuadRuleFactory<TQuadRule> ruleFactory, int order, TDomain domain, IIntegrationMetric scaling)
            where TDomain : ExecutionMask //
        {
            using (var tr = new FuncTrace()) {
                CompositeQuadRule<TQuadRule> compositeRule = new CompositeQuadRule<TQuadRule>(scaling, domain.GridData);
                compositeRule.IntegrationMetric = scaling;
                // BEWARE: This check may cause nasty trouble in parallel runs
                // where a domain is only present on some domains and has thus been
                // removed by Björn
                //if (domain.NoOfItemsLocally == 0) {
                //    return compositeRule;
                //}
#if TEST
                tr.Info("Checkpoint CC.1");
                tr.Info("This is: " + ruleFactory.GetType());
#endif
                IEnumerable<IChunkRulePair<TQuadRule>> ruleSet;
                using(new BlockTrace("Rule_Compilation_" + ruleFactory.GetType().Name, tr)) {
                    ruleSet = ruleFactory.GetQuadRuleSet(domain, order);
                    foreach(var rule in ruleSet) {
                        rule.Rule.Origin = ruleFactory;
                        if(rule.Rule.OrderOfPrecision < order) {
                            throw new ArithmeticException($"Requested quadrature rule of degree {order}, but rule reports order {rule.Rule.OrderOfPrecision}.");
                        }                       
                    }
                }

                var nodes = new List<NodeSet>();
                var nodesMap = new Dictionary<MultidimensionalArray, int>();

                var weights = new List<MultidimensionalArray>();
                var weightsMap = new Dictionary<MultidimensionalArray, int>();
#if TEST
                tr.Info("Checkpoint CC.2");
#endif

                foreach (var chunkRulePair in ruleSet) {
                    Chunk chunk = chunkRulePair.Chunk;
                    TQuadRule rule = chunkRulePair.Rule;

                    Debug.Assert(rule.Nodes.IsLocked, "Error in quadrature rule creation: factory delivered some rule non-locked node set.");

                    int iNode;
                    if (!nodesMap.TryGetValue(rule.Nodes, out iNode)) {
                        // nodes must be added
                        nodes.Add(rule.Nodes);
                        iNode = nodes.Count - 1;
                        nodesMap.Add(rule.Nodes, iNode);
                    }
#if TEST
                    tr.Info("Checkpoint CC.3");
#endif
                    int iWeight;
                    if (!weightsMap.TryGetValue(rule.Weights, out iWeight)) {
                        // weights must be added
                        weights.Add(rule.Weights);
                        iWeight = weights.Count - 1;
                        weightsMap.Add(rule.Weights, iWeight);
                    }
#if TEST
                    tr.Info("Checkpoint CC.4");
#endif
                    // Make sure arrays are not only equal but identical (i.e.,
                    // reference equals)
                    rule.Nodes = nodes[iNode];
                    rule.Weights = weights[iWeight];
#if TEST
                    tr.Info("Checkpoint CC.5");
#endif
                    compositeRule.chunkRulePairs.Add(
                        new ChunkRulePair<TQuadRule>(chunk, rule));
                }

                return compositeRule;
            }
        }

        /// <summary>
        /// Merges two quadrature rules into one. If for one quadrature item
        /// <paramref name="A"/> and <paramref name="B"/> both define quadrature
        /// rules, it will always be the one from <paramref name="B"/> that ends
        /// up in the result. 
        /// </summary>
        /// <param name="A">
        /// The inferior quadrature rule (i.e., is dominated by
        /// <paramref name="B"/> in case of conflicts)
        /// </param>
        /// <param name="B">
        /// The dominant quadrature rule
        /// </param>
        /// <returns>
        /// A merged quadrature rule
        /// </returns>
        public static CompositeQuadRule<TQuadRule> Merge(
            CompositeQuadRule<TQuadRule> A, CompositeQuadRule<TQuadRule> B) {
            if(!object.ReferenceEquals(A.GridData, B.GridData))
                throw new ArgumentException("grid mismatch");

            CompositeQuadRule<TQuadRule> result = new CompositeQuadRule<TQuadRule>(B.IntegrationMetric, B.GridData);
            result.chunkRulePairs = ZipHelper2(A.chunkRulePairs, B.chunkRulePairs);

            var nodes = new List<NodeSet>();
            var nodesMap = new Dictionary<MultidimensionalArray, int>();

            var weights = new List<MultidimensionalArray>();
            var weightsMap = new Dictionary<MultidimensionalArray, int>();

            // Make sure arrays are not only equal but also identical (i.e., in
            // terms of reference equals)
            foreach (var pair in result.chunkRulePairs) {
                int iNode;
                if (!nodesMap.TryGetValue(pair.Rule.Nodes, out iNode)) {
                    nodes.Add(pair.Rule.Nodes);
                    iNode = nodes.Count - 1;
                    nodesMap.Add(pair.Rule.Nodes, iNode);
                }

                int iWeight;
                if (!weightsMap.TryGetValue(pair.Rule.Weights, out iWeight)) {
                    weights.Add(pair.Rule.Weights);
                    iWeight = weights.Count - 1;
                    weightsMap.Add(pair.Rule.Weights, iWeight);
                }

                pair.Rule.Nodes = nodes[iNode];
                pair.Rule.Weights = weights[iWeight];
            }

            return result;
        }

        static List<IChunkRulePair<TQuadRule>> ZipHelper2(List<IChunkRulePair<TQuadRule>> A, List<IChunkRulePair<TQuadRule>> B) {

            // spezialfaelle
            // -------------
            if (A.Count == 0)
                return B;
            if (B.Count == 0)
                return A;

            List<IChunkRulePair<TQuadRule>> output = new List<IChunkRulePair<TQuadRule>>();

            var eA = A.GetEnumerator();
            var eB = B.GetEnumerator();

            bool runA = eA.MoveNext();
            bool runB = eB.MoveNext();
            ChunkRulePair<TQuadRule> cA = null;
            ChunkRulePair<TQuadRule> cB = null;
            if (runA)
                cA = new ChunkRulePair<TQuadRule>(eA.Current);
            if (runB) {
                cB = new ChunkRulePair<TQuadRule>(eB.Current);
            }

            // solong in boade lischtn wos drinnen isch,
            // miassn ma se zommentian...
            // -----------------------------------------
            while (runA && runB) {
                if (cA._Chunk.JE <= cB._Chunk.i0) {
                    // Fall 1
                    output.Add(cA);
                    runA = eA.MoveNext();
                    if (runA)
                        cA = new ChunkRulePair<TQuadRule>(eA.Current);
                    continue;
                }

                if (cB._Chunk.JE <= cA._Chunk.i0) {
                    output.Add(cB);
                    runB = eB.MoveNext();
                    if (runB) {
                        cB = new ChunkRulePair<TQuadRule>(eB.Current);
                    }
                    continue;
                }


                if (cB.Chunk.i0 <= cA.Chunk.i0) {

                    // iatz kimmt zersch't amol cB!
                    output.Add(cB);

                    if (cB._Chunk.JE < cA._Chunk.JE) {
                        // ein stückl von cA bleibt ibrig.
                        // -> aufkalten firn negscht'n loop
                        cA._Chunk.Len = cA._Chunk.JE - cB._Chunk.JE;
                        cA._Chunk.i0 = cB._Chunk.JE;
                        Debug.Assert(cA.Chunk.Len > 0);
                    } else {
                        while (runA && (cB._Chunk.JE >= cA._Chunk.JE)) {
                            // cA ist vollständig von cB überdeckt
                            // aktuelles cA wegschmeissen, weiter zum negscht'n
                            runA = eA.MoveNext();
                            if (runA)
                                cA = new ChunkRulePair<TQuadRule>(eA.Current);
                        }

                        if (runA && cA._Chunk.i0 <= cB._Chunk.JE) {
                            // vo dem cA wo ma jetz ham, bleibt a stückl übrig.
                            // -> aufkalten fir negscht'n loop
                            cA._Chunk.Len = cA._Chunk.JE - cB._Chunk.JE;
                            cA._Chunk.i0 = cB._Chunk.JE;
                            Debug.Assert(cA.Chunk.Len > 0);
                        }
                    }

                    // weiter zum negscht'n cB
                    runB = eB.MoveNext();
                    if (runB) {
                        cB = new ChunkRulePair<TQuadRule>(eB.Current);
                    }
                    continue;
                }

                if (cB.Chunk.i0 > cA.Chunk.i0) {
                    // vorne isch no a stückl cA!
                    Debug.Assert(cA._Chunk.JE > cB._Chunk.i0);
                    Debug.Assert(cB._Chunk.i0 < cA._Chunk.JE);
                    var cX = new ChunkRulePair<TQuadRule>(cA);
                    cX._Chunk.Len = cB._Chunk.i0 - cA._Chunk.i0;
                    Debug.Assert(cX.Chunk.Len > 0);
                    output.Add(cX);

                    cA._Chunk.Len = cA._Chunk.JE - cX._Chunk.JE;
                    Debug.Assert(cA.Chunk.Len > 0); // müsste eigentlich schon durch "Fall 1" (siehe oben) abgefangen werden.
                    cA._Chunk.i0 = cX._Chunk.JE;
                    continue;
                }


                Debug.Assert(false); // wenn mir do her kemmen, nor isch eppes gonz tschelawenget.
            }

            // izant isch nur mehr in A oder in B wos drinnen,
            // oder in koaner vu boaden.
            // -----------------------------------------------

            while (runB) {
                output.Add(cB);
                runB = eB.MoveNext();
                if (runB) {
                    cB = new ChunkRulePair<TQuadRule>(eB.Current);
                }
            }

            while (runA) {
                output.Add(cA);
                runA = eA.MoveNext();
                if (runA)
                    cA = new ChunkRulePair<TQuadRule>(eA.Current);
            }

            // return
            // ------

            return output;

        }

    }
}
