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
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Tracing;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using System.Threading;

namespace BoSSS.Foundation.Quadrature {

    public static class Quadrature_Bulksize {
        /// <summary>
        /// Number of floats (#items * #nodes * #integralg per cell) done at maximum
        /// </summary>
        public static int CHUNK_DATA_LIMIT = 12*1024*1024;
    }

    /// <summary>
    /// baseclass for vectorized quadrature
    /// </summary>
    public abstract class Quadrature<TQuadRule, TDomain> : IQuadrature
        where TQuadRule : QuadRule
        where TDomain : ExecutionMask {

        

        /// <summary>
        /// whether the integration should be performed in the reference or the physical coordinate system
        /// </summary>
        public CoordinateSystem CoordinateSystem {
            get;
            private set;
        }

        /// <summary>
        /// The total number of integrals per cell;
        /// </summary>
        protected int m_TotalNoOfIntegralsPerItem;

        /// <summary>
        /// see <see cref="IntegralCompDim"/>
        /// </summary>
        private int[] m_IntegralsComponent;

        /// <summary>
        /// Integrand dimension (the integrand can be a tensor, and this are its dimensions).
        /// </summary>
        protected int[] IntegralCompDim {
            get {
                return m_IntegralsComponent;
            }
        }

        private TQuadRule m_CurrentRule;

        /// <summary>
        /// the quadrature rule used for the current evaluation, only 
        /// available during the call to <see cref="Evaluate"/>
        /// </summary>
        public TQuadRule CurrentRule {
            get {
                if (m_CurrentRule == null)
                    throw new InvalidOperationException();
                return m_CurrentRule;
            }
        }

        /// <summary>
        /// For the currently used quadrature rule (see <see cref="CurrentRule"/>),
        /// the index of the reference element for which the rule is valid.
        /// </summary>
        public abstract int CurrentRuleRefElementIndex {
            get;
        }


        /// <summary>
        /// The grid on which this quadrature object operates on.
        /// </summary>
        protected IGridData gridData;

        /// <summary>
        /// The grid on which this quadrature object operates on.
        /// </summary>
        public IGridData GridDat {
            get {
                return gridData;
            }
        }

        
        /// <summary>
        /// results of evaluation<br/>
        /// 1st index: quadrature item (cell, edge)<br/>
        /// 2nd index: quadrature node <br/>
        /// 3rd to (<see cref="IntegralCompDim"/>.Length + 2)-th index: integral components
        /// </summary>
        protected MultidimensionalArray m_EvalResults;

        /// <summary>
        /// results of quadrature<br/>
        /// 1st index: quadrature item (cell, edge)<br/>
        /// 2nd to (<see cref="IntegralCompDim"/>.Length + 1)-th index: integral components
        /// </summary>
        protected MultidimensionalArray m_QuadResults;

        /// <summary>
        /// shallow copy of <see cref="m_QuadResults"/>, where the last (<see cref="IntegralCompDim"/>.Length) 
        /// indices/dimensions ale collapsed into one.
        /// </summary>
        protected MultidimensionalArray m_QuadResultsCollapsed;

        /// <summary>
        /// shallow copy of <see cref="m_EvalResults"/>, where the last (<see cref="IntegralCompDim"/>.Length) 
        /// indices/dimensions ale collapsed into one.
        /// </summary>
        protected MultidimensionalArray m_EvalResultsCollapsed;

        /// <summary>
        /// Standard constructor
        /// </summary>
        /// <param name="noOfIntegralsPerCell">
        /// The number of integrals to be computed
        /// </param>
        /// <param name="context">
        /// see <see cref="gridData"/>.
        /// </param>
        /// <param name="rule">
        /// quadrature domain and rules.
        /// </param>
        /// <param name="cs">integrate in physical or reference coordinates?</param>
        protected Quadrature(int[] noOfIntegralsPerCell, IGridData context, ICompositeQuadRule<TQuadRule> rule, CoordinateSystem cs) {
            m_TotalNoOfIntegralsPerItem = 1;
            foreach (var no in noOfIntegralsPerCell)
                m_TotalNoOfIntegralsPerItem *= no;
            m_IntegralsComponent = noOfIntegralsPerCell;

#if DEBUG
            if (rule.Count() > 0) {
                int J = rule.Max(crp => crp.Chunk.JE);
                System.Collections.BitArray ChunkTest = new System.Collections.BitArray(J);
                foreach (var chunk in rule) {
                    int IE = chunk.Chunk.JE;
                    for (int i = chunk.Chunk.i0; i < IE; i++) {
                        if (ChunkTest[i])
                            throw new ArgumentException("More than one quadrature rule defined for integration item " + i + ".");
                        ChunkTest[i] = true;

                    }
                }
            }
#endif


            gridData = context;
            m_compositeRule = rule;
            CoordinateSystem = cs;
        }

        /// <summary>
        /// Alternative constructor that allows for the usage of a scheme
        /// instead of a compiled quadrature rule.
        /// </summary>
        /// <param name="noOfIntegralsPerCell"></param>
        /// <param name="context"></param>
        /// <param name="scheme"></param>
        /// <param name="order"></param>
        /// <param name="cs">integrate in physical or reference coordinates?</param>
        protected Quadrature(
            int[] noOfIntegralsPerCell, Grid.Classic.GridData context, IQuadratureScheme<TQuadRule, TDomain> scheme, int order, CoordinateSystem cs)
            : this(noOfIntegralsPerCell, context, 
            scheme.Compile(context, order), cs) {
        }

        /// <summary>
        /// Quadrature domain (which cells/edges)
        /// and quadrature rules.
        /// </summary>
        protected ICompositeQuadRule<TQuadRule> m_compositeRule;

        /// <summary>
        /// 2nd phase of quadrature: allocation of memory for 
        /// the <see cref="Evaluate"/>-method;
        /// Called whenever the node set or the number of cells per evaluation is changed;
        /// </summary>
        /// <param name="NoOfItems">number of edges or cells to integrate</param>
        /// <param name="ruleNodes">quadrature rule nodes</param>
        protected virtual void AllocateBuffers(int NoOfItems, NodeSet ruleNodes) {
        }

        /// <summary>
        /// 2nd phase of quadrature: allocation of memory for 
        /// the <see cref="Evaluate"/>-method;
        /// Called whenever the node set or the number of cells per evaluation is changed;
        /// </summary>
        /// <param name="NoOfItems">number of edges or cells to integrate</param>
        /// <param name="ruleNodes">quadrature rule nodes</param>
        protected virtual void AllocateBuffersInternal(int NoOfItems, NodeSet ruleNodes) {
            Debug.Assert(ruleNodes.Dimension == 2);
            if(this.m_ExEvaluate == null) {
                if(m_EvalResults == null)
                    m_EvalResults = new MultidimensionalArray(2 + m_IntegralsComponent.Length);
                m_EvalResults.Allocate(((new int[] { NoOfItems, ruleNodes.GetLength(0) }).Concat(m_IntegralsComponent)).ToArray());
            }

            if(m_QuadResults == null)
                m_QuadResults = new MultidimensionalArray(1 + m_IntegralsComponent.Length);
            m_QuadResults.Allocate(((new int[] { NoOfItems }).Concat(m_IntegralsComponent)).ToArray());

            if(this.m_ExEvaluate == null) {
                m_EvalResultsCollapsed = m_EvalResults.ResizeShallow(
                    (m_EvalResults.Lengths.Take(2).Concat(new int[] { m_TotalNoOfIntegralsPerItem })).ToArray());
                m_QuadResultsCollapsed = m_QuadResults.ResizeShallow(
                    (m_QuadResults.Lengths.Take(1).Concat(new int[] { m_TotalNoOfIntegralsPerItem })).ToArray());
            }
            
            if(m_AllocateBuffers != null)
                m_AllocateBuffers(NoOfItems, ruleNodes);
        }


       
        /// <summary>
        /// 1st Phase of quadrature: preparation of nodes;
        /// This method will be called whenever a new 
        /// node set family was locked, i.e. whenever the set of quadrature nodes is changed.
        /// </summary>
        protected virtual void QuadNodesChanged(NodeSet newNodes) {
            if (m_quadNodesChanged != null)
                m_quadNodesChanged(newNodes);
        }

        /// <summary>
        /// 3rd phase of quadrature: vectorized evaluation of the integrand.
        /// Override this method to implement the integrand;
        /// </summary>
        /// <param name="i0">local index of first cell or edge</param>
        /// <param name="Length">number of cells or edges to process</param>
        /// <param name="QuadNodes">
        /// Quadrature nodes in reference coordinates.
        /// </param>
        /// <param name="EvalResult">
        /// On exit, the result of the integrand evaluation.
        /// Implementers can expect a cleared array, i.e. all entries are 0.0.
        ///  - 1st index: local cell or edge index minus <paramref name="i0"/>;
        ///  - 2nd index: Node Index;<br/>
        ///  - 3rd to (<see cref="IntegralCompDim"/>.Length + 2)-th index: integral component;
        /// </param>
        protected virtual void Evaluate(int i0, int Length, TQuadRule rule, MultidimensionalArray EvalResult) {
            m_Evaluate(i0, Length, rule, EvalResult);
        }

        /// <summary>
        /// 4th (and final) phase of quadrature: save the results somewhere;
        /// Override this method to write the results of the quadrature from temporary buffers
        /// to wherever they belong
        /// </summary>
        /// <param name="i0">local index of first cell or edge</param>
        /// <param name="Length">number of cells or edges to process</param>
        /// <param name="ResultsOfIntegration">
        /// 1st index: local cell or edge index minus <paramref name="i0"/>;
        /// 2nd to (<see cref="IntegralCompDim"/>.Length + 1)-th index:  integral component;
        /// </param>
        protected virtual void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
            m_SaveIntegrationResults(i0, Length, ResultsOfIntegration);
        }

        /// <summary>
        /// Additional timers, which may be derived in derived classes, that are added to the performance analysis.
        /// </summary>
        public Stopwatch[] CustomTimers {
            get {
                return m_CustomTimers;
            }
            set {
                m_CustomTimers = value;
            }
        }

        Stopwatch[] m_CustomTimers = new Stopwatch[0];

        /// <summary>
        /// Names of the custom timers (see <see cref="CustomTimers"/>), for reference reasons in the log files.
        /// Index correlates with <see cref="CustomTimers"/>;
        /// </summary>
        public string[] CustomTimers_Names {
            get {
                return m_CustomTimers_Names;
            }
            set {
                m_CustomTimers_Names = value;
            }
        }

        string[] m_CustomTimers_Names = new string[0];

        /// <summary>
        /// If one wants to establish a call tree among the custom timers, one can use this.
        /// </summary>
        public int[] CustomTimers_RootPointer {
            get {
                return m_CustomTimers_RootPointer;
            }
            set {
                m_CustomTimers_RootPointer = value;
            }
        }

        int[] m_CustomTimers_RootPointer = new int[0];

        /// <summary>
        /// performs an integration over cells or edges in the composite rule provided to the constructor;
        /// </summary>
        public virtual void Execute() {
            using (var tr = new FuncTrace()) {
                // Init
                // ====

                // timers
                Stopwatch stpwSaveIntRes = new Stopwatch();
                stpwSaveIntRes.Reset();
                Stopwatch stpwQuad = new Stopwatch();
                stpwQuad.Reset();
                Stopwatch stpwEval = new Stopwatch();
                stpwEval.Reset();
                Stopwatch stpwAlloc = new Stopwatch();
                stpwAlloc.Reset();
                Stopwatch stpwNdSet = new Stopwatch();
                stpwNdSet.Reset();

                for (int i = CustomTimers.Length - 1; i >= 0; i--)
                    CustomTimers[i].Reset();

                // check input ...
                IGridData grd = gridData;

                // do quadrature
                // =============
                MultidimensionalArray lastQuadRuleNodes = null;
                int oldBulksize = -1;
                int oldNoOfNodes = -1;
                int Bulkcnt = 0;
                foreach (var chunkRulePair in m_compositeRule) {
                    Chunk chunk = chunkRulePair.Chunk;
                    m_CurrentRule = chunkRulePair.Rule;

                    //// init node set
                    stpwNdSet.Start();
                    if (!object.ReferenceEquals(m_CurrentRule.Nodes, lastQuadRuleNodes)) {
                        QuadNodesChanged(m_CurrentRule.Nodes);
                    }
                    stpwNdSet.Stop();

                    // define bulk size
                    int NoOfNodes = m_CurrentRule.Nodes.GetLength(0);
                    int ItemSize = m_TotalNoOfIntegralsPerItem * NoOfNodes;
                    if(ItemSize <= 0)
                        continue;
                    int MaxChunkLength =  Quadrature_Bulksize.CHUNK_DATA_LIMIT / ItemSize;
                    if(MaxChunkLength < 1)
                        MaxChunkLength = 1;


                    int j = chunk.i0;
                    int ChunkLength = MaxChunkLength;
                    int ChunkEnd = chunk.i0 + chunk.Len;

                    
                    while (j < ChunkEnd) {
                        Bulkcnt++;

                        // limit bulksize 
                        if ((j + ChunkLength) > ChunkEnd) {
                            ChunkLength -= (j + ChunkLength - ChunkEnd);
                        }
                       

                        // DEBUG check
#if DEBUG
                        CheckQuadratureChunk(j, ChunkLength, CurrentRuleRefElementIndex);
#endif
                        
                        // reallocate buffers if bulksize was changed
                        stpwAlloc.Start();
                        if (ChunkLength != oldBulksize || m_CurrentRule.NoOfNodes != oldNoOfNodes) {
                            AllocateBuffersInternal(ChunkLength, m_CurrentRule.Nodes);
                            AllocateBuffers(ChunkLength, m_CurrentRule.Nodes);
                            oldBulksize = ChunkLength;
                            oldNoOfNodes = m_CurrentRule.NoOfNodes;
                        }
                        stpwAlloc.Stop();

                        
                        if(this.m_ExEvaluate == null) {

                            // evaluation of integrand
                            // =======================
                            stpwEval.Start();
                            m_EvalResults.Clear();
                            this.Evaluate(j, ChunkLength, this.CurrentRule, m_EvalResults);
                            stpwEval.Stop();

                            // quadrature
                            // ==========
                            stpwQuad.Start();
                            DoQuadrature(m_CurrentRule, j, ChunkLength);
                            stpwQuad.Stop();
                        } else {
                            Debug.Assert(this.m_Evaluate == null);

                            // evaluation of integrand
                            // =======================
                            stpwEval.Start();
                            this.m_ExEvaluate(j, ChunkLength, this.CurrentRule, m_QuadResults);
                            stpwEval.Stop();
                        }

                        // save results
                        // ============
                        stpwSaveIntRes.Start();
                        SaveIntegrationResults(j, ChunkLength, m_QuadResults);
                        stpwSaveIntRes.Stop();

                        // inc
                        j += ChunkLength;
                    }


                    lastQuadRuleNodes = m_CurrentRule.Nodes;
                }
                m_CurrentRule = null;

                //tr.Info("Quadrature performed in " + Bulkcnt + " chunk(s).");
                //Console.WriteLine("Quadrature performed in " + Bulkcnt + " chunk(s).");
               
                // finalize
                // ========
                
                {

                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    //
                    // note that StopWatch.Elapsed.Ticks != StopWatch.ElapsedTicks
                    //
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    var mcrEval = tr.LogDummyblock(stpwEval.Elapsed.Ticks, "integrand_evaluation");  
                    tr.LogDummyblock(stpwQuad.Elapsed.Ticks, "quadrature");
                    tr.LogDummyblock(stpwSaveIntRes.Elapsed.Ticks, "saving_results");
                    tr.LogDummyblock(stpwNdSet.Elapsed.Ticks, "node_set_management");
                    tr.LogDummyblock(stpwAlloc.Elapsed.Ticks, "buffer_allocation");

                    Debug.Assert(m_CustomTimers.Length == m_CustomTimers_Names.Length);
                    Debug.Assert(m_CustomTimers.Length == m_CustomTimers_RootPointer.Length);
                    MethodCallRecord[] mcrS = new MethodCallRecord[CustomTimers.Length];

                    for (int iTimer = 0; iTimer < mcrS.Length; iTimer++) {
                        int pt = m_CustomTimers_RootPointer[iTimer];
                        MethodCallRecord OwnerMcr = pt >= 0 ? mcrS[pt] : mcrEval;
                        mcrS[iTimer] = OwnerMcr.AddSubCall(CustomTimers_Names[iTimer], m_CustomTimers[iTimer].Elapsed.Ticks);
                    }

                    
                }
            }
        }

        /// <summary>
        /// performs the quadrature
        /// </summary>
        /// <param name="quadRule"></param>
        /// <param name="j0"></param>
        /// <param name="_Bulksize">number of quadrature items</param>
        protected virtual void DoQuadrature(TQuadRule quadRule, int j0, int _Bulksize) {
            var currentRuleWeights = quadRule.Weights;

            switch (CoordinateSystem) {

                case CoordinateSystem.Physical: {

                    int JE = j0 + _Bulksize;
                    int j = j0;
                    while(j < JE) {

                        // determine sub-chunk
                        bool Linear; int Bulksize;
                        NextPart(out Linear, out Bulksize, j, _Bulksize - (j - j0));
                        
                        unsafe {
                            int K = m_TotalNoOfIntegralsPerItem;
                            int N = currentRuleWeights.GetLength(0);
                            int NxK = N * K;

                            fixed (double* _pWeights = currentRuleWeights.Storage,
                                           pEvalResults = m_EvalResultsCollapsed.Storage,
                                           pQuadResults = m_QuadResultsCollapsed.Storage) {
                                double* pQr = pQuadResults;
                                double* pWeights = _pWeights + currentRuleWeights.Index(0);


                                if (Linear) {
                                    // codepath for integration of linear elements
                                    // +++++++++++++++++++++++++++++++++++++++++++

                                    var scalings = GetScalingsForLinearElements(j0, Bulksize);
                                    Debug.Assert(scalings.Dimension == 1);
                                    Debug.Assert(scalings.GetLength(0) == Bulksize);
                                    Debug.Assert(currentRuleWeights.IsContinious);
                                    Debug.Assert(currentRuleWeights.Dimension == 1);



                                    int iSc = scalings.Index(0);

                                    // loop over quadrature items (cells or edges)...
                                    for (int jj = 0; jj < Bulksize; jj++) {
                                        double* pEvalResuls_jj = pEvalResults + jj * NxK;
                                        double sc = scalings[jj];

                                        for (int k = 0; k < K; k++) { // loop over integral components
                                            double r = 0.0;

                                            for (int n = 0; n < N; n++) { // loop over quadrature nodes (summation)
                                                r += pEvalResuls_jj[n * K + k] * pWeights[n];
                                            }

                                            *pQr = r * sc;
                                            pQr++;
                                        }
                                    }



                                } else {
                                    // codepath for integration of nonlinear/curved elements
                                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                                    var scalings = GetScalingsForNonlinElements(j0, Bulksize);
                                    Debug.Assert(scalings.Dimension == 2);
                                    Debug.Assert(scalings.GetLength(0) == Bulksize);
                                    Debug.Assert(scalings.GetLength(1) == currentRuleWeights.GetLength(0));
                                    Debug.Assert(currentRuleWeights.IsContinious);
                                    Debug.Assert(currentRuleWeights.Dimension == 1);

                                    // loop over quadrature items (cells or edges)...
                                    for (int jj = 0; jj < Bulksize; jj++) {
                                        double* pEvalResuls_jj = pEvalResults + jj * NxK;


                                        for (int k = 0; k < K; k++) { // loop over integral components
                                            double r = 0.0;

                                            for (int n = 0; n < N; n++) { // loop over quadrature nodes (summation)
                                                r += pEvalResuls_jj[n * K + k]*pWeights[n]*scalings[jj, n];
                                            }

                                            *pQr = r;
                                            pQr++;
                                        }
                                    }
                                }
                            }
                        }

                        j += Bulksize;
                    }


                    break;
                }

                case CoordinateSystem.Reference: {
                    Debug.Assert(currentRuleWeights.IsContinious);
                    Debug.Assert(currentRuleWeights.Dimension == 1);

                    unsafe {
                        int K = m_TotalNoOfIntegralsPerItem;
                        int N = currentRuleWeights.GetLength(0);
                        int NxK = N * K;


                        fixed (double* _pWeights = currentRuleWeights.Storage,
                                       pEvalResults = m_EvalResultsCollapsed.Storage,
                                       pQuadResults = m_QuadResultsCollapsed.Storage) {
                            double* pQr = pQuadResults;
                            double* pWeights = _pWeights + currentRuleWeights.Index(0);

                            // loop over quadrature items (cells or edges)...
                            for (int jj = 0; jj < _Bulksize; jj++) {
                                double* pEvalResuls_jj = pEvalResults + jj * NxK;


                                for (int k = 0; k < K; k++) { // loop over integral components
                                    double r = 0.0;

                                    for (int n = 0; n < N; n++) { // loop over quadrature nodes (summation)
                                        r += pEvalResuls_jj[n * K + k] * pWeights[n];
                                    }

                                    *pQr = r;
                                    pQr++;
                                }
                            }
                        }
                    }
                    break;
                }

                default:
                    throw new NotImplementedException("not known: " + CoordinateSystem);
            }
        }

        /// <summary>
        /// various checks for the current quadrature check - only called in DEBUG mode;
        /// </summary>
        protected abstract void CheckQuadratureChunk(int j0, int Len, int iKRef);
        

        /// <summary>
        /// Implement this method by returning the correct scaling factors for
        /// elements with a linear mapping of the domain of integration (e.g., cell volumes).
        /// </summary>
        /// <returns>
        /// A list of scaling factors
        /// </returns>
        protected abstract MultidimensionalArray GetScalingsForLinearElements(int i0, int L);

        /// <summary>
        /// Returns Scaling metrics for
        /// elements with a non-linear mapping of the domain of integration (e.g., cell volumes).
        /// </summary>
        /// <param name="i0"></param>
        /// <param name="L"></param>
        /// <returns></returns>
        protected abstract MultidimensionalArray GetScalingsForNonlinElements(int i0, int L);

        /// <summary>
        /// identifies the chunk of elements (edges or cells), starting from <paramref name="i0"/>
        /// which are all either linear or all non-linear.
        /// </summary>
        /// <param name="Linear">
        /// on exit, true for linear elements, otherwise false
        /// </param>
        /// <param name="NoOfElm">
        /// on exit, the number of elements following <paramref name="i0"/> which are all either linear or all non-linear. 
        /// </param>
        /// <param name="i0"></param>
        /// <param name="Len">
        /// limit for the maximum chunk size, i.e. <paramref name="NoOfElm"/> will be bound by <paramref name="Len"/>.
        /// </param>
        protected abstract void NextPart(out bool Linear, out int NoOfElm, int i0, int Len);
                
        /// <summary>
        /// used by <see cref="CellQuadrature.GetQuadrature"/>,
        /// <see cref="EdgeQuadrature.GetQuadrature"/> and
        /// <see cref="CellBoundaryQuadrature{T}.GetQuadrature"/>
        /// </summary>
        public delegate void Del_Evaluate(int i0, int Length, TQuadRule rule, MultidimensionalArray EvalResult);

        


        /// <summary>
        /// used by <see cref="CellQuadrature.GetQuadrature"/>,
        /// <see cref="EdgeQuadrature.GetQuadrature"/> and
        /// <see cref="CellBoundaryQuadrature{T}.GetQuadrature"/>
        /// </summary>
        public delegate void Del_SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration);

        /// <summary>
        /// used by <see cref="CellQuadrature.GetQuadrature"/>,
        /// <see cref="EdgeQuadrature.GetQuadrature"/> and
        /// <see cref="CellBoundaryQuadrature{T}.GetQuadrature"/>
        /// </summary>
        public delegate void Del_AllocateBuffers(int NoOfItems, MultidimensionalArray ruleNodes);

        /// <summary>
        /// used by <see cref="CellQuadrature.GetQuadrature"/>,
        /// <see cref="EdgeQuadrature.GetQuadrature"/> and
        /// <see cref="CellBoundaryQuadrature{T}.GetQuadrature"/>
        /// </summary>
        public delegate void Del_QuadNodesChanged(NodeSet newNodes);

        internal Del_Evaluate m_Evaluate;
        internal Del_Evaluate m_ExEvaluate; // expert evaluate/ the user is responsible for multiplying with quad weigths
        internal Del_SaveIntegrationResults m_SaveIntegrationResults;
        internal Del_AllocateBuffers m_AllocateBuffers;
        internal Del_QuadNodesChanged m_quadNodesChanged;
    }
}
 
