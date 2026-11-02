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
using BoSSS.Foundation.Quadrature.FluxQuadCommon;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using ilPSP;
using System.Threading;

namespace BoSSS.Foundation.Quadrature.Linear {

    class LECEdgeQuadrature2<M, V>
        where M : IMutableMatrix
        where V : IList<double> {

        public LECEdgeQuadrature2(DifferentialOperator op) {
            Operator = op;

           
        }

        DifferentialOperator Operator;


        /// <summary>
        /// number of codomain variables
        /// </summary>
        int GAMMA {
            get {
                return Operator.CodomainVar.Count;
            }
        }

        /// <summary>
        /// number of domain variables
        /// </summary>
        int DELTA {
            get {
                return Operator.DomainVar.Count;
            }
        }

        protected void AllocateBuffers(int NoOfItems, MultidimensionalArray rule, int iThread, int NumThreads) {
            m_ThreadLocals[iThread].AllocateBuffers(NoOfItems, rule, iThread, NumThreads);
        }

        protected void OnCloneForThreadParallelization(IQuadrature q, int iThread, int NumThreads) {
            if (m_ThreadLocals == null) {
                m_ThreadLocals = new ThreadLocals[NumThreads];
            }

            m_ThreadLocals[iThread] = new ThreadLocals(iThread, this, q);
        }


        ThreadLocals[] m_ThreadLocals;

        class ThreadLocals {
            public ThreadLocals(int iThread, LECEdgeQuadrature2<M, V> owner, IQuadrature q) {
                m_owner = owner;
                m_iThread = iThread;

                var op = m_owner.Operator;
                m_Edgeform_UxV = EquationComponentArgMapping<IEdgeForm_UxV>.GetArgMapping(op, true,
                    eq => ((eq.BoundaryEdgeTerms & TermActivationFlags.UxV) != 0) || ((eq.InnerEdgeTerms & TermActivationFlags.UxV) != 0),
                    eq => (eq is IEdgeForm) ? new LinearEdgeFormVectorizer((IEdgeForm)eq) : null);
                m_Edgeform_GradUxV = EquationComponentArgMapping<IEdgeform_GradUxV>.GetArgMapping(op, true,
                    eq => ((eq.BoundaryEdgeTerms & TermActivationFlags.GradUxV) != 0) || ((eq.InnerEdgeTerms & TermActivationFlags.GradUxV) != 0),
                    eq => (eq is IEdgeForm) ? new LinearEdgeFormVectorizer((IEdgeForm)eq) : null);
                m_Edgeform_UxGradV = EquationComponentArgMapping<IEdgeform_UxGradV>.GetArgMapping(op, true,
                    eq => ((eq.BoundaryEdgeTerms & TermActivationFlags.UxGradV) != 0) || ((eq.InnerEdgeTerms & TermActivationFlags.UxGradV) != 0),
                    eq => (eq is IEdgeForm) ? new LinearEdgeFormVectorizer((IEdgeForm)eq) : null);
                m_Edgeform_GradUxGradV = EquationComponentArgMapping<IEdgeform_GradUxGradV>.GetArgMapping(op, true,
                    eq => ((eq.BoundaryEdgeTerms & TermActivationFlags.GradUxGradV) != 0) || ((eq.InnerEdgeTerms & TermActivationFlags.GradUxGradV) != 0),
                    eq => (eq is IEdgeForm) ? new LinearEdgeFormVectorizer((IEdgeForm)eq) : null);
                m_EdgeSourceV = EquationComponentArgMapping<IEdgeSource_V>.GetArgMapping(op, true,
                    eq => ((eq.BoundaryEdgeTerms & TermActivationFlags.V) != 0) || ((eq.InnerEdgeTerms & TermActivationFlags.V) != 0),
                    eq => (eq is IEdgeForm) ? new LinearEdgeFormVectorizer((IEdgeForm)eq) : null);
                m_EdgeSourceGradV = EquationComponentArgMapping<IEdgeSource_GradV>.GetArgMapping(op, true,
                    eq => ((eq.BoundaryEdgeTerms & TermActivationFlags.GradV) != 0) || ((eq.InnerEdgeTerms & TermActivationFlags.GradV) != 0),
                    eq => (eq is IEdgeForm) ? new LinearEdgeFormVectorizer((IEdgeForm)eq) : null);

                if (q == null)
                    return;

                q.CustomTimers = new Stopwatch[] { new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch() };
                q.CustomTimers_Names = new string[] { "Flux-Eval", "Basis-Eval", "LoopsLECedge", "ParametersAndNormals", "Flux-Trafo", "Flux-Lock-Wait" };
                q.CustomTimers_RootPointer = new int[q.CustomTimers_Names.Length];
                ArrayTools.SetAll(q.CustomTimers_RootPointer, -1);
                this.Flux_Eval = q.CustomTimers[0];
                this.Basis_Eval = q.CustomTimers[1];
                this.Loops = q.CustomTimers[2];
                this.ParametersAndNormals = q.CustomTimers[3];
                this.FluxTrafo = q.CustomTimers[4];
                this.FluxLockWait = q.CustomTimers[5];

                Debug.Assert(Array.IndexOf(q.CustomTimers_Names, "Flux-Eval") == 0);
                this.m_Edgeform_UxV_Watches = this.m_Edgeform_UxV.InitStopWatches(0, q);
                this.m_Edgeform_UxGradV_Watches = this.m_Edgeform_UxGradV.InitStopWatches(0, q);
                this.m_Edgeform_GradUxV_Watches = this.m_Edgeform_GradUxV.InitStopWatches(0, q);
                this.m_Edgeform_GradUxGradV_Watches = this.m_Edgeform_GradUxGradV.InitStopWatches(0, q);

                this.m_EdgeSourceV_Watches = this.m_EdgeSourceV.InitStopWatches(0, q);
                this.m_EdgeSourceGradV_Watches = this.m_EdgeSourceGradV.InitStopWatches(0, q);
            }

            int m_iThread;

            LECEdgeQuadrature2<M, V> m_owner;


            EquationComponentArgMapping<IEdgeForm_UxV>[] m_Edgeform_UxV;
            EquationComponentArgMapping<IEdgeform_GradUxV>[] m_Edgeform_GradUxV;
            EquationComponentArgMapping<IEdgeform_UxGradV>[] m_Edgeform_UxGradV;
            EquationComponentArgMapping<IEdgeform_GradUxGradV>[] m_Edgeform_GradUxGradV;
            EquationComponentArgMapping<IEdgeSource_V>[] m_EdgeSourceV;
            EquationComponentArgMapping<IEdgeSource_GradV>[] m_EdgeSourceGradV;

            /// <summary>
            /// true, if this integrator is responsible for any component
            /// </summary>
            public bool IsNonEmpty {
                get {
                    return m_Edgeform_UxV.IsNonEmpty() ||
                        m_Edgeform_GradUxV.IsNonEmpty() ||
                        m_Edgeform_UxGradV.IsNonEmpty() ||
                        m_Edgeform_GradUxGradV.IsNonEmpty() ||
                        m_EdgeSourceV.IsNonEmpty() ||
                        m_EdgeSourceGradV.IsNonEmpty();
                }
            }

            /// <summary>
            /// true, if integration of <see cref="m_Matrix"/> is required.
            /// </summary>
            public bool LinearRequired {
                get {
                    if (m_owner.m_Matrix == null)
                        return false;
                    for (int gamma = m_owner.GAMMA - 1; gamma >= 0; gamma--) {
                        if (m_Edgeform_UxV[gamma].m_AllComponentsOfMyType.Length > 0)
                            return true;
                        if (m_Edgeform_UxGradV[gamma].m_AllComponentsOfMyType.Length > 0)
                            return true;
                        if (m_Edgeform_GradUxV[gamma].m_AllComponentsOfMyType.Length > 0)
                            return true;
                        if (m_Edgeform_GradUxGradV[gamma].m_AllComponentsOfMyType.Length > 0)
                            return true;
                    }
                    return false;
                }
            }

            /// <summary>
            /// true, if integration of <see cref="m_AffineVector"/> is required.
            /// </summary>
            public bool AffineRequired {
                get {
                    if (m_owner.m_AffineVector == null)
                        return false;
                    for (int gamma = m_owner.GAMMA - 1; gamma >= 0; gamma--) {
                        if (m_EdgeSourceV[gamma].m_AllComponentsOfMyType.Length > 0)
                            return true;
                        if (m_EdgeSourceGradV[gamma].m_AllComponentsOfMyType.Length > 0)
                            return true;
                    }
                    return false;
                }
            }

            MultidimensionalArray ValuesXquadWgt = new MultidimensionalArray(3); // test function values   times quadrature weight
            MultidimensionalArray GradsXquadWgt = new MultidimensionalArray(4);  // test function gradient times quadrature weight



            public void AllocateBuffers(int NoOfItems, MultidimensionalArray rule, int iThread, int NumThreads) {
                int NoOfNodes = rule.GetLength(0);
                int D = m_owner.m_GridDat.SpatialDimension;



                // parameters
                // ==========
                {
                    if (m_ParameterFieldsValues_IN == null) {
                        m_ParameterFieldsValues_IN = m_owner.m_ParameterFields.Select(f => new MultidimensionalArray(2)).ToArray();
                    }
                    if (m_ParameterFieldsValues_OUT == null) {
                        m_ParameterFieldsValues_OUT = m_owner.m_ParameterFields.Select(f => new MultidimensionalArray(2)).ToArray();
                    }
                    foreach (var pfv in m_ParameterFieldsValues_IN) {
                        if (pfv.GetLength(0) != NoOfItems || pfv.GetLength(1) != NoOfNodes)
                            pfv.Allocate(NoOfItems, NoOfNodes);
                    }
                    foreach (var pfv in m_ParameterFieldsValues_OUT) {
                        if (pfv.GetLength(0) != NoOfItems || pfv.GetLength(1) != NoOfNodes)
                            pfv.Allocate(NoOfItems, NoOfNodes);
                    }
                }

                // component and summation buffers
                // ===============================
                {
                    AllocateComponentBuffer(ref m_UxVComponentBuffer, this.m_Edgeform_UxV, NoOfItems, NoOfNodes, D, 5, false);
                    AllocateComponentBuffer(ref m_UxGradVComponentBuffer, this.m_Edgeform_UxGradV, NoOfItems, NoOfNodes, D, 6, false);
                    AllocateComponentBuffer(ref m_GradUxVComponentBuffer, this.m_Edgeform_GradUxV, NoOfItems, NoOfNodes, D, 6, false);
                    AllocateComponentBuffer(ref m_GradUxGradVComponentBuffer, this.m_Edgeform_GradUxGradV, NoOfItems, NoOfNodes, D, 7, false);
                    AllocateComponentBuffer(ref m_VComponentBuffer, this.m_EdgeSourceV, NoOfItems, NoOfNodes, D, 3, true);
                    AllocateComponentBuffer(ref m_GradVComponentBuffer, this.m_EdgeSourceGradV, NoOfItems, NoOfNodes, D, 4, true);

                    AllocateSumBuffer(ref m_UxVSumBuffer, this.m_Edgeform_UxV, NoOfItems, NoOfNodes, D, 2, false);
                    AllocateSumBuffer(ref m_GradUxVSumBuffer, this.m_Edgeform_GradUxV, NoOfItems, NoOfNodes, D, 3, false);
                    AllocateSumBuffer(ref m_UxGradVSumBuffer, this.m_Edgeform_UxGradV, NoOfItems, NoOfNodes, D, 3, false);
                    AllocateSumBuffer(ref m_GradUxGradVSumBuffer, this.m_Edgeform_GradUxGradV, NoOfItems, NoOfNodes, D, 4, false);
                    AllocateSumBuffer(ref m_VSumBuffer, this.m_EdgeSourceV, NoOfItems, NoOfNodes, D, 2, true);
                    AllocateSumBuffer(ref m_GradVSumBuffer, this.m_EdgeSourceGradV, NoOfItems, NoOfNodes, D, 3, true);

                    m_Trf_UxVSumBuffer = m_UxVSumBuffer;
                    AllocateSumBuffer(ref m_Trf_GradUxVSumBuffer, this.m_Edgeform_GradUxV, NoOfItems, NoOfNodes, D, 3, false);
                    AllocateSumBuffer(ref m_Trf_UxGradVSumBuffer, this.m_Edgeform_UxGradV, NoOfItems, NoOfNodes, D, 3, false);
                    AllocateSumBuffer(ref m_Trf_GradUxGradVSumBuffer, this.m_Edgeform_GradUxGradV, NoOfItems, NoOfNodes, D, 4, false);
                    m_Trf_VSumBuffer = m_VSumBuffer;
                    AllocateSumBuffer(ref m_Trf_GradVSumBuffer, this.m_EdgeSourceGradV, NoOfItems, NoOfNodes, D, 3, true);
                }
            }

            Stopwatch[][] m_Edgeform_UxV_Watches;
            Stopwatch[][] m_Edgeform_GradUxV_Watches;
            Stopwatch[][] m_Edgeform_UxGradV_Watches;
            Stopwatch[][] m_Edgeform_GradUxGradV_Watches;
            Stopwatch[][] m_EdgeSourceV_Watches;
            Stopwatch[][] m_EdgeSourceGradV_Watches;

            Stopwatch Flux_Eval;
            Stopwatch Basis_Eval;
            Stopwatch Loops;
            Stopwatch ParametersAndNormals;
            Stopwatch FluxTrafo;
            Stopwatch FluxLockWait;


            // 1st index (staggered): test variable/component variable
            // 2nd index (staggered): equation component
            MultidimensionalArray[][] m_UxVComponentBuffer;
            MultidimensionalArray[][] m_UxGradVComponentBuffer;
            MultidimensionalArray[][] m_GradUxVComponentBuffer;
            MultidimensionalArray[][] m_GradUxGradVComponentBuffer;
            MultidimensionalArray[][] m_VComponentBuffer;
            MultidimensionalArray[][] m_GradVComponentBuffer;

            // 1st index: test/codomain variable
            // 2nd index: trial/domain component
            // 3rd index: IN/OUT for test/codomain
            // 4th index: IN/OUT for trial/domain
            MultidimensionalArray[,,,] m_UxVSumBuffer;
            MultidimensionalArray[,,,] m_UxGradVSumBuffer;
            MultidimensionalArray[,,,] m_GradUxVSumBuffer;
            MultidimensionalArray[,,,] m_GradUxGradVSumBuffer;

            // 1st index: test/codomain variable
            // 2nd index: trial/domain component
            // 3rd index: IN/OUT for test/codomain
            // 4th index: always 0
            MultidimensionalArray[,,,] m_VSumBuffer;
            MultidimensionalArray[,,,] m_GradVSumBuffer;

            // 1st index: test/codomain variable
            // 2nd index: trial/domain component
            // 3rd index: IN/OUT for test/codomain
            // 4th index: IN/OUT for trial/domain
            MultidimensionalArray[,,,] m_Trf_UxVSumBuffer;
            MultidimensionalArray[,,,] m_Trf_UxGradVSumBuffer;
            MultidimensionalArray[,,,] m_Trf_GradUxVSumBuffer;
            MultidimensionalArray[,,,] m_Trf_GradUxGradVSumBuffer;

            // 1st index: test/codomain variable
            // 2nd index: trial/domain component
            // 3rd index: IN/OUT for test/codomain
            // 4th index: always 0
            MultidimensionalArray[,,,] m_Trf_VSumBuffer;
            MultidimensionalArray[,,,] m_Trf_GradVSumBuffer;

            private void AllocateSumBuffer<Ee>(ref MultidimensionalArray[,,,] SumBuffer, EquationComponentArgMapping<Ee>[] FormMapping, int NoOfItems, int NoOfNodes, int D, int ArrayDim, bool affine)
                where Ee : IEquationComponent {
                if (SumBuffer == null) {
                    SumBuffer = new MultidimensionalArray[m_owner.GAMMA, affine ? 1 : m_owner.DELTA, 2, affine ? 1 : 2];
                }

                int[] LenToAlloc = new int[ArrayDim];
                LenToAlloc[0] = NoOfItems;
                LenToAlloc[1] = NoOfNodes;
                for (int k = 2; k < ArrayDim; k++) {
                    LenToAlloc[k] = D;
                }

                bool[] UsedMarker = new bool[m_owner.DELTA];
                for (int gamma = m_owner.GAMMA - 1; gamma >= 0; gamma--) {
                    var ecq = FormMapping[gamma];

                    UsedMarker.Clear();
                    for (int i = 0; i < ecq.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components
                        var comp = ecq.m_AllComponentsOfMyType[i];
                        int NoOfArgs = comp.ArgumentOrdering.Count;

                        for (int j = 0; j < NoOfArgs; j++) {
                            int targ = ecq.AllToSub[i, j];
                            UsedMarker[targ] = true;
                        }
                    }

                    if (!affine) {
                        for (int delta = m_owner.DELTA - 1; delta >= 0; delta--) {
                            for (int cr = 0; cr < 2; cr++) {

                                for (int cc = 0; cc < 2; cc++) {
                                    if (UsedMarker[delta]) {
                                        if (SumBuffer[gamma, delta, cr, cc] == null) {
                                            SumBuffer[gamma, delta, cr, cc] = new MultidimensionalArray(ArrayDim);
                                        }

                                        var mda = SumBuffer[gamma, delta, cr, cc];

                                        if (mda.GetLength(0) != NoOfItems || mda.GetLength(1) != NoOfNodes) {
                                            mda.Allocate(LenToAlloc);
                                        }
                                    }
                                }
                                //} else {

                                //}
                            }
                        }
                    } else {
                        if (ecq.m_AllComponentsOfMyType.Length > 0) {
                            for (int cr = 0; cr < 2; cr++) {
                                SumBuffer[gamma, 0, cr, 0] = new MultidimensionalArray(ArrayDim);

                                var mda = SumBuffer[gamma, 0, cr, 0];

                                if (mda.GetLength(0) != NoOfItems || mda.GetLength(1) != NoOfNodes) {
                                    mda.Allocate(LenToAlloc);
                                }
                            }
                        }
                    }
                }
            }

            private void AllocateComponentBuffer<Ee>(ref MultidimensionalArray[][] CompBuffer, EquationComponentArgMapping<Ee>[] FormMapping, int NoOfItems, int NoOfNodes, int D, int ArrayDim, bool affine)
                where Ee : IEquationComponent {

                int[] LenToAlloc = new int[ArrayDim];
                if (!affine) {
                    LenToAlloc[0] = NoOfItems;
                    LenToAlloc[1] = NoOfNodes;
                    LenToAlloc[2] = 2;
                    LenToAlloc[3] = 2;
                    LenToAlloc[4] = -1;
                    for (int k = 5; k < ArrayDim; k++)
                        LenToAlloc[k] = D;
                } else {
                    LenToAlloc[0] = NoOfItems;
                    LenToAlloc[1] = NoOfNodes;
                    LenToAlloc[2] = 2;
                    for (int k = 3; k < ArrayDim; k++)
                        LenToAlloc[k] = D;
                }

                if (CompBuffer == null) {
                    CompBuffer = FormMapping.Select(E => E.m_AllComponentsOfMyType.Select(C => new MultidimensionalArray(LenToAlloc.Length)).ToArray()).ToArray();
                }

                for (int gamma = m_owner.GAMMA - 1; gamma >= 0; gamma--) {
                    var ecq = FormMapping[gamma];

                    for (int i = 0; i < ecq.m_AllComponentsOfMyType.Length; i++) {  // loop over equation components

                        var comp = ecq.m_AllComponentsOfMyType[i];
                        int C = comp.ArgumentOrdering.Count();
                        if (!affine)
                            LenToAlloc[4] = C;

                        var mda = CompBuffer[gamma][i];
                        int Dim = mda.Dimension;
                        if (mda.GetLength(0) != NoOfItems || mda.GetLength(1) != NoOfNodes) {
                            mda.Allocate(LenToAlloc);
                        }
                    }
                }
            }

            /// <summary>
            /// Quadrature nodes in global coordinates
            /// - 1st index: edge
            /// - 2nd index: quadrature node
            /// - 3rd index: spatial direction
            /// </summary>
            MultidimensionalArray GlobalNodes = null;

            /// <summary>
            /// Guess what this is!
            /// </summary>
            MultidimensionalArray NormalBuffer;

            /// <summary>
            /// result of parameter field evaluation, for the IN-cell
            /// index: parameter field (correlates with <see cref="m_ParameterFields"/> 
            /// </summary>
            MultidimensionalArray[] m_ParameterFieldsValues_IN;

            /// <summary>
            /// result of parameter field evaluation, for the IN-cell
            /// index: parameter field (correlates with <see cref="m_ParameterFields"/> 
            /// </summary>
            MultidimensionalArray[] m_ParameterFieldsValues_OUT;

            struct MyChunk {
                public int i0;
                public int Len;
                //true for interior edges, false for a boundary edge
                public bool Type;
            }


            public void EvaluateEx(int i0, int Length, QuadRule QR, IIntegrationMetric metric, MultidimensionalArray EvalResult, int iThread, int NumOfThreads) {

                

                NodeSet NS = QR.Nodes;
                var gDat = m_owner.m_GridDat;
                bool AffineLinearEdge = gDat.iGeomEdges.IsEdgeAffineLinear(i0);
                int[,] Edge2Cell = gDat.iGeomEdges.CellIndices;
                int[,] TrafoIndx = gDat.iGeomEdges.Edge2CellTrafoIndex;
                int D = gDat.SpatialDimension;
                int NoOfTrafos = gDat.iGeomEdges.Edge2CellTrafos.Count;
                int GAMMA = m_owner.GAMMA;
                int DELTA = m_owner.DELTA;
#if DEBUG
                for (int i = 0; i < Length; i++) {
                    Debug.Assert(AffineLinearEdge == gDat.iGeomEdges.IsEdgeAffineLinear(i + i0));
                }
#endif

                // evaluate Normals, transform Nodes and Parameters
                // ================================================
                this.ParametersAndNormals.Start();
                MultidimensionalArray sqrtGram = null; // integral transformation metric
                this.GlobalNodes = gDat.GlobalNodes.GetValue_EdgeSV(NS, i0, Length);
                NormalBuffer = gDat.iGeomEdges.NormalsCache.GetNormals_Edge(NS, i0, Length);
                if (AffineLinearEdge && !metric.AlwaysUsePerNodeScaling)
                    //sqrtGram = gDat.iGeomEdges.SqrtGramian.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 });
                    sqrtGram = metric.GetScalingsForLinearElements(gDat, QR, i0, Length);
                else
                    //sqrtGram = gDat.iGeomEdges.NormalsCache.GetIntegrationMetric(NS, i0, Length);
                    sqrtGram = metric.GetScalingsForNonlinElements(gDat, QR, i0, Length);
                EvaluateParameters(i0, Length, NS);
                this.ParametersAndNormals.Stop();


                // evaluate forms
                // ==============
                #region FLUX_EVAL
                
                bool MustLock = this.m_owner.Operator.FluxesAreNOTMultithreadSafe;
                bool bLinearRequired, bAffineRequired;

                try {
                    this.Flux_Eval.Start();
                    if ( MustLock) {
                        FluxLockWait.Start();
                        Monitor.Enter(this.m_owner);
                        FluxLockWait.Stop();
                    }
                    {
                        List<MyChunk> Edges = new List<MyChunk>();
                        // sweep over edges and determine which are interior and which are boundary
                        {
                            var CellIdx = gDat.iGeomEdges.CellIndices;
                            int l = 0;
                            while ( l < Length ) {

                                // true for an interior edge, false for a boundary edge
                                MyChunk CurChunk = default(MyChunk);
                                CurChunk.Type = CellIdx[l + i0, 1] >= 0;
                                CurChunk.i0 = l + i0;
                                l++;
                                CurChunk.Len++;

                                while ( l < Length ) {
                                    bool NextType = CellIdx[l + i0, 1] >= 0;

                                    if ( NextType != CurChunk.Type )
                                        break;

                                    l++;
                                    CurChunk.Len++;
                                }

                                Edges.Add(CurChunk);
                            }
                        }

                        EdgeFormParams efp = default(EdgeFormParams);
                        efp.GridDat = gDat;
                        efp.time = m_owner.m_time;

                        bLinearRequired = LinearRequired;
                        bAffineRequired = AffineRequired;

                        // for rotationally periodic edges
                        if ( bLinearRequired )
                            DropStaleRotationSlots(Length, NS.NoOfNodes); // required for every chunk
                        var rotPeriodic = PeriodicEdgeChunk.TryCreate(this, i0, Length, gDat, NS.NoOfNodes, D); // null if none in this chunk
                        //

                        using ( rotPeriodic?.RotateOutParametersIntoInFrame() ) {
                            if ( bLinearRequired ) {
                                EvalNSumForm(ref efp, i0, m_Edgeform_UxV, Edges, m_UxVComponentBuffer, m_UxVSumBuffer, false, NS.NoOfNodes, D, m_Edgeform_UxV_Watches, MustLock,
                                    (E, mda) => E.InternalEdge_UxV(ref efp, mda),
                                    (E, mda) => E.BoundaryEdge_UxV(ref efp, mda));
                                EvalNSumForm(ref efp, i0, m_Edgeform_GradUxV, Edges, m_GradUxVComponentBuffer, m_GradUxVSumBuffer, false, NS.NoOfNodes, D, m_Edgeform_GradUxV_Watches, MustLock,
                                    (E, mda) => E.InternalEdge_GradUxV(ref efp, mda),
                                    (E, mda) => E.BoundaryEdge_GradUxV(ref efp, mda));
                                EvalNSumForm(ref efp, i0, m_Edgeform_UxGradV, Edges, m_UxGradVComponentBuffer, m_UxGradVSumBuffer, false, NS.NoOfNodes, D, m_Edgeform_UxGradV_Watches, MustLock,
                                    (E, mda) => E.InternalEdge_UxGradV(ref efp, mda),
                                    (E, mda) => E.BoundaryEdge_UxGradV(ref efp, mda));
                                EvalNSumForm(ref efp, i0, m_Edgeform_GradUxGradV, Edges, m_GradUxGradVComponentBuffer, m_GradUxGradVSumBuffer, false, NS.NoOfNodes, D, m_Edgeform_GradUxGradV_Watches, MustLock,
                                    (E, mda) => E.InternalEdge_GradUxGradV(ref efp, mda),
                                    (E, mda) => E.BoundaryEdge_GradUxGradV(ref efp, mda));
                            }
                            if ( bAffineRequired ) {
                                EvalNSumForm(ref efp, i0, m_EdgeSourceV, Edges, m_VComponentBuffer, m_VSumBuffer, true, NS.NoOfNodes, D, m_EdgeSourceV_Watches, MustLock,
                                    (E, mda) => E.InternalEdge_V(ref efp, mda),
                                    (E, mda) => E.BoundaryEdge_V(ref efp, mda));
                                EvalNSumForm(ref efp, i0, m_EdgeSourceGradV, Edges, m_GradVComponentBuffer, m_GradVSumBuffer, true, NS.NoOfNodes, D, m_EdgeSourceGradV_Watches, MustLock,
                                    (E, mda) => E.InternalEdge_GradV(ref efp, mda),
                                    (E, mda) => E.BoundaryEdge_GradV(ref efp, mda));
                            }
                        }
                        // for rotationally periodic edges
                        rotPeriodic?.ApplyRotation(bLinearRequired, bAffineRequired, efp, Edges, MustLock); // OUT-cell blocks in OUT's frame
                    }
                } finally {
                    this.Flux_Eval.Stop();
                    if ( MustLock )
                        Monitor.Exit(this.m_owner);
                }
                #endregion

                // evaluate test and trial basis functions
                // =======================================
                #region BASIS_EVAL
                this.Basis_Eval.Start();
                // 1st index
                MultidimensionalArray[] V_Xquadwgt = new MultidimensionalArray[GAMMA]; //.....test functions v times quadrature weigth
                MultidimensionalArray[] U = new MultidimensionalArray[DELTA]; //..............trial functions u
                MultidimensionalArray[] GradV_Xquadwgt = new MultidimensionalArray[GAMMA]; //.gradient of test functions v times quadrature weigth
                MultidimensionalArray[] GradU = new MultidimensionalArray[DELTA]; //..........gradient of trial functions u
                int jCellMin = int.MaxValue, jCellMax = int.MinValue;
                int NZ3 = -1, NZ4 = -1, MR = -1, NR = -1; // dimensions for intermediate result buffers
                bool JacobiRequired;
                int maxDeg;
                {
                    bool[] TestReq = new bool[GAMMA];
                    bool[] TestGradientreq = new bool[GAMMA];
                    bool[] TrialReq = new bool[DELTA];
                    bool[] TrialGradientreq = new bool[DELTA];

                    // maximum required degree for basis values/basis gradient values...
                    int MaxDeg_V = -1, MaxDeg_VGrad = -1;
                    int MaxDeg_U = -1, MaxDeg_UGrad = -1;
                    int MaxN_V = -1, MaxN_VGrad = -1;
                    int MaxN_U = -1, MaxN_UGrad = -1;

                    // determine what we really need
                    if (bLinearRequired || bAffineRequired) {
                        m_owner.TestFunctionRequired(TestReq, m_UxVSumBuffer, ref MaxDeg_V, ref MaxN_V);
                        m_owner.TestFunctionRequired(TestReq, m_GradUxVSumBuffer, ref MaxDeg_V, ref MaxN_V);
                        m_owner.TestFunctionRequired(TestGradientreq, m_UxGradVSumBuffer, ref MaxDeg_VGrad, ref MaxN_VGrad);
                        m_owner.TestFunctionRequired(TestGradientreq, m_GradUxGradVSumBuffer, ref MaxDeg_VGrad, ref MaxN_VGrad);
                        m_owner.TestFunctionRequired(TestReq, m_VSumBuffer, ref MaxDeg_V, ref MaxN_V);
                        m_owner.TestFunctionRequired(TestGradientreq, m_GradVSumBuffer, ref MaxDeg_VGrad, ref MaxN_VGrad);
                    }

                    m_owner.TrialFunctionRequired(TrialReq, m_UxVSumBuffer, bLinearRequired, ref MaxDeg_U, ref MaxN_U);
                    m_owner.TrialFunctionRequired(TrialReq, m_UxGradVSumBuffer, bLinearRequired, ref MaxDeg_U, ref MaxN_U);
                    m_owner.TrialFunctionRequired(TrialGradientreq, m_GradUxVSumBuffer, bLinearRequired, ref MaxDeg_UGrad, ref MaxN_UGrad);
                    m_owner.TrialFunctionRequired(TrialGradientreq, m_GradUxGradVSumBuffer, bLinearRequired, ref MaxDeg_UGrad, ref MaxN_UGrad);

                    // evaluate Basis
                    MultidimensionalArray Values = null, Grad = null;
                    int maxDeg_Values = Math.Max(MaxDeg_V, MaxDeg_U);
                    int maxDeg_Grad = Math.Max(MaxDeg_VGrad, MaxDeg_UGrad);
                    if (maxDeg_Values >= 0)
                        Values = gDat.ChefBasis.EdgeEval.GetValues(QR.Nodes, i0, Length, maxDeg_Values);
                    if (maxDeg_Grad >= 0)
                        Grad = gDat.ChefBasis.EdgeGradientEval.GetValues(QR.Nodes, i0, Length, maxDeg_Grad);

                    // multiply test functions with quadrature weigths
                    if ((MaxN_V > 0) && (ValuesXquadWgt.GetLength(0) != NoOfTrafos || ValuesXquadWgt.GetLength(1) != MaxN_V || ValuesXquadWgt.GetLength(2) != QR.NoOfNodes)) {
                        Debug.Assert(Values.GetLength(0) == NoOfTrafos);
                        ValuesXquadWgt.Allocate(NoOfTrafos, QR.NoOfNodes, MaxN_V);
                    }
                    if ((MaxN_VGrad > 0) && (GradsXquadWgt.GetLength(0) != NoOfTrafos || GradsXquadWgt.GetLength(1) != MaxN_V || GradsXquadWgt.GetLength(2) != QR.NoOfNodes)) {
                        Debug.Assert(Grad.GetLength(0) == NoOfTrafos);
                        GradsXquadWgt.Allocate(NoOfTrafos, QR.NoOfNodes, MaxN_VGrad, D);
                    }

                    bool[] marker_TstFuncXwgt = new bool[NoOfTrafos];
                    int[] I0 = new int[3], IE = new int[3], _I0 = new int[4], _IE = new int[4];
                    for (int i = 0; i < Length; i++) {
                        int iEdge = i + i0;

                        for (int ii = 0; ii < 2; ii++) {
                            int jCell = Edge2Cell[iEdge, ii];
                            if (jCell >= 0) {
                                jCellMin = Math.Min(jCell, jCellMin);
                                jCellMax = Math.Max(jCell, jCellMax);

                                int iTrfo = TrafoIndx[iEdge, ii];
                                if (!marker_TstFuncXwgt[iTrfo]) {
                                    if (MaxN_V > 0) {
                                        I0[0] = iTrfo; IE[0] = iTrfo - 1;
                                        I0[1] = 0; IE[1] = QR.NoOfNodes - 1;
                                        I0[2] = 0; IE[2] = MaxN_V - 1;
                                        ValuesXquadWgt.ExtractSubArrayShallow(I0, IE).Multiply(1.0, QR.Weights, Values.ExtractSubArrayShallow(I0, IE), 0.0, ref mp_kn_k_kn);
                                    }
                                    if (MaxN_VGrad > 0) {
                                        _I0[0] = iTrfo; _IE[0] = iTrfo - 1;
                                        _I0[1] = 0; _IE[1] = QR.NoOfNodes - 1;
                                        _I0[2] = 0; _IE[2] = MaxN_VGrad - 1;
                                        _I0[3] = 0; _IE[3] = D - 1;
                                        GradsXquadWgt.ExtractSubArrayShallow(_I0, _IE).Multiply(1.0, QR.Weights, Grad.ExtractSubArrayShallow(_I0, _IE), 0.0, ref mp_knd_k_knd);
                                    }
                                    marker_TstFuncXwgt[iTrfo] = true;
                                }
                            }
                        }
                    }

                    // extract data which fits da components
                    for (int delta = 0; delta < DELTA; delta++) {
                        if (TrialReq[delta])
                            U[delta] = Values.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Values.GetLength(0) - 1, QR.NoOfNodes - 1, m_owner.m_ColL[delta] - 1 });
                        if (TrialGradientreq[delta])
                            GradU[delta] = Grad.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Grad.GetLength(0) - 1, QR.NoOfNodes - 1, m_owner.m_ColL[delta] - 1, D - 1 });
                    }
                    for (int gamma = 0; gamma < GAMMA; gamma++) {
                        if (TestReq[gamma])
                            V_Xquadwgt[gamma] = ValuesXquadWgt.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Values.GetLength(0) - 1, QR.NoOfNodes - 1, m_owner.m_RowL[gamma] - 1 });
                        if (TestGradientreq[gamma])
                            GradV_Xquadwgt[gamma] = GradsXquadWgt.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Grad.GetLength(0) - 1, QR.NoOfNodes - 1, m_owner.m_RowL[gamma] - 1, D - 1 });
                    }

                    // blav
                    NZ3 = Math.Max(MaxN_VGrad, Math.Max(MaxN_V, MaxN_U));
                    NZ4 = MaxN_UGrad;
                    MR = Math.Max(MaxN_V, MaxN_VGrad);
                    NR = Math.Max(MaxN_U, MaxN_UGrad);
                    JacobiRequired = (MaxN_VGrad >= 0) || (MaxDeg_UGrad >= 0);
                    maxDeg = Math.Max(MaxDeg_V, Math.Max(MaxDeg_VGrad, Math.Max(MaxDeg_U, MaxDeg_UGrad)));
                }
                this.Basis_Eval.Stop();
                #endregion

                // transform fluxes
                // ================
                #region FLUX_TRAFO
                this.FluxTrafo.Start();
                {
                    MultidimensionalArray invJacobi = null;
                    MultidimensionalArray[] invJacobiInOt = null;
                    if (JacobiRequired) {
                        if (AffineLinearEdge ) {
                            invJacobi = gDat.iGeomCells.InverseTransformation;
                        } else {
                            invJacobiInOt = gDat.InverseJacobian.GetValue_EdgeDV(QR.Nodes, i0, Length).TupleToArray();
                        }
                    }

                    unsafe {
                        fixed (int* pEdge2Cell = Edge2Cell) {

                            for (int gamma = 0; gamma < GAMMA; gamma++) {

                                // linear part (matrix)
                                if (bLinearRequired) {
                                    for (int delta = 0; delta < DELTA; delta++) {
                                        for (int cr = 0; cr < 2; cr++) { // IN/OUT row/test function loop (variables V_in, V_out);
                                            MultidimensionalArray invJacobi_V = JacobiRequired ? (AffineLinearEdge ? invJacobi : invJacobiInOt[cr]) : null;

                                            for (int cc = 0; cc < 2; cc++) { // IN/OUT column/trial function loop (variables U_in, U_out);
                                                MultidimensionalArray invJacobi_U = JacobiRequired ? (AffineLinearEdge ? invJacobi : invJacobiInOt[cc]) : null;

                                                if (m_UxVSumBuffer[gamma, delta, cr, cc] != null) {
                                                    Debug.Assert(object.ReferenceEquals(m_Trf_UxVSumBuffer[gamma, delta, cr, cc], m_UxVSumBuffer[gamma, delta, cr, cc]));
                                                    if (!AffineLinearEdge || metric.AlwaysUsePerNodeScaling)
                                                        m_Trf_UxVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, sqrtGram, m_Trf_UxVSumBuffer[gamma, delta, cr, cc], 0.0, "ik", "ik", "ik");
                                                }
                                                if (m_UxGradVSumBuffer[gamma, delta, cr, cc] != null) {

                                                    // test function gradient ( \/V ) transform:
                                                    if (AffineLinearEdge) {
                                                        // linear edge
                                                        m_Trf_UxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, invJacobi, m_UxGradVSumBuffer[gamma, delta, cr, cc], 0.0, ref mp_ika_Tiad_ikd,
                                                            pEdge2Cell, pEdge2Cell,
                                                            trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);

                                                    } else {
                                                        // curved edge
                                                        m_Trf_UxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, invJacobi_V, m_UxGradVSumBuffer[gamma, delta, cr, cc], 0.0, ref mp_ika_ikad_ikd);
                                                    }

                                                    if (!AffineLinearEdge || metric.AlwaysUsePerNodeScaling)
                                                        m_Trf_UxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, sqrtGram, m_Trf_UxGradVSumBuffer[gamma, delta, cr, cc], 0.0, "ikd", "ik", "ikd");

                                                }
                                                if (m_GradUxVSumBuffer[gamma, delta, cr, cc] != null) {

                                                    // trial function gradient ( \/U ) transform:
                                                    if (AffineLinearEdge) {
                                                        m_Trf_GradUxVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, m_GradUxVSumBuffer[gamma, delta, cr, cc], invJacobi, 0.0, ref mp_ikb_ikd_Tibd,
                                                            pEdge2Cell, pEdge2Cell,
                                                            trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    } else {
                                                        m_Trf_GradUxVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, m_GradUxVSumBuffer[gamma, delta, cr, cc], invJacobi_U, 0.0, ref mp_ikb_ikd_ikbd);
                                                    }

                                                    if (!AffineLinearEdge || metric.AlwaysUsePerNodeScaling)
                                                        m_Trf_GradUxVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, sqrtGram, m_Trf_GradUxVSumBuffer[gamma, delta, cr, cc], 0.0, "ikd", "ik", "ikd");

                                                }
                                                if (m_GradUxGradVSumBuffer[gamma, delta, cr, cc] != null) {

                                                    // test function gradient ( \/V ) transform:
                                                    if (AffineLinearEdge) {
                                                        m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, invJacobi, m_GradUxGradVSumBuffer[gamma, delta, cr, cc], 0.0, ref mp_ikae_Tiad_ikde,
                                                            pEdge2Cell, pEdge2Cell,
                                                            trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);

                                                    } else {
                                                        m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, invJacobi_V, m_GradUxGradVSumBuffer[gamma, delta, cr, cc], 0.0, ref mp_ikae_ikad_ikde);
                                                    }

                                                    // swap the buffers for the second transform ...
                                                    MultidimensionalArray tmp = m_GradUxGradVSumBuffer[gamma, delta, cr, cc];
                                                    m_GradUxGradVSumBuffer[gamma, delta, cr, cc] = m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc];
                                                    m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc] = tmp;

                                                    // trial function gradient ( \/U ) transform:
                                                    if (AffineLinearEdge) {
                                                        m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, m_GradUxGradVSumBuffer[gamma, delta, cr, cc], invJacobi, 0.0, ref mp_ikab_ikae_Tibe,
                                                            pEdge2Cell, pEdge2Cell,
                                                            trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    } else {
                                                        m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, m_GradUxGradVSumBuffer[gamma, delta, cr, cc], invJacobi_U, 0.0, ref mp_ikab_ikae_ikbe);
                                                    }


                                                    // apply integral metric (only for curved edges)
                                                    if (!AffineLinearEdge || metric.AlwaysUsePerNodeScaling)
                                                        m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc].Multiply(1.0, sqrtGram, m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc], 0.0, "ikde", "ik", "ikde");
                                                }
                                            }
                                        }
                                    }
                                }

                                // affine part
                                if (bAffineRequired) {
                                    for (int cr = 0; cr < 2; cr++) {
                                        MultidimensionalArray invJacobi_V = JacobiRequired ? (AffineLinearEdge ? invJacobi : invJacobiInOt[cr]) : null;

                                        if (m_VSumBuffer[gamma, 0, cr, 0] != null) {
                                            Debug.Assert(object.ReferenceEquals(m_Trf_VSumBuffer[gamma, 0, cr, 0], m_VSumBuffer[gamma, 0, cr, 0]));
                                            if (!AffineLinearEdge || metric.AlwaysUsePerNodeScaling)
                                                m_Trf_VSumBuffer[gamma, 0, cr, 0].Multiply(1.0, sqrtGram, m_Trf_VSumBuffer[gamma, 0, cr, 0], 0.0, "ik", "ik", "ik");
                                        }

                                        if (m_GradVSumBuffer[gamma, 0, cr, 0] != null) {
                                            // test function gradient ( \/V ) transform:
                                            if (AffineLinearEdge && !metric.AlwaysUsePerNodeScaling) {
                                                m_Trf_GradVSumBuffer[gamma, 0, cr, 0].Multiply(1.0, invJacobi, m_GradVSumBuffer[gamma, 0, cr, 0], 0.0, ref mp_ika_Tiad_ikd,
                                                    pEdge2Cell, pEdge2Cell,
                                                    trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);

                                            } else {
                                                m_Trf_GradVSumBuffer[gamma, 0, cr, 0].Multiply(1.0, invJacobi_V, m_GradVSumBuffer[gamma, 0, cr, 0], 0.0, ref mp_ika_ikad_ikd);
                                            }

                                            if (!AffineLinearEdge || metric.AlwaysUsePerNodeScaling)
                                                m_Trf_GradVSumBuffer[gamma, 0, cr, 0].Multiply(1.0, sqrtGram, m_Trf_GradVSumBuffer[gamma, 0, cr, 0], 0.0, "ikd", "ik", "ikd");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                this.FluxTrafo.Stop();
                #endregion

                // do the summation
                // ================
                #region QUADLOOPS
                this.Loops.Start();
                unsafe {
                    fixed (int* pTrafoIndx = TrafoIndx, pEdge2Cell = Edge2Cell) {
                        MultidimensionalArray _R = null, _Z = null, _Q = null;
                        MultidimensionalArray z3Buf = null, z4Buf = null, RLbuf = null, RAbuf = null;
                        int iz3Buf = 0, iz4buf = 0, iRLbuf = 0, iRAbuf = 0;

                        MultidimensionalArray trafo = null, basisScale = null;
                        if (AffineLinearEdge) {
                            basisScale = gDat.ChefBasis.Scaling;
                        } else {
                            if (maxDeg >= 0)
                                trafo = gDat.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(jCellMin, jCellMax - jCellMin + 1, maxDeg);
                        }


                        int I0Row = 0;
                        //var Edge2Cell = this.m_GridDat.Edges.CellIndices;
                        //int Jup = this.m_GridDat.Cells.NoOfLocalUpdatedCells;

                        for(int gamma = 0; gamma < GAMMA; gamma++) { // loop over codomain variables V
                            int I0Col = 0;

                            // linear part (matrix)
                            // --------------------
                            if(bLinearRequired) {
                                for(int delta = 0; delta < DELTA; delta++) { // loop over domain variables U
                                    for(int cr = 0; cr < 2; cr++) { // in and out for test functions V
                                        for(int cc = 0; cc < 2; cc++) { // in and out for trial functions U
                                            bool bAny = (m_UxVSumBuffer[gamma, delta, cr, cc] != null)
                                                || (m_GradUxGradVSumBuffer[gamma, delta, cr, cc] != null)
                                                || (m_UxGradVSumBuffer[gamma, delta, cr, cc] != null)
                                                || (m_GradUxVSumBuffer[gamma, delta, cr, cc] != null);
                                            if(bAny && m_owner.m_RowL[gamma] > 0 && m_owner.m_ColL[delta] > 0) {
                                                GetQRbufferLinear(Length, EvalResult, cr, cc, AffineLinearEdge && !metric.AlwaysUsePerNodeScaling, MR, NR, ref iRLbuf, ref RLbuf, I0Row, gamma, I0Col, delta, out _R, out _Q);

                                                //var res = EvalResult.ExtractSubArrayShallow(new int[] { 0, 0, cr, cc, I0Row, I0Col }, new int[] { Length - 1, NS.NoOfNodes - 1, cr - 1, cc - 1, I0Row + m_RowL[gamma] - 1, I0Col + m_ColL[delta] - 1 });

                                                double cF = 0;

                                                if(m_UxVSumBuffer[gamma, delta, cr, cc] != null) {
                                                    GetZbuffer(Length, QR.NoOfNodes, -1, NZ3, ref iz3Buf, ref z3Buf, out _Z, m_owner.m_ColL[delta]);

                                                    _Z.Multiply(1.0, m_Trf_UxVSumBuffer[gamma, delta, cr, cc], U[delta], 0.0, ref mp_ikn_ik_Tikn,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    _R.Multiply(1.0, V_Xquadwgt[gamma], _Z, cF, ref mp_imn_Tikm_ikn,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                                    cF = 1.0;
                                                }
                                                if(m_GradUxGradVSumBuffer[gamma, delta, cr, cc] != null) {
                                                    GetZbuffer(Length, QR.NoOfNodes, D, NZ4, ref iz4buf, ref z4Buf, out _Z, m_owner.m_ColL[delta]);

                                                    _Z.Multiply(1.0, m_Trf_GradUxGradVSumBuffer[gamma, delta, cr, cc], GradU[delta], 0.0, ref mp_ikna_ikab_Tiknb,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    _R.Multiply(1.0, GradV_Xquadwgt[gamma], _Z, cF, ref mp_imn_Tikma_ikna,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                                    cF = 1.0;
                                                }
                                                if(m_UxGradVSumBuffer[gamma, delta, cr, cc] != null) {
                                                    GetZbuffer(Length, QR.NoOfNodes, -1, NZ3, ref iz3Buf, ref z3Buf, out _Z, m_owner.m_RowL[gamma]);

                                                    _Z.Multiply(1.0, m_Trf_UxGradVSumBuffer[gamma, delta, cr, cc], GradV_Xquadwgt[gamma], 0.0, ref mp_ikm_ika_Tikma,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cr), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    _R.Multiply(1.0, U[delta], _Z, cF, ref mp_imn_Tikn_ikm,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: (2 * i0 + cc), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                                    cF = 1.0;
                                                }
                                                if(m_GradUxVSumBuffer[gamma, delta, cr, cc] != null) {
                                                    GetZbuffer(Length, QR.NoOfNodes, -1, NZ3, ref iz3Buf, ref z3Buf, out _Z, m_owner.m_ColL[delta]);

                                                    _Z.Multiply(1.0, m_Trf_GradUxVSumBuffer[gamma, delta, cr, cc], GradU[delta], 0.0, ref mp_ikn_ikb_Tiknb,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    _R.Multiply(1.0, V_Xquadwgt[gamma], _Z, cF, ref mp_imn_Tikm_ikn,
                                                        pTrafoIndx, pTrafoIndx,
                                                        trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                                    cF = 1.0;

                                                }

                                                if(AffineLinearEdge) {
                                                    Debug.Assert(object.ReferenceEquals(_R.Storage, EvalResult.Storage));

                                                    // apply basis-transformation `basisScale` for affine-linear cells (a **scalar** per element), only test functions
                                                    _R.Multiply(1.0, _R, basisScale, 0.0, ref mp_imn_imn_Ti,
                                                        pEdge2Cell, pEdge2Cell,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cr), trfCycle_B: 2, trfPostOffset_B: 0);
                                                    _R.Multiply(1.0, _R, basisScale, 0.0, ref mp_imn_imn_Ti,
                                                        pEdge2Cell, pEdge2Cell,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: 0);

                                                    if(!metric.AlwaysUsePerNodeScaling) {
                                                        // apply `sqrtGram` (integration metric **per element**), in-place
                                                        _R.Multiply(1.0, _R, sqrtGram, 0.0, ref mp_imn_imn_i);
                                                    } else {
                                                        // `sqrtGram` (integration metric **per node**) has already been applied 
                                                        // noop, `_R` already contains values of (integrand * metric)
                                                    }

                                                } else {
                                                    Debug.Assert(!object.ReferenceEquals(EvalResult.Storage, _Q.Storage));
                                                    Debug.Assert(object.ReferenceEquals(EvalResult.Storage, _R.Storage));

                                                    // apply basis-transformation `trafo` for curved cells (**a matrix per element**), test and trial functions
                                                    // `sqrtGram` (integration metric **per node**) has already been applied
                                                    _Q.Multiply(1.0, _R, ExtractTrafo(trafo, m_owner.m_ColL[delta]), 0.0, ref mp_ibn_iba_Tian,
                                                    pEdge2Cell, pEdge2Cell,
                                                        trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cc), trfCycle_B: 2, trfPostOffset_B: -jCellMin);
                                                    _R.Multiply(1.0, ExtractTrafo(trafo, m_owner.m_RowL[gamma]), _Q, 0.0, ref mp_imn_Tibm_ibn,
                                                        pEdge2Cell, pEdge2Cell,
                                                        trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: -jCellMin, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);

                                                }
                                            }
                                        }
                                    }

                                    I0Col += m_owner.m_ColL[delta];
                                }
                            }

                            // affine vector
                            // -------------
                            if(bAffineRequired) {

                                for(int cr = 0; cr < 2; cr++) {
                                    var bAny = (m_VSumBuffer[gamma, 0, cr, 0] != null) || (m_GradVSumBuffer[gamma, 0, cr, 0] != null);
                                    if(bAny && m_owner.m_RowL[gamma] > 0) {
                                        GetRQbufferAffine(Length, EvalResult, cr, AffineLinearEdge, MR, ref iRAbuf, ref RAbuf, I0Row, gamma, I0Col, out _R, out _Q);

                                        double cF = 0;
                                        if(m_VSumBuffer[gamma, 0, cr, 0] != null) {

                                            _R.Multiply(1.0, m_Trf_VSumBuffer[gamma, 0, cr, 0], V_Xquadwgt[gamma], cF, ref mp_im_ik_Tkm,
                                                pTrafoIndx, pTrafoIndx,
                                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cr), trfCycle_B: 2, trfPostOffset_B: 0);
                                            cF = 1.0;
                                        }
                                        if(m_GradVSumBuffer[gamma, 0, cr, 0] != null) {

                                            _R.Multiply(1.0, m_Trf_GradVSumBuffer[gamma, 0, cr, 0], GradV_Xquadwgt[gamma], cF, ref mp_im_ik_Tikma,
                                                pTrafoIndx, pTrafoIndx,
                                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cr), trfCycle_B: 2, trfPostOffset_B: 0);

                                            cF = 1.0;
                                        }

                                        if(AffineLinearEdge) {
                                            Debug.Assert(object.ReferenceEquals(EvalResult.Storage, _R.Storage));

                                            // apply basis-transformation `basisScale` for affine-linear edges (a **scalar** per element), only test functions
                                            _R.Multiply(1.0, _R, basisScale, 0.0, ref mp_im_im_Ti,
                                                pEdge2Cell, pEdge2Cell,
                                                trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + cr), trfCycle_B: 2, trfPostOffset_B: 0);

                                            if(!metric.AlwaysUsePerNodeScaling) {
                                                // apply `sqrtGram` (integration metric **per element**), in-place
                                                _R.Multiply(1.0, _R, sqrtGram, 0.0, ref mp_im_im_i);
                                            } else {
                                                // `sqrtGram` (integration metric **per node**) has already been applied 
                                                // noop, `_R` already contains values of (integrand * metric)
                                            }

                                        } else {
                                            Debug.Assert(object.ReferenceEquals(EvalResult.Storage, _Q.Storage));

                                            // apply basis-transformation `trafo` for curved edges (a **matrix** per element), only test functions
                                            // `sqrtGram` (integration metric **per node**) has already been applied 
                                            Debug.Assert(object.ReferenceEquals(EvalResult.Storage, _Q.Storage));
                                            Debug.Assert(!object.ReferenceEquals(EvalResult.Storage, _R.Storage));
                                            _Q.Multiply(1.0, ExtractTrafo(trafo, m_owner.m_RowL[gamma]), _R, 0.0, ref mp_jn_Tjmn_jm,
                                                pEdge2Cell, pEdge2Cell,
                                                trfPreOffset_A: (2 * i0 + cr), trfCycle_A: 2, trfPostOffset_A: -jCellMin, trfPostOffset_B: 0, trfCycle_B: 0, trfPreOffset_B: 0);


                                        }
                                    }
                                }

                            }

                            // inc
                            // ---
                            I0Row += m_owner.m_RowL[gamma];
                        }

                        // free temp mem.
                        if (z3Buf != null)
                            TempBuffer.FreeTempBuffer(iz3Buf);
                        if (z4Buf != null)
                            TempBuffer.FreeTempBuffer(iz4buf);
                        if (RLbuf != null)
                            TempBuffer.FreeTempBuffer(iRLbuf);
                        if (RAbuf != null)
                            TempBuffer.FreeTempBuffer(iRAbuf);
                    }

                }
                this.Loops.Stop();
                

                #endregion

            }

            #region helper functions for rotated periodic boundaries
            /// <summary>
            /// Deep copy of a sum buffer (null entries stay null).
            /// </summary>
            private static MultidimensionalArray[, , ,] CloneSumBufferStructure(MultidimensionalArray[, , ,] src) {
                int n0 = src.GetLength(0), n1 = src.GetLength(1), n2 = src.GetLength(2), n3 = src.GetLength(3);
                var dst = new MultidimensionalArray[n0, n1, n2, n3];
                for (int a = 0; a < n0; a++)
                    for (int b = 0; b < n1; b++)
                        for (int c = 0; c < n2; c++)
                            for (int d = 0; d < n3; d++)
                                if (src[a, b, c, d] != null)
                                    dst[a, b, c, d] = src[a, b, c, d].CloneAs();
                return dst;
            }

            /// <summary>
            /// Copies the cr=1 (OUT-cell equation) blocks from <paramref name="scratch"/> into <paramref name="main"/>,
            /// restricted to the edges in <paramref name="periodicIdx"/>.
            /// </summary>
            /// <param name="scratch">Result of the rotated-normal re-evaluation.</param>
            /// <param name="main">Main (native-normal) buffer to splice the cr=1 block into.</param>
            /// <param name="GAMMA">Number of codomain variables.</param>
            /// <param name="NoOfNodes">Number of quadrature nodes per edge.</param>
            /// <param name="periodicIdx">Indices of the periodic edges in the current chunk.</param>
            private void SpliceCr1Block(MultidimensionalArray[, , ,] scratch, MultidimensionalArray[, , ,] main, int GAMMA, int NoOfNodes, List<int> periodicIdx) {
                int DELTA = scratch.GetLength(1);
                int CC = scratch.GetLength(3); // 2 for the bilinear (UxV etc.) buffers, 1 for the affine (EdgeSourceV etc.) buffers
                for (int gamma = 0; gamma < GAMMA; gamma++) {
                    for (int delta = 0; delta < DELTA; delta++) {
                        for (int cc = 0; cc < CC; cc++) {
                            var src = scratch[gamma, delta, 1, cc];
                            var dst = main[gamma, delta, 1, cc];
                            if (src == null || dst == null)
                                continue; // nothing to copy
                            int Dim = src.Dimension;
                            int extra1 = Dim >= 3 ? src.GetLength(2) : 1;
                            int extra2 = Dim >= 4 ? src.GetLength(3) : 1;
                            int[] idx = new int[Dim];
                            for (int k = 0; k < periodicIdx.Count; k++) {
                                idx[0] = periodicIdx[k];
                                for (int node = 0; node < NoOfNodes; node++) {
                                    idx[1] = node;
                                    for (int e1 = 0; e1 < extra1; e1++) {
                                        if (Dim >= 3) idx[2] = e1;
                                        for (int e2 = 0; e2 < extra2; e2++) {
                                            if (Dim >= 4) idx[3] = e2;
                                            dst[idx] = src[idx];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// Sets all entries of <paramref name="SumBuffer"/> to null which do not match the current chunk size (<paramref name="Length"/>, <paramref name="NoOfNodes"/>).
            /// </summary>
            private static void InvalidateStaleVectorTupleSlots(MultidimensionalArray[, , ,] SumBuffer, int Length, int NoOfNodes) {
                if (SumBuffer == null)
                    return;
                int n0 = SumBuffer.GetLength(0), n1 = SumBuffer.GetLength(1), n2 = SumBuffer.GetLength(2), n3 = SumBuffer.GetLength(3);
                for (int a = 0; a < n0; a++)
                    for (int b = 0; b < n1; b++)
                        for (int c = 0; c < n2; c++)
                            for (int d = 0; d < n3; d++) {
                                var mda = SumBuffer[a, b, c, d];
                                if (mda != null && (mda.GetLength(0) != Length || mda.GetLength(1) != NoOfNodes))
                                    SumBuffer[a, b, c, d] = null;
                            }
            }

            private static bool LengthsMatch(MultidimensionalArray a, int[] wantLengths) {
                if (a.Dimension != wantLengths.Length)
                    return false;
                for (int k = 0; k < wantLengths.Length; k++)
                    if (a.GetLength(k) != wantLengths[k])
                        return false;
                return true;
            }

            /// <summary>
            /// Rotates the (cr,cc) block of all vector-tuples of domain variables by <paramref name="Mtx"/>, on edge <paramref name="i"/>.
            /// Missing or wrongly sized entries of <paramref name="SumBuffer"/> and <paramref name="TrfBuffer"/> are (re-)allocated.
            /// </summary>
            private void RotateCrossBlock(MultidimensionalArray[, , ,] SumBuffer, MultidimensionalArray[, , ,] TrfBuffer, int[][] VectorTuples, int GAMMA, int Length, int NoOfNodes, MultidimensionalArray Mtx, int i, int cr, int cc) {
                bool aliased = object.ReferenceEquals(SumBuffer, TrfBuffer);

                foreach (int[] tuple in VectorTuples) {
                    int Dloc = tuple.Length;

                    for (int gamma = 0; gamma < GAMMA; gamma++) {
                        MultidimensionalArray template = null;
                        for (int d = 0; d < Dloc; d++) {
                            var b = SumBuffer[gamma, tuple[d], cr, cc];
                            if (b != null && b.GetLength(0) == Length && b.GetLength(1) == NoOfNodes) {
                                template = b; // a correctly sized slot
                                break;
                            }
                        }
                        if (template == null) {
                            // fallback: any non-null slot, only used for the shape of the extra axes
                            for (int d = 0; d < Dloc; d++) {
                                var b = SumBuffer[gamma, tuple[d], cr, cc];
                                if (b != null) {
                                    template = b;
                                    break;
                                }
                            }
                        }
                        if (template == null)
                            continue; // nothing to rotate for this (gamma,cr,cc)

                        int[] wantLengths = template.Lengths;
                        wantLengths[0] = Length;
                        wantLengths[1] = NoOfNodes;

                        var bufs = new MultidimensionalArray[Dloc];
                        for (int d = 0; d < Dloc; d++) {
                            var b = SumBuffer[gamma, tuple[d], cr, cc];
                            if (b == null || !LengthsMatch(b, wantLengths)) {
                                b = MultidimensionalArray.Create(wantLengths);
                                SumBuffer[gamma, tuple[d], cr, cc] = b;
                            }
                            bufs[d] = b;

                            if (!aliased && (TrfBuffer[gamma, tuple[d], cr, cc] == null || !LengthsMatch(TrfBuffer[gamma, tuple[d], cr, cc], wantLengths))) {
                                TrfBuffer[gamma, tuple[d], cr, cc] = MultidimensionalArray.Create(wantLengths);
                            }
                        }

                        int Dim = template.Dimension;
                        int extra1 = Dim >= 3 ? template.GetLength(2) : 1;
                        int extra2 = Dim >= 4 ? template.GetLength(3) : 1;

                        double[] vOld = new double[Dloc];
                        int[] idx = new int[Dim];
                        idx[0] = i;
                        for (int k = 0; k < NoOfNodes; k++) {
                            idx[1] = k;
                            for (int e1 = 0; e1 < extra1; e1++) {
                                if (Dim >= 3) idx[2] = e1;
                                for (int e2 = 0; e2 < extra2; e2++) {
                                    if (Dim >= 4) idx[3] = e2;

                                    for (int d = 0; d < Dloc; d++)
                                        vOld[d] = bufs[d][idx];

                                    for (int d = 0; d < Dloc; d++) {
                                        double acc = 0.0;
                                        for (int e = 0; e < Dloc; e++)
                                            acc += Mtx[d, e] * vOld[e];
                                        bufs[d][idx] = acc;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// Rotates the derivative-direction axis (at position <paramref name="gradAxisPos"/>) of the (cr,cc) block
            /// of all variables by <paramref name="Mtx"/>, on edge <paramref name="i"/>; used for the GradU axis of the GradUxV and GradUxGradV buffers.
            /// </summary>
            private void RotateGradientAxis(MultidimensionalArray[, , ,] SumBuffer, int GammaCount, int DeltaCount, int NoOfNodes, int gradAxisPos, MultidimensionalArray Mtx, int i, int cr, int cc) {
                for (int gamma = 0; gamma < GammaCount; gamma++) {
                    for (int delta = 0; delta < DeltaCount; delta++) {
                        var b = SumBuffer[gamma, delta, cr, cc];
                        if (b == null)
                            continue;

                        int Dim = b.Dimension;
                        int Dloc = b.GetLength(gradAxisPos);
                        int otherExtraPos = (Dim >= 4) ? (gradAxisPos == 2 ? 3 : 2) : -1;
                        int otherExtraLen = otherExtraPos >= 0 ? b.GetLength(otherExtraPos) : 1;

                        double[] vOld = new double[Dloc];
                        int[] idx = new int[Dim];
                        idx[0] = i;
                        for (int k = 0; k < NoOfNodes; k++) {
                            idx[1] = k;
                            for (int eo = 0; eo < otherExtraLen; eo++) {
                                if (otherExtraPos >= 0) idx[otherExtraPos] = eo;

                                for (int bIdx = 0; bIdx < Dloc; bIdx++) {
                                    idx[gradAxisPos] = bIdx;
                                    vOld[bIdx] = b[idx];
                                }
                                for (int bIdx = 0; bIdx < Dloc; bIdx++) {
                                    double acc = 0.0;
                                    for (int e = 0; e < Dloc; e++)
                                        acc += Mtx[bIdx, e] * vOld[e];
                                    idx[gradAxisPos] = bIdx;
                                    b[idx] = acc;
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// Groups the linearization-point gradient parameters "&lt;Comp&gt;_lin_d[&lt;dir&gt;]" into tuples of fixed direction and varying component.
            /// </summary>
            /// <remarks>
            /// Together with the same-component tuples found by <see cref="BoSSS.Foundation.Quadrature.PeriodicBoundaryUtils.GetVectorFieldIndices"/>,
            /// rotating both groupings yields the tensor rotation M*T*M^T of the gradient.
            /// </remarks>
            private int[][] GetLinDCrossComponentTuples(int D) {
                var names = m_owner.Operator.ParameterVar;
                const string marker = "_lin_d[";
                var byDirection = new Dictionary<int, List<int>>();
                for (int i = 0; i < names.Count; i++) {
                    string n = names[i];
                    int b = n.LastIndexOf(marker, StringComparison.Ordinal);
                    if (b < 0)
                        continue;
                    int closeBracket = n.IndexOf(']', b);
                    if (closeBracket < 0)
                        continue;
                    string dirDigits = n.Substring(b + marker.Length, closeBracket - (b + marker.Length));
                    if (!int.TryParse(dirDigits, System.Globalization.NumberStyles.Integer, System.Globalization.CultureInfo.InvariantCulture, out int dir))
                        continue;
                    if (!byDirection.TryGetValue(dir, out var list)) {
                        list = new List<int>();
                        byDirection[dir] = list;
                    }
                    list.Add(i);
                }

                var result = new List<int[]>();
                foreach (var kv in byDirection) {
                    var origIdx = kv.Value;
                    var baseNames = origIdx.Select(i => {
                        string n = names[i];
                        int b = n.LastIndexOf(marker, StringComparison.Ordinal);
                        return n.Substring(0, b);
                    }).ToArray();
                    var subTuples = BoSSS.Foundation.Quadrature.PeriodicBoundaryUtils.GetVectorFieldIndices(baseNames, D);
                    foreach (var st in subTuples)
                        result.Add(st.Select(localIdx => origIdx[localIdx]).ToArray());
                }
                return result.ToArray();
            }

            /// <summary>
            /// Vector tuples of the parameters, including the cross-component tuples of the "_lin_d" gradients.
            /// </summary>
            private int[][] GetParameterVectorTuples(int D) {
                var tuples = BoSSS.Foundation.Quadrature.PeriodicBoundaryUtils.GetVectorFieldIndices(m_owner.Operator.ParameterVar, D).ToArray();
                var crossTuples = GetLinDCrossComponentTuples(D);
                if (crossTuples.Length > 0)
                    tuples = tuples.Concat(crossTuples).ToArray();
                return tuples;
            }

            /// <summary>
            /// Rotates all vector tuples of <paramref name="paramValues"/> on the edges <paramref name="periodicIdx"/> by <paramref name="Mtxs"/>;
            /// returns the original arrays, for <see cref="RestoreParameters"/>.
            /// </summary>
            private static MultidimensionalArray[] RotateParameters(MultidimensionalArray[] paramValues, int[][] tuples, List<int> periodicIdx, List<MultidimensionalArray> Mtxs, int NoOfNodes) {
                var original = (MultidimensionalArray[])paramValues.Clone();
                foreach (var tuple in tuples)
                    RotateParamSlotsOnPeriodicEdges(paramValues, tuple, periodicIdx, Mtxs, NoOfNodes);
                return original;
            }

            /// <summary>
            /// Undoes <see cref="RotateParameters"/>.
            /// </summary>
            private static void RestoreParameters(MultidimensionalArray[] paramValues, MultidimensionalArray[] original) {
                Array.Copy(original, paramValues, original.Length);
            }

            /// <summary>
            /// Rotates the parameter values of one vector tuple on the edges <paramref name="periodicIdx"/> by the matrices <paramref name="Mtxs"/>.
            /// </summary>
            /// <remarks>
            /// Works on clones, which are assigned to <paramref name="paramValues"/>; the original arrays may be shared with domain-variable values and must not be modified.
            /// </remarks>
            private static void RotateParamSlotsOnPeriodicEdges(MultidimensionalArray[] paramValues, int[] tuple, List<int> periodicIdx, List<MultidimensionalArray> Mtxs, int NoOfNodes) {
                int Dloc = tuple.Length;
                var bufs = new MultidimensionalArray[Dloc];
                for (int d = 0; d < Dloc; d++) {
                    bufs[d] = paramValues[tuple[d]].CloneAs();
                    paramValues[tuple[d]] = bufs[d];
                }

                double[] vOld = new double[Dloc];
                for (int k = 0; k < periodicIdx.Count; k++) {
                    int i = periodicIdx[k];
                    var M = Mtxs[k];
                    for (int node = 0; node < NoOfNodes; node++) {
                        for (int d = 0; d < Dloc; d++)
                            vOld[d] = bufs[d][i, node];
                        for (int d = 0; d < Dloc; d++) {
                            double acc = 0.0;
                            for (int e = 0; e < Dloc; e++)
                                acc += M[d, e] * vOld[e];
                            bufs[d][i, node] = acc;
                        }
                    }
                }
            }

            /// <summary>
            /// Returns the periodic edges in the current chunk and, for each, the rotation matrix from IN's frame to OUT's frame.
            /// </summary>
            private void GetPeriodicEdgesAndMatrices(int i0, int Length, IGridData grid, out List<int> periodicIdx, out List<MultidimensionalArray> periodicMtx) {
                var EdgeTags = grid.iGeomEdges.EdgeTags;
                var gridDataConcrete = grid as BoSSS.Foundation.Grid.Classic.GridData;
                periodicIdx = new List<int>();
                periodicMtx = new List<MultidimensionalArray>();
                for (int i = 0; i < Length; i++) {
                    var eT = EdgeTags[i + i0];
                    if (eT < BoSSS.Foundation.Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG)
                        continue;
                    var Trafo = gridDataConcrete != null
                        ? gridDataConcrete.Edges.GetPeriodicTrafo(i + i0, true)
                        : grid.Grid.PeriodicTrafo[eT - BoSSS.Foundation.Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG];
                    periodicIdx.Add(i);
                    periodicMtx.Add(Trafo.Matrix);
                }
            }

            /// <summary>
            /// Rotated normals for the OUT-cell equation on the periodic edges <paramref name="periodicIdx"/> (other edges unchanged).
            /// </summary>
            private MultidimensionalArray GetRotatedNormals(List<int> periodicIdx, List<MultidimensionalArray> periodicMtx, int NoOfNodes) {
                var RotatedNormals = NormalBuffer.CloneAs();
                int D = NormalBuffer.GetLength(2);
                for (int k = 0; k < periodicIdx.Count; k++) {
                    int i = periodicIdx[k];
                    var Mtx = periodicMtx[k];
                    for (int node = 0; node < NoOfNodes; node++) {
                        for (int d = 0; d < D; d++) {
                            double acc = 0.0;
                            for (int e = 0; e < D; e++)
                                acc += Mtx[d, e] * NormalBuffer[i, node, e];
                            RotatedNormals[i, node, d] = acc;
                        }
                    }
                }
                return RotatedNormals;
            }

            /// <summary>
            /// Runs <paramref name="evaluate"/> in OUT's frame, i.e. with rotated normals and the IN-side parameters rotated by M;
            /// the OUT-side parameters stay native.
            /// </summary>
            private void EvaluateInOutFrame(List<int> periodicIdx, List<MultidimensionalArray> periodicMtx, int NoOfNodes, int D, Action evaluate) {
                var RotatedNormals = GetRotatedNormals(periodicIdx, periodicMtx, NoOfNodes);
                var savedNormals = this.NormalBuffer;
                var originalParamsIn = RotateParameters(m_ParameterFieldsValues_IN, GetParameterVectorTuples(D), periodicIdx, periodicMtx, NoOfNodes);
                this.NormalBuffer = RotatedNormals;
                try {
                    evaluate();
                } finally {
                    this.NormalBuffer = savedNormals;
                    RestoreParameters(m_ParameterFieldsValues_IN, originalParamsIn);
                }
            }

            /// <summary>
            /// Rotationally periodic edges: re-computes the cr=1 (OUT-cell equation) blocks in OUT's frame (Step 1)
            /// and rotates the IN/OUT cross blocks of vector-valued domain variables (Step 2).
            /// </summary>
            private void ApplyPeriodicRotation(int i0, int Length, List<int> periodicIdx, List<MultidimensionalArray> periodicMtx, EdgeFormParams efp, List<MyChunk> Edges, bool MustLock, int NoOfNodes, int D) {
                int GAMMA = m_owner.GAMMA;

                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Step 1: re-evaluate all forms in OUT's frame
                // and copy the cr=1 (OUT-cell equation) blocks into the main buffers
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                {
                    var scratchUxV = CloneSumBufferStructure(m_UxVSumBuffer);
                    var scratchGradUxV = CloneSumBufferStructure(m_GradUxVSumBuffer);
                    var scratchUxGradV = CloneSumBufferStructure(m_UxGradVSumBuffer);
                    var scratchGradUxGradV = CloneSumBufferStructure(m_GradUxGradVSumBuffer);

                    EvaluateInOutFrame(periodicIdx, periodicMtx, NoOfNodes, D, delegate () {
                        EvalNSumForm(ref efp, i0, m_Edgeform_UxV, Edges, m_UxVComponentBuffer, scratchUxV, false, NoOfNodes, D, m_Edgeform_UxV_Watches, MustLock,
                            (E, mda) => E.InternalEdge_UxV(ref efp, mda),
                            (E, mda) => E.BoundaryEdge_UxV(ref efp, mda));
                        EvalNSumForm(ref efp, i0, m_Edgeform_GradUxV, Edges, m_GradUxVComponentBuffer, scratchGradUxV, false, NoOfNodes, D, m_Edgeform_GradUxV_Watches, MustLock,
                            (E, mda) => E.InternalEdge_GradUxV(ref efp, mda),
                            (E, mda) => E.BoundaryEdge_GradUxV(ref efp, mda));
                        EvalNSumForm(ref efp, i0, m_Edgeform_UxGradV, Edges, m_UxGradVComponentBuffer, scratchUxGradV, false, NoOfNodes, D, m_Edgeform_UxGradV_Watches, MustLock,
                            (E, mda) => E.InternalEdge_UxGradV(ref efp, mda),
                            (E, mda) => E.BoundaryEdge_UxGradV(ref efp, mda));
                        EvalNSumForm(ref efp, i0, m_Edgeform_GradUxGradV, Edges, m_GradUxGradVComponentBuffer, scratchGradUxGradV, false, NoOfNodes, D, m_Edgeform_GradUxGradV_Watches, MustLock,
                            (E, mda) => E.InternalEdge_GradUxGradV(ref efp, mda),
                            (E, mda) => E.BoundaryEdge_GradUxGradV(ref efp, mda));
                    });

                    SpliceCr1Block(scratchUxV, m_UxVSumBuffer, GAMMA, NoOfNodes, periodicIdx);
                    SpliceCr1Block(scratchGradUxV, m_GradUxVSumBuffer, GAMMA, NoOfNodes, periodicIdx);
                    SpliceCr1Block(scratchUxGradV, m_UxGradVSumBuffer, GAMMA, NoOfNodes, periodicIdx);
                    SpliceCr1Block(scratchGradUxGradV, m_GradUxGradVSumBuffer, GAMMA, NoOfNodes, periodicIdx);
                }

                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Step 2: rotate the cross blocks of vector-valued domain variables:
                // (cr=0,cc=1) by M, (cr=1,cc=0) by M^-1; for GradUxV and GradUxGradV,
                // also the GradU direction axis (see RotateGradientAxis)
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                int[][] VectorTuples = m_owner.Operator.VectorFieldIndices.ToArray();
                if (VectorTuples.Length <= 0)
                    return;

                for (int k = 0; k < periodicIdx.Count; k++) {
                    int i = periodicIdx[k];
                    var M = periodicMtx[k];
                    var Minv = M.CloneAs();
                    Minv.InvertInPlace();

                    RotateCrossBlock(m_UxVSumBuffer, m_Trf_UxVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, M, i, 0, 1);
                    RotateCrossBlock(m_UxVSumBuffer, m_Trf_UxVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, Minv, i, 1, 0);

                    RotateCrossBlock(m_GradUxVSumBuffer, m_Trf_GradUxVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, M, i, 0, 1);
                    RotateCrossBlock(m_GradUxVSumBuffer, m_Trf_GradUxVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, Minv, i, 1, 0);
                    // GradU direction axis (position 2)
                    RotateGradientAxis(m_GradUxVSumBuffer, GAMMA, m_owner.DELTA, NoOfNodes, 2, M, i, 0, 1);
                    RotateGradientAxis(m_GradUxVSumBuffer, GAMMA, m_owner.DELTA, NoOfNodes, 2, Minv, i, 1, 0);

                    RotateCrossBlock(m_UxGradVSumBuffer, m_Trf_UxGradVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, M, i, 0, 1);
                    RotateCrossBlock(m_UxGradVSumBuffer, m_Trf_UxGradVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, Minv, i, 1, 0);
                    // no RotateGradientAxis: the extra axis is the GradV direction

                    RotateCrossBlock(m_GradUxGradVSumBuffer, m_Trf_GradUxGradVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, M, i, 0, 1);
                    RotateCrossBlock(m_GradUxGradVSumBuffer, m_Trf_GradUxGradVSumBuffer, VectorTuples, GAMMA, Length, NoOfNodes, Minv, i, 1, 0);
                    // GradU direction axis (position 3; position 2 is GradV)
                    RotateGradientAxis(m_GradUxGradVSumBuffer, GAMMA, m_owner.DELTA, NoOfNodes, 3, M, i, 0, 1);
                    RotateGradientAxis(m_GradUxGradVSumBuffer, GAMMA, m_owner.DELTA, NoOfNodes, 3, Minv, i, 1, 0);
                }
            }

            /// <summary>
            /// Affine counterpart of Step 1 in <see cref="ApplyPeriodicRotation"/>: re-computes the cr=1 (OUT-cell equation) part
            /// of <see cref="m_VSumBuffer"/>/<see cref="m_GradVSumBuffer"/> with rotated normal and IN-side parameters rotated into OUT's frame.
            /// </summary>
            private void ApplyPeriodicRotationAffine(int i0, List<int> periodicIdx, List<MultidimensionalArray> periodicMtx, EdgeFormParams efp, List<MyChunk> Edges, bool MustLock, int NoOfNodes, int D) {
                int GAMMA = m_owner.GAMMA;

                var scratchV = CloneSumBufferStructure(m_VSumBuffer);
                var scratchGradV = CloneSumBufferStructure(m_GradVSumBuffer);

                EvaluateInOutFrame(periodicIdx, periodicMtx, NoOfNodes, D, delegate () {
                    EvalNSumForm(ref efp, i0, m_EdgeSourceV, Edges, m_VComponentBuffer, scratchV, true, NoOfNodes, D, m_EdgeSourceV_Watches, MustLock,
                        (E, mda) => E.InternalEdge_V(ref efp, mda),
                        (E, mda) => E.BoundaryEdge_V(ref efp, mda));
                    EvalNSumForm(ref efp, i0, m_EdgeSourceGradV, Edges, m_GradVComponentBuffer, scratchGradV, true, NoOfNodes, D, m_EdgeSourceGradV_Watches, MustLock,
                        (E, mda) => E.InternalEdge_GradV(ref efp, mda),
                        (E, mda) => E.BoundaryEdge_GradV(ref efp, mda));
                });

                SpliceCr1Block(scratchV, m_VSumBuffer, GAMMA, NoOfNodes, periodicIdx);
                SpliceCr1Block(scratchGradV, m_GradVSumBuffer, GAMMA, NoOfNodes, periodicIdx);
            }

            /// <summary>
            /// Drops blocks created by <see cref="RotateCrossBlock"/> for a chunk of different size.
            /// </summary>
            private void DropStaleRotationSlots(int Length, int NoOfNodes) {
                InvalidateStaleVectorTupleSlots(m_UxVSumBuffer, Length, NoOfNodes);
                InvalidateStaleVectorTupleSlots(m_GradUxVSumBuffer, Length, NoOfNodes);
                InvalidateStaleVectorTupleSlots(m_Trf_GradUxVSumBuffer, Length, NoOfNodes);
                InvalidateStaleVectorTupleSlots(m_UxGradVSumBuffer, Length, NoOfNodes);
                InvalidateStaleVectorTupleSlots(m_Trf_UxGradVSumBuffer, Length, NoOfNodes);
                InvalidateStaleVectorTupleSlots(m_GradUxGradVSumBuffer, Length, NoOfNodes);
                InvalidateStaleVectorTupleSlots(m_Trf_GradUxGradVSumBuffer, Length, NoOfNodes);
            }

            /// <summary>
            /// Rotationally periodic edges of one chunk.
            /// </summary>
            private sealed class PeriodicEdgeChunk {

                readonly ThreadLocals m_tl;
                readonly int m_i0, m_Length, m_NoOfNodes, m_D;
                readonly List<int> m_Idx;
                readonly List<MultidimensionalArray> m_Mtx;

                PeriodicEdgeChunk(ThreadLocals tl, int i0, int Length, int NoOfNodes, int D, List<int> Idx, List<MultidimensionalArray> Mtx) {
                    m_tl = tl;
                    m_i0 = i0;
                    m_Length = Length;
                    m_NoOfNodes = NoOfNodes;
                    m_D = D;
                    m_Idx = Idx;
                    m_Mtx = Mtx;
                }

                /// <summary>
                /// null, if the chunk contains no rotationally periodic edges.
                /// </summary>
                public static PeriodicEdgeChunk TryCreate(ThreadLocals tl, int i0, int Length, IGridData grid, int NoOfNodes, int D) {
                    tl.GetPeriodicEdgesAndMatrices(i0, Length, grid, out var Idx, out var Mtx);
                    if (Idx.Count <= 0)
                        return null;
                    return new PeriodicEdgeChunk(tl, i0, Length, NoOfNodes, D, Idx, Mtx);
                }

                /// <summary>
                /// Rotates the OUT-side parameters into IN's frame, until the returned object is disposed.
                /// </summary>
                public IDisposable RotateOutParametersIntoInFrame() {
                    var MtxInv = m_Mtx.Select(M => { var Minv = M.CloneAs(); Minv.InvertInPlace(); return Minv; }).ToList();
                    var paramValues = m_tl.m_ParameterFieldsValues_OUT;
                    var original = RotateParameters(paramValues, m_tl.GetParameterVectorTuples(m_D), m_Idx, MtxInv, m_NoOfNodes);
                    return new ParameterRestore(paramValues, original);
                }

                /// <summary>
                /// See <see cref="ApplyPeriodicRotation"/> and <see cref="ApplyPeriodicRotationAffine"/>.
                /// </summary>
                public void ApplyRotation(bool LinearRequired, bool AffineRequired, EdgeFormParams efp, List<MyChunk> Edges, bool MustLock) {
                    if (LinearRequired)
                        m_tl.ApplyPeriodicRotation(m_i0, m_Length, m_Idx, m_Mtx, efp, Edges, MustLock, m_NoOfNodes, m_D);
                    if (AffineRequired)
                        m_tl.ApplyPeriodicRotationAffine(m_i0, m_Idx, m_Mtx, efp, Edges, MustLock, m_NoOfNodes, m_D);
                }

                sealed class ParameterRestore : IDisposable {
                    readonly MultidimensionalArray[] m_paramValues, m_original;

                    public ParameterRestore(MultidimensionalArray[] paramValues, MultidimensionalArray[] original) {
                        m_paramValues = paramValues;
                        m_original = original;
                    }

                    public void Dispose() {
                        RestoreParameters(m_paramValues, m_original);
                    }
                }
            }
            #endregion

            private MultidimensionalArray ExtractTrafo(MultidimensionalArray trafo, int N) {
                if (trafo.GetLength(1) == N)
                    return trafo;
                else
                    return trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { trafo.GetLength(0) - 1, N - 1, N - 1 });
            }

            private void GetQRbufferLinear(int Length, MultidimensionalArray QuadResult, int cr, int cc,
                bool Affine,
                int MR, int NR,
                ref int iQlBuf, ref MultidimensionalArray Ql,
                int I0Row, int gamma, int I0Col, int delta,
                out MultidimensionalArray _R, out MultidimensionalArray _Q) {
                _R = QuadResult.ExtractSubArrayShallow(new int[] { 0, cr, cc, I0Row, I0Col }, new int[] { Length - 1, cr - 1, cc - 1, I0Row + m_owner.m_RowL[gamma] - 1, I0Col + m_owner.m_ColL[delta] - 1 });


                if (Affine) {
                    _Q = null;
                } else {
                    if (Ql == null)
                        Ql = TempBuffer.GetTempMultidimensionalarray(out iQlBuf, Length, MR, NR);

                    if (m_owner.m_RowL[gamma] != MR || m_owner.m_ColL[delta] != NR)
                        _Q = Ql.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, m_owner.m_RowL[gamma] - 1, m_owner.m_ColL[delta] - 1 });
                    else
                        _Q = Ql;

                }
            }

            private void GetRQbufferAffine(int Length, MultidimensionalArray QuadResult, int cr, bool Affine, int MR, ref int iQaBuf, ref MultidimensionalArray Qa, int I0Row, int gamma, int I0Col, out MultidimensionalArray _R, out MultidimensionalArray _Q) {
                if (Affine) {
                    _R = QuadResult.ExtractSubArrayShallow(new int[] { 0, cr, 0, I0Row, I0Col }, new int[] { Length - 1, cr - 1, -1, I0Row + m_owner.m_RowL[gamma] - 1, I0Col - 1 });
                    _Q = null;
                } else {
                    if (Qa == null)
                        Qa = TempBuffer.GetTempMultidimensionalarray(out iQaBuf, Length, MR);

                    if (m_owner.m_RowL[gamma] != MR)
                        _R = Qa.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Length - 1, m_owner.m_RowL[gamma] - 1 });
                    else
                        _R = Qa;

                    _Q = QuadResult.ExtractSubArrayShallow(new int[] { 0, cr, 0, I0Row, I0Col }, new int[] { Length - 1, cr - 1, -1, I0Row + m_owner.m_RowL[gamma] - 1, I0Col - 1 });
                }
            }

            private void GetZbuffer(int Length, int NoOfNodes, int D, int NZ, ref int iBufZ, ref MultidimensionalArray Z, out MultidimensionalArray _Z, int N) {
                Debug.Assert(N >= 0);
                Debug.Assert(NZ >= 0);
                Debug.Assert(N <= NZ);

                if (Z == null) {
                    if (D < 0) {
                        Z = TempBuffer.GetTempMultidimensionalarray(out iBufZ, Length, NoOfNodes, NZ);
                    } else {
                        Z = TempBuffer.GetTempMultidimensionalarray(out iBufZ, Length, NoOfNodes, NZ, D);
                    }
                }

                if (Z.GetLength(2) != N)
                    if (D < 0) {
                        Debug.Assert(Z.Dimension == 3);
                        _Z = Z.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, N - 1 });
                    } else {
                        Debug.Assert(Z.Dimension == 4);
                        _Z = Z.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, N - 1, D - 1 });
                    }
                else
                    _Z = Z;
            }


            void EvalNSumForm<EE>(ref EdgeFormParams efp, int i0,
                EquationComponentArgMapping<EE>[] Comps, List<MyChunk> subChunkx,
                MultidimensionalArray[][] CompBuffer, MultidimensionalArray[,,,] SumBuffer,
                bool Source, int NoOfNodes, int D,
                Stopwatch[][] watches, bool ParrentLocked,
                Action<EE, MultidimensionalArray> InnerEdgeForm,
                Action<EE, MultidimensionalArray> BoundarydgeForm)
                where EE : IEquationComponent {
                int DELTA = m_owner.DELTA;

                Debug.Assert(SumBuffer.GetLength(3) == (Source ? 1 : 2));

                //byte[] EdgeTags = m_owner.m_GridDat.iGeomEdges.EdgeTags;

                for (int __gamma = m_owner.GAMMA - 1; __gamma >= 0; __gamma--) {
                    int gamma = (__gamma + this.m_iThread) % m_owner.GAMMA; // run through loop differently in each thread to reduce locking

                    var ecq = Comps[gamma];

                    for (int delta = 0; delta < SumBuffer.GetLength(1); delta++) {
                        for (int cr = 0; cr < 2; cr++) {
                            for (int cc = 0; cc < SumBuffer.GetLength(3); cc++) {
                                if (SumBuffer[gamma, delta, cr, cc] != null)
                                    SumBuffer[gamma, delta, cr, cc].Clear();
                            }
                        }
                    }

                    Stopwatch[] watches_gamma = watches[gamma];

                    for (int _i = 0; _i < ecq.m_AllComponentsOfMyType.Length; _i++) {  // loop over equation components
                        int i = (_i + this.m_iThread) % ecq.m_AllComponentsOfMyType.Length; // run through loop differently in each thread to reduce locking


                        EE comp = ecq.m_AllComponentsOfMyType[i];
                        var bLck = ParrentLocked ? null : ecq.m_LockObjects[i];
                        int NoOfArgs = ecq.NoOfArguments[i];
                        var CompBuf_gamma_i = CompBuffer[gamma][i];
                        CompBuf_gamma_i.Clear();

                        int[] I0 = new int[CompBuf_gamma_i.Dimension];
                        var IE = CompBuf_gamma_i.Lengths;
                        for (int ll = 1; ll < IE.Length; ll++)
                            IE[ll]--;
                        Debug.Assert(NoOfArgs == comp.ArgumentOrdering.Count);
                        int NoOfParams = ecq.NoOfParameters[i];
                        Debug.Assert(NoOfParams == ((comp.ParameterOrdering != null) ? comp.ParameterOrdering.Count : 0));

                        // map parameters
                        efp.ParameterVars_IN = new MultidimensionalArray[NoOfParams];
                        efp.ParameterVars_OUT = new MultidimensionalArray[NoOfParams];

                        for (int ii = 0; ii < subChunkx.Count; ii++) { // loop over chunks of inner edges/boundary edges
                            var chk = subChunkx[ii];

                            efp.e0 = chk.i0;
                            efp.Len = chk.Len;
                            efp.Normals = this.NormalBuffer.ExtractSubArrayShallow(new int[] { chk.i0 - i0, 0, 0 }, new int[] { chk.i0 + chk.Len - 1 - i0, NoOfNodes - 1, D - 1 });
                            efp.Nodes = this.GlobalNodes.ExtractSubArrayShallow(new int[] { chk.i0 - i0, 0, 0 }, new int[] { chk.i0 + chk.Len - 1 - i0, NoOfNodes - 1, D - 1 });

                            I0[0] = chk.i0 - i0;
                            IE[0] = chk.Len + I0[0] - 1;
                            for (int c = 0; c < NoOfParams; c++) {
                                int targ = ecq.AllToSub[i, c + NoOfArgs] - DELTA;
                                Debug.Assert(targ >= 0);
                                efp.ParameterVars_IN[c] = m_ParameterFieldsValues_IN[targ].ExtractSubArrayShallow(new int[] { I0[0], I0[1] }, new int[] { IE[0], IE[1] });
                                efp.ParameterVars_OUT[c] = m_ParameterFieldsValues_OUT[targ].ExtractSubArrayShallow(new int[] { I0[0], I0[1] }, new int[] { IE[0], IE[1] });
                            }


                            if (chk.Type) {
                                // inner edge
                                //comp.InternalEdge(ref efp, CompBuf_gamma_i.ExtractSubArrayShallow(I0, IE));
                                watches_gamma[i].Start();
                                if(bLck != null) {
                                    lock (bLck) {
                                        InnerEdgeForm(comp, CompBuf_gamma_i.ExtractSubArrayShallow(I0, IE));
                                    }
                                } else {
                                    InnerEdgeForm(comp, CompBuf_gamma_i.ExtractSubArrayShallow(I0, IE));
                                }
                                watches_gamma[i].Stop();
                            } else {
                                // boundary edge

                                if (!Source) {
                                    Debug.Assert(IE[2] == 1);
                                    Debug.Assert(IE[3] == 1);
                                    IE[2] = -1;
                                    IE[3] = -1;
                                } else {
                                    Debug.Assert(IE[2] == 1);
                                    IE[2] = -1;
                                }

                                //comp.BoundaryEdge(ref efp, null, CompBuf_gamma_i.ExtractSubArrayShallow(I0,IE));
                                watches_gamma[i].Start();
                                if (bLck != null) {
                                    lock (bLck) {
                                        BoundarydgeForm(comp, CompBuf_gamma_i.ExtractSubArrayShallow(I0, IE));
                                    }
                                } else {
                                    BoundarydgeForm(comp, CompBuf_gamma_i.ExtractSubArrayShallow(I0, IE));
                                }
                                watches_gamma[i].Stop();

                                if (!Source) {
                                    IE[2] = 1;
                                    IE[3] = 1;
                                } else {
                                    IE[2] = 1;
                                }
                            }
                        }

                        // sum up
                        if (SumBuffer.GetLength(3) == 2) {
                            // branch for bilinear forms -> contribute to operator matrix
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            int[] Sel = new int[CompBuf_gamma_i.Dimension];
                            Sel.SetAll(-1);

                            for (int c = 0; c < NoOfArgs; c++) {
                                int targ = ecq.AllToSub[i, c];
                                for (int cr = 0; cr < 2; cr++) {
                                    for (int cc = 0; cc < 2; cc++) {
                                        Sel[2] = cr;
                                        Sel[3] = cc;
                                        Sel[4] = c;
                                        SumBuffer[gamma, targ, cr, cc].Acc(1.0, CompBuffer[gamma][i].ExtractSubArrayShallow(Sel));
                                    }
                                }
                            }
                        } else {
                            // branch for linear forms -> they contribute to affine vector
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            Debug.Assert(SumBuffer.GetLength(1) == 1);
                            Debug.Assert(SumBuffer.GetLength(3) == 1);
                            int[] Sel = new int[CompBuf_gamma_i.Dimension];
                            Sel.SetAll(-1);

                            for (int cr = 0; cr < 2; cr++) {
                                Sel[2] = cr;
                                SumBuffer[gamma, 0, cr, 0].Acc(1.0, CompBuffer[gamma][i].ExtractSubArrayShallow(Sel));
                            }
                        }
                    }
                }
            }

            void EdgeEvalSomething(int i0, int Length, int NoOfNodes, MultidimensionalArray[] ResultBuffer, Action<int, int, MultidimensionalArray> eva) {
                var Edges = m_owner.m_GridDat.iGeomEdges;
                int D = m_owner.m_GridDat.SpatialDimension;
                int NP = m_owner.m_ParameterFields.Length;
                int RD = ResultBuffer[0].Dimension;
                int[] I0 = new int[RD];
                int[] IE = new int[RD];
                for (int rd = 2; rd < RD; rd++) {
                    IE[rd] = ResultBuffer[0].GetLength(rd) - 1;
                }

                for (int i = 0; i < Length; i++) {
                    int jEdge = i + i0;
                    int jCell1, edge1, jCell2, edge2;
                    jCell1 = Edges.CellIndices[jEdge, 0];
                    edge1 = Edges.Edge2CellTrafoIndex[jEdge, 0];
                    jCell2 = Edges.CellIndices[jEdge, 1];
                    edge2 = Edges.Edge2CellTrafoIndex[jEdge, 1];

                    Debug.Assert(NormalBuffer.GetLength(0) == Length);
                    Debug.Assert(NormalBuffer.GetLength(1) == NoOfNodes);
                    Debug.Assert(NormalBuffer.GetLength(2) == D);

                    I0[0] = i;
                    IE[0] = i;
                    IE[1] = NoOfNodes - 1;

                    var bufIn = ResultBuffer[0].ExtractSubArrayShallow(I0, IE);
                    var bufOt = ResultBuffer[1].ExtractSubArrayShallow(I0, IE);
                    eva(jCell1, edge1, bufIn);

                    if (jCell2 >= 0) {
                        eva(jCell2, edge2, bufOt);
                    } else {
                        bufOt.Clear();
                    }
                }
            }

            void EvaluateParameters(int i0, int Length, NodeSet NS) {
                //var Edges = this.m_GridDat.Edges;
                //int D = this.m_GridDat.SpatialDimension;
                int NP = m_owner.m_ParameterFields.Length;

                for (int np = 0; np < NP; np++) {
                    EdgeEvalSomething(i0, Length, NS.NoOfNodes, new MultidimensionalArray[] { m_ParameterFieldsValues_IN[np], m_ParameterFieldsValues_OUT[np] },
                        delegate (int jCell, int nodeSetIdx, MultidimensionalArray buf) {
                            if (m_owner.m_ParameterFields[np] != null)
                                m_owner.m_ParameterFields[np].Evaluate(jCell, 1, NS.GetVolumeNodeSet(m_owner.m_GridDat, nodeSetIdx, true), buf);
                            else
                                buf.Clear();
                        });
                }
            }
        }



        void EvaluateEx(int i0, int Length, QuadRule qr, IIntegrationMetric metric, MultidimensionalArray QuadResult, int iThread, int NumThreads) {
            m_ThreadLocals[iThread].EvaluateEx(i0, Length, qr, metric, QuadResult, iThread, NumThreads);
        }




        /// <summary>
        /// Execution of quadrature
        /// </summary>
        /// <param name="domNrule"></param>
        /// <param name="RowMap"></param>
        /// <param name="ParamsMap"></param>
        /// <param name="ColMap">
        /// Domain variables resp. trial variables resp. matrix column variables
        /// </param>
        /// <param name="Matrix">
        /// output: matrix components will be accumulated here
        /// </param>
        /// <param name="AffineVector">
        /// output: affine components will be accumulated here
        /// </param>
        /// <param name="time"></param>
        public void Execute(ICompositeQuadRule<QuadRule> domNrule,
            UnsetteledCoordinateMapping RowMap, IList<DGField> ParamsMap, UnsetteledCoordinateMapping ColMap,
            M Matrix, V AffineVector, double time) {
            if (RowMap.BasisS.Count != GAMMA)
                throw new ArgumentException("Mismatch in number of codomain (rew. row-variables, resp. test-variables) variables.", "RowMap");
            if (ColMap.BasisS.Count != DELTA)
                throw new ArgumentException("Mismatch in number of domain (rew. column-variables, resp. trial-variables) variables.", "ColMap");
            var temp = new ThreadLocals(0, this, null);
            if (temp.IsNonEmpty == false)
                return;

            m_GridDat = RowMap.GridDat;

            if (!object.ReferenceEquals(m_GridDat, ColMap.GridDat))
                throw new ArgumentException();

            m_ParameterFields = ParamsMap == null ? new DGField[0] : ParamsMap.ToArray();

            m_RowL = RowMap.BasisS.Select(b => b.Length).ToArray();
            m_ColL = ColMap.BasisS.Select(b => b.Length).ToArray();
            m_TestFunctions = RowMap.BasisS.ToArray();
            m_TrialFunctions = ColMap.BasisS.ToArray();
            m_RowMap = RowMap;
            m_ColMap = ColMap;

            m_time = time;
            m_Matrix = Matrix;
            m_AffineVector = AffineVector;

            if(m_Matrix != null) {
                if (Matrix.RowPartitioning.LocalLength != RowMap.LocalLength)
                    throw new ArgumentException("Matrix mismatch: local number of rows", "Matrix");
                if (Matrix.ColPartitioning.LocalLength != ColMap.LocalLength)
                    throw new ArgumentException("Matrix mismatch: number of columns", "Matrix");
            }
            if(m_AffineVector != null) {
                if(m_AffineVector.Count != RowMap.LocalLength) {
                    throw new ArgumentException("affine vector mismatch: local number of rows", "AffineVector");
                }
            }

            var q = EdgeQuadrature.GetQuadrature2(new int[] { 2, 2, m_RowL.Sum(), (temp.LinearRequired ? m_ColL.Sum() : 0) + (temp.AffineRequired ? 1 : 0) }, // the additional column carries the affine part
                m_GridDat, domNrule,
                this.EvaluateEx,
                this.SaveIntegrationResults,
                _OnCloneForThreadParallelization: this.OnCloneForThreadParallelization,
                _AllocateBuffers: this.AllocateBuffers);

            q.ExecuteParallel = true;
            q.Execute();
        }

        double m_time;

        

        UnsetteledCoordinateMapping m_RowMap;
        UnsetteledCoordinateMapping m_ColMap;

        /// <summary>
        /// parameter fields
        /// </summary>
        DGField[] m_ParameterFields;

        /// <summary>
        /// length of codomain basis per cell
        /// </summary>
        int[] m_RowL;

        /// <summary>
        /// 
        /// </summary>
        Basis[] m_TestFunctions;

        /// <summary>
        /// length of domain basis per cell
        /// </summary>
        int[] m_ColL;

        /// <summary>
        /// 
        /// </summary>
        Basis[] m_TrialFunctions;



        IGridData m_GridDat;

        /// <summary>
        /// Matrix where the result is saved
        /// </summary>
        M m_Matrix;

        /// <summary>
        /// Affine vector where the result is saved
        /// </summary>
        V m_AffineVector;

       
    

        static MultidimensionalArray.MultiplyProgram mp_im_ik_Tikma = MultidimensionalArray.MultiplyProgram.Compile("im", "ika", "T(i)kma", true);
        static MultidimensionalArray.MultiplyProgram mp_im_ik_Tkm = MultidimensionalArray.MultiplyProgram.Compile("im", "ik", "T(i)km", true);
        
        static MultidimensionalArray.MultiplyProgram mp_imn_imn_Ti = MultidimensionalArray.MultiplyProgram.Compile("imn", "imn", "T(i)", true);
        static MultidimensionalArray.MultiplyProgram mp_imn_imn_i = MultidimensionalArray.MultiplyProgram.Compile("imn", "imn", "i");

        static MultidimensionalArray.MultiplyProgram mp_im_im_Ti = MultidimensionalArray.MultiplyProgram.Compile("im", "im", "T(i)", true);
        static MultidimensionalArray.MultiplyProgram mp_im_im_i = MultidimensionalArray.MultiplyProgram.Compile("im", "im", "i");

        static MultidimensionalArray.MultiplyProgram mp_ikn_ik_Tikn = MultidimensionalArray.MultiplyProgram.Compile("ikn", "ik", "T(i)kn", true);
        static MultidimensionalArray.MultiplyProgram mp_imn_Tikm_ikn = MultidimensionalArray.MultiplyProgram.Compile("imn", "T(i)km", "ikn", true);

        static MultidimensionalArray.MultiplyProgram mp_ikna_ikab_Tiknb = MultidimensionalArray.MultiplyProgram.Compile("ikna", "ikab", "T(i)knb", true);
        static MultidimensionalArray.MultiplyProgram mp_imn_Tikma_ikna = MultidimensionalArray.MultiplyProgram.Compile("imn", "T(i)kma", "ikna", true);

        static MultidimensionalArray.MultiplyProgram mp_ikm_ika_Tikma = MultidimensionalArray.MultiplyProgram.Compile("ikm", "ika", "T(i)kma", true);
        static MultidimensionalArray.MultiplyProgram mp_imn_Tikn_ikm = MultidimensionalArray.MultiplyProgram.Compile("imn", "T(i)kn", "ikm", true);

        static MultidimensionalArray.MultiplyProgram mp_ikn_ikb_Tiknb = MultidimensionalArray.MultiplyProgram.Compile("ikn", "ikb", "T(i)knb", true);
        //static MultidimensionalArray.MultiplyProgram mp_imn_Tikm_ikn = MultidimensionalArray.MultiplyProgram.Compile("imn", "T(i)km", "ikn", true);

        static MultidimensionalArray.MultiplyProgram mp_ikb_ikd_Tibd = MultidimensionalArray.MultiplyProgram.Compile("ikb", "ikd", "T(i)bd", true);
        static MultidimensionalArray.MultiplyProgram mp_ikb_ikd_ikbd = MultidimensionalArray.MultiplyProgram.Compile("ikb", "ikd", "ikbd");

        static MultidimensionalArray.MultiplyProgram mp_ika_Tiad_ikd = MultidimensionalArray.MultiplyProgram.Compile("ika", "T(i)ad", "ikd", true);
        static MultidimensionalArray.MultiplyProgram mp_ika_ikad_ikd = MultidimensionalArray.MultiplyProgram.Compile("ika", "ikad", "ikd");

        static MultidimensionalArray.MultiplyProgram mp_ikae_Tiad_ikde = MultidimensionalArray.MultiplyProgram.Compile("ikae", "T(i)ad", "ikde", true);
        static MultidimensionalArray.MultiplyProgram mp_ikae_ikad_ikde = MultidimensionalArray.MultiplyProgram.Compile("ikae", "ikad", "ikde");
        static MultidimensionalArray.MultiplyProgram mp_ikab_ikae_Tibe = MultidimensionalArray.MultiplyProgram.Compile("ikab", "ikae", "T(i)be", true);
        static MultidimensionalArray.MultiplyProgram mp_ikab_ikae_ikbe = MultidimensionalArray.MultiplyProgram.Compile("ikab", "ikae", "ikbe");
        
        static MultidimensionalArray.MultiplyProgram mp_kn_k_kn = MultidimensionalArray.MultiplyProgram.Compile("kn", "k", "kn");
        static MultidimensionalArray.MultiplyProgram mp_knd_k_knd = MultidimensionalArray.MultiplyProgram.Compile("knd", "k", "knd");

        static MultidimensionalArray.MultiplyProgram mp_ibn_iba_Tian = MultidimensionalArray.MultiplyProgram.Compile("ibn", "iba", "T(i)an", true);
        static MultidimensionalArray.MultiplyProgram mp_imn_Tibm_ibn = MultidimensionalArray.MultiplyProgram.Compile("imn", "T(i)bm", "ibn", true);

        static MultidimensionalArray.MultiplyProgram mp_jn_Tjmn_jm = MultidimensionalArray.MultiplyProgram.Compile("jn", "T(j)mn", "jm", true);


        // true if any of the (cr,cc) blocks is present; affine buffers only have the (0,0) entry
        static bool AnyBlockNonNull(MultidimensionalArray[, , ,] SumBuf, int i, int j) {
            if(SumBuf.GetLength(3) <= 1)
                return SumBuf[i, j, 0, 0] != null; // affine/source buffer
            return SumBuf[i, j, 0, 0] != null || SumBuf[i, j, 0, 1] != null || SumBuf[i, j, 1, 0] != null || SumBuf[i, j, 1, 1] != null;
        }

        void TestFunctionRequired(bool[] Req, MultidimensionalArray[, , ,] SumBuf, ref int maxDeg, ref int maxN) {
            Debug.Assert(Req.Length == SumBuf.GetLength(0));
            
                for(int i = SumBuf.GetLength(0) - 1; i >= 0; i--) {
                    for(int j = SumBuf.GetLength(1) - 1; j >= 0; j--) {
                        if(AnyBlockNonNull(SumBuf, i, j)) {
                            maxDeg = Math.Max(maxDeg, this.m_TestFunctions[i].Degree);
                            maxN = Math.Max(maxN, this.m_RowL[i]);
                            Req[i] = true;
                            break;
                        }
                    }
                }
            
        }

        void TrialFunctionRequired(bool[] Req, MultidimensionalArray[, , ,] SumBuf, bool bLinreq, ref int maxDeg, ref int maxN) {
            Debug.Assert(Req.Length == SumBuf.GetLength(1));
            if(bLinreq) {
                for(int i = SumBuf.GetLength(1) - 1; i >= 0; i--) {
                    for(int j = SumBuf.GetLength(0) - 1; j >= 0; j--) {
                        if(AnyBlockNonNull(SumBuf, j, i)) {
                            maxDeg = Math.Max(maxDeg, this.m_TrialFunctions[i].Degree);
                            maxN = Math.Max(maxN, this.m_ColL[i]);
                            Req[i] = true;
                            break;
                        }
                    }
                }
            } else {
                Req.SetAll(false);
            }
        }



        /// <summary>
        /// scaling factor for accumulation
        /// </summary>
        internal double m_alpha = 1.0;

        protected void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
            var Edge2Cell = this.m_GridDat.iGeomEdges.LogicalCellIndices;
            int M = m_RowMap.NoOfCoordinatesPerCell;
            int N = m_ColMap.NoOfCoordinatesPerCell;
            int Jup = this.m_GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            double a = m_alpha;

            bool bLinearRequired = m_ThreadLocals[0].LinearRequired;
            bool bAffineRequired = m_ThreadLocals[0].AffineRequired;

            int offset = bLinearRequired ? N : 0;

            for (int i = 0; i < Length; i++) {
                int jEdge = i + i0;
                int CC = Edge2Cell[jEdge, 1] >= 0 ? 2 : 1;

                for (int cr = 0; cr < CC; cr++) {
                    int jCell_cr = Edge2Cell[jEdge, cr];

                    // Matrix part
                    if (bLinearRequired && jCell_cr < Jup) {
                        for (int cc = 0; cc < CC; cc++) {
                            int jCell_cc = Edge2Cell[jEdge, cc];


                            long m0 = this.m_RowMap.GlobalUnique1stCoordinate(jCell_cr);
                            long n0 = this.m_ColMap.GlobalUnique1stCoordinate(jCell_cc);


                            Debug.Assert(ResultsOfIntegration.GetLength(3) == M);
                            Debug.Assert(ResultsOfIntegration.GetLength(4) == (bAffineRequired ? N + 1 : N));

                            var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(
                                new int[] { i, cr, cc, 0, 0 },
                                new int[] { i - 1, cr - 1, cc - 1, M - 1, N - 1 });

                            //for (int m = 0; m < M; m++) {
                            //    for (int n = 0; n < N; n++) {
                            //        //m_Matrix[m0 + m, n0 + n] += BlockRes[m, n];
                            //    }
                            //}
                            m_Matrix.AccBlock(m0, n0, a, BlockRes);
                        }
                    }

                    // Affine offset part
                    if (bAffineRequired && jCell_cr < Jup) {
                        int m0 = this.m_RowMap.LocalUnique1stCoordinate(jCell_cr);

                        
                        var BlockRes = ResultsOfIntegration.ExtractSubArrayShallow(
                            new int[] { i, cr, 0, 0, offset },
                            new int[] { i - 1, cr - 1, 0 - 1, M - 1, offset - 1 });

                        for (int m = 0; m < M; m++)
                            m_AffineVector[m0 + m] += BlockRes[m]*a;
                    }
                }
            }
        }
    }
}
