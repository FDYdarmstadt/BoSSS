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

using System.Collections.Generic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;
using BoSSS.Platform;
using System.Diagnostics;
using System;
using ilPSP.Utils;
using System.Collections;
using System.Linq;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using System.Threading;

namespace BoSSS.Foundation.Quadrature.NonLin {

   
    /// <summary>
    /// edge quadrature of nonlinear equation components
    /// </summary>
    internal class NECQuadratureEdge : NECQuadratureCommon {

        /// <summary>
        /// In subgrid mode, i.e. with restricted evaluation, this denotes
        /// which cells are within the subgrid, usually acquired by 
        /// <see cref="CellMask.GetBitMaskWithExternal()"/>
        /// </summary>
        public BitArray SubGridCellsMarker;

        /// <summary>
        /// See <see cref="SubGridBoundaryModes"/>
        /// </summary>
        public SubGridBoundaryModes SubGridBoundaryTreatment;

        public double[] m_outputBndEdge;



        /// <summary>
        /// ctor.
        /// </summary>
        public NECQuadratureEdge(IGridData context,
                                 DifferentialOperator DiffOp,
                                 IEnumerable<DGField> __DomainFields,
                                 IEnumerable<DGField> __ParameterFields,
                                 UnsetteledCoordinateMapping CodomainMapping,
                                 ICompositeQuadRule<QuadRule> domNrule) :
            base(context, DiffOp, __DomainFields, __ParameterFields, CodomainMapping) {


            // ----
            // init 
            // ----
            var temp = new ThreadLocalsEdg(0, this, null);
            base.IsNonEmpty = temp.IsNonEmpty;

            // determine for which fields evaluation is required.
            m_GradientRequired = new bool[base.m_DomainFields.Length];
            m_ValueRequired = new bool[base.m_DomainAndParamFields.Length]; // note that m_DomainFields may also contain concatenated parameters.
            m_MeanValueRequired = new bool[base.m_DomainFields.Length];

            temp.m_EdgeForm_V.DetermineReqFields(m_GradientRequired,
                comp => (((comp.BoundaryEdgeTerms | comp.InnerEdgeTerms) & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0));
            temp.m_EdgeForm_GradV.DetermineReqFields(m_GradientRequired,
                comp => (((comp.BoundaryEdgeTerms | comp.InnerEdgeTerms) & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0));
            temp.m_EdgeForm_V.DetermineReqFields(m_ValueRequired,
                comp => (((comp.BoundaryEdgeTerms | comp.InnerEdgeTerms) & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0));
            temp.m_EdgeForm_GradV.DetermineReqFields(m_ValueRequired,
                comp => (((comp.BoundaryEdgeTerms | comp.InnerEdgeTerms) & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0));

            temp.m_NonlinFluxes.DetermineReqFields(m_ValueRequired, comp => true);
            temp.m_NonlinFluxesEx.DetermineReqFields(m_ValueRequired, comp => true);
            temp.m_NonlinFluxesEx.DetermineReqFields(m_MeanValueRequired, comp => true);

            Debug.Assert(base.m_DomainFields.Length <= m_DomainAndParamFields.Length); // note that 'm_DomainFields' may contain concatenated parameter fields
            for (int i = base.m_DomainFields.Length; i < m_DomainAndParamFields.Length; i++) {
                m_ValueRequired[i] = true; // parameters are always required!
            }

            VectorComponentIndices = PeriodicBoundaryUtils.GetVectorFieldIndices(Operator.DomainVar.Cat(Operator.ParameterVar), GridDat.SpatialDimension).ToArray();

            // --------------------------------
            // basis
            // --------------------------------

            int Gamma = base.m_DomainFields.Length;

            for (int i = 0; i < m_CodomainBasisS.Length; i++) {
                if (temp.m_EdgeForm_V[i].m_AllComponentsOfMyType.Length > 0
                    || temp.m_NonlinFluxes[i].m_AllComponentsOfMyType.Length > 0
                    || temp.m_NonlinFluxesEx[i].m_AllComponentsOfMyType.Length > 0) {
                    if (maxTestBasis == null || maxTestBasis.Degree < m_CodomainBasisS[i].Degree) {
                        maxTestBasis = m_CodomainBasisS[i];
                    }
                }
                if (temp.m_EdgeForm_GradV[i].m_AllComponentsOfMyType.Length > 0) {
                    if (maxTestGradientBasis == null || maxTestGradientBasis.Degree < m_CodomainBasisS[i].Degree) {
                        maxTestGradientBasis = m_CodomainBasisS[i];
                    }
                }
            }


            // -----------------
            // quadrature object
            // -----------------

            var _Quad = EdgeQuadrature.GetQuadrature2(new int[] { CodomainMapping.BasisS.Sum(x => x.Length), 2 }, context, domNrule,
                this.EvaluateEx,
                this.SaveIntegrationResults,
                this.AllocateBuffers,
                _OnCloneForThreadParallelization:OnCloneForThreadParallelization);
            _Quad.ExecuteParallel = true;
            m_Quad = _Quad;
        }

      



        /// <summary>
        /// - Describes which tuples of the domain and parameter variables form vector fields
        /// - Required mainly for non-parallel periodic edges;
        /// </summary>
        int[][] VectorComponentIndices;

        readonly Basis maxTestBasis = null;
        readonly Basis maxTestGradientBasis = null;

        /// <summary>
        /// true, if the evaluation of the a domain variable is required. 
        /// index: correlates with domain variables
        /// </summary>
        readonly bool[] m_ValueRequired;

        /// <summary>
        /// true, if the evaluation of the a domain variable mean value is required. 
        /// index: correlates with domain variables
        /// </summary>
        readonly bool[] m_MeanValueRequired;

        /// <summary>
        /// true, if the evaluation of the gradient of a domain variable is required. 
        /// index: correlates with domain variables
        /// </summary>
        readonly bool[] m_GradientRequired;


        class ThreadLocalsEdg : ThreadLocals {

            readonly new NECQuadratureEdge m_owner;

            public ThreadLocalsEdg(int iThread, NECQuadratureEdge _owner, IQuadrature q) : base(iThread, _owner, q) {
                m_owner = _owner;


                // ------------------------
                // sort equation components
                // ------------------------

                m_EdgeForm_V = EquationComponentArgMapping<INonlinEdgeForm_V>.GetArgMapping(m_owner.Operator, true,
                    comp => (((comp.BoundaryEdgeTerms | comp.InnerEdgeTerms) & (TermActivationFlags.V | TermActivationFlags.UxV | TermActivationFlags.GradUxV)) != 0),
                    eq => (eq is IEdgeForm ? new NonlinEdgeFormVectorizer((IEdgeForm)eq) : null));

                m_EdgeForm_GradV = EquationComponentArgMapping<INonlinEdgeForm_GradV>.GetArgMapping(m_owner.Operator, true,
                    comp => (((comp.BoundaryEdgeTerms | comp.InnerEdgeTerms) & (TermActivationFlags.GradV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxGradV)) != 0),
                    eq => (eq is IEdgeForm ? new NonlinEdgeFormVectorizer((IEdgeForm)eq) : null));


                if (q == null)
                    return;

                // ---------
                // profiling
                // ---------

                var _CustomTimers = new Stopwatch[] { new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch() };
                var _CustomTimers_Names = new string[] { "Flux-Eval", "Basis-Eval", "Field-Eval", "Loops", "ParametersAndNormals", "Flux-Trafo" };
                Flux_Eval = _CustomTimers[0];
                Field_Eval = _CustomTimers[2];
                Basis_Eval = _CustomTimers[1];
                Loops = _CustomTimers[3];
                Flux_Trafo = _CustomTimers[5];
                ParametersAndNormals = _CustomTimers[4];
                q.CustomTimers = _CustomTimers;
                q.CustomTimers_Names = _CustomTimers_Names;
                q.CustomTimers_RootPointer = new int[_CustomTimers_Names.Length];
                ArrayTools.SetAll(q.CustomTimers_RootPointer, -1);

                this.m_EdgeForm_V_Watches = this.m_EdgeForm_V.InitStopWatches(0, q);
                this.m_EdgeForm_GradV_Watches = this.m_EdgeForm_GradV.InitStopWatches(0, q);
                base.m_NonlinFluxesWatches = base.m_NonlinFluxes.InitStopWatches(0, q);
                base.m_NonlinFluxesExWatches = base.m_NonlinFluxesEx.InitStopWatches(0, q);

                int Gamma = m_owner.m_CodomainBasisS.Length;

                for (int i = 0; i < Gamma; i++) {
                    var b = m_owner.m_CodomainBasisS[i];

                    if (this.m_EdgeForm_V[i].m_AllComponentsOfMyType.Length > 0) {
                        if (m_owner.maxTestBasis == null)
                            throw new NotSupportedException();
                        if (!b.IsSubBasis(m_owner.maxTestBasis))
                            throw new NotSupportedException();
                    }



                    if (this.m_EdgeForm_GradV[i].m_AllComponentsOfMyType.Length > 0) {
                        if (m_owner.maxTestGradientBasis == null)
                            throw new NotSupportedException();
                        if (!b.IsSubBasis(m_owner.maxTestGradientBasis))
                            throw new NotSupportedException();
                        //foreach (var b in base.m_CodomainBasisS) {
                        //    if (!b.IsSubBasis(maxTestGradientBasis))
                        //        throw new NotSupportedException();
                        //}

                    }
                }


                // ---------------------
                // alloc multidim arrays
                // ---------------------

                m_FluxValuesIN = new MultidimensionalArray[Gamma];
                m_FluxValuesOT = new MultidimensionalArray[Gamma];
                m_GradientFluxValuesIN = new MultidimensionalArray[Gamma];
                m_GradientFluxValuesOT = new MultidimensionalArray[Gamma];
                m_GradientFluxValuesINtrf = new MultidimensionalArray[Gamma];
                m_GradientFluxValuesOTtrf = new MultidimensionalArray[Gamma];
                for (int i = 0; i < Gamma; i++) {
                    if (this.m_EdgeForm_V[i].m_AllComponentsOfMyType.Length > 0
                        || this.m_NonlinFluxes[i].m_AllComponentsOfMyType.Length > 0
                        || this.m_NonlinFluxesEx[i].m_AllComponentsOfMyType.Length > 0) {
                        m_FluxValuesIN[i] = new MultidimensionalArray(2);
                        m_FluxValuesOT[i] = new MultidimensionalArray(2);
                    }

                    if (this.m_EdgeForm_GradV[i].m_AllComponentsOfMyType.Length > 0) {
                        m_GradientFluxValuesIN[i] = new MultidimensionalArray(3);
                        m_GradientFluxValuesOT[i] = new MultidimensionalArray(3);
                        m_GradientFluxValuesINtrf[i] = new MultidimensionalArray(3);
                        m_GradientFluxValuesOTtrf[i] = new MultidimensionalArray(3);
                    }
                }
                m_QuadResultIN = new MultidimensionalArray[Gamma];
                m_QuadResultOT = new MultidimensionalArray[Gamma];

                // --------------------------------------------
                // storage for field values at quadrature nodes
                // --------------------------------------------

                m_FieldValuesIN = new MultidimensionalArray[m_owner.m_DomainAndParamFields.Length];
                m_FieldValuesOT = new MultidimensionalArray[m_owner.m_DomainAndParamFields.Length];
                m_MeanFieldValuesIN = new MultidimensionalArray[m_owner.m_DomainFields.Length];
                m_MeanFieldValuesOT = new MultidimensionalArray[m_owner.m_DomainFields.Length];
                m_FieldGradientIN = new MultidimensionalArray[m_owner.m_DomainFields.Length];
                m_FieldGradientOT = new MultidimensionalArray[m_owner.m_DomainFields.Length];

                for (int i = 0; i < m_owner.m_DomainAndParamFields.Length; i++) {
                    if (m_owner.m_ValueRequired[i]) {
                        m_FieldValuesIN[i] = new MultidimensionalArray(2);
                        m_FieldValuesOT[i] = new MultidimensionalArray(2);
                    }

                    if (i < m_owner.m_MeanValueRequired.Length && m_owner.m_MeanValueRequired[i]) {
                        m_MeanFieldValuesIN[i] = new MultidimensionalArray(1);
                        m_MeanFieldValuesOT[i] = new MultidimensionalArray(1);
                    }

                    if (i < m_owner.m_GradientRequired.Length && m_owner.m_GradientRequired[i]) {
                        m_FieldGradientIN[i] = new MultidimensionalArray(3);
                        m_FieldGradientOT[i] = new MultidimensionalArray(3);
                    }
                }

                

            }


            Stopwatch Flux_Eval;
            Stopwatch Flux_Trafo;
            Stopwatch Field_Eval;
            Stopwatch Basis_Eval;
            Stopwatch Loops;
            Stopwatch ParametersAndNormals;
            Stopwatch[][] m_EdgeForm_V_Watches;
            Stopwatch[][] m_EdgeForm_GradV_Watches;

            /// <summary>
            /// array index: codomain variable
            /// </summary>
            internal EquationComponentArgMapping<INonlinEdgeForm_GradV>[] m_EdgeForm_GradV;

            /// <summary>
            /// array index: codomain variable
            /// </summary>
            internal EquationComponentArgMapping<INonlinEdgeForm_V>[] m_EdgeForm_V;

            /// <summary>
            /// true, if this integrator is responsible for any component
            /// </summary>
            override public bool IsNonEmpty {
                get {
                    return
                        base.IsNonEmpty ||
                        m_EdgeForm_GradV.IsNonEmpty() ||
                        m_EdgeForm_V.IsNonEmpty();
                }
            }


            /// <summary>
            /// storage for field values at quadrature nodes in the first (in) cell which bounds to the edge;
            /// </summary>
            /// <remarks>
            /// array index: field index, where field order is defined by <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>;
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index: local edge index, with some offset;</item>
            ///   <item> 2nd index: node index;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FieldValuesIN;

            /// <summary>
            /// storage for field values at quadrature nodes in the second (out) cell which bounds to the edge;
            /// </summary>
            /// <remarks>
            /// array index: field index, where field order is defined by <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>;
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index: local edge index, with some offset;</item>
            ///   <item> 2nd index: node index;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FieldValuesOT;

            /// <summary>
            /// storage for field gradients at quadrature nodes in the first (in) cell which bounds to the edge;
            /// </summary>
            /// <remarks>
            /// array index: field index, where field order is defined by <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>;
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index: local edge index, with some offset;</item>
            ///   <item>2nd index: node index;</item>
            ///   <item>3rd index: spatial direction of derivative</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FieldGradientIN;

            /// <summary>
            /// storage for field gradients at quadrature nodes in the second (out) cell which bounds to the edge;
            /// </summary>
            /// <remarks>
            /// array index: field index, where field order is defined by <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>;
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index: local edge index, with some offset;</item>
            ///   <item> 2nd index: node index;</item>
            ///   <item>3rd index: spatial direction of derivative</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FieldGradientOT;

            /// <summary>
            /// storage for mean field values at quadrature nodes in the second cell which bounds to the edge;
            /// array index: field index, where field order is defined by <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>
            /// if one array entry is null, no mean value evaluation of the corresponding field is needed,
            /// because ther is no <see cref="INonlinearFluxEx"/>-object which requires it.
            /// 1st index: local edge index, with some offset;
            /// of course, the mean value is equal on all quadrature nodes on an edge, so there is no index 
            /// for the quadrature node;
            /// </summary>
            MultidimensionalArray[] m_MeanFieldValuesIN;

            /// <summary>
            /// storage for mean field values at quadrature nodes in the second cell which bounds to the edge;
            /// array index: field index, where field order is defined by <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>
            /// if one array entry is null, no mean value evaluation of the corresponding field is needed,
            /// because there is no <see cref="INonlinearFluxEx"/>-object which requires it.
            /// 1st index: local edge index, with some offset;
            /// of course, the mean value is equal on all quadrature nodes on an edge, so there is no index 
            /// for the quadrature node;
            /// </summary>
            MultidimensionalArray[] m_MeanFieldValuesOT;

            /// <summary>
            /// values of Riemann fluxes at quadrature nodes, for the IN-cell (in order to support arbitrary forms, IN and OUT flux may differ)
            /// </summary>
            /// <remarks>
            /// Array index: codomain variable index: (see 
            /// <see cref="DifferentialOperator.CodomainVar"/>,
            /// <see cref="NECQuadratureCommon.Operator"/>);
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index/array index: codomain variable index</item>
            ///   <item>2nd index: local edge index, with some offset;</item>
            ///   <item>3rd index: node index;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FluxValuesIN;

            /// <summary>
            /// values of Riemann fluxes at quadrature nodes, for the IN-cell (in order to support arbitrary forms, IN and OUT flux may differ)
            /// </summary>
            /// <remarks>
            /// Array index: codomain variable index: (see 
            /// <see cref="DifferentialOperator.CodomainVar"/>,
            /// <see cref="NECQuadratureCommon.Operator"/>);
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index/array index: codomain variable index</item>
            ///   <item>2nd index: local edge index, with some offset;</item>
            ///   <item>3rd index: node index;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FluxValuesOT;


            /// <summary>
            /// values of 'Gradient fluxes' at quadrature nodes, for the IN-cell (in order to support arbitrary forms, IN and OUT flux may differ)
            /// </summary>
            /// <remarks>
            /// Array index: codomain variable index: (see 
            /// <see cref="DifferentialOperator.CodomainVar"/>,
            /// <see cref="NECQuadratureCommon.Operator"/>);
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index/array index: codomain variable index</item>
            ///   <item>2nd index: local edge index, with some offset;</item>
            ///   <item>3rd index: node index;</item>
            ///   <item>4th index: spatial direction;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_GradientFluxValuesIN;

            /// <summary>
            /// values of 'Gradient fluxes' at quadrature nodes, for the IN-cell (in order to support arbitary forms, IN and OUT flux may differ)
            /// </summary>
            /// <remarks>
            /// Array index: codomain variable index: (see 
            /// <see cref="DifferentialOperator.CodomainVar"/>,
            /// <see cref="NECQuadratureCommon.Operator"/>);
            /// <br/>
            /// For each <see cref="MultidimensionalArray"/>:
            /// <list type="bullet">
            ///   <item>1st index/array index: codomain variable index</item>
            ///   <item>2nd index: local edge index, with some offset;</item>
            ///   <item>3rd index: node index;</item>
            ///   <item>4th index: spatial direction;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_GradientFluxValuesOT;

            MultidimensionalArray[] m_GradientFluxValuesINtrf;
            MultidimensionalArray[] m_GradientFluxValuesOTtrf;

            MultidimensionalArray[] m_QuadResultIN;
            MultidimensionalArray[] m_QuadResultOT;

            int[] m_InnerEdgesI0s;
            int[] m_InnerEdgesIEs;


       

            public void AllocateBuffers(int NoOfItems, MultidimensionalArray rule, int iThread, int NumThreads) {

                int NoOfNodes = rule.GetLength(0);
                int D = this.m_owner.GridDat.SpatialDimension;

                // ----------------------
                // array for field values
                // ----------------------

                for (int f = 0; f < m_FieldValuesIN.Length; f++) {
                    if (m_FieldValuesIN[f] != null) {
                        m_FieldValuesIN[f].Allocate(NoOfItems, NoOfNodes);
                        m_FieldValuesOT[f].Allocate(NoOfItems, NoOfNodes);
                    }
                }


                // ---------------------------
                // array for mean field values
                // ---------------------------


                for (int f = 0; f < m_MeanFieldValuesIN.Length; f++) {
                    if (m_MeanFieldValuesIN[f] != null) {
                        m_MeanFieldValuesIN[f].Allocate(NoOfItems);
                        m_MeanFieldValuesOT[f].Allocate(NoOfItems);
                    }
                }

                // -------------------------------
                // array for field gradient values
                // -------------------------------


                for (int f = 0; f < m_FieldGradientIN.Length; f++) {
                    if (m_FieldGradientIN[f] != null) {
                        m_FieldGradientIN[f].Allocate(NoOfItems, NoOfNodes, D);
                        m_FieldGradientOT[f].Allocate(NoOfItems, NoOfNodes, D);
                    }
                }



                // ---------------------
                // array for flux values
                // ---------------------

                for (int i = 0; i < m_owner.m_CodomainMapping.BasisS.Count; i++) {
                    if (m_FluxValuesIN[i] != null) {
                        m_FluxValuesIN[i].Allocate(NoOfItems, NoOfNodes);
                        m_FluxValuesOT[i].Allocate(NoOfItems, NoOfNodes);
                    }
                }

                for (int i = 0; i < m_owner.m_CodomainMapping.BasisS.Count; i++) {
                    if (m_GradientFluxValuesIN[i] != null) {
                        m_GradientFluxValuesIN[i].Allocate(NoOfItems, NoOfNodes, D);
                        m_GradientFluxValuesOT[i].Allocate(NoOfItems, NoOfNodes, D);
                        m_GradientFluxValuesINtrf[i].Allocate(NoOfItems, NoOfNodes, D);
                        m_GradientFluxValuesOTtrf[i].Allocate(NoOfItems, NoOfNodes, D);
                    }
                }
            }

            /// <summary>
            /// Integrand evaluation.
            /// </summary>
            public void EvaluateEx(int i0, int Length, QuadRule QR, MultidimensionalArray QuadResult, int iFred, int NoOfFreds) {

                NodeSet qrNodes = QR.Nodes;
                IGridData grid = m_owner.GridDat;
                int D = grid.SpatialDimension;
                int NoOfFields = m_owner.m_DomainAndParamFields.Length;
                int NoOfEquations = m_owner.m_CodomainBasisS.Length;
                int NoOfNodes = qrNodes.NoOfNodes;
                bool affine = grid.iGeomEdges.IsEdgeAffineLinear(i0);
                int iKrefEdge = grid.iGeomEdges.GetRefElementIndex(i0);


#if DEBUG
            for(int i = 1; i < Length; i++) {
                int iedge = i + i0;
                Debug.Assert(iKrefEdge == grid.iGeomEdges.GetRefElementIndex(iedge));
                Debug.Assert(affine == grid.iGeomEdges.IsEdgeAffineLinear(iedge));
            }
#endif

                if (!affine) {
                    for (int gamma = 0; gamma < NoOfEquations; gamma++) {
                        if (m_QuadResultIN[gamma] == null)
                            m_QuadResultIN[gamma] = new MultidimensionalArray(2);
                        if (m_QuadResultOT[gamma] == null)
                            m_QuadResultOT[gamma] = new MultidimensionalArray(2);

                        int Ne = m_owner.m_CodomainBasisS[gamma].Length;

                        m_QuadResultIN[gamma].Allocate(Length, Ne);
                        m_QuadResultOT[gamma].Allocate(Length, Ne);
                    }
                }
                // =============================================================
                // Evaluate basis functions and their gradients (test functions)
                // =============================================================
                #region EVAL_BASIS
                this.Basis_Eval.Start();

                MultidimensionalArray TstFuncXwgt = null;         // value of test functions times quadrature weight
                MultidimensionalArray TstFuncGradXwgt = null; // value of test function gradient times quadrature weight
                int jCellMin = int.MaxValue, jCellMax = int.MinValue;
                int NoOfSec = -1;
                {
                    // array for test functions
                    bool[] marker_TstFuncXwgt = null;
                    MultidimensionalArray TstFunc = null;         // value of test functions
                    MultidimensionalArray TstFuncGrad = null; // value of test function gradient

                    int Nmax = m_owner.maxTestBasis.Length;


                    // evaluate test functions and test function gradients
                    // ----------------------------------------------------
                    if (m_owner.maxTestBasis != null) {
                        TstFunc = grid.ChefBasis.EdgeEval.GetValues(qrNodes, i0, Length, m_owner.maxTestBasis.Degree);
                        Debug.Assert(TstFunc.Dimension == 3);
                        Debug.Assert(TstFunc.GetLength(2) >= Nmax);
                        if (TstFunc.GetLength(2) > Nmax) {
                            int[] I0 = new int[] { 0, 0, 0 };
                            int[] IE = new int[] { TstFunc.GetLength(0) - 1, qrNodes.NoOfNodes - 1, Nmax - 1 };
                            TstFunc = TstFunc.ExtractSubArrayShallow(I0, IE);
                        }


                        TstFuncXwgt = MultidimensionalArray.Create(TstFunc.Lengths);
                        marker_TstFuncXwgt = new bool[TstFunc.GetLength(0)];
                    }

                    if (m_owner.maxTestGradientBasis != null) {
                        TstFuncGrad = grid.ChefBasis.EdgeGradientEval.GetValues(qrNodes, i0, Length, m_owner.maxTestBasis.Degree);
                        Debug.Assert(TstFuncGrad.Dimension == 4);
                        Debug.Assert(TstFuncGrad.GetLength(2) >= Nmax);
                        if (TstFuncGrad.GetLength(2) > Nmax) {
                            int[] I0 = new int[] { 0, 0, 0, 0 };
                            int[] IE = new int[] { TstFuncGrad.GetLength(0) - 1, qrNodes.NoOfNodes - 1, Nmax - 1, D - 1 };
                            TstFuncGrad = TstFuncGrad.ExtractSubArrayShallow(I0, IE);
                        }


                        TstFuncGradXwgt = MultidimensionalArray.Create(TstFuncGrad.Lengths);
                        marker_TstFuncXwgt = new bool[TstFunc.GetLength(0)];
                    }

                    // sweep over all edges to process in order to check, ...
                    // ------------------------------------------------------

                    if (m_InnerEdgesI0s == null) {
                        m_InnerEdgesI0s = new int[Length];
                        m_InnerEdgesIEs = new int[Length];
                    }

                    if (m_InnerEdgesI0s.Length < Length) {
                        m_InnerEdgesI0s = new int[Length];
                        m_InnerEdgesIEs = new int[Length];
                    }

                    bool InInnerEsec = false;

                    int[,] E2C = grid.iGeomEdges.CellIndices;
                    int[,] TrIdx = grid.iGeomEdges.Edge2CellTrafoIndex;
                    for (int jEdge = 0; jEdge < Length; jEdge++) {
                        int jCell1, jCell2, edge1, edge2;
                        jCell1 = E2C[jEdge + i0, 0];
                        jCell2 = E2C[jEdge + i0, 1];
                        edge1 = TrIdx[jEdge + i0, 0];
                        edge2 = TrIdx[jEdge + i0, 1];
                        Debug.Assert(m_owner.SubGridCellsMarker == null || (m_owner.SubGridCellsMarker[jCell1] || (jCell2 >= 0 && m_owner.SubGridCellsMarker[jCell2])));
                        bool SubGridBnd = m_owner.SubGridCellsMarker != null && (jCell2 >= 0) && (m_owner.SubGridCellsMarker[jCell1] != m_owner.SubGridCellsMarker[jCell2]);

                        // find minimal and maximal index of invoked cells
                        {
                            jCellMin = Math.Min(jCellMin, jCell1);
                            jCellMax = Math.Max(jCellMax, jCell1);
                        }
                        if (jCell2 >= 0) {
                            jCellMin = Math.Min(jCellMin, jCell2);
                            jCellMax = Math.Max(jCellMax, jCell2);
                        }

                        // multiply basis functions with quadrature weights
                        if (marker_TstFuncXwgt != null) {
                            if (!marker_TstFuncXwgt[edge1]) {
                                if (TstFuncXwgt != null) {
                                    int[] I0 = new int[] { edge1, 0, 0 };
                                    int[] IE = new int[] { edge1 - 1, NoOfNodes - 1, TstFunc.GetLength(2) - 1 };
                                    TstFuncXwgt.ExtractSubArrayShallow(I0, IE).Multiply(1.0, QR.Weights, TstFunc.ExtractSubArrayShallow(I0, IE), 0.0, ref mp_kn_k_kn);
                                }

                                if (TstFuncGradXwgt != null) {
                                    int[] I0 = new int[] { edge1, 0, 0, 0 };
                                    int[] IE = new int[] { edge1 - 1, NoOfNodes - 1, TstFunc.GetLength(2) - 1, D - 1 };
                                    TstFuncGradXwgt.ExtractSubArrayShallow(I0, IE).Multiply(1.0, QR.Weights, TstFuncGrad.ExtractSubArrayShallow(I0, IE), 0.0, ref mp_knd_k_knd);
                                    marker_TstFuncXwgt[edge1] = true;


                                }

                                marker_TstFuncXwgt[edge1] = true;
                            }
                            if (jCell2 >= 0 && !marker_TstFuncXwgt[edge2]) {
                                if (TstFuncXwgt != null) {
                                    int[] I0 = new int[] { edge2, 0, 0 };
                                    int[] IE = new int[] { edge2 - 1, NoOfNodes - 1, TstFunc.GetLength(2) - 1 };
                                    TstFuncXwgt.ExtractSubArrayShallow(I0, IE).Multiply(1.0, QR.Weights, TstFunc.ExtractSubArrayShallow(I0, IE), 0.0, ref mp_kn_k_kn);
                                }

                                if (TstFuncGradXwgt != null) {
                                    int[] I0 = new int[] { edge2, 0, 0, 0 };
                                    int[] IE = new int[] { edge2 - 1, NoOfNodes - 1, TstFunc.GetLength(2) - 1, D - 1 };
                                    TstFuncGradXwgt.ExtractSubArrayShallow(I0, IE).Multiply(1.0, QR.Weights, TstFuncGrad.ExtractSubArrayShallow(I0, IE), 0.0, ref mp_knd_k_knd);
                                }

                                marker_TstFuncXwgt[edge2] = true;
                            }

                        }

                        // sections to switch between boundary and inner edge fluxes
                        if (jCell2 >= 0 && !SubGridBnd) {
                            if (InInnerEsec) {
                                m_InnerEdgesIEs[NoOfSec]++;
                            } else {
                                InInnerEsec = true;
                                NoOfSec++;
                                m_InnerEdgesI0s[NoOfSec] = jEdge;
                                m_InnerEdgesIEs[NoOfSec] = jEdge;
                            }
                        } else {

                            if (InInnerEsec) {
                                InInnerEsec = false;
                            }
                        }

                        //// evaluate all fields ...  
                        //for(int f = 0; f < NoOfFields; f++) {
                        //    if(m_DomainFields[f] != null) {

                        //        if(f < m_MeanFieldValuesIN.Length && m_MeanFieldValuesIN[f] != null) {
                        //            Debug.Assert((m_MeanFieldValuesIN[f] == null) || (m_MeanFieldValuesOT[f] != null));

                        //            m_DomainFields[f].EvaluateMean(jCell1, 1, m_MeanFieldValuesIN[f], jEdge, 0.0);

                        //            if(jCell2 >= 0) {
                        //                m_DomainFields[f].EvaluateMean(jCell2, 1, m_MeanFieldValuesOT[f], jEdge, 0.0);
                        //            }
                        //        }


                        //        if(f < m_FieldGradientIN.Length && m_FieldGradientIN[f] != null) {
                        //            Debug.Assert((m_FieldGradientIN[f] == null) || (m_FieldGradientOT[f] != null));

                        //            m_DomainFields[f].EvaluateGradient(jCell1, 1, qrNodes.GetVolumeNodeSet(grid, edge1), m_FieldGradientIN[f], jEdge, 0.0);

                        //            if(jCell2 >= 0) {
                        //                m_DomainFields[f].EvaluateGradient(jCell2, 1, qrNodes.GetVolumeNodeSet(grid, edge2), m_FieldGradientOT[f], jEdge, 0.0);
                        //            }
                        //        }
                        //    }
                        //}


                    }

                    NoOfSec++;
                }

                this.Basis_Eval.Stop();



                #endregion

                // ===================================================================
                // Evaluate all fields/transform nodes to global space/and the normals
                // ===================================================================
                #region EVAL_ALL_FIELDS


                this.ParametersAndNormals.Start();

                // evaluate normals and quadrature transformation metric
                MultidimensionalArray NormalsGlobalCoords = grid.iGeomEdges.NormalsCache.GetNormals_Edge(qrNodes, i0, Length);
                MultidimensionalArray QuadScalings;
                if (affine) {
                    QuadScalings = grid.iGeomEdges.SqrtGramian.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 });
                } else {
                    QuadScalings = grid.iGeomEdges.NormalsCache.GetIntegrationMetric(qrNodes, i0, Length);
                }

                // nodes in global coordinates
                MultidimensionalArray NodesGlobalCoords = grid.GlobalNodes.GetValue_EdgeSV(qrNodes, i0, Length);

                this.ParametersAndNormals.Stop();


                this.Field_Eval.Start();

                for (int f = 0; f < NoOfFields; f++) {
                    if (m_owner.m_DomainAndParamFields[f] != null) {
                        m_owner.m_DomainAndParamFields[f].EvaluateEdge(i0, Length, qrNodes,
                            m_FieldValuesIN[f], m_FieldValuesOT[f],
                            f < m_MeanFieldValuesIN.Length ? m_MeanFieldValuesIN[f] : null, f < m_MeanFieldValuesOT.Length ? m_MeanFieldValuesOT[f] : null,
                            f < m_FieldGradientIN.Length ? m_FieldGradientIN[f] : null, f < m_FieldGradientOT.Length ? m_FieldGradientOT[f] : null,
                            0, 0.0);
                    }
                }

                NormalsGlobalCoords = PeriodicTransform(i0, Length, grid, D, NoOfNodes, NormalsGlobalCoords);


#if DEBUG
                for(int f = 0; f < NoOfFields; f++) {
                    if(m_owner.m_DomainAndParamFields[f] == null) {

                        if(m_FieldValuesIN[f] != null) {
                            Debug.Assert(m_FieldValuesIN[f].L2Norm() == 0.0);
                            Debug.Assert(m_FieldValuesOT[f].L2Norm() == 0.0);
                        }


                        if(f < m_MeanFieldValuesIN.Length && m_MeanFieldValuesIN[f] != null) {
                            Debug.Assert(m_MeanFieldValuesIN[f].L2Norm() == 0.0);
                            Debug.Assert(m_MeanFieldValuesOT[f].L2Norm() == 0.0);
                        }


                        if(f < m_FieldGradientIN.Length && m_FieldGradientIN[f] != null) {
                            Debug.Assert(m_FieldGradientIN[f].L2Norm() == 0.0);
                            Debug.Assert(m_FieldGradientOT[f].L2Norm() == 0.0);
                        }

                    }
                }
#endif



                this.Field_Eval.Stop();

                #endregion


                // =======================
                // Evaluate Flux functions
                // =======================
                #region FLUX_EVAL

                bool MustLock = this.m_owner.Operator.FluxesAreNOTMultithreadSafe;

                this.Flux_Eval.Start();

                if (MustLock)
                    Monitor.Enter(m_owner);

                {
                    // loop over all equations ...
                    for (int _e = 0; _e < NoOfEquations; _e++) {
                        int e = (_e + this.m_iThread) % NoOfEquations; // shuffling in threads to reduce locking

                        if (m_FluxValuesIN[e] != null) {
                            m_FluxValuesIN[e].Clear();
                            m_FluxValuesOT[e].Clear();
                        }

                        MultidimensionalArray FluxValuesIN = m_FluxValuesIN[e];
                        MultidimensionalArray FluxValuesOT = m_FluxValuesOT[e];

                        // -------------------------------
                        // All INonlinearFlux - Components
                        // -------------------------------

                        EvalFlux(m_NonlinFluxes[e], i0, Length, grid, NoOfSec, true, false, base.m_NonlinFluxesWatches[e],
                            delegate (INonlinearFlux nonlinFlx, int _jEdge, int _IndexOffset, int _L, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] UinMean, MultidimensionalArray[] UoutMean, MultidimensionalArray[] UinGrad, MultidimensionalArray[] UoutGrad) {
                                nonlinFlx.InnerEdgeFlux(m_owner.m_Time, _jEdge,
                                  NodesGlobalCoords,
                                  NormalsGlobalCoords,
                                  Uin, Uout,
                                  _IndexOffset, _L,
                                  FluxValuesIN);
                            },
                            delegate (INonlinearFlux nonlinFlx, int _jEdge, int _IndexOffset, int _L, int _EdgeTagsOffset, bool flipNormal, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] UinMean, MultidimensionalArray[] UinGrad) {
                                nonlinFlx.BorderEdgeFlux(m_owner.m_Time, _jEdge,
                                                         NodesGlobalCoords,
                                                         NormalsGlobalCoords, flipNormal,
                                                         grid.iGeomEdges.EdgeTags, _EdgeTagsOffset,
                                                         Uin,
                                                         _IndexOffset, _L,
                                                         FluxValuesIN);
                            });
#if DEBUG
                        if (FluxValuesIN != null)
                            FluxValuesIN.CheckForNanOrInf(true, true, true);
#endif

                        // ---------------------------------
                        // All INonlinearFluxEx - Components
                        // ---------------------------------
                        EvalFlux(m_NonlinFluxesEx[e], i0, Length, grid, NoOfSec, true, false, base.m_NonlinFluxesExWatches[e],
                            delegate (INonlinearFluxEx nonlinFlx, int _jEdge, int _IndexOffset, int _L, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] UinMean, MultidimensionalArray[] UoutMean, MultidimensionalArray[] UinGrad, MultidimensionalArray[] UoutGrad) {
                                nonlinFlx.InnerEdgeFlux(m_owner.m_Time, _jEdge,
                                                        NodesGlobalCoords,
                                                        NormalsGlobalCoords,
                                                        Uin, Uout, UinMean, UoutMean,
                                                        _IndexOffset, _L,
                                                        FluxValuesIN);
                            },
                            delegate (INonlinearFluxEx nonlinFlx, int _jEdge, int _IndexOffset, int _L, int _EdgeTagsOffset, bool flipNormal, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] UinMean, MultidimensionalArray[] UinGrad) {
                                nonlinFlx.BorderEdgeFlux(m_owner.m_Time, _jEdge,
                                                         NodesGlobalCoords,
                                                         NormalsGlobalCoords, flipNormal,
                                                         grid.iGeomEdges.EdgeTags, _EdgeTagsOffset,
                                                         Uin, UinMean,
                                                         _IndexOffset, _L,
                                                         FluxValuesIN);
                            });
#if DEBUG
                        if (FluxValuesIN != null)
                            FluxValuesIN.CheckForNanOrInf(true, true, true);
#endif

                        // set out-cell flux
                        if (FluxValuesIN != null)
                            FluxValuesOT.Acc(-1.0, FluxValuesIN);


                        // ----------------------------------
                        // All INonlinEdgeform_V - components
                        // ----------------------------------

                        //bla 1
                        EvalFlux(m_EdgeForm_V[e], i0, Length, grid, NoOfSec, false, true, this.m_EdgeForm_V_Watches[e],
                            delegate (INonlinEdgeForm_V edgeform, int _jEdge, int _IndexOffset, int _L, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] UinMean, MultidimensionalArray[] UoutMean, MultidimensionalArray[] UinGrad, MultidimensionalArray[] UoutGrad) {
                                EdgeFormParams efp;
                                efp.GridDat = m_owner.GridDat;
                                efp.Len = _L;
                                efp.e0 = _jEdge;
                                efp.time = m_owner.m_Time;
                                Debug.Assert(NoOfNodes == NormalsGlobalCoords.GetLength(1));
                                Debug.Assert(NoOfNodes == NodesGlobalCoords.GetLength(1));

                                int[] I0vec = new int[] { _IndexOffset, 0, 0 };
                                int[] IEvec = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1, D - 1 };
                                int[] I0scl = new int[] { _IndexOffset, 0 };
                                int[] IEscl = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1 };

                                efp.Normals = NormalsGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);
                                efp.Nodes = NodesGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);

                                var _FluxValuesIN = FluxValuesIN.ExtractSubArrayShallow(I0scl, IEscl);
                                var _FluxValuesOT = FluxValuesOT.ExtractSubArrayShallow(I0scl, IEscl);

                                MultidimensionalArray[] _Uin = new MultidimensionalArray[Uin.Length], _Uout = new MultidimensionalArray[Uout.Length];
                                for (int i = 0; i < Uin.Length; i++) {
                                    _Uin[i] = Uin[i] != null ? Uin[i].ExtractSubArrayShallow(I0scl, IEscl) : null;
                                    _Uout[i] = Uout[i] != null ? Uout[i].ExtractSubArrayShallow(I0scl, IEscl) : null;
                                }
                                MultidimensionalArray[] _UinGrad = new MultidimensionalArray[UinGrad.Length], _UoutGrad = new MultidimensionalArray[UoutGrad.Length];
                                for (int i = 0; i < UinGrad.Length; i++) {
                                    _UinGrad[i] = UinGrad[i] != null ? UinGrad[i].ExtractSubArrayShallow(I0vec, IEvec) : null;
                                    _UoutGrad[i] = UoutGrad[i] != null ? UoutGrad[i].ExtractSubArrayShallow(I0vec, IEvec) : null;
                                }
                                efp.ParameterVars_IN = _Uin.GetSubVector(NoArgs, NoParams);
                                efp.ParameterVars_OUT = _Uout.GetSubVector(NoArgs, NoParams);

                                edgeform.NonlinInternalEdge_V(ref efp, _Uin.GetSubVector(0, NoArgs), _Uout.GetSubVector(0, NoArgs), _UinGrad, _UoutGrad, _FluxValuesIN, _FluxValuesOT);
                            },
                            delegate (INonlinEdgeForm_V nonlinFlx, int _jEdge, int _IndexOffset, int _L, int _EdgeTagsOffset, bool flipNormal, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] UinMean, MultidimensionalArray[] UinGrad) {
                                EdgeFormParams efp;
                                efp.GridDat = m_owner.GridDat;
                                efp.Len = _L;
                                efp.e0 = _jEdge;
                                efp.time = m_owner.m_Time;
                                Debug.Assert(NoOfNodes == NormalsGlobalCoords.GetLength(1));
                                Debug.Assert(NoOfNodes == NodesGlobalCoords.GetLength(1));

                                int[] I0vec = new int[] { _IndexOffset, 0, 0 };
                                int[] IEvec = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1, D - 1 };
                                int[] I0scl = new int[] { _IndexOffset, 0 };
                                int[] IEscl = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1 };

                                efp.Normals = NormalsGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);
                                efp.Nodes = NodesGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);

                                var _FluxValuesIN = (flipNormal ? FluxValuesOT : FluxValuesIN).ExtractSubArrayShallow(I0scl, IEscl);

                                MultidimensionalArray[] _Uin = new MultidimensionalArray[Uin.Length];
                                for (int i = 0; i < Uin.Length; i++) {
                                    _Uin[i] = Uin[i] != null ? Uin[i].ExtractSubArrayShallow(I0scl, IEscl) : null;
                                }
                                MultidimensionalArray[] _UinGrad = new MultidimensionalArray[UinGrad.Length];
                                for (int i = 0; i < UinGrad.Length; i++) {
                                    _UinGrad[i] = UinGrad[i] != null ? UinGrad[i].ExtractSubArrayShallow(I0vec, IEvec) : null;
                                }
                                efp.ParameterVars_IN = _Uin.GetSubVector(NoArgs, NoParams);
                                efp.ParameterVars_OUT = null;

                                if (flipNormal)
                                    efp.Normals.Scale(-1.0);
                                nonlinFlx.NonlinBoundaryEdge_V(ref efp, _Uin.GetSubVector(0, NoArgs), _UinGrad, _FluxValuesIN);
                                if (flipNormal)
                                    efp.Normals.Scale(-1.0);
                            });
#if DEBUG
                        if (FluxValuesIN != null)
                            FluxValuesIN.CheckForNanOrInf(true, true, true);
                        if (FluxValuesOT != null)
                            FluxValuesOT.CheckForNanOrInf(true, true, true);
#endif

                        // --------------------------------------
                        // All INonlinEdgeform_GradV - components
                        // --------------------------------------

                        if (m_GradientFluxValuesIN[e] != null) {
                            m_GradientFluxValuesIN[e].Clear();
                            m_GradientFluxValuesOT[e].Clear();
                        }


                        EvalFlux(m_EdgeForm_GradV[e], i0, Length, grid, NoOfSec, false, true, this.m_EdgeForm_GradV_Watches[e],
                            delegate (INonlinEdgeForm_GradV edgeform, int _jEdge, int _IndexOffset, int _L, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] UinMean, MultidimensionalArray[] UoutMean, MultidimensionalArray[] UinGrad, MultidimensionalArray[] UoutGrad) {
                                EdgeFormParams efp;
                                efp.GridDat = m_owner.GridDat;
                                efp.Len = _L;
                                efp.e0 = _jEdge;
                                efp.time = m_owner.m_Time;

                                Debug.Assert(NoOfNodes == NormalsGlobalCoords.GetLength(1));
                                Debug.Assert(NoOfNodes == NodesGlobalCoords.GetLength(1));
                                int[] I0vec = new int[] { _IndexOffset, 0, 0 };
                                int[] IEvec = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1, D - 1 };
                                int[] I0scl = new int[] { _IndexOffset, 0 };
                                int[] IEscl = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1 };

                                efp.Normals = NormalsGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);
                                efp.Nodes = NodesGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);

                                var _GradFluxIN = m_GradientFluxValuesIN[e].ExtractSubArrayShallow(I0vec, IEvec);
                                var _GradFluxOT = m_GradientFluxValuesOT[e].ExtractSubArrayShallow(I0vec, IEvec);

                                MultidimensionalArray[] _Uin = new MultidimensionalArray[Uin.Length], _Uout = new MultidimensionalArray[Uout.Length];
                                for (int i = 0; i < Uin.Length; i++) {
                                    _Uin[i] = Uin[i] != null ? Uin[i].ExtractSubArrayShallow(I0scl, IEscl) : null;
                                    _Uout[i] = Uout[i] != null ? Uout[i].ExtractSubArrayShallow(I0scl, IEscl) : null;
                                }
                                MultidimensionalArray[] _UinGrad = new MultidimensionalArray[UinGrad.Length], _UoutGrad = new MultidimensionalArray[UoutGrad.Length];
                                for (int i = 0; i < UinGrad.Length; i++) {
                                    _UinGrad[i] = UinGrad[i] != null ? UinGrad[i].ExtractSubArrayShallow(I0vec, IEvec) : null;
                                    _UoutGrad[i] = UoutGrad[i] != null ? UoutGrad[i].ExtractSubArrayShallow(I0vec, IEvec) : null;
                                }
                                efp.ParameterVars_IN = _Uin.GetSubVector(NoArgs, NoParams);
                                efp.ParameterVars_OUT = _Uout.GetSubVector(NoArgs, NoParams);

                                edgeform.NonlinInternalEdge_GradV(ref efp, _Uin.GetSubVector(0, NoArgs), _Uout.GetSubVector(0, NoArgs), _UinGrad, _UoutGrad, _GradFluxIN, _GradFluxOT);
                            },
                            delegate (INonlinEdgeForm_GradV nonlinFlx, int _jEdge, int _IndexOffset, int _L, int _EdgeTagsOffset, bool flipNormal, int NoArgs, int NoParams, MultidimensionalArray[] Uin, MultidimensionalArray[] UinMean, MultidimensionalArray[] UinGrad) {
                                EdgeFormParams efp;
                                efp.GridDat = m_owner.GridDat;
                                efp.Len = _L;
                                efp.e0 = _jEdge;
                                efp.time = m_owner.m_Time;
                                Debug.Assert(NoOfNodes == NormalsGlobalCoords.GetLength(1));
                                Debug.Assert(NoOfNodes == NodesGlobalCoords.GetLength(1));
                                int[] I0vec = new int[] { _IndexOffset, 0, 0 };
                                int[] IEvec = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1, D - 1 };
                                int[] I0scl = new int[] { _IndexOffset, 0 };
                                int[] IEscl = new int[] { _IndexOffset + _L - 1, NoOfNodes - 1 };

                                efp.Normals = NormalsGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);
                                efp.Nodes = NodesGlobalCoords.ExtractSubArrayShallow(I0vec, IEvec);

                                var _GradFluxIN = (flipNormal ? m_GradientFluxValuesOT[e] : m_GradientFluxValuesIN[e]).ExtractSubArrayShallow(I0vec, IEvec);

                                MultidimensionalArray[] _Uin = new MultidimensionalArray[Uin.Length];
                                for (int i = 0; i < Uin.Length; i++) {
                                    _Uin[i] = Uin[i] != null ? Uin[i].ExtractSubArrayShallow(I0scl, IEscl) : null;
                                }
                                MultidimensionalArray[] _UinGrad = new MultidimensionalArray[UinGrad.Length];
                                for (int i = 0; i < UinGrad.Length; i++) {
                                    _UinGrad[i] = UinGrad[i] != null ? UinGrad[i].ExtractSubArrayShallow(I0vec, IEvec) : null;
                                }
                                efp.ParameterVars_IN = _Uin.GetSubVector(NoArgs, NoParams);
                                efp.ParameterVars_OUT = null;

                                if (flipNormal)
                                    efp.Normals.Scale(-1.0);
                                nonlinFlx.NonlinBoundaryEdge_GradV(ref efp, _Uin.GetSubVector(0, NoArgs), _UinGrad, _GradFluxIN);
                                if (flipNormal)
                                    efp.Normals.Scale(-1.0);
                            });
#if DEBUG
                        if (m_GradientFluxValuesIN[e] != null)
                            m_GradientFluxValuesIN[e].CheckForNanOrInf(true, true, true);
                        if (m_GradientFluxValuesOT[e] != null)
                            m_GradientFluxValuesOT[e].CheckForNanOrInf(true, true, true);
#endif
                    }

                }
                this.Flux_Eval.Stop();

                #endregion

                // ===================
                // Flux Transformation
                // ===================
                #region FLUX_TRAFO
                this.Flux_Trafo.Start();


                // multiply fluxes with Jacobi determinant (integral transformation metric):
                for (int i = 0; i < NoOfEquations; i++) {
                    if (m_FluxValuesIN[i] != null) {
                        m_FluxValuesIN[i].Multiply(1.0, m_FluxValuesIN[i], QuadScalings, 0.0, "jk", "jk", affine ? "j" : "jk");
                        m_FluxValuesOT[i].Multiply(1.0, m_FluxValuesOT[i], QuadScalings, 0.0, "jk", "jk", affine ? "j" : "jk");
                    }

                    if (m_GradientFluxValuesIN[i] != null) {
                        m_GradientFluxValuesIN[i].Multiply(1.0, m_GradientFluxValuesIN[i], QuadScalings, 0.0, "jkd", "jkd", affine ? "j" : "jk");
                        m_GradientFluxValuesOT[i].Multiply(1.0, m_GradientFluxValuesOT[i], QuadScalings, 0.0, "jkd", "jkd", affine ? "j" : "jk");
                    }
                }

                // transform gradient fluxes -- less costly than transforming basis
                if (m_owner.maxTestGradientBasis != null) {
#if DEBUG
                    for (int i = 0; i < NoOfEquations; i++) {
                        Debug.Assert((m_GradientFluxValuesINtrf[i] != null) == (m_GradientFluxValuesOTtrf[i] != null));
                        Debug.Assert((m_GradientFluxValuesIN[i] != null) == (m_GradientFluxValuesINtrf[i] != null));
                        Debug.Assert((m_GradientFluxValuesOT[i] != null) == (m_GradientFluxValuesOTtrf[i] != null));
                    }
#endif

                    if (affine) {
                        MultidimensionalArray invJacobi = grid.iGeomCells.InverseTransformation;
                        int[,] Edge2Cell = grid.iGeomEdges.CellIndices;

                        unsafe {
                            fixed (int* pEdge2Cell = Edge2Cell) {
                                for (int i = 0; i < NoOfEquations; i++) {

                                    if (m_GradientFluxValuesINtrf[i] != null) {
                                        m_GradientFluxValuesINtrf[i].Multiply(1.0, m_GradientFluxValuesIN[i], invJacobi, 0.0, ref mp_jke_jkd_Tjed,
                                            pEdge2Cell, pEdge2Cell,
                                            trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + 0), trfCycle_B: 2, trfPostOffset_B: 0);
                                        m_GradientFluxValuesOTtrf[i].Multiply(1.0, m_GradientFluxValuesOT[i], invJacobi, 0.0, ref mp_jke_jkd_Tjed,
                                            pEdge2Cell, pEdge2Cell,
                                            trfPreOffset_A: 0, trfCycle_A: 0, trfPostOffset_A: 0, trfPreOffset_B: (2 * i0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);
                                    }
                                }
                            }
                        }
                    } else {
                        MultidimensionalArray invJacobiIN, invJacobiOT;
                        var TiJ = grid.InverseJacobian.GetValue_EdgeDV(qrNodes, i0, Length);
                        invJacobiIN = TiJ.Item1;
                        invJacobiOT = TiJ.Item2;

                        for (int i = 0; i < NoOfEquations; i++) {
                            if (m_GradientFluxValuesINtrf[i] != null) {
                                m_GradientFluxValuesINtrf[i].Multiply(1.0, m_GradientFluxValuesIN[i], invJacobiIN, 0.0, "jke", "jkd", "jked");
                                m_GradientFluxValuesOTtrf[i].Multiply(1.0, m_GradientFluxValuesOT[i], invJacobiOT, 0.0, "jke", "jkd", "jked");
                            }
                        }

                    }
                }

                if (MustLock)
                    Monitor.Exit(m_owner);

                this.Flux_Trafo.Stop();
                #endregion


                // ==========================================
                // multiply with test functions / save result
                // ==========================================
                #region QUADLOOPS
                this.Loops.Start();
                {
                    MultidimensionalArray OrthoTrf = null; // to transform back to ONB on physical space...
                    if (affine) {
                        OrthoTrf = null;
                    } else {
                        int MaxDegree = Math.Max(m_owner.maxTestBasis != null ? m_owner.maxTestBasis.Degree : 0, m_owner.maxTestGradientBasis != null ? m_owner.maxTestGradientBasis.Degree : 0);
                        OrthoTrf = grid.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(jCellMin, jCellMax - jCellMin + 1, MaxDegree); // sollte irgendwann sowieso
                                                                                                                                       //                                                                                                                für alle Zellen vorliegen, also kein Stress
                                                                                                                                       //                                                                                                                falls jCellMin und jCellMax weit auseinander
                    }


                    int N0 = 0;
                    for (int gamma = 0; gamma < NoOfEquations; gamma++) {
                        Basis testBasis = m_owner.m_CodomainBasisS[gamma];
                        int Ne = testBasis.Length;

                        MultidimensionalArray QuadResultIN, QuadResultOT;
                        if (affine) {
                            QuadResultIN = QuadResult.ExtractSubArrayShallow(new int[] { 0, N0, 0 }, new int[] { Length - 1, N0 + Ne - 1, -11 });
                            QuadResultOT = QuadResult.ExtractSubArrayShallow(new int[] { 0, N0, 1 }, new int[] { Length - 1, N0 + Ne - 1, -11 });
                        } else {
                            QuadResultIN = this.m_QuadResultIN[gamma];
                            QuadResultOT = this.m_QuadResultOT[gamma];
                        }

                        double cF = 0.0;
                        if (m_owner.maxTestBasis != null) {

                            int Nmax = TstFuncXwgt.GetLength(2);
                            MultidimensionalArray _TestFunction;
                            if (Nmax > Ne) {
                                _TestFunction = TstFuncXwgt.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { TstFuncXwgt.GetLength(0) - 1, NoOfNodes - 1, Ne - 1 });
                            } else {
                                _TestFunction = TstFuncXwgt;
                            }
                            if (m_FluxValuesIN[gamma] != null) {
                                int[,] trfIdx = grid.iGeomEdges.Edge2CellTrafoIndex;
                                unsafe {
                                    fixed (int* pTrfIdx = trfIdx) {
                                        // QuadResultIN[j,n] = sum_{k}  _TestFunction[T(j),k,n]*m_FluxValuesIN[gamma][j,k] 
                                        //    where j: edge index,
                                        //          n: DG mode index/test function index
                                        //          k: quadrature node index
                                        //       T(j) = trfIdx[i0 + j, 0]
                                        QuadResultIN.Multiply(1.0, _TestFunction, m_FluxValuesIN[gamma], cF, ref mp_jn_Tjkn_jk,
                                            pTrfIdx, pTrfIdx,
                                            trfPreOffset_A: (2 * i0), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);

                                        // analog, aber mit T(j) = trfIdx[i0 + j, 1]
                                        QuadResultOT.Multiply(1.0, _TestFunction, m_FluxValuesOT[gamma], cF, ref mp_jn_Tjkn_jk,
                                            pTrfIdx, pTrfIdx,
                                            trfPreOffset_A: (2 * i0 + 1), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                    }

                                }
                                cF = 1;
                            }

                        }

                        if (m_owner.maxTestGradientBasis != null) {
                            int Nmax = TstFuncGradXwgt.GetLength(2);
                            MultidimensionalArray _TstFuncGradXwgt;
                            if (Nmax > Ne) {
                                _TstFuncGradXwgt = TstFuncGradXwgt.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { TstFuncGradXwgt.GetLength(0) - 1, NoOfNodes - 1, Ne - 1, D - 1 });
                            } else {
                                _TstFuncGradXwgt = TstFuncGradXwgt;
                            }

                            if (m_GradientFluxValuesINtrf[gamma] != null) {
                                int[,] trfIdx = grid.iGeomEdges.Edge2CellTrafoIndex;
                                unsafe {
                                    fixed (int* pTrfIdx = trfIdx) {
                                        // QuadResultIN[j,n] = sum_{k,d}  _TestFunctionGradient[T(j),k,n,d]*m_GradientFluxValuesINtrf[gamma][j,k,d] 
                                        //   ansonsten wie oben
                                        QuadResultIN.Multiply(1.0, _TstFuncGradXwgt, m_GradientFluxValuesINtrf[gamma], cF, ref mp_jn_Tjknd_jkd,
                                            pTrfIdx, pTrfIdx,
                                            trfPreOffset_A: (2 * i0 + 0), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                        QuadResultOT.Multiply(1.0, _TstFuncGradXwgt, m_GradientFluxValuesOTtrf[gamma], cF, ref mp_jn_Tjknd_jkd,
                                            pTrfIdx, pTrfIdx,
                                            trfPreOffset_A: (2 * i0 + 1), trfCycle_A: 2, trfPostOffset_A: 0, trfPreOffset_B: 0, trfCycle_B: 0, trfPostOffset_B: 0);
                                    }
                                }
                                cF = 1;
                            }


                        }

                        if (cF == 0) {
                            QuadResultIN.Clear();
                            QuadResultOT.Clear();
                        }


                        // final transformations after quadrature ...

                        if (affine) {
                            // nop: done for all components at once, see bolow


                        } else {

                            MultidimensionalArray QuadResultINtrf, QuadResultOTtrf;
                            QuadResultINtrf = QuadResult.ExtractSubArrayShallow(new int[] { 0, N0, 0 }, new int[] { Length - 1, N0 + Ne - 1, -1 });
                            QuadResultOTtrf = QuadResult.ExtractSubArrayShallow(new int[] { 0, N0, 1 }, new int[] { Length - 1, N0 + Ne - 1, -1 });


                            MultidimensionalArray _OrthoTrfIN, _OrthoTrfOT;
                            if (OrthoTrf.GetLength(1) == Ne) {
                                _OrthoTrfIN = OrthoTrf;
                                _OrthoTrfOT = _OrthoTrfIN;
                            } else {
                                _OrthoTrfIN = OrthoTrf.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { OrthoTrf.GetLength(0) - 1, Ne - 1, Ne - 1 });
                                _OrthoTrfOT = _OrthoTrfIN;
                            }

                            int[,] Edge2Cell = grid.iGeomEdges.CellIndices;

                            unsafe {
                                fixed (int* pEdge2Cell = Edge2Cell) {

                                    QuadResultINtrf.Multiply(1.0, _OrthoTrfIN, QuadResultIN, 0.0, ref mp_jn_Tjmn_jm,
                                        pEdge2Cell, pEdge2Cell,
                                        trfPreOffset_A: (2 * i0), trfCycle_A: 2, trfPostOffset_A: -jCellMin, trfPostOffset_B: 0, trfCycle_B: 0, trfPreOffset_B: 0);
                                    QuadResultOTtrf.Multiply(1.0, _OrthoTrfOT, QuadResultOT, 0.0, ref mp_jn_Tjmn_jm,
                                        pEdge2Cell, pEdge2Cell,
                                        trfPreOffset_A: (2 * i0 + 1), trfCycle_A: 2, trfPostOffset_A: -jCellMin, trfPostOffset_B: 0, trfCycle_B: 0, trfPreOffset_B: 0);
                                }
                            }
                        }

                        N0 += Ne;
                    }

                    if (affine) {
                        // apply scaling of test functions in physcal domain

                        MultidimensionalArray QuadResultINtrf, QuadResultOTtrf;
                        QuadResultINtrf = QuadResult.ExtractSubArrayShallow(-1, -1, 0);
                        QuadResultOTtrf = QuadResult.ExtractSubArrayShallow(-1, -1, 1);

                        MultidimensionalArray scaling = grid.ChefBasis.Scaling;

                        int[,] Edge2Cell = grid.iGeomEdges.CellIndices;

                        unsafe {
                            fixed (int* pEdge2Cell = Edge2Cell) {
                                // QuadResultINtrf[i,n] = QuadResultINtrf[i,n]*scaling[T(i)], where T(i) = Edge2Cell[i,0]
                                QuadResultINtrf.Multiply(1.0, QuadResultINtrf, scaling, 0.0, ref mp_in_in_Ti,
                                    pEdge2Cell, pEdge2Cell,
                                    trfPostOffset_A: 0, trfCycle_A: 0, trfPreOffset_A: 0, trfPreOffset_B: (2 * i0), trfCycle_B: 2, trfPostOffset_B: 0);

                                // QuadResultINtrf[i,n] = QuadResultINtrf[i,n]*scaling[T(i)], where T(i) = Edge2Cell[i,1]
                                QuadResultOTtrf.Multiply(1.0, QuadResultOTtrf, scaling, 0.0, ref mp_in_in_Ti,
                                    pEdge2Cell, pEdge2Cell,
                                    trfPostOffset_A: 0, trfCycle_A: 0, trfPreOffset_A: 0, trfPreOffset_B: (2 * i0 + 1), trfCycle_B: 2, trfPostOffset_B: 0);
                            }
                        }

                    } else {
                        // nop
                    }

                }
                this.Loops.Stop();
                #endregion


#if DEBUG
            QuadResult.CheckForNanOrInf(true, true, true);
#endif
            }

            private MultidimensionalArray PeriodicTransform(int i0, int Length, IGridData grid, int D, int NoOfNodes, MultidimensionalArray NormalsGlobalCoords) {
                if (m_owner._PeriodicVectorTrafo != PeriodicVectorTrafo.nix) {
                    // ++++++++++++++++++++++++++++
                    // Periodic transformation
                    // ++++++++++++++++++++++++++++

                    // If the grid contains some periodic boundaries which are not parallel (e.g. some cake-pie-subsection of a rotational domain)
                    // periodicity required additional transformations/rotations of vectors for **both** sides of the periodic edge;
                    // Furthermore, these rotations are different (inverse) for the in- and the out-edge, 
                    // therefore the contribution to the out-cell is computed in a second pass, by this integrator.


                    var EdgeTags = grid.iGeomEdges.EdgeTags;
                    int[][] VectorTuples = m_owner.VectorComponentIndices;
                    List<int> ScalarIndices = new List<int>(m_owner.m_DomainAndParamFields.Length.ForLoop(i => i));
                    foreach (int[] tt in VectorTuples) {
                        foreach (int iVar in tt) {
                            ScalarIndices.Remove(iVar);
                        }
                    }




                    for (int i = 0; i < Length; i++) {
                        int iEdge = i + i0;
                        var eT = EdgeTags[iEdge];
                        if (eT >= Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG) {

                            //(new CellMask(GridDat, Chunk.GetSingleElementChunk(GridDat.iLogicalEdges.CellIndices[iEdge, 0]))).SaveToTextFile("PerInCell.csv", false);
                            //(new CellMask(GridDat, Chunk.GetSingleElementChunk(GridDat.iLogicalEdges.CellIndices[iEdge, 1]))).SaveToTextFile("PerOtCell.csv", false);


                            var Trafo = grid.Grid.PeriodicTrafo[eT - Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG];


                            void TransformVector(MultidimensionalArray Mtx, MultidimensionalArray[] Vals) {
                                Vector Uorg = new Vector(D);
                                foreach (int[] VecPtr in VectorTuples) {
                                    Debug.Assert(VecPtr.Length == D);

                                    for (int k = 0; k < NoOfNodes; k++) {
                                        if (VecPtr[0] >= Vals.Length || Vals[VecPtr[0]] == null)
                                            continue;

                                        for (int d = 0; d < D; d++) {
                                            Uorg[d] = Vals[VecPtr[d]][i, k];
                                        }

                                        Vector Udest = Mtx.MtxVecMul(Uorg);

                                        for (int d = 0; d < D; d++) {
                                            Vals[VecPtr[d]][i, k] = Udest[d];
                                        }
                                    }
                                }
                            }

                            void TransfromGradientTensor(MultidimensionalArray Mtx, MultidimensionalArray[] Gradients) {

                                var MtxT = Mtx.TransposeTo();
                                var GT = MultidimensionalArray.Create(D, D);

                                foreach (int[] VecPtr in VectorTuples) {
                                    Debug.Assert(VecPtr.Length == D);

                                    /*
                                    Vector Uorg = new Vector(D);
                                    for(int k = 0; k < NoOfNodes; k++) {
                                        if(VecPtr[0] >= Gradients.Length || Gradients[VecPtr[0]] == null)
                                            continue;

                                        for(int e = 0; e < D; e++) {
                                            for(int d = 0; d < D; d++) {
                                                Uorg[d] = Gradients[VecPtr[d]][i, k, e];
                                            }

                                            Vector Udest = Mtx.MtxVecMul(Uorg);

                                            for(int d = 0; d < D; d++) {
                                                Gradients[VecPtr[d]][i, k, e] = Udest[d];
                                            }
                                        }
                                    }
                                    */

                                    for (int k = 0; k < NoOfNodes; k++) {
                                        if (VecPtr[0] >= Gradients.Length || Gradients[VecPtr[0]] == null)
                                            continue;

                                        for (int e = 0; e < D; e++) {
                                            for (int d = 0; d < D; d++) {
                                                GT[e, d] = Gradients[VecPtr[e]][i, k, d];
                                            }
                                        }

                                        var GTD = Mtx.GEMM(GT, MtxT);


                                        for (int e = 0; e < D; e++) {
                                            for (int d = 0; d < D; d++) {
                                                Gradients[VecPtr[e]][i, k, d] = GTD[e, d];
                                            }
                                        }
                                    }
                                }
                            }

                            void TransfromGradientVector(MultidimensionalArray Mtx, MultidimensionalArray[] Gradients) {

                                //var MtxT = Mtx.TransposeTo();
                                //var GT = MultidimensionalArray.Create(D, D);

                                foreach (int iVar in ScalarIndices) {
                                    if (iVar >= Gradients.Length)
                                        continue;
                                    MultidimensionalArray Gradient = Gradients[iVar];
                                    if (Gradient == null)
                                        continue;

                                    Vector Uorg = new Vector(D);
                                    for (int k = 0; k < NoOfNodes; k++) {


                                        for (int d = 0; d < D; d++) {
                                            Uorg[d] = Gradient[i, k, d];
                                        }

                                        Vector Udest = Mtx.MtxVecMul(Uorg);

                                        for (int d = 0; d < D; d++) {
                                            Gradient[i, k, d] = Udest[d];
                                        }

                                    }

                                }
                            }


                            if (m_owner._PeriodicVectorTrafo == PeriodicVectorTrafo.fwd) {
                                var MatrixFW = Trafo.Matrix.TransposeTo(); // transform from IN to OT
                                MatrixFW.InvertInPlace();

                                TransformVector(MatrixFW, m_FieldValuesIN);
                                TransformVector(MatrixFW, m_MeanFieldValuesIN);
                                TransfromGradientTensor(MatrixFW, m_FieldGradientIN);
                                TransfromGradientVector(MatrixFW, m_FieldGradientIN);

                                for (int k = 0; k < NoOfNodes; k++) {
                                    for (int d = 0; d < D; d++) {
                                        //m_FieldGradientIN[0][i, k, d] = m_FieldGradientOT[0][i, k, d];
                                    }
                                }


                                // transform normal
                                Vector N = new Vector(D);
                                NormalsGlobalCoords = NormalsGlobalCoords.CloneAs(); // Don't mess with cached normals.
                                for (int k = 0; k < NoOfNodes; k++) {
                                    for (int d = 0; d < D; d++) {
                                        N[d] = NormalsGlobalCoords[i, k, d];
                                    }

                                    Vector ND = MatrixFW.MtxVecMul(N);

                                    for (int d = 0; d < D; d++) {
                                        NormalsGlobalCoords[i, k, d] = ND[d];
                                    }
                                }
                            } else if (m_owner._PeriodicVectorTrafo == PeriodicVectorTrafo.bck) {
                                var MatrixBK = Trafo.Matrix.TransposeTo(); // transform from OT to IN

                                TransformVector(MatrixBK, m_FieldValuesOT);
                                TransformVector(MatrixBK, m_MeanFieldValuesOT);
                                TransfromGradientTensor(MatrixBK, m_FieldGradientOT);
                                TransfromGradientVector(MatrixBK, m_FieldGradientOT);

                                for (int k = 0; k < NoOfNodes; k++) {
                                    for (int d = 0; d < D; d++) {
                                        //m_FieldGradientOT[0][i, k, d] = m_FieldGradientIN[0][i, k, d];
                                    }
                                }

                            } else {
                                throw new NotImplementedException();
                            }
                        }
                    }
                }

                return NormalsGlobalCoords;
            }


            private void EvalFlux<T>(EquationComponentArgMapping<T> components, int i0, int Length, IGridData grid, int NoOfSec,
                bool MapAlsoMean, bool MapAlsoGradient, Stopwatch[] timers,
                Action<T, int, int, int, int, int, MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[]> CallInner,
                Action<T, int, int, int, int, bool, int, int, MultidimensionalArray[], MultidimensionalArray[], MultidimensionalArray[]> CallBorder)
                where T : IEquationComponent //
            {
                int L = components.m_AllComponentsOfMyType.Length;
                if (L > 0) {
                    Debug.Assert(timers.Length == components.m_AllComponentsOfMyType.Length);

                    // sum up all fluxes - inner edges
                    for (int __iComp = 0; __iComp < L; __iComp++) {
                        int iComp = (__iComp + this.m_iThread) % L; // shuffling in threads to reduce locking

                        object blck = components.m_LockObjects[iComp];

                        timers[iComp].Start();
                        T nonlinFlx = components.m_AllComponentsOfMyType[iComp];
                        int NoArgs = components.NoOfArguments[iComp];
                        int NoParams = components.NoOfParameters[iComp];

                        for (int l = 0; l < NoOfSec; l++) {
                            int IndexOffset = m_InnerEdgesI0s[l];
                            int __Len = m_InnerEdgesIEs[l] - m_InnerEdgesI0s[l] + 1;
                            int jEdge = m_InnerEdgesI0s[l] + i0;

                            //TempTimers.intItems += __Len;
                            //TempTimers.intWatch.Start();

                            if(blck != null)
                                Monitor.Enter(blck);

                            CallInner(nonlinFlx, jEdge, IndexOffset, __Len, NoArgs, NoParams,
                                components.MapArguments(m_FieldValuesIN, nonlinFlx, false),
                                components.MapArguments(m_FieldValuesOT, nonlinFlx, false),
                                MapAlsoMean ? components.MapArguments(m_MeanFieldValuesIN, nonlinFlx, true) : null,
                                MapAlsoMean ? components.MapArguments(m_MeanFieldValuesOT, nonlinFlx, true) : null,
                                MapAlsoGradient ? components.MapArguments(m_FieldGradientIN, nonlinFlx, true) : null,
                                MapAlsoGradient ? components.MapArguments(m_FieldGradientOT, nonlinFlx, true) : null);

                            if (blck != null)
                                Monitor.Exit(blck);

                            //TempTimers.intWatch.Stop();

                        }
                        timers[iComp].Stop();
                    }


                    // sum up all fluxes - border edges
                    for (int __iComp = 0; __iComp < L; __iComp++) {
                        int iComp = (__iComp + this.m_iThread) % L; // shuffling in threads to reduce locking
                        object blck = components.m_LockObjects[iComp];

                        timers[iComp].Start();
                        T nonlinFlx = components.m_AllComponentsOfMyType[iComp];
                        int NoArgs = components.NoOfArguments[iComp];
                        int NoParams = components.NoOfParameters[iComp];

                        //foreach (T nonlinFlx in components.m_AllComponentsOfMyType) {

                        if (blck != null)
                            Monitor.Enter(blck);

                        int iSec = 0;
                        for (int IndexOffset = 0; IndexOffset < Length; IndexOffset++) {
                            if (iSec < NoOfSec && IndexOffset == m_InnerEdgesI0s[iSec]) {
                                IndexOffset = m_InnerEdgesIEs[iSec];
                                iSec++;
                                continue;
                            }

                            int jCell1 = grid.iGeomEdges.CellIndices[IndexOffset + i0, 0];
                            int jCell2 = grid.iGeomEdges.CellIndices[IndexOffset + i0, 1];

                            bool DomainBnd = jCell2 < 0;
                            bool SubGrdBnd = (m_owner.SubGridCellsMarker != null) && (!DomainBnd) && (m_owner.SubGridCellsMarker[jCell1] != m_owner.SubGridCellsMarker[jCell2]);
                            Debug.Assert(DomainBnd || SubGrdBnd);

                            int jEdge = IndexOffset + i0;

                            int __Len = 1;
                            for (; IndexOffset + __Len < Length; __Len++) {
                                int _jCell1 = grid.iGeomEdges.CellIndices[IndexOffset + i0 + __Len, 0];
                                int _jCell2 = grid.iGeomEdges.CellIndices[IndexOffset + i0 + __Len, 1];
                                bool _DomainBnd = _jCell2 < 0;
                                bool _SubGrdBnd = (m_owner.SubGridCellsMarker != null) && (!_DomainBnd) && (m_owner.SubGridCellsMarker[_jCell1] != m_owner.SubGridCellsMarker[_jCell2]);

                                if (_DomainBnd != DomainBnd)
                                    break;
                                if (_SubGrdBnd != SubGrdBnd)
                                    break;
                            }
#if DEBUG
                            for(int jEdgeTest = jEdge; jEdgeTest < jEdge + __Len; jEdgeTest++) {
                                int _jCell1 = grid.iGeomEdges.CellIndices[jEdgeTest, 0];
                                int _jCell2 = grid.iGeomEdges.CellIndices[jEdgeTest, 1];
                                bool _DomainBnd = _jCell2 < 0;
                                bool _SubGrdBnd = (m_owner.SubGridCellsMarker != null) && (!_DomainBnd) && (m_owner.SubGridCellsMarker[_jCell1] != m_owner.SubGridCellsMarker[_jCell2]);
                                Debug.Assert(_DomainBnd || _SubGrdBnd);
                                Debug.Assert(_DomainBnd == DomainBnd);
                                Debug.Assert(_SubGrdBnd == SubGrdBnd);
                            }

#endif
                            if (DomainBnd) {
                                // an edge on the domain boundary
                                // ++++++++++++++++++++++++++++++

                                //TempTimers.bndItems += __Len;
                                //TempTimers.bndWatch.Start();



                                // Vektorisierung für Rand-Flussfunktionen im Moment ungenutzt
                                CallBorder(nonlinFlx, jEdge, IndexOffset, __Len, jEdge, false, NoArgs, NoParams,
                                    components.MapArguments(m_FieldValuesIN, nonlinFlx, false),
                                    MapAlsoMean ? components.MapArguments(m_MeanFieldValuesIN, nonlinFlx, true) : null,
                                    MapAlsoGradient ? components.MapArguments(m_FieldGradientIN, nonlinFlx, true) : null);


                                //TempTimers.bndWatch.Stop();

                            } else {
                                // an internal edge, but on the boundary of the subgrid
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++

                                __Len = 1; // turn off vectorization

                                Debug.Assert(jCell1 >= 0 && jCell2 >= 0); // must be an internal edge !
                                Debug.Assert(SubGrdBnd);  // must be a subgrid boundary
                                Debug.Assert(m_owner.SubGridCellsMarker[jCell1] != m_owner.SubGridCellsMarker[jCell2]);

                                bool Cell1In = m_owner.SubGridCellsMarker[jCell1];

                                switch (m_owner.SubGridBoundaryTreatment) {
                                    case SubGridBoundaryModes.BoundaryEdge:
                                        MultidimensionalArray[] FieldVals, FieldValsMean, FieldGrad;
                                        if (Cell1In) {
                                            FieldVals = components.MapArguments(m_FieldValuesIN, nonlinFlx, false);
                                            FieldValsMean = MapAlsoMean ? components.MapArguments(m_MeanFieldValuesIN, nonlinFlx, true) : null;
                                            FieldGrad = MapAlsoGradient ? components.MapArguments(m_FieldGradientIN, nonlinFlx, true) : null;
                                        } else {
                                            FieldVals = components.MapArguments(m_FieldValuesOT, nonlinFlx, false);
                                            FieldValsMean = MapAlsoMean ? components.MapArguments(m_MeanFieldValuesOT, nonlinFlx, true) : null;
                                            FieldGrad = MapAlsoGradient ? components.MapArguments(m_FieldGradientOT, nonlinFlx, true) : null;
                                        }
                                        CallBorder(nonlinFlx, jEdge, IndexOffset, __Len, jEdge, !Cell1In, NoArgs, NoParams,
                                            FieldVals,
                                            FieldValsMean,
                                            FieldGrad);
                                        break;

                                    case SubGridBoundaryModes.InnerEdge:
                                    case SubGridBoundaryModes.InnerEdgeLTS:
                                        CallInner(nonlinFlx, jEdge, IndexOffset, __Len, NoArgs, NoParams,
                                            components.MapArguments(m_FieldValuesIN, nonlinFlx, false),
                                            components.MapArguments(m_FieldValuesOT, nonlinFlx, false),
                                            MapAlsoMean ? components.MapArguments(m_MeanFieldValuesIN, nonlinFlx, true) : null,
                                            MapAlsoMean ? components.MapArguments(m_MeanFieldValuesOT, nonlinFlx, true) : null,
                                            MapAlsoGradient ? components.MapArguments(m_FieldGradientIN, nonlinFlx, true) : null,
                                            MapAlsoGradient ? components.MapArguments(m_FieldGradientOT, nonlinFlx, true) : null);
                                        break;

                                    case SubGridBoundaryModes.OpenBoundary:
                                        if (Cell1In) {
                                            FieldVals = components.MapArguments(m_FieldValuesIN, nonlinFlx, false);
                                            FieldValsMean = MapAlsoMean ? components.MapArguments(m_MeanFieldValuesIN, nonlinFlx, true) : null;
                                            FieldGrad = MapAlsoGradient ? components.MapArguments(m_FieldGradientIN, nonlinFlx, true) : null;
                                        } else {
                                            FieldVals = components.MapArguments(m_FieldValuesOT, nonlinFlx, false);
                                            FieldValsMean = MapAlsoMean ? components.MapArguments(m_MeanFieldValuesOT, nonlinFlx, true) : null;
                                            FieldGrad = MapAlsoGradient ? components.MapArguments(m_FieldGradientOT, nonlinFlx, true) : null;
                                        }

                                        CallInner(nonlinFlx, jEdge, IndexOffset, __Len, NoArgs, NoParams,
                                            FieldVals,
                                            FieldVals,
                                            FieldValsMean,
                                            FieldValsMean,
                                            FieldGrad,
                                            FieldGrad);
                                        break;

                                    default:
                                        throw new NotImplementedException();
                                }
                            }


                            IndexOffset += (__Len - 1);
                        }
                        timers[iComp].Stop();

                        if (blck != null)
                            Monitor.Exit(blck);
                    }
                }
            }

        }


        protected void AllocateBuffers(int NoOfItems, MultidimensionalArray rule, int iThread, int NumThreads) {
            m_ThreadLocals[iThread].AllocateBuffers(NoOfItems, rule, iThread, NumThreads);
        }


        protected void OnCloneForThreadParallelization(IQuadrature q, int iThread, int NumThreads) {
            if (m_ThreadLocals == null) {
                m_ThreadLocals = new ThreadLocalsEdg[NumThreads];
            }

            m_ThreadLocals[iThread] = new ThreadLocalsEdg(iThread, this, q);
        }


        ThreadLocalsEdg[] m_ThreadLocals;


        void EvaluateEx(int i0, int Length, QuadRule qr, MultidimensionalArray QuadResult, int iThread, int NumThreads) {
            m_ThreadLocals[iThread].EvaluateEx(i0, Length, qr, QuadResult, iThread, NumThreads);
        }

        /// <summary>
        /// Stores the result of the quadrature.
        /// </summary>
        protected void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
            int NoOfFields = m_CodomainBasisS.Length;
            IGridData grid = GridDat;
            int nupdate = grid.iLogicalCells.NoOfLocalUpdatedCells;

            double alpha = m_alpha;

            BitArray cellMarker = this.SubGridCellsMarker;

            int[,] E2C = grid.iGeomEdges.LogicalCellIndices;
            byte[] tags = grid.iGeomEdges.EdgeTags;

            for (int jEdge = 0; jEdge < Length; jEdge++) {

                int jCell1, jCell2;
                jCell1 = E2C[jEdge + i0, 0];
                jCell2 = E2C[jEdge + i0, 1];

                //bool touchCell1 = true;
                bool touchCell2 = (jCell2 >= 0 && jCell2 < nupdate);
                bool touchCell1 = (cellMarker == null) ? true : cellMarker[jCell1];

                //bool touchCell2 = (jCell2 >= 0 && jCell2 < nupdate);
                if (cellMarker != null && touchCell2) {
                    touchCell2 = cellMarker[jCell2];
                }


                if (tags[jEdge + i0] >= Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG) {
                    if (_PeriodicVectorTrafo == PeriodicVectorTrafo.fwd) {
                        touchCell1 = false;
                    } else if (_PeriodicVectorTrafo == PeriodicVectorTrafo.bck) {
                        touchCell2 = false;
                    }
                }

                // Only active in case of Local timestepping:
                // We save more edgeFluxes than necessary, e.g. across possible MPI-borders.
                // We filter the relevant edgeFluxes later! Here not possible, because we 
                // only have the cellMarker of the own subGrid 
                bool srgdBndEdge = false;
                int side = 0;
                if (m_outputBndEdge != null) {
                    if (touchCell1 && !touchCell2 && jCell2 >=0) {
                        srgdBndEdge = true;
                        side = 1;
                    } else if (!touchCell1 && jCell2 >=0) {
                        srgdBndEdge = true;
                        side = 0;
                    }
                }

                for (int f = 0; f < NoOfFields; f++) {
                    //Field fld = m_CodomainBasisS[f];
                    int f_offset = m_MyMap[f];
                    int mE = m_NoOfTestFunctions[f];

                    int i0in = 0;
                    if (touchCell1)
                        i0in = m_CodomainMapping.LocalUniqueCoordinateIndex(f, jCell1, 0);
                    int i0ot = 0;
                    if (touchCell2)
                        i0ot = m_CodomainMapping.LocalUniqueCoordinateIndex(f, jCell2, 0);

                    for (int m = 0; m < mE; m++) {
                        int idx = f_offset + m;

                        if (touchCell1) {
                            Debug.Assert(m_CodomainMapping.LocalUniqueCoordinateIndex(f, jCell1, m) == i0in);
                            m_Output[i0in] += ResultsOfIntegration[jEdge, idx, 0] * alpha;
                        }

                        if (touchCell2) {
                            Debug.Assert(m_CodomainMapping.LocalUniqueCoordinateIndex(f, jCell2, m) == i0ot);
                            m_Output[i0ot] += ResultsOfIntegration[jEdge, idx, 1] * alpha;
                        }

                        if (srgdBndEdge && m == 0) {
                            m_outputBndEdge[NoOfFields * (jEdge + i0) + f] += ResultsOfIntegration[jEdge, idx, side] * alpha;
                        }

                        i0in++;
                        i0ot++;
                    }
                }
            }
        }



        static MultidimensionalArray.MultiplyProgram mp_in_in_Ti = MultidimensionalArray.MultiplyProgram.Compile("in", "in", "T(i)", true);
        static MultidimensionalArray.MultiplyProgram mp_kn_k_kn = MultidimensionalArray.MultiplyProgram.Compile("kn", "k", "kn");
        static MultidimensionalArray.MultiplyProgram mp_knd_k_knd = MultidimensionalArray.MultiplyProgram.Compile("knd", "k", "knd");
        static MultidimensionalArray.MultiplyProgram mp_jke_jkd_Tjed = MultidimensionalArray.MultiplyProgram.Compile("jke", "jkd", "T(j)ed", true);
        static MultidimensionalArray.MultiplyProgram mp_jn_Tjkn_jk = MultidimensionalArray.MultiplyProgram.Compile("jn", "T(j)kn", "jk", true);
        static MultidimensionalArray.MultiplyProgram mp_jn_Tjknd_jkd = MultidimensionalArray.MultiplyProgram.Compile("jn", "T(j)knd", "jkd", true);
        static MultidimensionalArray.MultiplyProgram mp_jn_Tjmn_jm = MultidimensionalArray.MultiplyProgram.Compile("jn", "T(j)mn", "jm", true);

        /// <summary>
        /// <see cref="PeriodicVectorTrafo"/>
        /// </summary>
        public PeriodicVectorTrafo _PeriodicVectorTrafo = PeriodicVectorTrafo.nix;

        /// <summary>
        /// Feature to support non-parallel periodic boundaries, e.g. cake-pies of rotational symmetrical domains.
        /// These must be integrated in two passes, one for the IN and one for the OUT-cell of the respective periodic edges.
        /// In both passes, a different (inverse) transformation is applied onto vector and normal fields.
        /// </summary>
        public enum PeriodicVectorTrafo {

            /// <summary>
            /// Periodic vector transform turned off.
            /// Can also be used for parallel periodic boundaries which do **not** require vector transformations.
            /// </summary>
            nix = 0,

            /// <summary>
            /// writes to the OUT-cell, no contribution to IN-cell
            /// </summary>
            fwd = 1,

            /// <summary>
            /// writes to the IN-cell, no contribution to OUT-cell
            /// </summary>
            bck = 2
        }
    }
}
