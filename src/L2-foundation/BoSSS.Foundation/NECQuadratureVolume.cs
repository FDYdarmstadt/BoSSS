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
using ilPSP.Tracing;
using System.Diagnostics;
using System.Linq;
using System;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using System.Configuration;
using System.Threading;

namespace BoSSS.Foundation.Quadrature.NonLin {

    //public static class Arsch {
    //    static public bool ShutTheFuckUp = false;
    //}

    /// <summary>
    /// edge quadrature of nonlinear equation components
    /// </summary>
    internal class NECQuadratureVolume : NECQuadratureCommon {



        /// <summary>
        /// ctor.
        /// </summary>
        public NECQuadratureVolume(IGridData context,
                                   DifferentialOperator DiffOp,
                                   IEnumerable<DGField> _DomainFields,
                                   IEnumerable<DGField> _ParameterFields,
                                   UnsetteledCoordinateMapping CodomainMapping,
                                   ICompositeQuadRule<QuadRule> domNrule)
            : base(context, DiffOp, _DomainFields, _ParameterFields, CodomainMapping) {

            // ----
            // init 
            // ----
            var temp = new ThreadLocalsVol(0, this, null);
            base.IsNonEmpty = temp.IsNonEmpty;

            int Gamma = base.m_DomainFields.Length;

            Debug.Assert(m_DomainAndParamFields.Length >= Gamma);
            m_ValueRequired = new bool[m_DomainAndParamFields.Length];
            m_GradientRequired = new bool[Gamma];

            // base.m_DomainFields may also contain parameter fields:
            for (int i = Gamma; i < base.m_DomainAndParamFields.Length; i++) {
                m_ValueRequired[i] = true;
            }

            temp.m_NonlinFormV.DetermineReqFields(m_GradientRequired,
                comp => ((comp.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0));
            temp.m_NonlinFormGradV.DetermineReqFields(m_GradientRequired,
                comp => ((comp.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0));
            temp.m_NonlinFormV.DetermineReqFields(m_ValueRequired,
                comp => ((comp.VolTerms & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0));
            temp.m_NonlinFormGradV.DetermineReqFields(m_ValueRequired,
                comp => ((comp.VolTerms & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0));
            temp.m_NonlinSources.DetermineReqFields(m_ValueRequired, comp => true);
            temp.m_NonlinFluxes.DetermineReqFields(m_ValueRequired, comp => true);
            temp.m_NonlinFluxesEx.DetermineReqFields(m_ValueRequired, comp => true);

            // -------------------
            // basis
            // -------------------

            
            for (int i = 0; i < m_CodomainBasisS.Length; i++) {

                if (temp.m_NonlinFluxes[i].m_AllComponentsOfMyType.Length > 0
                    || temp.m_NonlinFluxesEx[i].m_AllComponentsOfMyType.Length > 0
                    || temp.m_NonlinFormGradV[i].m_AllComponentsOfMyType.Length > 0) {

                    Basis GradBasis = m_CodomainBasisS[i];
                    if (m_MaxCodBasis_Gradient == null || m_MaxCodBasis_Gradient.Degree < GradBasis.Degree)
                        m_MaxCodBasis_Gradient = GradBasis;
                }
            }

            for (int i = 0; i < m_CodomainBasisS.Length; i++) {
                if (temp.m_NonlinSources[i].m_AllComponentsOfMyType.Length > 0
                    || temp.m_NonlinFormV[i].m_AllComponentsOfMyType.Length > 0) {

                    Basis ValBasis = m_CodomainBasisS[i];
                    if (m_MaxCodBasis == null || m_MaxCodBasis.Degree < ValBasis.Degree)
                        m_MaxCodBasis = ValBasis;
                }
            }



            // -----------------
            // quadrature object
            // -----------------

            var _Quad = CellQuadrature.GetQuadrature2(new int[] { CodomainMapping.NoOfCoordinatesPerCell }, context, domNrule,
                this.EvaluateEx,
                this.SaveIntegrationResults,
                this.AllocateBuffers,
                _OnCloneForThreadParallelization: this.OnCloneForThreadParallelization);
            _Quad.ExecuteParallel = true;
            m_Quad = _Quad;
        }




        /// <summary>
        /// Maximal codomain/testfunction basis for which the gradient values are required.
        /// </summary>
        readonly Basis m_MaxCodBasis_Gradient;

        /// <summary>
        /// Maximal codomain/testfunction basis for which the values are required.
        /// </summary>
        readonly Basis m_MaxCodBasis;

        /// <summary>
        /// true, if the evaluation of the a domain variable is required. 
        /// index: correlates with domain variables
        /// </summary>
        readonly bool[] m_ValueRequired;

        /// <summary>
        /// true, if the evaluation of the gradient of a domain variable is required. 
        /// index: correlates with domain variables
        /// </summary>
        readonly bool[] m_GradientRequired;

        class ThreadLocalsVol : ThreadLocals {

            readonly new NECQuadratureVolume m_owner;


            public ThreadLocalsVol(int iThread, NECQuadratureVolume _owner, IQuadrature q) : base(iThread, _owner, q) {
                m_owner = _owner;

                // ------------------------
                // sort equation components
                // ------------------------

                m_NonlinSources = EquationComponentArgMapping<INonlinearSource>.GetArgMapping(m_owner.Operator, true);
                m_NonlinFormV = EquationComponentArgMapping<INonlinVolumeForm_V>.GetArgMapping(m_owner.Operator, true,
                    eq => ((eq.VolTerms & (TermActivationFlags.V |TermActivationFlags.UxV | TermActivationFlags.GradUxV)) != 0),
                    eq => (eq is IVolumeForm ? new NonlinVolumeFormVectorizer((IVolumeForm)eq) : null));
                m_NonlinFormGradV = EquationComponentArgMapping<INonlinVolumeForm_GradV>.GetArgMapping(m_owner.Operator, true,
                    eq => ((eq.VolTerms & (TermActivationFlags.UxGradV | TermActivationFlags.GradV | TermActivationFlags.GradUxGradV)) != 0),
                    eq => (eq is IVolumeForm ? new NonlinVolumeFormVectorizer((IVolumeForm)eq) : null));



                if (q == null)
                    return;

                int Gamma = m_owner.m_DomainFields.Length;

                // ---------
                // profiling
                // ---------

                var _CustomTimers = new Stopwatch[] { new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch(), new Stopwatch() };
                var _CustomTimers_Names = new string[] { "Flux-Eval", "Basis-Eval", "Field-Eval", "Loops", "ParametersAndNormals", "Flux-Trafo" };
                Flux_Eval = _CustomTimers[0];
                Flux_Trafo = _CustomTimers[5];
                Field_Eval = _CustomTimers[2];
                Basis_Eval = _CustomTimers[1];
                Loops = _CustomTimers[3];
                ParametersAndNormals = _CustomTimers[4];
                q.CustomTimers = _CustomTimers;
                q.CustomTimers_Names = _CustomTimers_Names;
                q.CustomTimers_RootPointer = new int[_CustomTimers_Names.Length];
                ArrayTools.SetAll(q.CustomTimers_RootPointer, -1);

                this.m_NonlinSources_watch = this.m_NonlinSources.InitStopWatches(0, q);
                this.m_NonlinFormV_watch = this.m_NonlinFormV.InitStopWatches(0, q);
                this.m_NonlinFormGradV_watch = this.m_NonlinFormGradV.InitStopWatches(0, q);
                base.m_NonlinFluxesWatches = base.m_NonlinFluxes.InitStopWatches(0, q);
                base.m_NonlinFluxesExWatches = base.m_NonlinFluxesEx.InitStopWatches(0, q);

                // ---------------------
                // alloc multidim arrays
                // ---------------------

                m_FluxValues = new MultidimensionalArray[m_owner.m_CodomainBasisS.Length];
                m_FluxValuesTrf = new MultidimensionalArray[m_owner.m_CodomainBasisS.Length];
                for (int i = 0; i < m_FluxValues.Length; i++) {

                    if (m_NonlinFluxes[i].m_AllComponentsOfMyType.Length > 0 || m_NonlinFluxesEx[i].m_AllComponentsOfMyType.Length > 0 || m_NonlinFormGradV[i].m_AllComponentsOfMyType.Length > 0) {
                        m_FluxValues[i] = new MultidimensionalArray(3);
                        m_FluxValuesTrf[i] = new MultidimensionalArray(3);


                    }
                }

                m_SourceValues = new MultidimensionalArray[m_owner.m_CodomainBasisS.Length];
                for (int i = 0; i < m_SourceValues.Length; i++) {
                    if (m_NonlinSources[i].m_AllComponentsOfMyType.Length > 0 || m_NonlinFormV[i].m_AllComponentsOfMyType.Length > 0) {
                        m_SourceValues[i] = new MultidimensionalArray(2);


                    }
                }


                m_FieldValues = new MultidimensionalArray[m_owner.m_DomainAndParamFields.Length];
                m_FieldGradients = new MultidimensionalArray[Gamma];
                for (int i = 0; i < m_owner.m_DomainAndParamFields.Length; i++) {
                    if (m_owner.m_ValueRequired[i])
                        m_FieldValues[i] = new MultidimensionalArray(2);
                    if (i < Gamma && m_owner.m_GradientRequired[i])
                        m_FieldGradients[i] = new MultidimensionalArray(3);
                }


                m_TestFuncWeighted = new MultidimensionalArray(2);
                m_TestFuncGradWeighted = new MultidimensionalArray(3);
            }

            Stopwatch Flux_Eval;
            Stopwatch Flux_Trafo;
            Stopwatch Field_Eval;
            Stopwatch Basis_Eval;
            Stopwatch Loops;
            Stopwatch ParametersAndNormals;
            Stopwatch[][] m_NonlinSources_watch;
            Stopwatch[][] m_NonlinFormV_watch;
            Stopwatch[][] m_NonlinFormGradV_watch;


            /// <summary>
            /// array index: equation index/codomain variable index
            /// </summary>
            internal EquationComponentArgMapping<INonlinearSource>[] m_NonlinSources;

            /// <summary>
            /// array index: equation index/codomain variable index
            /// </summary>
            internal EquationComponentArgMapping<INonlinVolumeForm_V>[] m_NonlinFormV;

            /// <summary>
            /// array index: equation index/codomain variable index
            /// </summary>
            internal EquationComponentArgMapping<INonlinVolumeForm_GradV>[] m_NonlinFormGradV;

            /// <summary>
            /// true, if this integrator is responsible for any component
            /// </summary>
            override public bool IsNonEmpty {
                get {
                    return
                        base.IsNonEmpty ||
                        m_NonlinSources.IsNonEmpty() ||
                        m_NonlinFormV.IsNonEmpty() ||
                        m_NonlinFormGradV.IsNonEmpty();
                }
            }


            /// <summary>
            /// values of fields in the domain (of the operator);
            /// </summary>
            /// <remarks>
            /// index: corresponds with <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>
            /// </remarks>
            MultidimensionalArray[] m_FieldValues;


            /// <summary>
            /// values of fields in the domain (of the operator);
            /// </summary>
            /// <remarks>
            /// index: corresponds with <see cref="NECQuadratureCommon.m_DomainAndParamFields"/>
            /// </remarks>
            MultidimensionalArray[] m_FieldGradients;


            /// <summary>
            /// values of Riemann flux functions
            /// </summary>
            /// <remarks>
            /// index: codomain variable index;
            /// for each entry:
            /// <list type="bullet">
            /// <item>1st index: local cell index, with some offset;</item>
            /// <item>2nd index: node index;</item>
            /// <item>3rd index: spatial dimension index;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FluxValues;

            /// <summary>
            /// values of Riemann flux functions, transformed to reference coordinate system.
            /// </summary>
            /// <remarks>
            /// index: codomain variable index;
            /// for each entry:
            /// <list type="bullet">
            /// <item>1st index: local cell index, with some offset;</item>
            /// <item>2nd index: node index;</item>
            /// <item>3rd index: spatial dimension index;</item>
            /// </list>
            /// </remarks>
            MultidimensionalArray[] m_FluxValuesTrf;

            /// <summary>
            /// index: equation;
            /// an entry is null, if the corresponding equation contains no sources 
            /// (<see cref="m_NonlinSources"/>);
            /// for each entry:
            /// 1st index: local edge index, with some offset;
            /// 2nd index: node index;
            /// </summary>
            MultidimensionalArray[] m_SourceValues;

            public void AllocateBuffers(int NoOfItems, MultidimensionalArray rule, int iThread, int NumThreads) {

                int D = this.m_owner.GridDat.SpatialDimension;
                int NoOfNodes = rule.GetLength(0);

                // ----------------------
                // array for field values
                // ----------------------
                Debug.Assert(m_owner.m_DomainAndParamFields.Length == m_FieldValues.Length);
                Debug.Assert(m_owner.m_DomainAndParamFields.Length >= m_FieldGradients.Length);
                for (int f = 0; f < m_owner.m_DomainAndParamFields.Length; f++) {
                    if (m_FieldValues[f] != null)
                        m_FieldValues[f].Allocate(NoOfItems, NoOfNodes);
                    if (f < m_FieldGradients.Length && m_FieldGradients[f] != null)
                        m_FieldGradients[f].Allocate(NoOfItems, NoOfNodes, D);
                }

                // ---------------------
                // array for flux values
                // ---------------------
                for (int i = 0; i < m_FluxValues.Length; i++) {
                    Debug.Assert((m_FluxValues[i] != null) == (m_FluxValuesTrf[i] != null));
                    if (m_FluxValues[i] != null) {
                        m_FluxValues[i].Allocate(NoOfItems, NoOfNodes, D);
                        m_FluxValuesTrf[i].Allocate(NoOfItems, NoOfNodes, D);
                    }
                }

                // -----------------------
                // array for source values
                // -----------------------
                for (int i = 0; i < m_SourceValues.Length; i++) {
                    if (m_SourceValues[i] != null)
                        m_SourceValues[i].Allocate(NoOfItems, NoOfNodes);
                }
            }




            /// <summary>
            /// 
            /// </summary>
            public void EvaluateEx(int i0, int Length, QuadRule QR, MultidimensionalArray QuadResult, int iFred, int NoOfFreds) {
                NodeSet NodesUntransformed = QR.Nodes;
                IGridData grid = m_owner.GridDat;
                int D = grid.SpatialDimension;
                int NoOfEquations = m_owner.m_CodomainBasisS.Length;
                int NoOfNodes = NodesUntransformed.NoOfNodes;
                bool isAffine = grid.iGeomCells.IsCellAffineLinear(i0);
                int[] geom2log = grid.iGeomCells.GeomCell2LogicalCell;
#if DEBUG
                for (int i = 1; i < Length; i++) {
                    Debug.Assert(grid.iGeomCells.IsCellAffineLinear(i + i0) == isAffine);

                    if (geom2log == null) {
                        for (int e = 0; e < m_owner.m_CodomainBasisS.Length; e++) {
                            Debug.Assert(m_owner.m_CodomainBasisS[e].GetLength(i + i0) == m_owner.m_CodomainBasisS[e].GetLength(i0));
                        }
                    } else {
                        for (int e = 0; e < m_owner.m_CodomainBasisS.Length; e++) {
                            Debug.Assert(m_owner.m_CodomainBasisS[e].GetLength(geom2log[i + i0]) == m_owner.m_CodomainBasisS[e].GetLength(geom2log[i0]));
                        }
                    }
                }
#endif

                // this is an EvaluateEx: we are also responsible for multiplying with quadrature weights and summing up!
                Debug.Assert(QuadResult.Dimension == 2);
                Debug.Assert(QuadResult.GetLength(0) == Length);
                Debug.Assert(QuadResult.GetLength(1) == m_owner.m_CodomainMapping.BasisS.Sum(basis => basis.Length));



                // ===================
                // Evaluate all fields
                // ===================

                this.Field_Eval.Start();
                {
                    for (int f = 0; f < m_owner.m_DomainAndParamFields.Length; f++) {
                        if (m_owner.m_ValueRequired[f]) {
                            Debug.Assert(m_FieldValues[f] != null);
                            if (m_owner.m_DomainAndParamFields[f] != null) {
                                m_owner.m_DomainAndParamFields[f].Evaluate(i0, Length, NodesUntransformed, m_FieldValues[f]);
                            } else {
                                // field is null => set to 0.0
                                m_FieldValues[f].Clear();
                            }
                        }
                    }

                    for (int f = 0; f < m_owner.m_GradientRequired.Length; f++) {
                        if (m_owner.m_GradientRequired[f]) {
                            Debug.Assert(m_FieldGradients[f] != null);
                            if (m_owner.m_DomainAndParamFields[f] != null) {
                                m_owner.m_DomainAndParamFields[f].EvaluateGradient(i0, Length, NodesUntransformed, m_FieldGradients[f]);
                            } else {
                                // field is null => set to 0.0
                                m_FieldValues[f].Clear();
                            }
                        }
                    }
                }
                this.Field_Eval.Stop();



                // =====================================
                // Transform Nodes to global coordinates
                // =====================================

                this.ParametersAndNormals.Start();
                MultidimensionalArray NodesGlobalCoords;
                {
                    NodesGlobalCoords = grid.GlobalNodes.GetValue_Cell(NodesUntransformed, i0, Length);
                }
                this.ParametersAndNormals.Stop();

                // =======================
                // Evaluate Flux functions
                // =======================


                bool[] RequireTestFunctionGradient = new bool[NoOfEquations];
                bool[] Cleared_m_FluxValues = new bool[NoOfEquations];

                this.Flux_Eval.Start();
                {
                    // loop over all equations ...
                    for (int _e = 0; _e < NoOfEquations; _e++) {
                        int e = (_e + this.m_iThread) % NoOfEquations; // shuffling in threads to reduce locking
                        //                Field fld = m_CodomainFields[e];

                        if (m_NonlinFluxes[e].m_AllComponentsOfMyType.Length + m_NonlinFluxesEx[e].m_AllComponentsOfMyType.Length > 0) {
                            m_FluxValues[e].Clear();
                            Cleared_m_FluxValues[e] = true;

                            RequireTestFunctionGradient[e] = true;

                            // sum up all INonlinearFlux - objects 
                            int jjj = 0;
                            int L = m_NonlinFluxes[e].m_AllComponentsOfMyType.Length;
                            for (int _i = 0; _i < L; _i++) {
                                int i = (_i + m_iThread)%L;  // shuffling in threads to reduce locking
                                INonlinearFlux nonlinFlx = m_NonlinFluxes[e].m_AllComponentsOfMyType[i];

                                bool mustLock = m_NonlinFluxes[e].m_LockObjects[i] != null;
                                if (mustLock)
                                    Monitor.Enter(m_NonlinFluxes[e].m_LockObjects[i]);

                                m_NonlinFluxesWatches[e][jjj].Start();
                                nonlinFlx.Flux(m_owner.m_Time,
                                               NodesGlobalCoords,
                                               m_NonlinFluxes[e].MapArguments(m_FieldValues, nonlinFlx),
                                               0, Length,
                                               m_FluxValues[e]);
                                m_NonlinFluxesWatches[e][jjj].Stop();
                                jjj++;

                                if (mustLock)
                                    Monitor.Exit(m_NonlinFluxes[e].m_LockObjects[i]);
                            }

                            // sum up all INonlinearFluxEx - objects
                            jjj = 0;
                            int L2 = m_NonlinFluxesEx[e].m_AllComponentsOfMyType.Length;
                            for (int _i = 0; _i < L2; _i++) {
                                int i = (_i + m_iThread)%L2;  // shuffling in threads to reduce locking
                                INonlinearFluxEx nonlinFlxEx = m_NonlinFluxesEx[e].m_AllComponentsOfMyType[i];
                                
                                bool mustLock = m_NonlinFluxesEx[e].m_LockObjects[i] != null;
                                if (mustLock)
                                    Monitor.Enter(m_NonlinFluxesEx[e].m_LockObjects[i]);

                                m_NonlinFluxesExWatches[e][jjj].Start();
                                nonlinFlxEx.Flux(m_owner.m_Time,
                                                 NodesGlobalCoords,
                                                 m_NonlinFluxesEx[e].MapArguments(m_FieldValues, nonlinFlxEx),
                                                 0, Length,
                                                 m_FluxValues[e],
                                                 i0);
                                m_NonlinFluxesExWatches[e][jjj].Stop();
                                jjj++;

                                if (mustLock)
                                    Monitor.Exit(m_NonlinFluxesEx[e].m_LockObjects[i]);
                            }

                            m_FluxValues[e].Scale(-1.0);
                        }
                    }
                }


                // =========================
                // Evaluate Source functions
                // =========================

                bool[] RequireTestfunction = new bool[NoOfEquations];
                bool[] Cleared_m_SourceValues = new bool[NoOfEquations];
                {
                    for (int _e = 0; _e < NoOfEquations; _e++) {
                        int e = (_e + this.m_iThread) % NoOfEquations; // shuffling in threads to reduce locking
                        //                Equation eq = m_Equations[e];
                        //                Field fld = eq.MyField;


                        if (m_NonlinSources[e].m_AllComponentsOfMyType.Length > 0) {
                            m_SourceValues[e].Clear();
                            Cleared_m_SourceValues[e] = true;

                            RequireTestfunction[e] = true;

                            // sum up all sources
                            int jjj = 0;
                            int L = m_NonlinSources[e].m_AllComponentsOfMyType.Length;
                            for (int _i = 0; _i < L; _i++) {
                                int i = (_i + m_iThread)%L;  // shuffling in threads to reduce locking
                                INonlinearSource nonlinSrc = m_NonlinSources[e].m_AllComponentsOfMyType[i];
                                
                                
                                bool mustLock = m_NonlinSources[e].m_LockObjects[i] != null;
                                if (mustLock)
                                    Monitor.Enter(m_NonlinSources[e].m_LockObjects[i]);

                                m_NonlinSources_watch[e][jjj].Start();
                                nonlinSrc.Source(m_owner.m_Time,
                                                 NodesGlobalCoords,
                                                 m_NonlinSources[e].MapArguments(m_FieldValues, nonlinSrc),
                                                 0, i0, Length,
                                                 m_SourceValues[e]);
                                m_NonlinSources_watch[e][jjj].Stop();
                                jjj++;

                                if (mustLock)
                                    Monitor.Exit(m_NonlinSources[e].m_LockObjects[i]);
                            }
                        }
                    }
                }
                // ==============
                // Evaluate Forms
                // ==============
                {
                    for (int _e = 0; _e < NoOfEquations; _e++) {
                        int e = (_e + this.m_iThread) % NoOfEquations; // shuffling in threads to reduce locking

                        if (m_NonlinFormV[e].m_AllComponentsOfMyType.Length > 0) {

                            int L = m_NonlinFormV[e].m_AllComponentsOfMyType.Length;
                            for (int __icomp = 0; __icomp < L; __icomp++) {
                                int icomp = (__icomp + m_iThread) % L; // shuffling in threads to reduce locking
                                INonlinVolumeForm_V nonlinform = m_NonlinFormV[e].m_AllComponentsOfMyType[icomp];

                                if ((nonlinform.VolTerms & (TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradUxV)) == 0) {
                                    continue;
                                } else {
                                    m_NonlinFormV_watch[e][icomp].Start();
                                    RequireTestfunction[e] = true;

                                    if (!Cleared_m_SourceValues[e]) {
                                        m_SourceValues[e].Clear();
                                        Cleared_m_SourceValues[e] = true;
                                    }

                                    VolumFormParams vfp;
                                    vfp.GridDat = m_owner.GridDat;
                                    vfp.j0 = i0;
                                    vfp.Len = Length;
                                    vfp.Xglobal = NodesGlobalCoords;
                                    vfp.time = m_owner.m_Time;
                                    int NoArgs = m_NonlinFormV[e].NoOfArguments[icomp];
                                    int NoParams = m_NonlinFormV[e].NoOfParameters[icomp];
                                    var MappedArgsAndParams = m_NonlinFormV[e].MapArguments(this.m_FieldValues, nonlinform);
                                    vfp.ParameterVars = MappedArgsAndParams.GetSubVector(NoArgs, NoParams);
                                    var MappedArgs = MappedArgsAndParams.GetSubVector(0, NoArgs);

                                    var MappedGradients = m_NonlinFormV[e].MapArguments(this.m_FieldGradients, nonlinform, true);

                                    if (m_NonlinFormV[e].m_LockObjects[icomp] != null)
                                        Monitor.Enter(m_NonlinFormV[e].m_LockObjects[icomp]);

                                    nonlinform.Form(ref vfp, MappedArgs, MappedGradients, this.m_SourceValues[e]);

                                    if (m_NonlinFormV[e].m_LockObjects[icomp] != null)
                                        Monitor.Exit(m_NonlinFormV[e].m_LockObjects[icomp]);

                                    m_NonlinFormV_watch[e][icomp].Stop();
                                }
                            }
                        }
                    }

                    for (int _e = 0; _e < NoOfEquations; _e++) {
                        int e = (_e + this.m_iThread) % NoOfEquations; // shuffling in threads to reduce locking
                        if (m_NonlinFormGradV[e].m_AllComponentsOfMyType.Length > 0) {

                            int L = m_NonlinFormGradV[e].m_AllComponentsOfMyType.Length;
                            for (int __icomp = 0; __icomp < L; __icomp++) {
                                int icomp = (__icomp + m_iThread) % L; // shuffling in threads to reduce locking
                                INonlinVolumeForm_GradV nonlinform = m_NonlinFormGradV[e].m_AllComponentsOfMyType[icomp];

                                if ((nonlinform.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV | TermActivationFlags.GradV)) == 0) {
                                    continue;
                                } else {
                                    m_NonlinFormGradV_watch[e][icomp].Start();
                                    RequireTestFunctionGradient[e] = true;

                                    if (!Cleared_m_FluxValues[e]) {
                                        this.m_FluxValues[e].Clear();
                                        Cleared_m_FluxValues[e] = true;
                                    }

                                    VolumFormParams vfp;
                                    vfp.GridDat = m_owner.GridDat;
                                    vfp.j0 = i0;
                                    vfp.Len = Length;
                                    vfp.Xglobal = NodesGlobalCoords;
                                    vfp.time = m_owner.m_Time;
                                    int NoArgs = m_NonlinFormGradV[e].NoOfArguments[icomp];
                                    int NoParams = m_NonlinFormGradV[e].NoOfParameters[icomp];
                                    var MappedArgsAndParams = m_NonlinFormGradV[e].MapArguments(this.m_FieldValues, nonlinform);
                                    vfp.ParameterVars = MappedArgsAndParams.GetSubVector(NoArgs, NoParams);
                                    var MappedArgs = MappedArgsAndParams.GetSubVector(0, NoArgs);

                                    var MappedGradients = m_NonlinFormGradV[e].MapArguments(this.m_FieldGradients, nonlinform, true);

                                    if (m_NonlinFormGradV[e].m_LockObjects[icomp] != null)
                                        Monitor.Enter(m_NonlinFormGradV[e].m_LockObjects[icomp]);

                                    nonlinform.Form(ref vfp, MappedArgs, MappedGradients, this.m_FluxValues[e]);

                                    if (m_NonlinFormGradV[e].m_LockObjects[icomp] != null)
                                        Monitor.Exit(m_NonlinFormGradV[e].m_LockObjects[icomp]);

                                    m_NonlinFormGradV_watch[e][icomp].Stop();
                                }
                            }
                        }
                    }
                }
                this.Flux_Eval.Stop();

                // ================
                // Transform fluxes
                // ================

                // its less work to multiply the fluxes by the inverse Jacobi,
                // than each test function gradient by inverse Jacobi.
                // Could be interpreted as transforming fluxes to Refelem, i think....

                this.Flux_Trafo.Start();

                MultidimensionalArray InverseJacobi = null, JacobiDet = null;
                for (int e = 0; e < NoOfEquations; e++) {
                    Debug.Assert((m_FluxValues[e] != null) == (m_FluxValuesTrf[e] != null));

                    if (m_FluxValues[e] != null) {
                        if (InverseJacobi == null) {
                            if (isAffine) {
                                InverseJacobi = grid.iGeomCells.InverseTransformation.ExtractSubArrayShallow(new int[] { i0, 0, 0 }, new int[] { i0 + Length - 1, D - 1, D - 1 });
                            } else {
                                InverseJacobi = grid.InverseJacobian.GetValue_Cell(QR.Nodes, i0, Length);
                                //InverseJacobi = MultidimensionalArray.Create(Length, NoOfNodes, D, D);
                            }
                        }

                        if (JacobiDet == null && !isAffine)
                            JacobiDet = grid.JacobianDeterminat.GetValue_Cell(QR.Nodes, i0, Length);

                        if (isAffine) {
                            m_FluxValuesTrf[e].Multiply(1.0, m_FluxValues[e], InverseJacobi, 0.0, "jke", "jkd", "jed");
                            // for affine-linear cells the multiplication with Jacobi determinant is done AFTER quadrature, since it is constant per cell.
                        } else {
                            m_FluxValuesTrf[e].Multiply(1.0, m_FluxValues[e], InverseJacobi, 0.0, "jke", "jkd", "jked");
                            m_FluxValuesTrf[e].Multiply(1.0, m_FluxValuesTrf[e], JacobiDet, 0.0, "jke", "jke", "jk"); // apply scaling with Jacobi determinant, for integral transformation
                        }
                    }

                    if (m_SourceValues[e] != null) {

                        if (JacobiDet == null && !isAffine)
                            JacobiDet = grid.JacobianDeterminat.GetValue_Cell(QR.Nodes, i0, Length);
                        //JacobiDet = MultidimensionalArray.Create(Length, NoOfNodes);

                        // apply scaling with Jacobi determinant, for integral transformation
                        if (isAffine) {
                            // nop: for affine-linear cells the multiplication with Jacobi determinant is done AFTER quadrature, since it is constant per cell.
                        } else {
                            m_SourceValues[e].Multiply(1.0, m_SourceValues[e], JacobiDet, 0.0, "jk", "jk", "jk");
                        }
                    }
                }

                this.Flux_Trafo.Stop();



                // =======================
                // evaluate test functions
                // =======================

                this.Basis_Eval.Start();

                if (m_owner.m_MaxCodBasis != null && QR.NoOfNodes != this.m_TestFuncWeighted.GetLength(0)) {
                    
                    this.m_TestFuncWeighted.Allocate(QR.NoOfNodes, m_owner.m_MaxCodBasis.GetLength(geom2log != null ? geom2log[i0] : i0));
                }
                if (m_owner.m_MaxCodBasis_Gradient != null && QR.NoOfNodes != this.m_TestFuncGradWeighted.GetLength(0)) {
                    this.m_TestFuncGradWeighted.Allocate(QR.NoOfNodes, m_owner.m_MaxCodBasis_Gradient.GetLength(geom2log != null ? geom2log[i0] : i0), D);
                }

                if (m_owner.m_MaxCodBasis != null) {
                    var testFunc = m_owner.m_MaxCodBasis.Evaluate(QR.Nodes);
                    m_TestFuncWeighted.Multiply(1.0, QR.Weights, testFunc, 0.0, "kn", "k", "kn");
                }

                if (m_owner.m_MaxCodBasis_Gradient != null) {
                    var testFuncGrad = m_owner.m_MaxCodBasis_Gradient.EvaluateGradient(QR.Nodes);
                    m_TestFuncGradWeighted.Multiply(1.0, QR.Weights, testFuncGrad, 0.0, "knd", "k", "knd");
                }
                this.Basis_Eval.Stop();

                // ==========================================
                // multiply with test functions / save result
                // ==========================================

                MultidimensionalArray OrthoTrf = null; // to transform back to ONB on physical space...
                int iBufOrthoTrf;
                if (isAffine) {
                    OrthoTrf = TempBuffer.GetTempMultidimensionalarray(out iBufOrthoTrf, Length);
                    OrthoTrf.Multiply(1.0,
                        grid.iGeomCells.JacobiDet.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 }),
                        grid.ChefBasis.Scaling.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + Length - 1 }),
                        0.0,
                        "j", "j", "j");
                } else {
                    int MaxDegree = Math.Max(m_owner.m_MaxCodBasis != null ? m_owner.m_MaxCodBasis.Degree : 0, m_owner.m_MaxCodBasis_Gradient != null ? m_owner.m_MaxCodBasis_Gradient.Degree : 0);
                    OrthoTrf = grid.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(i0, Length, MaxDegree);
                    iBufOrthoTrf = int.MinValue;
                }

                this.Loops.Start();
                int N0, N;
                N0 = 0;
                for (int e = 0; e < NoOfEquations; e++) { // loop over equations...

                    if (geom2log != null)
                        N = m_owner.m_CodomainBasisS[e].GetLength(geom2log[i0]);
                    else
                        N = m_owner.m_CodomainBasisS[e].GetLength(i0);

                    int iBuf;
                    MultidimensionalArray QuadResult_e = TempBuffer.GetTempMultidimensionalarray(out iBuf, Length, N);

                    // fluxes
                    // ------
                    if (RequireTestFunctionGradient[e]) {
                        MultidimensionalArray
                            Fluxes_e = m_FluxValuesTrf[e],
                            testFuncGrad_e = null;
                        if (m_TestFuncGradWeighted.GetLength(1) == N)
                            testFuncGrad_e = m_TestFuncGradWeighted;
                        else
                            testFuncGrad_e = m_TestFuncGradWeighted.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { NoOfNodes - 1, N - 1, D - 1 });

                        QuadResult_e.Multiply(1.0, Fluxes_e, testFuncGrad_e, 0.0, "jn", "jke", "kne");
                    }

                    // sources
                    // -------

                    if (RequireTestfunction[e]) {
                        MultidimensionalArray
                            SourceValues_e = m_SourceValues[e],
                            testFunc_e = null;
                        if (m_TestFuncWeighted.GetLength(1) == N)
                            testFunc_e = m_TestFuncWeighted;
                        else
                            testFunc_e = m_TestFuncWeighted.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NoOfNodes - 1, N - 1 });


                        QuadResult_e.Multiply(1.0, SourceValues_e, testFunc_e, 1.0, "jn", "jk", "kn");

                    }

                    // final transformations
                    // ---------------------
                    MultidimensionalArray trfQuadResult_e = QuadResult.ExtractSubArrayShallow(new int[] { 0, N0 }, new int[] { Length - 1, N0 + N - 1 });

                    if (isAffine) {
                        trfQuadResult_e.Multiply(1.0, OrthoTrf, QuadResult_e, 0.0, "jn", "j", "jn");

                    } else {
                        MultidimensionalArray _OrthoTrf;
                        if (OrthoTrf.GetLength(1) == N)
                            _OrthoTrf = OrthoTrf;
                        else
                            _OrthoTrf = OrthoTrf.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, N - 1, N - 1 });

                        trfQuadResult_e.Multiply(1.0, _OrthoTrf, QuadResult_e, 0.0, "jn", "jmn", "jm");
                    }

                    // next
                    // ----

                    TempBuffer.FreeTempBuffer(iBuf);
                    N0 += N;
                }

                if (isAffine)
                    TempBuffer.FreeTempBuffer(iBufOrthoTrf);

                this.Loops.Stop();

#if DEBUG
                QuadResult.CheckForNanOrInf(true, true, true);
#endif
            }

            /// <summary>
            /// Buffer to store values of test functions multiplied by quadrature weights.
            /// </summary>
            MultidimensionalArray m_TestFuncWeighted;

            /// <summary>
            /// Buffer to store values of test function gradients multiplied by quadrature weights.
            /// </summary>
            MultidimensionalArray m_TestFuncGradWeighted;




        }

        protected void AllocateBuffers(int NoOfItems, MultidimensionalArray rule, int iThread, int NumThreads) {
            m_ThreadLocals[iThread].AllocateBuffers(NoOfItems, rule, iThread, NumThreads);
        }


        protected void OnCloneForThreadParallelization(IQuadrature q, int iThread, int NumThreads) {
            if (m_ThreadLocals == null) {
                m_ThreadLocals = new ThreadLocalsVol[NumThreads];
            }

            m_ThreadLocals[iThread] = new ThreadLocalsVol(iThread, this, q);
        }


        ThreadLocalsVol[] m_ThreadLocals;


        void EvaluateEx(int i0, int Length, QuadRule qr, MultidimensionalArray QuadResult, int iThread, int NumThreads) {
            m_ThreadLocals[iThread].EvaluateEx(i0, Length, qr, QuadResult, iThread, NumThreads);
        }


        /// <summary>
        /// 
        /// </summary>
        protected void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
            int NoOfFields = m_CodomainBasisS.Length;
            IGridData grid = this.GridDat;
            int[] geom2log = grid.iGeomCells.GeomCell2LogicalCell;


            double alpha = m_alpha;

            for (int j = 0; j < Length; j++) {
                int jCell;
                if (geom2log == null) {
                    // standard grid - geometrical and logical cells are 1-to-1
                    jCell = j + i0;
                } else {
                    // aggregate cell grid - multiple geometrical cells map to a logical cell
                    jCell = geom2log[j + i0];
                }


                for (int f = 0; f < NoOfFields; f++) {
                    int mE = m_NoOfTestFunctions[f];
                    int f_offset = m_MyMap[f];
                    int idx0 = m_CodomainMapping.LocalUniqueCoordinateIndex(f, jCell, 0);

                    for (int m = 0; m < mE; m++) {
                        int idx = f_offset + m;
                        Debug.Assert(idx0 + m == m_CodomainMapping.LocalUniqueCoordinateIndex(f, jCell, m));
                        m_Output[idx0 + m] += ResultsOfIntegration[j, idx] * alpha;
                    }
                }
            }
        }
    }
}
