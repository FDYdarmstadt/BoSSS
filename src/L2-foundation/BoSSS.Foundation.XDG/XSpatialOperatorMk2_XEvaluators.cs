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
using System.Linq;

using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;

using static BoSSS.Foundation.DifferentialOperator;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// An operator which is specialized in XDG Fields, i.e.
    /// it can have components which couple the phases.
    /// Mk2: enables the definition of different equation components for each phase 
    /// </summary>
    partial class XDifferentialOperatorMk2 {

        /// <summary>
        /// Assembly of matrices for linear (or linearized) XDG operators
        /// </summary>
        public class XEvaluatorLinear : XEvaluatorBase, IEvaluatorLinear {

            /// <summary>
            /// ctor
            /// </summary>
            internal XEvaluatorLinear(XDifferentialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory) :
                base(ownr, lsTrk, DomainVarMap, null, ParameterMap, CodomainVarMap, TrackerHistory) //
            {
                using (var tr = new FuncTrace()) {
                    base.MPITtransceive = true;


                }
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld4Species) {
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);


                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, this.m_lsTrk, spcName, quadOrder, edgeScheme, cellScheme, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var BulkMtxBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap);
                BulkMtxBuilder.MPITtransceive = false;
                SpeciesBulkMtxBuilder.Add(SpeciesId, BulkMtxBuilder);
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params4Spc, DGField[] DomFld4Species) {
                Debug.Assert(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                
                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.GhostEdgesOperator, this.m_lsTrk, spcName, quadOrder, ghostEdgeScheme, nullvolumeScheme, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var GhostEdgeBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params4Spc, CodomFrame.FrameMap);
                GhostEdgeBuilder.MPITtransceive = false;
                SpeciesGhostEdgeBuilder.Add(SpeciesId, GhostEdgeBuilder);
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                Debug.Assert(m_Xowner.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.SurfaceElementOperator_Ls0, this.m_lsTrk, spcName, quadOrder, SurfaceElement_Edge, SurfaceElement_volume, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var SurfElmBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap);
                SurfElmBuilder.MPITtransceive = false;
                SpeciesSurfElmBuilder.Add(SpeciesId, SurfElmBuilder);

            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorContactLineSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld4Species) {
                Debug.Assert(m_Xowner.ContactLineOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.ContactLineOperator_Ls0, this.m_lsTrk, spcName, quadOrder, eqs, cqs, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var ContactLineBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params_4Species, CodomFrame.FrameMap);
                ContactLineBuilder.MPITtransceive = false;
                SpeciesContactLineBuilder.Add(SpeciesId, ContactLineBuilder);
            }

            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesBulkMtxBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesGhostEdgeBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesSurfElmBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesContactLineBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();

            /// <summary>
            /// <see cref="IEvaluatorLinear.ComputeAffine{V}(V)"/>
            /// </summary>
            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                ComputeMatrix_Internal(default(BlockMsrMatrix), AffineOffset, true, 1.0);
            }

            /// <summary>
            /// <see cref="IEvaluatorLinear.ComputeMatrix{M,V}"/>
            /// </summary>
            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset, double alpha = 1.0)
                where M : IMutableMatrixEx
                where V : IList<double> {
                ComputeMatrix_Internal(Matrix, AffineOffset, false, alpha);
            }

            /// <summary>
            /// fields for MPI exchange
            /// </summary>
            protected override DGField[] GetTrxFields() {
                return base.Parameters.ToArray();
            }

            /// <summary>
            /// computation of operator matrix, currently only two species are supported
            /// </summary>
            void ComputeMatrix_Internal<M, V>(
                M Matrix, V AffineOffset, bool OnlyAffine, double alpha)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {
                if(base.MPITtransceive == true)
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                using(var tr = new FuncTrace()) {
                    var lsTrk = base.m_lsTrk;
                    IGridData GridDat = lsTrk.GridDat;

                    #region Check Input Arguments

                    if(OnlyAffine == false) {
                        if(!Matrix.RowPartitioning.Equals(base.CodomainMapping))
                            throw new ArgumentException("wrong number of columns in matrix.", "Matrix");
                        if(!Matrix.ColPartition.Equals(base.DomainMapping))
                            throw new ArgumentException("wrong number of rows in matrix.", "Matrix");
                    }

                    #endregion

                    #region MPI exchange of parameter fields
                    // --------------------------------
                    Transceiver trx = base.m_TRX;
                    if(trx != null) {
                        trx.TransceiveStartImReturn();
                    }
                    #endregion



                    // build matrix, bulk
                    // ---------------------
                    //MsrMatrix BulkMatrix = null;
                    //double[] BulkAffineOffset = null;
                    using(new BlockTrace("bulk_integration", tr)) {

                        // create the frame matrices & vectors...
                        // this is an MPI-collective operation, so it must be executed before the program may take different branches...
                        SpeciesFrameMatrix<M>[] mtx_spc = new SpeciesFrameMatrix<M>[ReqSpecies.Length];
                        SpeciesFrameVector<V>[] vec_spc = new SpeciesFrameVector<V>[ReqSpecies.Length];
                        for(int i = 0; i < ReqSpecies.Length; i++) {
                            SpeciesId SpId = ReqSpecies[i];
                            mtx_spc[i] = new SpeciesFrameMatrix<M>(Matrix, this.SpeciesCodomFrame[SpId], this.SpeciesDomainFrame[SpId]);
                            vec_spc[i] = (AffineOffset != null) ?
                                    (new SpeciesFrameVector<V>(AffineOffset, this.SpeciesCodomFrame[SpId]))
                                    :
                                    null;
                        }

                        // do the Bulk integration...
                        foreach(var SpeciesId in ReqSpecies) {
                            int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);



                            SpeciesFrameMatrix<M> mtx = mtx_spc[iSpecies];
                            var _mtx = Matrix != null ? mtx : default(SpeciesFrameMatrix<M>);

                            SpeciesFrameVector<V> vec = vec_spc[iSpecies];

                            var SpeciesBuilders = DoEdge ? new[] { SpeciesBulkMtxBuilder, SpeciesGhostEdgeBuilder, SpeciesSurfElmBuilder, SpeciesContactLineBuilder } : new[] { SpeciesBulkMtxBuilder };
                            foreach (var SpeciesBuilder in SpeciesBuilders) {

                                if(SpeciesBuilder.ContainsKey(SpeciesId)) {

                                    var builder = SpeciesBuilder[SpeciesId];
                                    //builder.OperatorCoefficients = this.SpeciesOperatorCoefficients[SpeciesId];
                                    NotifySpecies((DifferentialOperator)(builder.Owner), this.m_lsTrk, SpeciesId);

                                    if(trx != null) {
                                        trx.TransceiveFinish();
                                        trx = null;
                                    }

                                    builder.time = base.time;

                                    if(OnlyAffine) {
                                        if (alpha != 1.0)
                                            throw new NotSupportedException();
                                        builder.ComputeAffine(vec);
                                    } else {
                                        builder.ComputeMatrix(_mtx, vec, alpha);
                                    }
                                }

                            }
                        }

                    }

                    // build matrix, coupling
                    ///////////////////


                    using(new BlockTrace("surface_integration", tr)) {
                        if (DoEdge) {
#if DEBUG
                            {
                                // test if the 'coupling rules' are synchronous among MPI processes - otherwise, deadlock!
                                int[] crAllCount = CouplingRules.Count.MPIAllGather();
                                for(int rnk = 0; rnk < crAllCount.Length; rnk++) {
                                    Debug.Assert(crAllCount[rnk] == CouplingRules.Count);
                                }
                            }
#endif

                            var allBuilders = new List<LECQuadratureLevelSet<M, V>>();
                            foreach(var tt in this.CouplingRules) {
                                int iLevSet = tt.Item1;
                                var SpeciesA = tt.Item2;
                                var SpeciesB = tt.Item3;
                                var rule = tt.Item4;
#if DEBUG
                                int[] all_iLs = iLevSet.MPIAllGather();
                                int[] allSpcA = SpeciesA.cntnt.MPIAllGather();
                                int[] allSpcB = SpeciesB.cntnt.MPIAllGather();
                                for(int rnk = 0; rnk < all_iLs.Length; rnk++) {
                                    Debug.Assert(all_iLs[rnk] == iLevSet);
                                    Debug.Assert(allSpcA[rnk] == SpeciesA.cntnt);
                                    Debug.Assert(allSpcB[rnk] == SpeciesB.cntnt);
                                }

#endif
                                
                                var MtxBuilder = new LECQuadratureLevelSet<M, V>(GridDat,
                                                                 m_Xowner,
                                                                 OnlyAffine ? default(M) : Matrix, AffineOffset,
                                                                 CodomainMapping, Parameters, DomainMapping,
                                                                 lsTrk, iLevSet, TrackerHistoryIndex, new Tuple<SpeciesId, SpeciesId>(SpeciesA, SpeciesB),
                                                                 rule);
                                allBuilders.Add(MtxBuilder);
                                
                                if(trx != null) {
                                    trx.TransceiveFinish();
                                    trx = null; // we only need to do comm once!
                                }
                            }

                            // Note: this kind of deferred execution
                            // (first, collecting all integrators in a list and second, executing them in a separate loop)
                            // should prevent waiting for unevenly balanced level sets

                            foreach(var MtxBuilder in allBuilders) {
                                MtxBuilder.time = time;
                                UpdateLevelSetCoefficients(MtxBuilder.m_LevSetIdx, MtxBuilder.SpeciesA, MtxBuilder.SpeciesB);
                                MtxBuilder.Execute();

#if DEBUG
                                if(Matrix != null && OnlyAffine == false)
                                    Matrix.CheckForNanOrInfM("Matrix assembly, after surface integration: ");
                                if(AffineOffset != null)
                                    GenericBlas.CheckForNanOrInfV(AffineOffset, messageprefix: "Affine vector assembly, after surface integration: ");
#endif

                            }
                        }

                        // allow all processes to catch up
                        // -------------------------------
                        if(trx != null) {
                            trx.TransceiveFinish();
                            trx = null;
                        }

                    }
                }                
            }

            

            /// <summary>
            /// nix
            /// </summary>
            protected override void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule) {
            }
        }


        /// <summary>
        /// Explicit evaluation of (nonlinear and linear) XDG operators
        /// </summary>
        public class XEvaluatorNonlin : XEvaluatorBase, IEvaluatorNonLin {



            internal XEvaluatorNonlin(XDifferentialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                CoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory) :
                base(ownr, lsTrk, DomainVarMap, DomainVarMap.Fields, ParameterMap, CodomainVarMap, TrackerHistory) //
            {
                this.DomainFields = DomainVarMap;
                base.MPITtransceive = true;

            }

            /// <summary>
            /// Returns domain fields and parameters.
            /// </summary>
            protected override DGField[] GetTrxFields() {
                return ArrayTools.Cat(DomainFields.Fields, (base.Parameters != null) ? base.Parameters : new DGField[0]);
            }


            public CoordinateMapping DomainFields {
                get;
                private set;
            }

            public void Evaluate<Tout>(double alpha, double beta, Tout output, double[] outputBndEdge = null) where Tout : IList<double> {
                if (base.MPITtransceive == true)
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                using (var tr = new FuncTrace()) {
                    var lsTrk = base.m_lsTrk;
                    IGridData GridDat = lsTrk.GridDat;


                    if (outputBndEdge != null)
                        throw new NotImplementedException("todo");

                    if (beta != 1.0)
                        output.ScaleV(beta);


                    #region MPI exchange of parameter fields
                    // --------------------------------
                    Transceiver trx = base.m_TRX;
                    if (trx != null) {
                        trx.TransceiveStartImReturn();
                    }
                    #endregion


                    // bulk
                    // ---------------------
                    using (new BlockTrace("bulk_integration", tr)) {

                        // create the frame matrices & vectors...
                        // this is an MPI-collective operation, so it must be executed before the program may take different branches...
                        SpeciesFrameVector<Tout>[] vec_spc = new SpeciesFrameVector<Tout>[ReqSpecies.Length];
                        for (int i = 0; i < ReqSpecies.Length; i++) {
                            SpeciesId SpId = ReqSpecies[i];
                            vec_spc[i] = (new SpeciesFrameVector<Tout>(output, this.SpeciesCodomFrame[SpId]));
                        }

                        // do the Bulk integration...
                        foreach (var SpeciesId in ReqSpecies) {
                            int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);
                            var vec = vec_spc[iSpecies];

                            var SpeciesBuilders = DoEdge ? new[] { SpeciesBulkEval, SpeciesGhostEval, SpeciesSurfElmEval, SpeciesContactLineEval } : new[] { SpeciesBulkEval };
                            foreach (var SpeciesEval in SpeciesBuilders) {

                                if (SpeciesEval.ContainsKey(SpeciesId)) {

                                    var eval = SpeciesEval[SpeciesId];
                                    //eval.OperatorCoefficients = this.SpeciesOperatorCoefficients[SpeciesId];
                                    NotifySpecies((DifferentialOperator)(eval.Owner), this.m_lsTrk, SpeciesId);

                                    if (trx != null) {
                                        trx.TransceiveFinish();
                                        trx = null;
                                    }

                                    eval.time = base.time;
                                    eval.Evaluate(alpha, 1.0, vec, null);
                                }

                            }
                        }

                    }

                    //  coupling
                    ///////////////////

                    using (new BlockTrace("surface_integration", tr)) {
                        if (DoEdge) {
                            // TODO: are the quadrature rules non-empty?
#if DEBUG
                            {
                                // test if the 'coupling rules' are synchronous among MPI processes - otherwise, deadlock!
                                int[] crAllCount = CouplingRules.Count.MPIAllGather();
                                for(int rnk = 0; rnk < crAllCount.Length; rnk++) {
                                    Debug.Assert(crAllCount[rnk] == CouplingRules.Count);
                                }
                            }
#endif
                            var necList = new List<NECQuadratureLevelSet<Tout>> ();

                            foreach(var tt in this.CouplingRules) {
                                int iLevSet = tt.Item1;
                                var SpeciesA = tt.Item2;
                                var SpeciesB = tt.Item3;
                                var rule = tt.Item4;
#if DEBUG
                                int[] all_iLs = iLevSet.MPIAllGather();
                                int[] allSpcA = SpeciesA.cntnt.MPIAllGather();
                                int[] allSpcB = SpeciesB.cntnt.MPIAllGather();
                                for(int rnk = 0; rnk < all_iLs.Length; rnk++) {
                                    Debug.Assert(all_iLs[rnk] == iLevSet);
                                    Debug.Assert(allSpcA[rnk] == SpeciesA.cntnt);
                                    Debug.Assert(allSpcB[rnk] == SpeciesB.cntnt);
                                }
#endif

                                
                                // constructor is a collective operation
                                var LsEval = new NECQuadratureLevelSet<Tout>(GridDat,
                                                                 m_Xowner,
                                                                 output,
                                                                 this.DomainFields.Fields, Parameters, base.CodomainMapping,
                                                                 lsTrk, iLevSet, TrackerHistoryIndex, new Tuple<SpeciesId, SpeciesId>(SpeciesA, SpeciesB),
                                                                 rule);
                                necList.Add(LsEval);
                                
                                if(trx != null) {
                                    trx.TransceiveFinish();
                                    trx = null;
                                }
                            }

                            // Note: this kind of deferred execution
                            // (first, collecting all integrators in a list and second, executing them in a separate loop)
                            // should prevent waiting for unevenly balanced level sets
                            foreach(var LsEval in necList) {
                                LsEval.time = time;
                                UpdateLevelSetCoefficients(LsEval.m_LevSetIdx, LsEval.SpeciesA, LsEval.SpeciesB);
                                LsEval.Execute();

#if DEBUG
                                GenericBlas.CheckForNanOrInfV(output);
#endif
                            }

                        }
                    }


                    // allow all processes to catch up
                    // -------------------------------
                    if (trx != null) {
                        trx.TransceiveFinish();
                        trx = null;
                    }


                }
            }


            /// <summary>
            /// nix
            /// </summary>
            protected override void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule) {

            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld4Species) {
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = base.m_Xowner.FilterSpeciesOperator(base.m_Xowner, m_lsTrk, spcName, quadOrder, edgeScheme, cellScheme, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var BulkEval = tempOp.GetEvaluatorEx(DomFld4Species, Params_4Species, CodomFrame.FrameMap);

                BulkEval.MPITtransceive = false;
                SpeciesBulkEval.Add(SpeciesId, BulkEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld) {
                Debug.Assert(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.GhostEdgesOperator, m_lsTrk, spcName, quadOrder, ghostEdgeScheme, nullvolumeScheme, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var GhostEdgeEval = tempOp.GetEvaluatorEx(DomFld, Params_4Species, CodomFrame.FrameMap);
                GhostEdgeEval.MPITtransceive = false;
                SpeciesGhostEval.Add(SpeciesId, GhostEdgeEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld) {
                Debug.Assert(m_Xowner.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.SurfaceElementOperator_Ls0, m_lsTrk, spcName, quadOrder, SurfaceElement_Edge, SurfaceElement_volume, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var SurfElmEval = tempOp.GetEvaluatorEx(DomFld, Params_4Species, CodomFrame.FrameMap);
                SurfElmEval.MPITtransceive = false;
                SpeciesSurfElmEval.Add(SpeciesId, SurfElmEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorContactLineSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld) {
                Debug.Assert(m_Xowner.ContactLineOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.ContactLineOperator_Ls0, m_lsTrk, spcName, quadOrder, SurfaceElement_Edge, SurfaceElement_volume, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var ContactLineEval = tempOp.GetEvaluatorEx(DomFld, Params_4Species, CodomFrame.FrameMap);
                ContactLineEval.MPITtransceive = false;
                SpeciesContactLineEval.Add(SpeciesId, ContactLineEval);
            }

            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesBulkEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesGhostEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesSurfElmEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesContactLineEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();

        }


        abstract public class XEvaluatorBase : IEvaluator {

            /// <summary>
            /// equal to <see cref="Owner"/>
            /// </summary>
            internal XDifferentialOperatorMk2 m_Owner;

            /// <summary>
            /// the operator used to construct this object
            /// </summary>
            public IDifferentialOperator Owner {
                get {
                    return m_Owner;
                }
            }


            Dictionary<SpeciesId, MultidimensionalArray> m_CellLengthScales = new Dictionary<SpeciesId, MultidimensionalArray>();

            /// <summary>
            /// Dirty hack to provide length scales for the operator evaluation/matrix assembly, <see cref="CoefficientSet.CellLengthScales"/>
            /// </summary>
            public Dictionary<SpeciesId, MultidimensionalArray> CellLengthScales {
                get {
                    return m_CellLengthScales;
                }
            }

            Dictionary<SpeciesId, MultidimensionalArray> m_EdgeLengthScales = new Dictionary<SpeciesId, MultidimensionalArray>();

            /// <summary>
            /// Dirty hack to provide length scales for the operator evaluation/matrix assembly, <see cref="CoefficientSet.EdgeLengthScales"/>
            /// </summary>
            public Dictionary<SpeciesId, MultidimensionalArray> EdgeLengthScales {
                get {
                    return m_EdgeLengthScales;
                }
            }

            /// <summary>
            /// Write quadrature rules to text file, for debugging
            /// </summary>
            static private bool ruleDiagnosis = false;

            /// <summary>
            /// ctor
            /// </summary>
            protected internal XEvaluatorBase(
                XDifferentialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int __TrackerHistoryIndex) //
            {
                using (var tr = new FuncTrace()) {
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                    if (DomainVarMap.NoOfVariables != ownr.DomainVar.Count) {
                        throw new ArgumentException("wrong number of domain variables provided.");
                    }
                    this.m_Parameters = new DGField[ownr.ParameterVar.Count];
                    if (CodomainVarMap.NoOfVariables != ownr.CodomainVar.Count) {
                        throw new ArgumentException("wrong number of codomain variables provided.");
                    }

                    if (!object.ReferenceEquals(DomainVarMap.GridDat, CodomainVarMap.GridDat))
                        throw new ArgumentException("Domain and Codomain map must be assigned to the same grid");

                    foreach (var f in Parameters) {
                        if (f != null) {
                            if (!object.ReferenceEquals(DomainVarMap.GridDat, f.GridDat))
                                throw new ArgumentException("Parameter fields, domain and codomain basis must be assigned to the same grid");
                        }
                    }


                    m_Owner = ownr;
                    CodomainMapping = CodomainVarMap;
                    DomainMapping = DomainVarMap;
                    m_Parameters = (ParameterMap != null) ? ParameterMap.ToArray() : new DGField[0];


                    if (!m_Owner.IsCommitted)
                        throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");

                    int quadOrder = ownr.GetOrderFromQuadOrderFunction(DomainMapping.BasisS, GetBasisS(ParameterMap), CodomainVarMap.BasisS);

                 
                    if (!object.ReferenceEquals(GridData, lsTrk.GridDat))
                        throw new ArgumentException("grid data mismatch");
                    m_lsTrk = lsTrk;
                    m_Xowner = ownr;
                    ReqSpecies = ownr.Species.Select(spcNmn => lsTrk.GetSpeciesId(spcNmn)).ToArray();

                    this.UsedQuadOrder = quadOrder;
                    this.TrackerHistoryIndex = __TrackerHistoryIndex;

                    var SchemeHelper = lsTrk.GetXDGSpaceMetrics(ReqSpecies, quadOrder, __TrackerHistoryIndex).XQuadSchemeHelper;
                    var TrackerRegions = lsTrk.RegionsHistory[__TrackerHistoryIndex];

                    tr.Info("XSpatialOperator.ComputeMatrixEx quad order: " + quadOrder);

                    // compile quadrature rules & create matrix builders for each species 
                    // ------------------------------------------------------------------
                    foreach (var SpeciesId in ReqSpecies) {
                        //if(lsTrk.GetSpeciesName(SpeciesId) != "A") {
                        //    Console.WriteLine("REM: !!!!!!!!!!!!!!!!!!!!!!!1   species integration deact");
                        //    continue;
                        //}

                        int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);

                        // parameters for species
                        // ----------------------
                        DGField[] Params_4Species = (from f in (Parameters ?? new DGField[0])
                                            select ((f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f)).ToArray<DGField>();
                        SpeciesParams.Add(SpeciesId, Params_4Species);

                        DGField[] DomFld_4Species;
                        if (DomainFields != null) {
                            DomFld_4Species = DomainFields.Select(f => (f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f).ToArray();
                        } else {
                            DomFld_4Species = null;
                        }

                        // species frames
                        // --------------

                        var CodomFrame = new FrameBase(TrackerRegions, SpeciesId, CodomainMapping, false);
                        var DomainFrame = new FrameBase(TrackerRegions, SpeciesId, DomainMapping, true);
                        SpeciesCodomFrame.Add(SpeciesId, CodomFrame);
                        SpeciesDomainFrame.Add(SpeciesId, DomainFrame);

                        // quadrature rules
                        // ----------------

                        if (m_Xowner.TotalNoOfComponents > 0) {
                            EdgeQuadratureScheme edgeScheme = m_Xowner.EdgeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                            CellQuadratureScheme cellScheme = m_Xowner.VolumeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                            

                            if(ruleDiagnosis) {
                                var edgeRule = edgeScheme.Compile(this.GridData, quadOrder);
                                var volRule = cellScheme.Compile(this.GridData, quadOrder);
                                //edgeRule.(GridData, $"Edge{iLevSet}-{lsTrk.GetSpeciesName(SpeciesA)}{lsTrk.GetSpeciesName(SpeciesB)}.csv");
                                edgeRule.ToTextFileEdge(GridData, $"Edge-{lsTrk.GetSpeciesName(SpeciesId)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}.csv");
                                edgeRule.SumOfWeightsToTextFileEdge(this.GridData, $"Edge-{lsTrk.GetSpeciesName(SpeciesId)}-MPI{this.GridData.MpiRank}.csv");

                                volRule.ToTextFileCell(GridData, $"Volume-{lsTrk.GetSpeciesName(SpeciesId)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}.csv");
                                volRule.SumOfWeightsToTextFileVolume(GridData, $"Volume-{lsTrk.GetSpeciesName(SpeciesId)}-MPI{this.GridData.MpiRank}.csv");
                            }


                            ctorSpeciesIntegrator(SpeciesId, quadOrder, cellScheme, edgeScheme, DomainFrame, CodomFrame, Params_4Species, DomFld_4Species);
                        }

                        if (m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0) {
                            CellQuadratureScheme nullvolumeScheme = new CellQuadratureScheme(false, CellMask.GetEmptyMask(GridData));
                            EdgeQuadratureScheme ghostEdgeScheme = m_Xowner.GhostEdgeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                            ctorGhostSpeciesIntegrator(SpeciesId, quadOrder, nullvolumeScheme, ghostEdgeScheme, DomainFrame, CodomFrame, Params_4Species, DomFld_4Species);
                        }

                        //Only for ls0 so far:
                        //Add species, if it is separated from another species by level set 0
                        //For species not separated by ls0, nothing happens
                        var levelSetSpecies = lsTrk.GetSpeciesSeparatedByLevSet(0);
                        if(levelSetSpecies.Contains(lsTrk.GetSpeciesName(SpeciesId))) {
                            if (m_Xowner.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0) {
                                EdgeQuadratureScheme SurfaceElement_Edge = m_Xowner.SurfaceElement_EdgeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                CellQuadratureScheme SurfaceElement_volume = m_Xowner.SurfaceElement_VolumeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                if (ruleDiagnosis) {
                                    SurfaceElement_volume.ToTextFileCell(GridData, quadOrder, $"surfaceElementOperator_volume_{lsTrk.GetSpeciesName(SpeciesId)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}.txt");
                                    SurfaceElement_Edge.ToTextFileEdge(GridData, quadOrder, $"surfaceElementOperator_edge_{lsTrk.GetSpeciesName(SpeciesId)}-MPI{this.GridData.MpiRank}.txt");
                                    SurfaceElement_volume.Compile(GridData, 0).SumOfWeightsToTextFileVolume(GridData, $"surfaceElementOperator_volume_{lsTrk.GetSpeciesName(SpeciesId)}-MPI{this.GridData.MpiRank}.txt");
                                }
                                ctorSurfaceElementSpeciesIntegrator(SpeciesId, quadOrder, SurfaceElement_volume, SurfaceElement_Edge, DomainFrame, CodomFrame, Params_4Species, DomFld_4Species);
                            }
                            if (m_Xowner.ContactLineOperator_Ls0.TotalNoOfComponents > 0) {
                                EdgeQuadratureScheme ContactLine_Edge = new EdgeQuadratureScheme(false, EdgeMask.GetEmptyMask(GridData));
                                CellQuadratureScheme ContactLine_Volume = m_Xowner.ContactLine_VolumeQuadratureSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                if (ruleDiagnosis) {
                                    ContactLine_Volume.ToTextFileCell(GridData, quadOrder, $"contactLineOperator_{lsTrk.GetSpeciesName(SpeciesId)}-MPI{this.GridData.MpiRank}.csv");
                                }
                                ctorContactLineSpeciesIntegrator(SpeciesId, quadOrder, ContactLine_Volume, ContactLine_Edge, DomainFrame, CodomFrame, Params_4Species, DomFld_4Species);
                            }
                        }
                    }

                    // coupling terms
                    // --------------

                    using (new BlockTrace("surface_integration", tr)) {
                        if (m_Xowner.ContainesComponentType(typeof(ILevelSetForm))) {


                            var AllSpc = lsTrk.SpeciesIdS;

                            // loop over all possible pairs of species
                            for (int iSpcA = 0; iSpcA < AllSpc.Count; iSpcA++) {
                                var SpeciesA = AllSpc[iSpcA];
                                var SpeciesADom = TrackerRegions.GetSpeciesMask(SpeciesA);
                                //if (SpeciesADom.NoOfItemsLocally <= 0)
                                //    continue;

                                int _iSpcA = Array.IndexOf(ReqSpecies, SpeciesA);

                                for (int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                                    var SpeciesB = AllSpc[iSpcB];

                                    int _iSpcB = Array.IndexOf(ReqSpecies, SpeciesB);
                                    if (_iSpcA < 0 && _iSpcB < 0)
                                        continue;

                                    var SpeciesBDom = TrackerRegions.GetSpeciesMask(SpeciesB);
                                    var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);

                                    // Checks removed since they can cause parallel problems
                                    //if (SpeciesBDom.NoOfItemsLocally <= 0)
                                    //    continue;
                                    //if (SpeciesCommonDom.NoOfItemsLocally <= 0)
                                    //    continue;


                                    // loop over level-sets
                                    int NoOfLs = lsTrk.LevelSets.Count;
                                    for (int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                                        if (SchemeHelper.SpeciesAreSeparatedByLevSet(iLevSet, SpeciesA, SpeciesB)) {
                                            var LsDom = TrackerRegions.GetCutCellMask4LevSet(iLevSet);
                                            var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                            Chunk c = IntegrationDom.FirstOrDefault();
                                            if (c.Len > 0) {
                                                Debug.Assert(IntegrationDom.IsEmptyOnRank == false);
                                                int jtest = c.i0;

                                                LevelsetCellSignCode csc = lsTrk.RegionsHistory[__TrackerHistoryIndex].GetCellSignCode(jtest);

                                                if (!(csc.GetSign(iLevSet) == LevelsetSign.Both))
                                                    throw new ApplicationException("Seem to perform level-set integration in a non-cut cell.");

                                                var cscNeg = csc; cscNeg.SetSign(iLevSet, LevelsetSign.Negative);
                                                var cscPos = csc; cscPos.SetSign(iLevSet, LevelsetSign.Positive);


                                                bool SpeciesA_inNeg = lsTrk.ContainesSpecies(SpeciesA, cscNeg);
                                                bool SpeciesA_inPos = lsTrk.ContainesSpecies(SpeciesA, cscPos);
                                                bool SpeciesB_inNeg = lsTrk.ContainesSpecies(SpeciesB, cscNeg);
                                                bool SpeciesB_inPos = lsTrk.ContainesSpecies(SpeciesB, cscPos);

                                                if (SpeciesA_inPos == SpeciesA_inNeg) {
                                                    throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesA)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                }
                                                if (SpeciesB_inPos == SpeciesB_inNeg) {
                                                    throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesB)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                }


                                                if (SpeciesA_inNeg && SpeciesB_inPos) {
                                                    // nothing to do

                                                    if (SpeciesA_inPos == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesA)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }
                                                    if (SpeciesB_inNeg == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesB)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }
                                                } else if (SpeciesA_inPos && SpeciesB_inNeg) {
                                                    // flip species

                                                    if (SpeciesA_inNeg == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesA)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }
                                                    if (SpeciesB_inPos == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesB)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }

                                                    SpeciesId tmp = SpeciesA;
                                                    SpeciesA = SpeciesB;
                                                    SpeciesB = tmp;


                                                } else {
                                                    throw new ArgumentException($"Unable to determine negative (aka. In, A) and positive (aka. Out, B) species out of {m_lsTrk.GetSpeciesName(SpeciesA)}, {m_lsTrk.GetSpeciesName(SpeciesA)} w.r.t. Level-Set {iLevSet}.");
                                                }
                                            }
                                            ICompositeQuadRule<QuadRule> rule;
                                            using (new BlockTrace("QuadRule-compilation", tr)) {
                                                CellQuadratureScheme SurfIntegration = SchemeHelper.GetLevelSetquadScheme(iLevSet, SpeciesA,  IntegrationDom);
                                                rule = SurfIntegration.Compile(GridData, quadOrder);

                                                if (ruleDiagnosis) {
                                                    rule.ToTextFileCell(GridData, $"Levset{iLevSet}-{lsTrk.GetSpeciesName(SpeciesA)}{lsTrk.GetSpeciesName(SpeciesB)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}.csv");
                                                    rule.SumOfWeightsToTextFileVolume(GridData, $"Levset{iLevSet}-{lsTrk.GetSpeciesName(SpeciesA)}{lsTrk.GetSpeciesName(SpeciesB)}-MPI{this.GridData.MpiRank}.csv");
                                                }
                                            }

                                            LECQuadratureLevelSet<IMutableMatrix, double[]>.TestNegativeAndPositiveSpecies(rule, m_lsTrk, __TrackerHistoryIndex, SpeciesA, SpeciesB, iLevSet);

                                            CouplingRules.Add((iLevSet, SpeciesA, SpeciesB, rule));
                                            ctorLevSetFormIntegrator(iLevSet, SpeciesA, SpeciesB, rule);
                                        }
                                    }
                                }
                            }
#if DEBUG
                            int[] crAllCount = CouplingRules.Count.MPIAllGather();
                            for(int rnk = 0; rnk < crAllCount.Length; rnk++) {
                                Debug.Assert(crAllCount[rnk] == CouplingRules.Count);
                            }
#endif
                        }
                    }

                }
           }

            /// <summary>
            /// 
            /// </summary>
            protected int UsedQuadOrder;

            /// <summary>
            /// Index into <see cref="LevelSetTracker.RegionsHistory"/> and similar stacks.
            /// </summary>
            protected int TrackerHistoryIndex;


            /// <summary>
            /// create integrator for bulk phase
            /// </summary>
            abstract protected void ctorSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld4Species);

            /// <summary>
            /// Create integrator for <see cref="XDifferentialOperatorMk2.GhostEdgesOperator"/>
            /// </summary>
            abstract protected void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld4Species);

            /// <summary>
            /// Create integrator for <see cref="XDifferentialOperatorMk2.SurfaceElementOperator_Ls0"/>
            /// </summary>
            abstract protected void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld4Species);

            /// <summary>
            /// Create Integrator for <see cref="XDifferentialOperatorMk2.ContactLineOperator_Ls0"/>
            /// </summary>
            abstract protected void ctorContactLineSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params_4Species, DGField[] DomFld4Species);

            /// <summary>
            /// Create integrator for <see cref="ILevelSetForm"/> components
            /// </summary>
            abstract protected void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule);

            /// <summary>
            /// all species to integrate, defined through <see cref="XDifferentialOperatorMk2.Species"/>
            /// </summary>
            protected SpeciesId[] ReqSpecies;

            /// <summary>
            /// %
            /// </summary>
            protected LevelSetTracker m_lsTrk;

            /// <summary>
            /// the spatial operator
            /// </summary>
            protected XDifferentialOperatorMk2 m_Xowner;

            /// <summary>
            /// Parameter fields for each species, for <see cref="XDGField"/>s the species shadow, see <see cref="XDGField.GetSpeciesShadowField(SpeciesId)"/>
            /// </summary>
            protected Dictionary<SpeciesId, DGField[]> SpeciesParams = new Dictionary<SpeciesId, DGField[]>();

            /// <summary>
            /// for each species, the frame of the co-domain
            /// </summary>
            protected Dictionary<SpeciesId, FrameBase> SpeciesCodomFrame = new Dictionary<SpeciesId, FrameBase>();

            /// <summary>
            /// for each species, the frame of the domain
            /// </summary>
            protected Dictionary<SpeciesId, FrameBase> SpeciesDomainFrame = new Dictionary<SpeciesId, FrameBase>();

            /// <summary>
            /// quadrature rules for each level set/species pair;
            /// - 1st item: level-set index
            /// - 2nd item: negative species/species A
            /// - 3rd item positive species/species B
            /// - 4th item respective quadrature rule
            /// </summary>
            protected List<(int LsIdx, SpeciesId spcA, SpeciesId spcB, ICompositeQuadRule<QuadRule> quadRule)> CouplingRules = new List<(int, SpeciesId, SpeciesId, ICompositeQuadRule<QuadRule>)>();

            /// <summary>
            /// 
            /// </summary>
            private void GetCoefficients(SpeciesId spcA, SpeciesId spcB, out CoefficientSet csA, out CoefficientSet csB) {
                
                void FF(SpeciesId spc, out CoefficientSet cs) {
                    if(this.ReqSpecies.Contains(spc)) {
                        cs = m_Xowner.OperatorCoefficientsProvider(m_lsTrk, spc, UsedQuadOrder, TrackerHistoryIndex, time);

                        // hackedihack: i hate this shit
                        this.CellLengthScales.TryGetValue(spc, out var cls); cs.CellLengthScales = cls;
                        this.EdgeLengthScales.TryGetValue(spc, out var els); cs.EdgeLengthScales = els;
                    } else {
                        cs = null;
                    }
                }
                
                FF(spcA, out csA);
                FF(spcB, out csB);
            }

            /// <summary>
            /// calls all <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/> methods
            /// </summary>
            protected void UpdateLevelSetCoefficients(int iLevSet, SpeciesId spcA, SpeciesId spcB) {
                GetCoefficients(spcA, spcB, out var csA, out var csB);


                int[] DomDGdeg = this.DomainMapping.BasisS.Select(b => b.Degree).ToArray();
                int[] CodDGdeg = this.CodomainMapping.BasisS.Select(b => b.Degree).ToArray();
                string[] DomNames = Owner.DomainVar.ToArray();
                string[] CodNames = Owner.CodomainVar.ToArray();


                Debug.Assert(CodNames.Length == CodDGdeg.Length);
                for (int iCod = 0; iCod < CodDGdeg.Length; iCod++) {
                    var comps = Owner.EquationComponents[CodNames[iCod]];
                    foreach (var c in comps) {
                        if (c is ILevelSetFormSetup cs) {


                            cs.Setup(this.m_lsTrk);
                        }

                        if(c is ILevelSetForm lsc) {
                            // test if the component is actually relevant
                            if(lsc.PositiveSpecies != this.m_lsTrk.GetSpeciesName(spcB))
                                continue;
                            if(lsc.NegativeSpecies != this.m_lsTrk.GetSpeciesName(spcA))
                                continue;
                            if(lsc.LevelSetIndex != iLevSet)
                                continue;
                        }

                        if (c is ILevelSetEquationComponentCoefficient ce) {
                            
                            int[] DomDGdeg_cd = new int[ce.ArgumentOrdering.Count];
                            for (int i = 0; i < DomDGdeg_cd.Length; i++) {
                                string domName = ce.ArgumentOrdering[i];
                                int idx = Array.IndexOf(DomNames, domName);
                                DomDGdeg_cd[i] = DomDGdeg[idx];
                            }

                            ce.CoefficientUpdate(csA, csB, DomDGdeg_cd, CodDGdeg[iCod]);
                        }
                    }
                }
            }

            /// <summary>
            /// calls all <see cref="IEquationComponentSpeciesNotification.SetParameter"/> methods
            /// </summary>
            protected static void NotifySpecies(DifferentialOperator Owner, LevelSetTracker lsTrk, SpeciesId id) {
                string sNmn = lsTrk.GetSpeciesName(id);

                string[] CodNames = Owner.CodomainVar.ToArray();
                for (int iCod = 0; iCod < CodNames.Length; iCod++) {
                    var comps = Owner.EquationComponents[CodNames[iCod]];
                    foreach (var c in comps) {
                        if (c is IEquationComponentSpeciesNotification) {
                            var ce = c as IEquationComponentSpeciesNotification;
                            ce.SetParameter(sNmn);
                        }
                    }
                }
            }

          
            /// <summary>
            /// 
            /// </summary>
            protected CellMask m_SubGrid_InCells;

            SubGridBoundaryModes m_SubGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge;

            /// <summary>
            /// State set by <see cref="ActivateSubgridBoundary"/>
            /// </summary>
            public SubGridBoundaryModes SubGridBoundaryTreatment {
                get {
                    return m_SubGridBoundaryTreatment;
                }
            }

            /// <summary>
            /// Restricts the evaluation of the operator to a specific cell mask.
            /// </summary>
            /// <param name="Mask">
            /// cell mask where the operator should be evaluated
            /// </param>
            /// <param name="subGridBoundaryTreatment">
            /// defines what is to be done at edges where one neighboring cell is part of the cell mask <paramref name="Mask"/>, 
            /// but the other neighboring cell is *not* part of the cell mask <paramref name="Mask"/>.
            /// </param>
            public void ActivateSubgridBoundary(CellMask Mask, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                if (!object.ReferenceEquals(Mask.GridData, this.GridData))
                    throw new ArgumentException("grid mismatch");
                if (Mask != null && Mask.MaskType != MaskType.Logical)
                    throw new ArgumentException("expecting logical mask");
                m_SubGrid_InCells = Mask;
                m_SubGridBoundaryTreatment = subGridBoundaryTreatment;
            }


            /// <summary>
            /// coordinate mapping for the codomain variables;
            /// </summary>
            public UnsetteledCoordinateMapping CodomainMapping {
                get;
                private set;
            }



            /// <summary>
            /// <see cref="Parameters"/>
            /// </summary>
            DGField[] m_Parameters;

            /// <summary>
            /// parameter mapping
            /// </summary>
            public IList<DGField> Parameters {
                get {
                    return m_Parameters;
                }
            }

            /// <summary>
            /// coordinate mapping for the domain variables;
            /// </summary>
            public UnsetteledCoordinateMapping DomainMapping {
                get;
                private set;
            }

            /// <summary>
            /// Grid, on which this evaluator operates on.
            /// </summary>
            public IGridData GridData {
                get {
                    Debug.Assert(object.ReferenceEquals(DomainMapping.GridDat, CodomainMapping.GridDat));
                    return DomainMapping.GridDat;
                }
            }




            double m_time = 0.0;

            /// <summary>
            /// Time passed e.g. to <see cref="CommonParams.time"/>, <see cref="CommonParamsBnd.time"/> and <see cref="CommonParamsVol.time"/>.
            /// </summary>
            public double time {
                get {
                    return m_time;
                }
                set {
                    m_time = value;
                }
            }

            /// <summary>
            /// Should return all DG fields which should be exchanged, see <see cref="m_TRX"/>.
            /// </summary>
            abstract protected DGField[] GetTrxFields();

            bool m_MPITtransceive = false;

            /// <summary>
            /// Turn MPI sending/receiving of parameters and domain fields on/off.
            /// </summary>
            public bool MPITtransceive {
                get {
                    //return m_TRX != null;
                    var RealTrxFields = this.GetTrxFields().Where(f => f != null).ToArray();

                    //Debug.Assert((m_MPITtransceive == (m_TRX != null)) || (RealTrxFields.Length <= 0));

                    return m_MPITtransceive;
                }
                set {
                    m_MPITtransceive = value;

                    //if((m_TRX != null) && (value == true)) {
                    //    ArrayTools.ListEquals(m_TRX.)
                    //}
                    var RealTrxFields = this.GetTrxFields().Where(f => f != null).ToArray();

                    if ((value == true) && (RealTrxFields.Length > 0)) {
                        // + + + + + + + + + +
                        // create transceiver
                        // + + + + + + + + + +

                        if (m_TRX == null)
                            m_TRX = new Transceiver(RealTrxFields);
                    } else {
                        // + + + + + + + + + + + + +
                        // no communication required.
                        // + + + + + + + + + + + + +
                        m_TRX = null;
                    }
                }
            }

            /// <summary>
            /// Transceiver for the fields within <see cref="DomainMapping"/>
            /// </summary>
            protected Transceiver m_TRX;
        }
    }
}
