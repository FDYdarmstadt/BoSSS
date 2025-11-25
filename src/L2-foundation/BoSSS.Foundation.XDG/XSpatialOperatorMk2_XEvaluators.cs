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

using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using NUnit.Framework;

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
                using(var tr = new FuncTrace()) {
                    base.MPITtransceive = true;


                }
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params, DGField[] DomFld4Species) {
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);
    
                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, RowSwitch, this.m_lsTrk, spcName, quadOrder, edgeScheme, cellScheme, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var BulkMtxBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap);
                BulkMtxBuilder.MPITtransceive = false;
                SpeciesBulkMtxBuilder.Add(SpeciesId, BulkMtxBuilder);
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params4Spc, DGField[] DomFld4Species) {
                Debug.Assert(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);
                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.GhostEdgesOperator, RowSwitch, this.m_lsTrk, spcName, quadOrder, ghostEdgeScheme, nullvolumeScheme, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var GhostEdgeBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params4Spc, CodomFrame.FrameMap);
                GhostEdgeBuilder.MPITtransceive = false;
                SpeciesGhostEdgeBuilder.Add(SpeciesId, GhostEdgeBuilder);
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params, DGField[] DomFld) {
                Debug.Assert(m_Xowner.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.SurfaceElementOperator_Ls0, RowSwitch, this.m_lsTrk, spcName, quadOrder, SurfaceElement_Edge, SurfaceElement_volume, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var SurfElmBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap);
                SurfElmBuilder.MPITtransceive = false;
                SpeciesSurfElmBuilder.Add(SpeciesId, SurfElmBuilder);

            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorContactLineSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld4Species) {
                Debug.Assert(m_Xowner.ContactLineOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.ContactLineOperator_Ls0, RowSwitch, this.m_lsTrk, spcName, quadOrder, eqs, cqs, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var ContactLineBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params_4Species, CodomFrame.FrameMap);
                ContactLineBuilder.MPITtransceive = false;
                SpeciesContactLineBuilder.Add(SpeciesId, ContactLineBuilder);
            }



            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorTraceDGBulkIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params, DGField[] DomFld) {

                // Here so far, all the implementations of traceDg are based on single spc. (Theoretically, we don't need it)
                // But in order to use XDG ready-to-use
                // features, i need to specify the spc name. I cannot pass null to the spc name otherwise the following MatrixBuilder will cause problem.
                // Maybe i need to refactor it one day....
                // var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, RowSwitch, this.m_lsTrk, null, quadOrder, eqs, cqs, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, RowSwitch, this.m_lsTrk, spcName, quadOrder, eqs, cqs, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);


                TraceBulkMtxBuilder = tempOp.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap);
                TraceBulkMtxBuilder.MPITtransceive = false;
            }


            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesBulkMtxBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesGhostEdgeBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesSurfElmBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesContactLineBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();

            IEvaluatorLinear TraceBulkMtxBuilder;

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


                    void BulkIntegrator(IEvaluatorLinear builder, SpeciesId spId, SpeciesFrameMatrix<M> _mtx, SpeciesFrameVector<V> _vec) {
                        if(spId != default(SpeciesId))
                            NotifySpecies((DifferentialOperator)(builder.Owner), this.m_lsTrk, spId);

                        if(trx != null) {
                            trx.TransceiveFinish();
                            trx = null;
                        }

                        builder.time = base.time;

                        if(OnlyAffine) {
                            if(alpha != 1.0)
                                throw new NotSupportedException();
                            builder.ComputeAffine(_vec);
                        } else {
                            builder.ComputeMatrix(_mtx, _vec, alpha);
                        }
                    }

                    using(new BlockTrace("bulk_integration", tr)) {

                        // create the frame matrices & vectors...
                        // this is an MPI-collective operation, so it must be executed before the program may take different branches...
                        SpeciesFrameMatrix<M>[] mtx_spc_incTgd = new SpeciesFrameMatrix<M>[ReqSpecies.Length];
                        SpeciesFrameVector<V>[] vec_spc_incTdg = new SpeciesFrameVector<V>[ReqSpecies.Length];
                        SpeciesFrameMatrix<M>[] mtx_spc_noTgd = new SpeciesFrameMatrix<M>[ReqSpecies.Length];
                        SpeciesFrameVector<V>[] vec_spc_noTdg = new SpeciesFrameVector<V>[ReqSpecies.Length];
                        for(int i = 0; i < ReqSpecies.Length; i++) {
                            SpeciesId SpId = ReqSpecies[i];
                            mtx_spc_incTgd[i] = new SpeciesFrameMatrix<M>(Matrix, this.SpeciesCodomFrame_WithTraceDg[SpId], this.SpeciesDomainFrame_WithTraceDg[SpId]);
                            vec_spc_incTdg[i] = (AffineOffset != null) ?
                                                (new SpeciesFrameVector<V>(AffineOffset, this.SpeciesCodomFrame_WithTraceDg[SpId]))
                                                :
                                                null;
                            mtx_spc_noTgd[i] = new SpeciesFrameMatrix<M>(Matrix, this.SpeciesCodomFrame_WithoutTraceDg[SpId], this.SpeciesDomainFrame_WithoutTraceDg[SpId]);
                            vec_spc_noTdg[i] = (AffineOffset != null) ?
                                               (new SpeciesFrameVector<V>(AffineOffset, this.SpeciesCodomFrame_WithoutTraceDg[SpId]))
                                               :
                                               null;
                        }

                        // do the Bulk integration...
                        foreach(var SpeciesId in ReqSpecies) {
                            int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);

                            SpeciesFrameMatrix<M> mtx_incTdg = mtx_spc_incTgd[iSpecies];
                            var _mtx_incTdg = Matrix != null ? mtx_incTdg : default(SpeciesFrameMatrix<M>);
                            SpeciesFrameVector<V> vec_incTdg = vec_spc_incTdg[iSpecies];

                            SpeciesFrameMatrix<M> mtx_noTdg = mtx_spc_noTgd[iSpecies];
                            var _mtx_noTdg = Matrix != null ? mtx_noTdg : default(SpeciesFrameMatrix<M>);
                            SpeciesFrameVector<V> vec_noTdg = vec_spc_noTdg[iSpecies];

                            var SpeciesBuilders = new[] { (1, SpeciesBulkMtxBuilder, mtx_noTdg, vec_noTdg),
                                                          (2, SpeciesGhostEdgeBuilder, mtx_noTdg, vec_noTdg),
                                                          (3, SpeciesSurfElmBuilder, mtx_incTdg, vec_incTdg),
                                                          (4, SpeciesContactLineBuilder, mtx_incTdg, vec_incTdg) };

                            foreach(var t4 in SpeciesBuilders) {
                                var SpeciesBuilder = t4.Item2;
                                var _mtx = t4.Item3;
                                var vec = t4.Item4;

                                if(onlyfordebugging_DoBulk == false && (t4.Item1 == 1 || t4.Item1 == 2)) {
                                    Console.WriteLine("skipping bulk and ghost " + t4.Item1);
                                    continue;
                                }

                                if(onlyfordebugging_DoContactline == false && t4.Item1 == 4) {
                                    Console.WriteLine("skipping contact line " + t4.Item1);
                                    continue;
                                }
                                if(onlyfordebugging_DoSurfaceelm == false && t4.Item1 == 3) {
                                    Console.WriteLine("skipping surface element " + t4.Item1);
                                    continue;
                                }

                                if(SpeciesBuilder.ContainsKey(SpeciesId)) {
                                    var builder = SpeciesBuilder[SpeciesId];
                                    BulkIntegrator(builder, SpeciesId, _mtx, vec);
                                }
                            }
                        }
                    }

                    
                    if(base.AnyTraceDG_Codomn || base.AnyTraceDG_Domain) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // integrate equation components in the normal operator, which are applied to trace DG
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if(onlyfordebugging_DoTraceDG) {
                            using(new BlockTrace("bulk_integration_tracedg", tr)) {
                                var mtx = new SpeciesFrameMatrix<M>(Matrix, this.TraceDgCodomFrame, this.TraceDgDomainFrame);
                                var vec = (AffineOffset != null) ?
                                          (new SpeciesFrameVector<V>(AffineOffset, this.TraceDgCodomFrame))
                                          :
                                          null;

                                BulkIntegrator(TraceBulkMtxBuilder, default(SpeciesId), mtx, vec);

                            }
                        } else {
                            Console.WriteLine("skipping TraceDG");
                        }
                    }



                    // build matrix, coupling
                    ///////////////////



                    using(var bt = new BlockTrace("surface_integration", tr)) {
                        bt.IntermediateReportOfChildCalls = true;
                        if(onlyfordebugging_DoSurface) {
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

                                int[] NoOfItems = rule.Count().MPIAllGather();
                                int[] all_iLs = iLevSet.MPIAllGather();
                                int[] allSpcA = SpeciesA.cntnt.MPIAllGather();
                                int[] allSpcB = SpeciesB.cntnt.MPIAllGather();
                                for(int rnk = 0; rnk < all_iLs.Length; rnk++) {
                                    Assert.IsTrue(all_iLs[rnk] == iLevSet, "Level-set index mismatch.");

                                    if(NoOfItems.Min() > 0) {
                                        Assert.IsTrue(allSpcA[rnk] == SpeciesA.cntnt, "Mismatch of species designators across MPI ranks for species A");
                                        Assert.IsTrue(allSpcB[rnk] == SpeciesB.cntnt, "Mismatch of species designators across MPI ranks for species B");
                                    }
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
                                MtxBuilder.ExecuteParallel = true;
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
        }

        /// <summary>
        /// Only for debugging; can be used to turn surface integration in spatial operators off.
        /// </summary>
        static bool onlyfordebugging_DoSurface = true;

        /// <summary>
        /// Only for debugging; can be used to turn bulk component integration in spatial operators off.
        /// </summary>
        static bool onlyfordebugging_DoBulk = true;

        /// <summary>
        /// Only for debugging; can be used to turn contact line contribution (<see cref="ContactLineOperator_Ls0"/>) in spatial operators off.
        /// </summary>
        static bool onlyfordebugging_DoContactline = true;

        /// <summary>
        /// Only for debugging; can be used to turn surface element contribution (<see cref="SurfaceElementOperator_Ls0"/>) in spatial operators off.
        /// </summary>
        static bool onlyfordebugging_DoSurfaceelm = true;

        /// <summary>
        /// Only for debugging; can be used to turn trace DG integration in spatial operators off.
        /// </summary>
        static bool onlyfordebugging_DoTraceDG = true;


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
                if(base.MPITtransceive == true)
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                using(var tr = new FuncTrace()) {
                    var lsTrk = base.m_lsTrk;
                    IGridData GridDat = lsTrk.GridDat;


                    if(outputBndEdge != null)
                        throw new NotImplementedException("todo");

                    if(beta != 1.0)
                        output.ScaleV(beta);


                    #region MPI exchange of parameter fields
                    // --------------------------------
                    Transceiver trx = base.m_TRX;
                    if(trx != null) {
                        trx.TransceiveStartImReturn();
                    }
                    #endregion


                    // bulk
                    // ---------------------

                    void BulkIntegrator(IEvaluatorNonLin eval, SpeciesId spId, SpeciesFrameVector<Tout> vec) {
                        if(spId != default)
                            NotifySpecies((DifferentialOperator)(eval.Owner), this.m_lsTrk, spId);

                        if(trx != null) {
                            trx.TransceiveFinish();
                            trx = null;
                        }

                        eval.time = base.time;
                        //var bkup = output.ToArray();
                        eval.Evaluate(alpha, 1.0, vec, null);
                    }

                    using(new BlockTrace("bulk_integration", tr)) {

                        // create the frame matrices & vectors...
                        // this is an MPI-collective operation, so it must be executed before the program may take different branches...
                        SpeciesFrameVector<Tout>[] vec_spc_incTgd = new SpeciesFrameVector<Tout>[ReqSpecies.Length];
                        SpeciesFrameVector<Tout>[] vec_spc_noTgd = new SpeciesFrameVector<Tout>[ReqSpecies.Length];
                        for(int i = 0; i < ReqSpecies.Length; i++) {
                            SpeciesId SpId = ReqSpecies[i];
                            vec_spc_incTgd[i] = (new SpeciesFrameVector<Tout>(output, this.SpeciesCodomFrame_WithTraceDg[SpId]));
                            vec_spc_noTgd[i] = (new SpeciesFrameVector<Tout>(output, this.SpeciesCodomFrame_WithoutTraceDg[SpId]));
                        }

                        // do the Bulk integration...
                        foreach(var SpeciesId in ReqSpecies) {
                            int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);
                            var vec_incTdg = vec_spc_incTgd[iSpecies];
                            var vec_noTdg = vec_spc_noTgd[iSpecies];

                            var SpeciesEvals = new[] {
                                    (1, SpeciesBulkEval, vec_noTdg),
                                    (2, SpeciesGhostEval, vec_noTdg),
                                    (3, SpeciesSurfElmEval, vec_incTdg),
                                    (4, SpeciesContactLineEval, vec_incTdg) };

                            foreach(var t3 in SpeciesEvals) {
                                var SpeciesEval = t3.Item2;
                                var vec = t3.Item3;

                                if(onlyfordebugging_DoBulk == false && (t3.Item1 == 1 || t3.Item1 == 2)) {
                                    Console.WriteLine("skipping bulk and ghost " + t3.Item1);
                                    continue;
                                }

                                if(onlyfordebugging_DoContactline == false && t3.Item1 == 4) {
                                    Console.WriteLine("skipping contact line " + t3.Item1);
                                    continue;
                                }
                                if(onlyfordebugging_DoSurfaceelm == false && t3.Item1 == 3) {
                                    Console.WriteLine("skipping surface element " + t3.Item1);
                                    continue;
                                }



                                if(SpeciesEval.ContainsKey(SpeciesId)) {
                                    var eval = SpeciesEval[SpeciesId];
                                    BulkIntegrator(eval, SpeciesId, vec);
                                }
                            }
                        }
                    }

                    if(base.AnyTraceDG_Codomn || base.AnyTraceDG_Domain) {
                        if(onlyfordebugging_DoTraceDG) {

                            using(new BlockTrace("bulk_integration_tracedg", tr)) {
                                var vec = new SpeciesFrameVector<Tout>(output, this.TraceDgCodomFrame);

                                BulkIntegrator(TraceBulkEval, default, vec);
                            }
                        } else {
                            Console.WriteLine("skipping TraceDG");
                        }
                    }

                    //  coupling
                    ///////////////////

                    using(new BlockTrace("surface_integration", tr)) {
                        if(onlyfordebugging_DoSurface) {
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
                            var necList = new List<NECQuadratureLevelSet<Tout>>();

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
                                    //Debug.WriteLine("Rank {0} : LS {1} SpcA {2} SpcB {3}", rnk, all_iLs[rnk], allSpcA[rnk], allSpcB[rnk]);
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
                                string prefix = $"coupling-{LsEval.m_LevSetIdx}-{lsTrk.GetSpeciesName(LsEval.SpeciesA)}{lsTrk.GetSpeciesName(LsEval.SpeciesB)}";

                                LsEval.time = time;
                                LsEval.ExecuteParallel = true;
                                UpdateLevelSetCoefficients(LsEval.m_LevSetIdx, LsEval.SpeciesA, LsEval.SpeciesB);
                                LsEval.Execute();

#if DEBUG
                                GenericBlas.CheckForNanOrInfV(output);
#endif
                            }

                        } else {
                            Console.WriteLine("skipping level-set surface integration");
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

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld4Species) {
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = base.m_Xowner.FilterSpeciesOperator(base.m_Xowner, RowSwitch, m_lsTrk, spcName, quadOrder, edgeScheme, cellScheme, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var BulkEval = tempOp.GetEvaluatorEx(DomFld4Species, Params_4Species, CodomFrame.FrameMap);

                BulkEval.MPITtransceive = false;
                SpeciesBulkEval.Add(SpeciesId, BulkEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld) {
                Debug.Assert(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.GhostEdgesOperator, RowSwitch, m_lsTrk, spcName, quadOrder, ghostEdgeScheme, nullvolumeScheme, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var GhostEdgeEval = tempOp.GetEvaluatorEx(DomFld, Params_4Species, CodomFrame.FrameMap);
                GhostEdgeEval.MPITtransceive = false;
                SpeciesGhostEval.Add(SpeciesId, GhostEdgeEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld) {
                Debug.Assert(m_Xowner.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.SurfaceElementOperator_Ls0, RowSwitch, m_lsTrk, spcName, quadOrder, SurfaceElement_Edge, SurfaceElement_volume, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var SurfElmEval = tempOp.GetEvaluatorEx(DomFld, Params_4Species, CodomFrame.FrameMap);
                SurfElmEval.MPITtransceive = false;
                SpeciesSurfElmEval.Add(SpeciesId, SurfElmEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorContactLineSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld) {
                Debug.Assert(m_Xowner.ContactLineOperator_Ls0.TotalNoOfComponents > 0);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                CellLengthScales.TryGetValue(SpeciesId, out var cls);
                EdgeLengthScales.TryGetValue(SpeciesId, out var els);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner.ContactLineOperator_Ls0, RowSwitch, m_lsTrk, spcName, quadOrder, SurfaceElement_Edge, SurfaceElement_volume, TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                var ContactLineEval = tempOp.GetEvaluatorEx(DomFld, Params_4Species, CodomFrame.FrameMap);
                ContactLineEval.MPITtransceive = false;
                SpeciesContactLineEval.Add(SpeciesId, ContactLineEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorTraceDGBulkIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params, DGField[] DomFld) {
                //Here so far, all the implementations of traceDg are based on single spc. (Theoretically we don't need it) But in order to use XDG ready-to-use
                // features, i need to specify the spc name. I cannot pass null to the spc name otherwise the following MatrixBuilder will cause problem.
                //Maybe i need to refactor it one day....
                //var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, RowSwitch, this.m_lsTrk, null, quadOrder, eqs, cqs, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);
                string spcName = m_lsTrk.GetSpeciesName(SpeciesId);

                var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, RowSwitch, this.m_lsTrk, spcName, quadOrder, eqs, cqs, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);
                //var tempOp = m_Xowner.FilterSpeciesOperator(m_Xowner, RowSwitch, this.m_lsTrk, null, quadOrder, eqs, cqs, this.TrackerHistoryIndex, CellLengthScales, EdgeLengthScales);

                TraceBulkEval = tempOp.GetEvaluatorEx(DomFld, Params, CodomFrame.FrameMap);
                TraceBulkEval.MPITtransceive = false;
            }

            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesBulkEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesGhostEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesSurfElmEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesContactLineEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            IEvaluatorNonLin TraceBulkEval;
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
            static bool onlyfordebugging_RuleDiagnosis = false;

            /// <summary>
            /// ctor
            /// </summary>
            protected internal XEvaluatorBase(
                XDifferentialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int __TrackerHistoryIndex) //
            {
                using(var tr = new FuncTrace()) {
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                    if(DomainVarMap.NoOfVariables != ownr.DomainVar.Count) {
                        throw new ArgumentException("wrong number of domain variables provided.");
                    }
                    this.m_Parameters = new DGField[ownr.ParameterVar.Count];
                    if(CodomainVarMap.NoOfVariables != ownr.CodomainVar.Count) {
                        throw new ArgumentException("wrong number of codomain variables provided.");
                    }

                    if(!object.ReferenceEquals(DomainVarMap.GridDat, CodomainVarMap.GridDat))
                        throw new ArgumentException("Domain and Codomain map must be assigned to the same grid");

                    foreach(var f in Parameters) {
                        if(f != null) {
                            if(!object.ReferenceEquals(DomainVarMap.GridDat, f.GridDat))
                                throw new ArgumentException("Parameter fields, domain and codomain basis must be assigned to the same grid");
                        }
                    }


                    m_Owner = ownr;
                    CodomainMapping = CodomainVarMap;
                    DomainMapping = DomainVarMap;
                    m_Parameters = (ParameterMap != null) ? ParameterMap.ToArray() : new DGField[0];


                    if(!m_Owner.IsCommitted)
                        throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");

                    int quadOrder = ownr.GetOrderFromQuadOrderFunction(DomainMapping.BasisS, GetBasisS(ParameterMap), CodomainVarMap.BasisS);


                    if(!object.ReferenceEquals(GridData, lsTrk.GridDat))
                        throw new ArgumentException("grid data mismatch");
                    m_lsTrk = lsTrk;
                    m_Xowner = ownr;
                    ReqSpecies = ownr.Species.Select(lsTrk.GetSpeciesId).ToArray();

                    this.UsedQuadOrder = quadOrder;
                    this.TrackerHistoryIndex = __TrackerHistoryIndex;

                    var SchemeHelper = lsTrk.GetXDGSpaceMetrics(ReqSpecies, quadOrder, __TrackerHistoryIndex).XQuadSchemeHelper;
                    var TrackerRegions = lsTrk.RegionsHistory[__TrackerHistoryIndex];

                    tr.Info("XSpatialOperator.ComputeMatrixEx quad order: " + quadOrder);

                    // compile quadrature rules & create matrix builders for each species 
                    // ------------------------------------------------------------------
                    {
                        foreach(var SpeciesId in ReqSpecies) {
                            int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);

                            // parameters for species
                            // ----------------------
                            DGField[] Params_4Species = (from f in (Parameters ?? new DGField[0])
                                                         select ((f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f)).ToArray<DGField>();
                            SpeciesParams.Add(SpeciesId, Params_4Species);

                            DGField[] DomFld_4Species;
                            if(DomainFields != null) {
                                DomFld_4Species = DomainFields.Select(f => (f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f).ToArray();
                            } else {
                                DomFld_4Species = null;
                            }

                            // species frames
                            // --------------

                            var CodomFrame_WithoutTraceDg = new FrameBase(TrackerRegions, SpeciesId, CodomainMapping, false, FrameBase_TraceDGhandling.Without_TraceDG);
                            var DomainFrame_WithoutTraceDg = new FrameBase(TrackerRegions, SpeciesId, DomainMapping, true, FrameBase_TraceDGhandling.Without_TraceDG);
                            SpeciesCodomFrame_WithoutTraceDg.Add(SpeciesId, CodomFrame_WithoutTraceDg);
                            SpeciesDomainFrame_WithoutTraceDg.Add(SpeciesId, DomainFrame_WithoutTraceDg);

                            var CodomFrame_WithTraceDg = new FrameBase(TrackerRegions, SpeciesId, CodomainMapping, false, FrameBase_TraceDGhandling.DG_XDG_and_TraceDG);
                            var DomainFrame_WithTraceDg = new FrameBase(TrackerRegions, SpeciesId, DomainMapping, true, FrameBase_TraceDGhandling.DG_XDG_and_TraceDG);
                            SpeciesCodomFrame_WithTraceDg.Add(SpeciesId, CodomFrame_WithTraceDg);
                            SpeciesDomainFrame_WithTraceDg.Add(SpeciesId, DomainFrame_WithTraceDg);

                            // quadrature rules
                            // ----------------

                            if(m_Xowner.TotalNoOfComponents > 0) {
                                EdgeQuadratureScheme edgeScheme = m_Xowner.EdgeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                CellQuadratureScheme cellScheme = m_Xowner.VolumeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);


                                if(onlyfordebugging_RuleDiagnosis) {
                                    var edgeRule = edgeScheme.Compile(this.GridData, quadOrder);
                                    var volRule = cellScheme.Compile(this.GridData, quadOrder);

                                    string suffix = $"{lsTrk.GetSpeciesName(SpeciesId)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}of{this.GridData.MpiSize}";
                                    edgeRule.SaveToTextFileEdge(GridData, $"Edge-{suffix}.csv");
                                    edgeRule.ToVtpFilesEdge(GridData, $"Edge-{suffix}");
                                    edgeRule.SumOfWeightsToTextFileEdge(this.GridData, $"WgtSumEdge-{suffix}.csv");

                                    volRule.SaveToTextFileCell(GridData, $"Volume-{suffix}.csv");
                                    volRule.ToVtpFilesCell(GridData, $"Volume-{suffix}");
                                    volRule.SumOfWeightsToTextFileVolume(GridData, $"WgtSumVolume-{suffix}.csv");

                                }


                                ctorSpeciesIntegrator(SpeciesId, quadOrder, cellScheme, edgeScheme, DomainFrame_WithoutTraceDg, CodomFrame_WithoutTraceDg,
                                    CodomFrame_WithoutTraceDg.FrameMap.BasisS.Select(b => b.MaximalLength > 0).ToArray(),
                                    Params_4Species, DomFld_4Species);
                            }

                            if(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0) {
                                CellQuadratureScheme nullvolumeScheme = new CellQuadratureScheme(false, CellMask.GetEmptyMask(GridData));
                                EdgeQuadratureScheme ghostEdgeScheme = m_Xowner.GhostEdgeQuadraturSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                bool[] rowSwitch = DomainFrame_WithoutTraceDg.FrameMap.BasisS.Select(b => b.MaximalLength < 0).ToArray();
                                ctorGhostSpeciesIntegrator(SpeciesId, quadOrder, nullvolumeScheme, ghostEdgeScheme, DomainFrame_WithoutTraceDg, CodomFrame_WithoutTraceDg,
                                    CodomFrame_WithoutTraceDg.FrameMap.BasisS.Select(b => b.MaximalLength > 0).ToArray(),
                                    Params_4Species, DomFld_4Species);
                            }

                            // Only for ls0 so far:
                            // Add species, if it is separated from another species by level set 0
                            // For species not separated by ls0, nothing happens
                            var levelSetSpecies = lsTrk.GetSpeciesSeparatedByLevSet(0);
                            if(levelSetSpecies.Contains(lsTrk.GetSpeciesName(SpeciesId))) {
                                if(m_Xowner.SurfaceElementOperator_Ls0.TotalNoOfComponents > 0) {
                                    foreach(var speciesB in levelSetSpecies) {
                                        var speciesBId = lsTrk.GetSpeciesId(speciesB);
                                        if(SpeciesId == speciesBId)
                                            continue;

                                        EdgeQuadratureScheme SurfaceElement_Edge = m_Xowner.SurfaceElement_EdgeQuadraturSchemeProvider(lsTrk, SpeciesId, speciesBId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                        CellQuadratureScheme SurfaceElement_volume = m_Xowner.SurfaceElement_VolumeQuadraturSchemeProvider(lsTrk, SpeciesId, speciesBId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                        if(onlyfordebugging_RuleDiagnosis) {

                                            //if(GridData.MpiRank == 1)
                                            //    Debugger.Launch();
                                            var coEdgRule = SurfaceElement_Edge.Compile(GridData, quadOrder);
                                            var coVolRole = SurfaceElement_volume.Compile(GridData, quadOrder);

                                            string suffix = $"{lsTrk.GetSpeciesName(SpeciesId)}{speciesB}-ls{0}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}of{this.GridData.MpiSize}";
                                            coVolRole.SaveToTextFileCell(GridData, $"surfaceElementOperator_volume_{suffix}.csv");
                                            coEdgRule.SaveToTextFileEdge(GridData, $"surfaceElementOperator_edge_{suffix}.csv");

                                            coVolRole.ToVtpFilesCell(GridData, $"surfaceElementOperator_volume_{suffix}");
                                            coEdgRule.ToVtpFilesEdge(GridData, $"surfaceElementOperator_edge_{suffix}");

                                            SurfaceElement_volume.Compile(GridData, 0).SumOfWeightsToTextFileVolume(GridData, $"wgtSumSurfaceElementOperator_volume_{suffix}.csv");
                                        }

                                        ctorSurfaceElementSpeciesIntegrator(SpeciesId, quadOrder, SurfaceElement_volume, SurfaceElement_Edge, DomainFrame_WithTraceDg, CodomFrame_WithTraceDg,
                                            CodomFrame_WithTraceDg.FrameMap.BasisS.Select(b => b.MaximalLength > 0).ToArray(),
                                            Params_4Species, DomFld_4Species);
                                    }
                                }
                            }

                            if(m_Xowner.ContactLineOperator_Ls0.TotalNoOfComponents > 0) {
                                EdgeQuadratureScheme ContactLine_Edge = new EdgeQuadratureScheme(false, EdgeMask.GetEmptyMask(GridData));
                                CellQuadratureScheme ContactLine_Volume = m_Xowner.ContactLine_VolumeQuadratureSchemeProvider(lsTrk, SpeciesId, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                if(onlyfordebugging_RuleDiagnosis) {
                                    string suffix = $"{lsTrk.GetSpeciesName(SpeciesId)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}of{this.GridData.MpiSize}";

                                    ContactLine_Volume.SaveToTextFileCell(GridData, quadOrder, $"contactLineOperator_{suffix}.csv");
                                    ContactLine_Volume.Compile(GridData, quadOrder).ToVtpFilesCell(GridData, $"contactLineOperator_{suffix}");
                                }
                                ctorContactLineSpeciesIntegrator(SpeciesId, quadOrder, ContactLine_Volume, ContactLine_Edge, DomainFrame_WithTraceDg, CodomFrame_WithTraceDg,
                                    CodomFrame_WithTraceDg.FrameMap.BasisS.Select(b => b.MaximalLength > 0).ToArray(),
                                    Params_4Species, DomFld_4Species);
                            }
                        }
                    }


                    if(AnyTraceDG_Codomn || AnyTraceDG_Domain) {

                        TraceDgCodomFrame = new FrameBase(TrackerRegions, default(SpeciesId), CodomainMapping, false, FrameBase_TraceDGhandling.Without_DG_and_XDG);
                        TraceDgDomainFrame = new FrameBase(TrackerRegions, default(SpeciesId), DomainMapping, true, FrameBase_TraceDGhandling.Without_DG_and_XDG);

                        var sgrd = lsTrk.Regions.GetCutCellSubGrid();
                        var cellScheme = new CellQuadratureScheme(UseDefaultFactories: true, domain: sgrd.VolumeMask);
                        var edgeScheme = new EdgeQuadratureScheme(UseDefaultFactories: true, domain: sgrd.InnerEdgesMask);
                        //foreach(var SpeciesId in ReqSpecies) {
                        //    ctorTraceDGBulkIntegrator(SpeciesId, quadOrder, cellScheme, edgeScheme, TraceDgDomainFrame, TraceDgCodomFrame,
                        //    TraceDgCodomFrame.FrameMap.BasisS.Select(b => b.MaximalLength > 0).ToArray(),
                        //    Parameters?.ToArray(), DomainFields?.ToArray());
                        //}

                        ctorTraceDGBulkIntegrator(ReqSpecies[0], quadOrder, cellScheme, edgeScheme, TraceDgDomainFrame, TraceDgCodomFrame,
                        TraceDgCodomFrame.FrameMap.BasisS.Select(b => b.MaximalLength > 0).ToArray(),
                        Parameters?.ToArray(), DomainFields?.ToArray());

                    }


                    // coupling terms
                    // --------------

                    using(new BlockTrace("surface_integration", tr)) {
                        if(m_Xowner.ContainesComponentType(typeof(ILevelSetForm))) {


                            var AllSpc = lsTrk.SpeciesIdS;

                            // loop over all possible pairs of species
                            for(int iSpcA = 0; iSpcA < AllSpc.Count; iSpcA++) {
                                var SpeciesA = AllSpc[iSpcA];
                                var SpeciesADom = TrackerRegions.GetSpeciesMask(SpeciesA);
                                //if (SpeciesADom.NoOfItemsLocally <= 0)
                                //    continue;

                                int _iSpcA = Array.IndexOf(ReqSpecies, SpeciesA);

                                for(int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                                    var SpeciesB = AllSpc[iSpcB];

                                    int _iSpcB = Array.IndexOf(ReqSpecies, SpeciesB);
                                    if(_iSpcA < 0 && _iSpcB < 0)
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
                                    for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {
                                        if(SchemeHelper.SpeciesAreSeparatedByLevSet(iLevSet, SpeciesA, SpeciesB)) {
                                            var LsDom = TrackerRegions.GetCutCellMask4LevSet(iLevSet);
                                            var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                            Chunk c = IntegrationDom.FirstOrDefault();
                                            if(c.Len > 0) {
                                                Debug.Assert(IntegrationDom.IsEmptyOnRank == false);
                                                int jtest = c.i0;

                                                LevelsetCellSignCode csc = lsTrk.RegionsHistory[__TrackerHistoryIndex].GetCellSignCode(jtest);

                                                if(!(csc.GetSign(iLevSet) == LevelsetSign.Both))
                                                    throw new ApplicationException("Seem to perform level-set integration in a non-cut cell.");

                                                var cscNeg = csc; cscNeg.SetSign(iLevSet, LevelsetSign.Negative);
                                                var cscPos = csc; cscPos.SetSign(iLevSet, LevelsetSign.Positive);


                                                bool SpeciesA_inNeg = lsTrk.ContainesSpecies(SpeciesA, cscNeg);
                                                bool SpeciesA_inPos = lsTrk.ContainesSpecies(SpeciesA, cscPos);
                                                bool SpeciesB_inNeg = lsTrk.ContainesSpecies(SpeciesB, cscNeg);
                                                bool SpeciesB_inPos = lsTrk.ContainesSpecies(SpeciesB, cscPos);
/*
                                                bool[] SpeciesA_inNeg_all = SpeciesA_inNeg.MPIAllGatherO();
                                                bool[] SpeciesA_inPos_all = SpeciesA_inPos.MPIAllGatherO();
                                                bool[] SpeciesB_inNeg_all = SpeciesB_inNeg.MPIAllGatherO();
                                                bool[] SpeciesB_inPos_all = SpeciesB_inPos.MPIAllGatherO();
                                                for(int rnk = 0; rnk < m_lsTrk.GridDat.MpiSize; rnk++) {
                                                    if(SpeciesA_inNeg_all[rnk] != SpeciesA_inNeg
                                                    || SpeciesA_inPos_all[rnk] != SpeciesA_inPos
                                                    || SpeciesB_inNeg_all[rnk] != SpeciesB_inNeg
                                                    || SpeciesB_inPos_all[rnk] != SpeciesB_inPos) {
                                                        WaitForDebugger();
                                                   }

                                                


                                                }
*/

                                                if(SpeciesA_inPos == SpeciesA_inNeg) {
                                                    throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesA)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                }
                                                if(SpeciesB_inPos == SpeciesB_inNeg) {
                                                    throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesB)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                }


                                                if(SpeciesA_inNeg && SpeciesB_inPos) {
                                                    // nothing to do

                                                    if(SpeciesA_inPos == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesA)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }
                                                    if(SpeciesB_inNeg == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesB)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }
                                                } else if(SpeciesA_inPos && SpeciesB_inNeg) {
                                                    // flip species

                                                    if(SpeciesA_inNeg == true) {
                                                        throw new ApplicationException($"Species {m_lsTrk.GetSpeciesName(SpeciesA)} seems to be present in negative and positive domain of Level-Set No. {iLevSet} - internal error or illegal Level-Set and species map..");
                                                    }
                                                    if(SpeciesB_inPos == true) {
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
                                            using(new BlockTrace("QuadRule-compilation", tr)) {
                                                //CellQuadratureScheme SurfIntegration = SchemeHelper.GetLevelSetquadScheme(iLevSet, SpeciesA, IntegrationDom);
                                                CellQuadratureScheme SurfIntegration = m_Xowner.CouplingQuadraturSchemeProvider(lsTrk, iLevSet, SpeciesA, SpeciesB, IntegrationDom, SchemeHelper, quadOrder, __TrackerHistoryIndex);
                                                rule = SurfIntegration.Compile(GridData, quadOrder);

                                                if(onlyfordebugging_RuleDiagnosis) {
                                                    var suffix = $"{iLevSet}-{lsTrk.GetSpeciesName(SpeciesA)}{lsTrk.GetSpeciesName(SpeciesB)}-{lsTrk.CutCellQuadratureType}-MPI{this.GridData.MpiRank}";
                                                    rule.SaveToTextFileCell(GridData, $"Levset{suffix}.csv");
                                                    rule.SumOfWeightsToTextFileVolume(GridData, $"SumOfWgtLevset{suffix}.csv");
                                                }
                                            }

                                            LECQuadratureLevelSet<IMutableMatrix, double[]>.TestNegativeAndPositiveSpecies(rule, m_lsTrk, __TrackerHistoryIndex, SpeciesA, SpeciesB, iLevSet);

                                            CouplingRules.Add((iLevSet, SpeciesA, SpeciesB, rule));
                                            /*{
                                                var tt = CouplingRules.Last();

                                                int __iLevSet = tt.LsIdx;
                                                var __SpeciesA = tt.spcA;
                                                var __SpeciesB = tt.spcB;
                                                int NoOfItems = tt.quadRule.Count();
                                                int[] all_iLs = __iLevSet.MPIAllGather();
                                                int[] allSpcA = __SpeciesA.cntnt.MPIAllGather();
                                                int[] allSpcB = __SpeciesB.cntnt.MPIAllGather();
                                                for(int rnk = 0; rnk < all_iLs.Length; rnk++) {
                                                    if(allSpcA[rnk] != __SpeciesA.cntnt)
                                                        WaitForDebugger();
                                                    Assert.IsTrue(all_iLs[rnk] == __iLevSet, "Level-set index mismatch.");
                                                    Assert.IsTrue(allSpcA[rnk] == __SpeciesA.cntnt, "Mismatch of species designators across MPI ranks for species A");
                                                    Assert.IsTrue(allSpcB[rnk] == __SpeciesB.cntnt, "Mismatch of species designators across MPI ranks for species B");
                                                }
                                            }*/

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

            bool AnyTraceDG(UnsetteledCoordinateMapping map) {
                foreach(var b in map.BasisS) {
                    if(b is TraceDGBasis)
                        return true;
                }
                return false;
            }
            bool AnyNonTrace(UnsetteledCoordinateMapping map) {
                foreach(var b in map.BasisS) {
                    if(b is TraceDGBasis)
                        continue;
                    return true;
                }
                return false;
            }

            protected bool AnyTraceDG_Domain => AnyTraceDG(DomainMapping);
            protected bool AnyTraceDG_Codomn => AnyTraceDG(DomainMapping);
            protected bool AnyNonTraceDG_Domain => AnyNonTrace(DomainMapping);
            protected bool AnyNonTraceDG_Codomn => AnyNonTrace(DomainMapping);


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
            abstract protected void ctorSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld4Species);
            
            /// <summary>
            /// create integrator for bulk components (stabilization components) for Trace DG
            /// </summary>
            //abstract protected void ctorTraceDGBulkIntegrator(int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params, DGField[] DomFld);
            abstract protected void ctorTraceDGBulkIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params, DGField[] DomFld);


            /// <summary>
            /// Create integrator for <see cref="XDifferentialOperatorMk2.GhostEdgesOperator"/>
            /// </summary>
            abstract protected void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld4Species);

            /// <summary>
            /// Create integrator for <see cref="XDifferentialOperatorMk2.SurfaceElementOperator_Ls0"/>
            /// </summary>
            abstract protected void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld4Species);

            /// <summary>
            /// Create Integrator for <see cref="XDifferentialOperatorMk2.ContactLineOperator_Ls0"/>
            /// </summary>
            abstract protected void ctorContactLineSpeciesIntegrator(SpeciesId SpeciesId, int quadOrder, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, bool[] RowSwitch, DGField[] Params_4Species, DGField[] DomFld4Species);


            /// <summary>
            /// Create integrator for <see cref="ILevelSetForm"/> components
            /// </summary>
            virtual protected void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule) { }

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
            /// for each species, the frame of the co-domain, no TraceDG
            /// </summary>
            protected Dictionary<SpeciesId, FrameBase> SpeciesCodomFrame_WithoutTraceDg = new Dictionary<SpeciesId, FrameBase>();

            /// <summary>
            /// for each species, the frame of the domain, no TraceDG 
            /// </summary>
            protected Dictionary<SpeciesId, FrameBase> SpeciesDomainFrame_WithoutTraceDg = new Dictionary<SpeciesId, FrameBase>();

            /// <summary>
            /// for each species, the frame of the co-domain, TraceDG included
            /// </summary>
            protected Dictionary<SpeciesId, FrameBase> SpeciesCodomFrame_WithTraceDg = new Dictionary<SpeciesId, FrameBase>();

            /// <summary>
            /// for each species, the frame of the domain, TraceDG included
            /// </summary>
            protected Dictionary<SpeciesId, FrameBase> SpeciesDomainFrame_WithTraceDg = new Dictionary<SpeciesId, FrameBase>();


            /// <summary>
            /// If TraceDG is used, the frame of the co-domain, without any DG or XDG components
            /// </summary>
            protected FrameBase TraceDgCodomFrame = null;

            /// <summary>
            /// If TraceDG is used,  the frame of the domain, without any DG or XDG components
            /// </summary>
            protected FrameBase TraceDgDomainFrame = null;

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
