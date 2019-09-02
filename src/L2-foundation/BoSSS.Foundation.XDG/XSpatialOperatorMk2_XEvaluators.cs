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

using static BoSSS.Foundation.SpatialOperator;


namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// An operator which is specialized in XDG Fields, i.e.
    /// it can have components which couple the phases.
    /// Mk2: enables the definition of different equation components for each phase 
    /// </summary>
    partial class XSpatialOperatorMk2 {

        /// <summary>
        /// Assembly of matrices for XDG operators
        /// </summary>
        public class XEvaluatorLinear : XEvaluatorBase, IXEvaluatorLinear {

            /// <summary>
            /// ctor
            /// </summary>
            internal XEvaluatorLinear(XSpatialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory,
                IDictionary<SpeciesId, QrSchemPair> __SpeciesSchemes) :
                base(ownr, lsTrk, DomainVarMap, null, ParameterMap, CodomainVarMap, TrackerHistory, __SpeciesSchemes) //
            {
                using (var tr = new FuncTrace()) {
                    base.MPITtransceive = true;


                }
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                var BulkMtxBuilder = base.m_Xowner.GetSpeciesMatrixBuilderBase(SpeciesId, DomainFrame.FrameMap, Params, CodomFrame.FrameMap,
                                edgeScheme, cellScheme);
                Debug.Assert(((EvaluatorBase)BulkMtxBuilder).order == base.order);
                BulkMtxBuilder.MPITtransceive = false;
                SpeciesBulkMtxBuilder.Add(SpeciesId, BulkMtxBuilder);
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                Debug.Assert(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0);

                var GhostEdgeBuilder = m_Xowner.GhostEdgesOperator.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap,
                    ghostEdgeScheme, nullvolumeScheme);
                Debug.Assert(((EvaluatorBase)GhostEdgeBuilder).order == base.order);
                GhostEdgeBuilder.MPITtransceive = false;
                SpeciesGhostEdgeBuilder.Add(SpeciesId, GhostEdgeBuilder);

            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                Debug.Assert(m_Xowner.SurfaceElementOperator.TotalNoOfComponents > 0);

                var SurfElmBuilder = m_Xowner.SurfaceElementOperator.GetMatrixBuilder(DomainFrame.FrameMap, Params, CodomFrame.FrameMap, SurfaceElement_Edge, SurfaceElement_volume);
                Debug.Assert(((EvaluatorBase)SurfElmBuilder).order == base.order);
                SurfElmBuilder.MPITtransceive = false;
                SpeciesSurfElmBuilder.Add(SpeciesId, SurfElmBuilder);

            }


            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesBulkMtxBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesGhostEdgeBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();
            Dictionary<SpeciesId, IEvaluatorLinear> SpeciesSurfElmBuilder = new Dictionary<SpeciesId, IEvaluatorLinear>();



            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                ComputeMatrix_Internal(default(BlockMsrMatrix), AffineOffset, true);
            }

            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset)
                where M : IMutableMatrixEx
                where V : IList<double> {
                ComputeMatrix_Internal(Matrix, AffineOffset, false);
            }

            protected override DGField[] GetTrxFields() {
                return base.Parameters.ToArray();
            }



            /// <summary>
            /// computation of operator matrix, currently only two species are supported
            /// </summary>
            void ComputeMatrix_Internal<M, V>(
                M Matrix, V AffineOffset, bool OnlyAffine)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {
                if (base.MPITtransceive == true)
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                using (var tr = new FuncTrace()) {
                    var lsTrk = base.m_lsTrk;
                    IGridData GridDat = lsTrk.GridDat;

                    #region Check Input Arguments
                    // --------------------------------
                    //if (!this.IsCommited)
                    //    throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");

                    //if (DomainMap.BasisS.Count != this.DomainVar.Count)
                    //    throw new ArgumentException("mismatch between specified domain variables and number of DG fields in domain mapping", "DomainMap");
                    //if (CodomainMap.BasisS.Count != this.CodomainVar.Count)
                    //    throw new ArgumentException("mismatch between specified codomain variables and number of DG fields in codomain mapping", "CodomainMap");

                    //if (this.ParameterVar.Count == 0) {
                    //    if (Parameters != null && Parameters.Count > 0)
                    //        throw new ArgumentException("mismatch between specified parameter variables and number of DG fields in parameter mapping", "Parameters");
                    //} else {
                    //    if (Parameters == null)
                    //        throw new ArgumentNullException("Parameters", "parameters must be specified");
                    //    if (Parameters.Count != this.ParameterVar.Count)
                    //        throw new ArgumentException("mismatch between specified parameter variables and number of DG fields in parameter mapping", "Parameters");
                    //}

                    if (OnlyAffine == false) {
                        if (!Matrix.RowPartitioning.Equals(base.CodomainMapping))
                            throw new ArgumentException("wrong number of columns in matrix.", "Matrix");
                        if (!Matrix.ColPartition.Equals(base.DomainMapping))
                            throw new ArgumentException("wrong number of rows in matrix.", "Matrix");
                    }


                    //if (!ReqSpecies.IsSubsetOf(agg.SpeciesList))
                    //    throw new ArgumentException("HMF mismatch");

                    //if (momentFittingVariant != agg.HMFvariant)
                    //    throw new ArgumentException("HMF mismatch");

                    #endregion

                    #region MPI exchange of parameter fields
                    // --------------------------------
                    Transceiver trx = base.m_TRX;
                    if (trx != null) {
                        trx.TransceiveStartImReturn();
                    }
                    #endregion

                    #region find quadrature instructions
                    // ----------------------------



                    #endregion


                    // build matrix, bulk
                    // ---------------------
                    //MsrMatrix BulkMatrix = null;
                    //double[] BulkAffineOffset = null;
                    using (new BlockTrace("bulk_integration", tr)) {

                        // create the frame matrices & vectors...
                        // this is an MPI-collective operation, so it must be executed before the program may take different branches...
                        SpeciesFrameMatrix<M>[] mtx_spc = new SpeciesFrameMatrix<M>[ReqSpecies.Length];
                        SpeciesFrameVector<V>[] vec_spc = new SpeciesFrameVector<V>[ReqSpecies.Length];
                        for (int i = 0; i < ReqSpecies.Length; i++) {
                            SpeciesId SpId = ReqSpecies[i];
                            mtx_spc[i] = new SpeciesFrameMatrix<M>(Matrix, this.SpeciesCodomFrame[SpId], this.SpeciesDomainFrame[SpId]);
                            vec_spc[i] = (AffineOffset != null) ?
                                    (new SpeciesFrameVector<V>(AffineOffset, this.SpeciesCodomFrame[SpId]))
                                    :
                                    null;
                        }
                        //// Create Masks before the Loops, so it doesn't affect MPI
                        //CellMask SubGridCellMask = null;
                        //EdgeMask SubGridEdgeMask = null;
                        //if (SubGrid != null) {
                        //    SubGridCellMask = SubGrid.VolumeMask;
                        //    /// I don't know why, but this seems to work:
                        //    SubGridEdgeMask = SubGrid.AllEdgesMask;
                        //    /// And this does not:
                        //    //SubGridEdgeMask = SubGrid.InnerEdgesMask;
                        //}

                        // do the Bulk integration...
                        foreach (var SpeciesId in ReqSpecies) {
                            int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);


                            //if(m_Xowner.OnIntegratingBulk != null)
                            //    m_Xowner.OnIntegratingBulk(lsTrk.GetSpeciesName(SpeciesId), SpeciesId);

                            SpeciesFrameMatrix<M> mtx = mtx_spc[iSpecies];
                            var _mtx = Matrix != null ? mtx : default(SpeciesFrameMatrix<M>);

                            SpeciesFrameVector<V> vec = vec_spc[iSpecies];


                            //DGField[] Params = 


                            //#if DEBUG
                            //                            // switch the diagnostic output on or off
                            //                            bool SubGridRuleDiagnosis = false;
                            //                            if (SubGrid == null && SubGridRuleDiagnosis == true) {
                            //                                Console.WriteLine("Warning SubGrid Rule Diagnosis is Switched on!");
                            //                            }
                            //                            if (SubGridRuleDiagnosis) {
                            //                                edgeRule.SumOfWeightsToTextFileEdge(GridDat, string.Format("C:\\tmp\\BoSSS_Diagnosis\\PhysEdge_{0}.csv", lsTrk.GetSpeciesName(SpeciesId)));
                            //                                volRule.SumOfWeightsToTextFileVolume(GridDat, string.Format("C:\\tmp\\BoSSS_Diagnosis\\PhysVol_{0}.csv", lsTrk.GetSpeciesName(SpeciesId)));
                            //                            }
                            //#endif

                            foreach (var SpeciesBuilder in new[] { SpeciesBulkMtxBuilder, SpeciesGhostEdgeBuilder, SpeciesSurfElmBuilder }) {

                                if (SpeciesBuilder.ContainsKey(SpeciesId)) {

                                    var builder = SpeciesBuilder[SpeciesId];
                                    builder.OperatorCoefficients = this.SpeciesOperatorCoefficients[SpeciesId];
                                    NotifySpecies(builder.Owner, this.m_lsTrk, SpeciesId);

                                    if (trx != null) {
                                        trx.TransceiveFinish();
                                        trx = null;
                                    }

                                    builder.time = base.time;

                                    if (OnlyAffine) {
                                        builder.ComputeAffine(vec);
                                    } else {
                                        builder.ComputeMatrix(_mtx, vec);
                                    }
                                }

                            }
                        }

                    }

                    // build matrix, coupling
                    ///////////////////

                    using (new BlockTrace("surface_integration", tr)) {
                        foreach (var tt in this.CouplingRules) {
                            int iLevSet = tt.Item1;
                            var SpeciesA = tt.Item2;
                            var SpeciesB = tt.Item3;
                            var rule = tt.Item4;



                            if (trx != null) {
                                trx.TransceiveFinish();
                                trx = null;
                            }

                            var MtxBuilder = new LECQuadratureLevelSet<M, V>(GridDat,
                                                             m_Xowner,
                                                             OnlyAffine ? default(M) : Matrix, AffineOffset,
                                                             CodomainMapping, Parameters, DomainMapping,
                                                             lsTrk, iLevSet, new Tuple<SpeciesId, SpeciesId>(SpeciesA, SpeciesB),
                                                             rule);
                            MtxBuilder.time = time;
                            this.SpeciesOperatorCoefficients.TryGetValue(SpeciesA, out var csA);
                            this.SpeciesOperatorCoefficients.TryGetValue(SpeciesB, out var csB);
                            UpdateLevelSetCoefficients(csA, csB);
                            MtxBuilder.Execute();

#if DEBUG
                            if (Matrix != null && OnlyAffine == false)
                                Matrix.CheckForNanOrInfM();
                            if (AffineOffset != null)
                                GenericBlas.CheckForNanOrInfV(AffineOffset);
#endif

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
        }


        public class XEvaluatorNonlin : XEvaluatorBase, IXEvaluatorNonLin {



            internal XEvaluatorNonlin(XSpatialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                CoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory,
                IDictionary<SpeciesId, QrSchemPair> __SpeciesSchemes) :
                base(ownr, lsTrk, DomainVarMap, DomainVarMap.Fields, ParameterMap, CodomainVarMap, TrackerHistory, __SpeciesSchemes) //
            {
                this.DomainFields = DomainVarMap;
                base.MPITtransceive = true;

            }

            ///// <summary>
            ///// Domain fields for each species, for <see cref="XDGField"/>s the species shadow, see <see cref="XDGField.GetSpeciesShadowField(SpeciesId)"/>
            ///// </summary>
            //protected Dictionary<SpeciesId, DGField[]> SpeciesDomFieleds = new Dictionary<SpeciesId, DGField[]>();


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
                    //MsrMatrix BulkMatrix = null;
                    //double[] BulkAffineOffset = null;
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

                            foreach (var SpeciesEval in new[] { SpeciesBulkEval, SpeciesGhostEval, SpeciesSurfElmEval }) {

                                if (SpeciesEval.ContainsKey(SpeciesId)) {

                                    var eval = SpeciesEval[SpeciesId];
                                    eval.OperatorCoefficients = this.SpeciesOperatorCoefficients[SpeciesId];
                                    NotifySpecies(eval.Owner, this.m_lsTrk, SpeciesId);

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


                    //if(OnIntegratingSurfaceElement != null)
                    //    OnIntegratingSurfaceElement(lsTrk.GetSpeciesName(SpeciesId), SpeciesId, InterfaceLengths[SpeciesId]);


                    //SurfaceElementOperator.ComputeMatrixEx(
                    //    mtx.ColMapping, Params, mtx.RowMapping,
                    //    _mtx, vec, OnlyAffine, time,
                    //    SurfaceElement_Edge, SurfaceElement_volume, null);


                    using (new BlockTrace("surface_integration", tr)) {
                        foreach (var tt in this.CouplingRules) {
                            int iLevSet = tt.Item1;
                            var SpeciesA = tt.Item2;
                            var SpeciesB = tt.Item3;
                            var rule = tt.Item4;

                            if (trx != null) {
                                trx.TransceiveFinish();
                                trx = null;
                            }

                            var LsEval = new NECQuadratureLevelSet<Tout>(GridDat,
                                                             m_Xowner,
                                                             output,
                                                             this.DomainFields.Fields, Parameters, base.CodomainMapping,
                                                             lsTrk, iLevSet, new Tuple<SpeciesId, SpeciesId>(SpeciesA, SpeciesB),
                                                             rule);
                            LsEval.time = time;
                            this.SpeciesOperatorCoefficients.TryGetValue(SpeciesA, out var csA);
                            this.SpeciesOperatorCoefficients.TryGetValue(SpeciesB, out var csB);
                            UpdateLevelSetCoefficients(csA, csB);
                            LsEval.Execute();

#if DEBUG
                            GenericBlas.CheckForNanOrInfV(output);
#endif
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



            protected override void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule) {

            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                var BulkEval = base.m_Xowner.GetSpeciesEvaluatorExBase(SpeciesId, DomFld, base.SpeciesParams[SpeciesId], CodomFrame.FrameMap,
                                edgeScheme, cellScheme);
                Debug.Assert(((EvaluatorBase)BulkEval).order == base.order);
                BulkEval.MPITtransceive = false;
                SpeciesBulkEval.Add(SpeciesId, BulkEval);
            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                Debug.Assert(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0);

                var GhostEdgeEval = m_Xowner.GhostEdgesOperator.GetEvaluatorEx(DomFld, Params, CodomFrame.FrameMap,
                    ghostEdgeScheme, nullvolumeScheme);
                Debug.Assert(((EvaluatorBase)GhostEdgeEval).order == base.order);
                GhostEdgeEval.MPITtransceive = false;
                SpeciesGhostEval.Add(SpeciesId, GhostEdgeEval);

            }

            /// <summary>
            /// creates an evaluator
            /// </summary>
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld) {
                Debug.Assert(m_Xowner.SurfaceElementOperator.TotalNoOfComponents > 0);

                var SurfElmEval = m_Xowner.SurfaceElementOperator.GetEvaluatorEx(DomFld, Params, CodomFrame.FrameMap, SurfaceElement_Edge, SurfaceElement_volume);
                Debug.Assert(((EvaluatorBase)SurfElmEval).order == base.order);
                SurfElmEval.MPITtransceive = false;
                SpeciesSurfElmEval.Add(SpeciesId, SurfElmEval);

            }


            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesBulkEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesGhostEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();
            Dictionary<SpeciesId, IEvaluatorNonLin> SpeciesSurfElmEval = new Dictionary<SpeciesId, IEvaluatorNonLin>();


        }


        abstract public class XEvaluatorBase : IXEvaluator {

            XSpatialOperatorMk2 m_Owner;

            /// <summary>
            /// the operator used to construct this object
            /// </summary>
            public XSpatialOperatorMk2 Owner {
                get {
                    return m_Owner;
                }
            }

            /// <summary>
            /// ctor
            /// </summary>
            protected internal XEvaluatorBase(
                XSpatialOperatorMk2 ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory,
                IDictionary<SpeciesId, QrSchemPair> __SpeciesSchemes) //
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


                    //IEnumerable<Basis> allBasis = DomainVarMap.BasisS;
                    //if (ParameterMap != null) {
                    //    allBasis = allBasis.Union(ParameterMap.Select(f => f.Basis));
                    //}
                    //allBasis = allBasis.Union(CodomainVarMap.BasisS);
                    //IGridData grdDat = allBasis.First().GridDat;
                    //foreach (var b in allBasis) {
                    //    if (!object.ReferenceEquals(grdDat, b.GridDat)) {
                    //        throw new ArgumentException("all fields (domain, parameter, codomain) must be defined on the same grid.");
                    //    }
                    //}

                    if (!m_Owner.IsCommited)
                        throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");

                    order = ownr.GetOrderFromQuadOrderFunction(DomainMapping, ParameterMap, CodomainVarMap);

                    m_OperatorCoefficients = new CoefficientSet() {
                        UserDefinedValues = new Dictionary<string, object>(),
                        GrdDat = this.GridData
                    };

                    if (this.GridData is Grid.Classic.GridData) {
                        m_OperatorCoefficients.CellLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.CellLengthScale;
                        m_OperatorCoefficients.EdgeLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge;

                    } else {
                        Console.WriteLine("Rem: still missing cell length scales");
                    }
                }

                using (var tr = new FuncTrace()) {
                    if (!object.ReferenceEquals(GridData, lsTrk.GridDat))
                        throw new ArgumentException("grid data mismatch");
                    m_lsTrk = lsTrk;
                    m_Xowner = ownr;
                    SpeciesSchemes = __SpeciesSchemes;
                    ReqSpecies = SpeciesSchemes.Keys.ToArray();


                    var SchemeHelper = lsTrk.GetXDGSpaceMetrics(ReqSpecies, order, TrackerHistory).XQuadSchemeHelper;
                    var TrackerRegions = lsTrk.RegionsHistory[TrackerHistory];

                    ((FixedOrder_SpatialOperator)(m_Xowner.GhostEdgesOperator)).m_Order = order;
                    ((FixedOrder_SpatialOperator)(m_Xowner.SurfaceElementOperator)).m_Order = order;
                    foreach (SpeciesId spcid in ReqSpecies) {
                        ((FixedOrder_SpatialOperator)m_Xowner.m_SpeciesOperator[spcid]).m_Order = order;
                    }
                    //((FixedOrder_SpatialOperator)(m_Xowner.SpeciesOperator)).m_Order = order;
                    tr.Info("XSpatialOperator.ComputeMatrixEx quad order: " + order);


                    // compile quadrature rules & create matrix builders for each species 
                    // ------------------------------------------------------------------
                    foreach (var SpeciesId in ReqSpecies) {
                        int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);

                        // parameters for species
                        // ----------------------
                        DGField[] Params = (from f in (Parameters ?? new DGField[0])
                                            select ((f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f)).ToArray<DGField>();
                        SpeciesParams.Add(SpeciesId, Params);

                        DGField[] DomFld;
                        if (DomainFields != null) {
                            DomFld = DomainFields.Select(f => (f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f).ToArray();
                        } else {
                            DomFld = null;
                        }

                        // species frames
                        // --------------

                        var CodomFrame = new FrameBase(TrackerRegions, SpeciesId, CodomainMapping, false);
                        var DomainFrame = new FrameBase(TrackerRegions, SpeciesId, DomainMapping, true);
                        SpeciesCodomFrame.Add(SpeciesId, CodomFrame);
                        SpeciesDomainFrame.Add(SpeciesId, DomainFrame);

                        // quadrature rules
                        // ----------------
                        EdgeQuadratureScheme edgeScheme;
                        CellQuadratureScheme cellScheme;
                        using (new BlockTrace("QuadRule-compilation", tr)) {

                            var qrSchemes = SpeciesSchemes[SpeciesId];

                            //bool AssembleOnFullGrid = (SubGrid == null);
                            if (qrSchemes.EdgeScheme == null) {
                                edgeScheme = SchemeHelper.GetEdgeQuadScheme(SpeciesId);// AssembleOnFullGrid, SubGridEdgeMask);
                            } else {
                                edgeScheme = qrSchemes.EdgeScheme;
                                //throw new NotSupportedException();
                            }
                            if (qrSchemes.CellScheme == null) {
                                cellScheme = SchemeHelper.GetVolumeQuadScheme(SpeciesId);//, AssembleOnFullGrid, SubGridCellMask);
                            } else {
                                cellScheme = qrSchemes.CellScheme;
                                //throw new NotSupportedException();
                            }
                        }

                        if (ruleDiagnosis) {
                            var edgeRule = edgeScheme.Compile(GridData, order);
                            var volRule = cellScheme.Compile(GridData, order);
                            edgeRule.SumOfWeightsToTextFileEdge(GridData, string.Format("PhysEdge_{0}.csv", lsTrk.GetSpeciesName(SpeciesId)));
                            volRule.SumOfWeightsToTextFileVolume(GridData, string.Format("Volume_{0}.csv", lsTrk.GetSpeciesName(SpeciesId)));
                        }


                        if (m_Xowner.TotalNoOfComponents > 0) {
                            ctorSpeciesIntegrator(SpeciesId, cellScheme, edgeScheme, DomainFrame, CodomFrame, Params, DomFld);
                        }

                        if (m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0) {
                            CellQuadratureScheme nullvolumeScheme = new CellQuadratureScheme(false, CellMask.GetEmptyMask(GridData));
                            EdgeQuadratureScheme ghostEdgeScheme = null;
                            ghostEdgeScheme = SchemeHelper.GetEdgeGhostScheme(SpeciesId);//, SubGridEdgeMask);
                            ctorGhostSpeciesIntegrator(SpeciesId, nullvolumeScheme, ghostEdgeScheme, DomainFrame, CodomFrame, Params, DomFld);
                        }

                        if (m_Xowner.SurfaceElementOperator.TotalNoOfComponents > 0) {
                            EdgeQuadratureScheme SurfaceElement_Edge;
                            CellQuadratureScheme SurfaceElement_volume;
                            SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(SpeciesId);
                            SurfaceElement_volume = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(SpeciesId);
                            ctorSurfaceElementSpeciesIntegrator(SpeciesId, SurfaceElement_volume, SurfaceElement_Edge, DomainFrame, CodomFrame, Params, DomFld);
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
                                var SpeciesADom = lsTrk.Regions.GetSpeciesMask(SpeciesA);
                                if (SpeciesADom.NoOfItemsLocally <= 0)
                                    continue;

                                int _iSpcA = Array.IndexOf(ReqSpecies, SpeciesA);

                                for (int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                                    var SpeciesB = AllSpc[iSpcB];

                                    int _iSpcB = Array.IndexOf(ReqSpecies, SpeciesB);
                                    if (_iSpcA < 0 && _iSpcB < 0)
                                        continue;

                                    var SpeciesBDom = lsTrk.Regions.GetSpeciesMask(SpeciesB);
                                    var SpeciesCommonDom = SpeciesADom.Intersect(SpeciesBDom);

                                    // Checks removed since they can cause parallel problems
                                    //if (SpeciesBDom.NoOfItemsLocally <= 0)
                                    //    continue;
                                    //if (SpeciesCommonDom.NoOfItemsLocally <= 0)
                                    //    continue;


                                    // loop over level-sets
                                    int NoOfLs = lsTrk.LevelSets.Count;
                                    for (int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {



                                        var LsDom = lsTrk.Regions.GetCutCellMask4LevSet(iLevSet);
                                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                        // Check removed since it can cause parallel problems
                                        //if (IntegrationDom.NoOfItemsLocally > 0) {

                                        ICompositeQuadRule<QuadRule> rule;
                                        using (new BlockTrace("QuadRule-compilation", tr)) {
                                            CellQuadratureScheme SurfIntegration = SchemeHelper.GetLevelSetquadScheme(iLevSet, IntegrationDom);
                                            rule = SurfIntegration.Compile(GridData, order);

                                            if (ruleDiagnosis) {
                                                rule.SumOfWeightsToTextFileVolume(GridData, string.Format("Levset_{0}.csv", iLevSet));
                                            }
                                        }

                                        CouplingRules.Add(new Tuple<int, SpeciesId, SpeciesId, ICompositeQuadRule<QuadRule>>(
                                            iLevSet, SpeciesA, SpeciesB, rule));
                                        ctorLevSetFormIntegrator(iLevSet, SpeciesA, SpeciesB, rule);
                                    }
                                }
                            }
                        }
                    }

                    // coeff kacke
                    // -----------

                    this.SpeciesOperatorCoefficients = new Dictionary<SpeciesId, CoefficientSet>();
                    foreach (var SpeciesId in ReqSpecies) {
                        this.SpeciesOperatorCoefficients.Add(SpeciesId,
                            new CoefficientSet() {
                                CellLengthScales = null, // ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.cj,
                                EdgeLengthScales = null, //((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge,
                                UserDefinedValues = new Dictionary<string, object>(),
                                GrdDat = this.GridData
                            });
                    }
                    // m_OperatorCoefficients = new CoefficientSet() {
                    //    CellLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.cj,
                    //    EdgeLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge,
                    //    UserDefinedValues = new Dictionary<string, object>()
                    //};
                }
            }

            CoefficientSet m_OperatorCoefficients;


            /// <summary>
            /// create integrator for bulk phase
            /// </summary>
            abstract protected void ctorSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld);

            /// <summary>
            /// Create integrator for <see cref="XSpatialOperator.GhostEdgesOperator"/>
            /// </summary>
            abstract protected void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld);

            /// <summary>
            /// Create integrator for <see cref="XSpatialOperator.SurfaceElementOperator"/>
            /// </summary>
            abstract protected void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params, DGField[] DomFld);

            /// <summary>
            /// Create integrator for <see cref="ILevelSetForm"/> components
            /// </summary>
            abstract protected void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule);

            /// <summary>
            /// all species to integrate
            /// </summary>
            protected SpeciesId[] ReqSpecies;

            /// <summary>
            /// %
            /// </summary>
            protected LevelSetTracker m_lsTrk;

            /// <summary>
            /// the spatial operator
            /// </summary>
            protected XSpatialOperatorMk2 m_Xowner;

            /// <summary>
            /// edge and volume scheme for each species
            /// </summary>
            protected IDictionary<SpeciesId, QrSchemPair> SpeciesSchemes;

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
            protected List<Tuple<int, SpeciesId, SpeciesId, ICompositeQuadRule<QuadRule>>> CouplingRules = new List<Tuple<int, SpeciesId, SpeciesId, ICompositeQuadRule<QuadRule>>>();



            static bool ruleDiagnosis = false;

            /// <summary>
            /// calls all <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate(CoefficientSet, CoefficientSet, int[], int)"/> methods
            /// </summary>
            protected void UpdateLevelSetCoefficients(CoefficientSet csA, CoefficientSet csB) {
                int[] DomDGdeg = this.DomainMapping.BasisS.Select(b => b.Degree).ToArray();
                int[] CodDGdeg = this.CodomainMapping.BasisS.Select(b => b.Degree).ToArray();
                string[] DomNames = Owner.DomainVar.ToArray();
                string[] CodNames = Owner.CodomainVar.ToArray();


                Debug.Assert(CodNames.Length == CodDGdeg.Length);
                for (int iCod = 0; iCod < CodDGdeg.Length; iCod++) {
                    var comps = Owner.EquationComponents[CodNames[iCod]];
                    foreach (var c in comps) {
                        if (c is ILevelSetEquationComponentCoefficient) {
                            var ce = c as ILevelSetEquationComponentCoefficient;

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
            /// calls all <see cref="IEquationComponentSpeciesNotification.SetParameter(string, SpeciesId)"/> methods
            /// </summary>
            protected static void NotifySpecies(SpatialOperator Owner, LevelSetTracker lsTrk, SpeciesId id) {
                string sNmn = lsTrk.GetSpeciesName(id);

                string[] CodNames = Owner.CodomainVar.ToArray();
                for (int iCod = 0; iCod < CodNames.Length; iCod++) {
                    var comps = Owner.EquationComponents[CodNames[iCod]];
                    foreach (var c in comps) {
                        if (c is IEquationComponentSpeciesNotification) {
                            var ce = c as IEquationComponentSpeciesNotification;
                            ce.SetParameter(sNmn, id);
                        }
                    }
                }
            }


            /// <summary>
            /// Operator coefficients for each species
            /// </summary>
            public Dictionary<SpeciesId, CoefficientSet> SpeciesOperatorCoefficients {
                get;
            }


            /// <summary>
            /// Sets the coefficients for all equation components of the operator which implement <see cref="IEquationComponentCoefficient"/>.
            /// </summary>
            protected void UpdateCoefficients() {
                int[] DomDGdeg = this.DomainMapping.BasisS.Select(b => b.Degree).ToArray();
                int[] CodDGdeg = this.CodomainMapping.BasisS.Select(b => b.Degree).ToArray();
                string[] DomNames = m_Owner.DomainVar.ToArray();
                string[] CodNames = m_Owner.CodomainVar.ToArray();


                Debug.Assert(CodNames.Length == CodDGdeg.Length);
                for (int iCod = 0; iCod < CodDGdeg.Length; iCod++) {
                    var comps = m_Owner.m_EquationComponents[CodNames[iCod]];
                    foreach (var c in comps) {
                        if (c is IEquationComponentCoefficient) {
                            var ce = c as IEquationComponentCoefficient;

                            int[] DomDGdeg_cd = new int[ce.ArgumentOrdering.Count];
                            for (int i = 0; i < DomDGdeg_cd.Length; i++) {
                                string domName = ce.ArgumentOrdering[i];
                                int idx = Array.IndexOf(DomNames, domName);
                                DomDGdeg_cd[i] = DomDGdeg[idx];
                            }


                            ce.CoefficientUpdate(m_OperatorCoefficients, DomDGdeg_cd, CodDGdeg[iCod]);
                        }
                    }
                }
            }


            /// <summary>
            /// Quadrature order to compile quadrature schemes
            /// </summary>
            public int order {
                get;
                private set;
            }

            SubGridBoundaryModes m_SubGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge;

            /// <summary>
            /// 
            /// </summary>
            protected CellMask m_SubGrid_InCells;


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



        /// <summary>
        /// The 'X' version of the Jacobian builder, computes the (approximate) Jacobian matrix of the spatial operator by finite differences.
        /// </summary>
        public class XFDJacobianBuilder : IXEvaluatorLinear {


            /// <summary>
            /// Not for direct user interaction
            /// </summary>
            internal protected XFDJacobianBuilder(IXEvaluatorNonLin __XEval, DelParameterUpdate __delParameterUpdate) {

                eps = 1.0;
                while (1.0 + eps > 1.0) {
                    eps = eps / 2;
                }
                eps = Math.Sqrt(eps);

                XEval = __XEval;
                XEval.MPITtransceive = true;
                DelParamUpdate = __delParameterUpdate;

                BuildGridColoring();
            }

            /// <summary>
            /// Approximately the square root of the machine epsilon.
            /// </summary>
            double eps;

            /// <summary>
            /// Epsilon used for the finite difference
            /// </summary>
            public double Eps {
                get {
                    return eps;
                }
                set {
                    if (value <= 0.0)
                        throw new ArgumentOutOfRangeException();
                    eps = value;
                }
            }

            /// <summary>
            /// Grid
            /// </summary>
            public IGridData GridData {
                get {
                    return XEval.GridData;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public UnsetteledCoordinateMapping CodomainMapping {
                get {
                    return XEval.CodomainMapping;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public UnsetteledCoordinateMapping DomainMapping {
                get {
                    return XEval.DomainMapping;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public IList<DGField> Parameters {
                get {
                    return XEval.Parameters;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public SubGridBoundaryModes SubGridBoundaryTreatment {
                get {
                    return XEval.SubGridBoundaryTreatment;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public double time {
                get { return XEval.time; }
                set { XEval.time = value; }
            }

            /// <summary>
            /// 
            /// </summary>
            public XSpatialOperatorMk2 Owner {
                get {
                    return XEval.Owner;
                }
            }

            /// <summary>
            /// can be dangerous to turn off
            /// </summary>
            public bool MPITtransceive {
                get => throw new NotImplementedException();
                set => throw new NotImplementedException();
            }


            /// <summary>
            /// Operator coefficients for each species
            /// </summary>
            public Dictionary<SpeciesId, CoefficientSet> SpeciesOperatorCoefficients {
                get {
                    return XEval.SpeciesOperatorCoefficients;
                }
            }


            IXEvaluatorNonLin XEval;


            DelParameterUpdate DelParamUpdate;

            /// <summary>
            /// - 1st index: enumeration of color lists
            /// - 2nd index: enumeration of cells in color list
            /// - content: local cell index in which an Epsilon-distortion is applied
            /// </summary>
            int[][] ColorLists;

            /// <summary>
            /// Cells which are distorted on an other processor, but influence the result on this processor
            /// - 1st index: correlates to 1st index of <see cref="ColorLists"/>
            /// - 2nd index: enumeration of cells
            /// - content: some external cell index; this determines the column index of the finite difference result into the Jacobian matrix
            /// </summary>
            int[][] ExternalColorLists;

            /// <summary>
            /// - 1st index: correlates to 1st index of <see cref="ColorLists"/>
            /// - 2nd index: correlates with 2nd index of <see cref="ExternalColorLists"/>
            /// - 3rd index: enumeration
            /// - content: a local cell index of some cell which is affected by an epsilon-distortion on some other processor; this determines the row index of the finite difference result into the Jacobian matrix
            /// </summary>
            int[][][] ExternalColorListsNeighbors;


            void BuildGridColoring() {
                var gDat = XEval.GridData;


                int[][] Neighs = gDat.iLogicalCells.CellNeighbours;
                int J = gDat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = gDat.iLogicalCells.Count;
                long[] GlidxExt = gDat.iParallel.GlobalIndicesExternalCells;
                var Gl2LocExt = gDat.iParallel.Global2LocalIdx;
                var CellPart = gDat.CellPartitioning;

                Random rnd = new Random(gDat.MpiRank + 1);

#if DEBUG
                for (int j = 0; j < J; j++) {
                    int[] CN = Neighs[j];
                    foreach (int jN in CN) {
                        Debug.Assert(jN >= 0);
                        Debug.Assert(jN != j);
                    }
                }
#endif
                /*
                this.ColorLists = new int[J][];
                this.ExternalColorLists = new int[J][];
                this.ExternalColorListsNeighbors = new int[J][][];
                for(int j = 0; j < J; j++) {
                    this.ColorLists[j] = new int[] { j };
                    this.ExternalColorLists[j] = new int[0];
                    this.ExternalColorListsNeighbors[j] = new int[0][];
                }

                return;
                */


                int[] LocalMarker = new int[JE]; //    marker for blocked in the current pass 
                int[] ExchangedMarker = new int[JE]; //  accumulation buffer for MPI exchange
                BitArray Colored = new BitArray(JE); // all cells which are already colored (in previous passes)
                BitArray ColoredPass = new BitArray(JE); // all cells which are colored in current pass
                int[] LocalColorCause = new int[JE];

                List<int> CellList = new List<int>();
                List<int[]> ColorListsTmp = new List<int[]>();
                List<int[]> ExternalColorListsTmp = new List<int[]>();
                List<int[][]> ExternalColorListsNeighborsTmp = new List<int[][]>();

                int myMarkerToken = gDat.MpiRank + 1;
                int bRun = 0xFFFFFF;
                int locColoredCells = 0;
                int DeadlockWatch = 0;
                while (bRun != 0) {

                    // find next color list
                    // ====================

                    Array.Clear(LocalMarker, 0, LocalMarker.Length);// LocalMarker.SetAll(int.MaxValue);
                    ColoredPass.SetAll(false);

                    for (int j = 0; j < J; j++) {
                        if (Colored[j] == true)
                            continue;
                        if (LocalMarker[j] != 0)
                            continue;

                        int[] Neighs_j = Neighs[j];
                        Debug.Assert(Neighs_j.Contains(j) == false, "Cell seems to be its own neighbor.");
                        bool cont = false;
                        foreach (int jn in Neighs_j) {
                            if (LocalMarker[jn] != 0) {
                                cont = true;
                                break;
                            }
                        }
                        if (cont)
                            continue;

                        // if we reached this point, we finally found a cell which we are allowed to add to the current color set.
                        ColoredPass[j] = true;
                        LocalMarker[j] = myMarkerToken;
                        LocalColorCause[j] = j;
                        foreach (int jn in Neighs_j) {
                            LocalMarker[jn] = -myMarkerToken; // mark neighbor cells with negative token!
                            LocalColorCause[jn] = j;
                        }
                    }

                    // fix parallel conflicts
                    // ======================

                    var LocalMarker_Bkup = LocalMarker.CloneAs();
                    var Removed = new List<int>();
                    int[] ExchangedMarker_Bkup = null;

                    if (gDat.MpiSize > 1) {
                        ExchangedMarker_Bkup = ExchangedMarker.CloneAs();

                        int GlobalConflicts = 999;
                        do {
                            Array.Copy(LocalMarker, ExchangedMarker, JE);
                            ExchangedMarker.MPIExchange(gDat);

                            int LocalConflicts = 0;

                            for (int je = J; je < JE; je++) {
                                if (LocalMarker[je] != 0 && ExchangedMarker[je] != 0) {
                                    Debug.Assert(LocalMarker[je] != ExchangedMarker[je]);
                                    LocalConflicts++;

                                    double rndVal = rnd.NextDouble();

                                    //// some parallel conflict detected: one of the two ranks has to yield

                                    //if (ExchangedMarker[je] > 0 && Math.Abs(ExchangedMarker[je]) > myMarkerToken) {
                                    //    // the other rank should yield
                                    //} else {
                                    //    // this rank has to yield
                                    if (rndVal >= 0.5) {
                                        int jToRemove = LocalColorCause[je];
                                        Debug.Assert(jToRemove < J);
                                        Debug.Assert(ColoredPass[jToRemove] == true);

                                        Removed.Add(jToRemove);

                                        ColoredPass[jToRemove] = false;
                                        LocalMarker[jToRemove] = 0;
                                        int[] Neighs_jToRemove = Neighs[jToRemove];
                                        foreach (int jn in Neighs_jToRemove) {
                                            LocalMarker[jn] = 0;
                                        }

                                        //}
                                    }
                                }
                            }

                            GlobalConflicts = LocalConflicts.MPISum();

                        } while (GlobalConflicts > 0);

#if DEBUG
                        Array.Copy(LocalMarker, ExchangedMarker, JE);
                        ExchangedMarker.MPIExchange(gDat);
                        for (int je = J; je < JE; je++) {
                            Debug.Assert(ExchangedMarker[je] == 0 || LocalMarker[je] == 0, "Error in conflict resolution.");
                        }
#endif
                    }



                    // remember recently found color list
                    // ==================================
                    CellList.Clear();
                    for (int j = 0; j < J; j++) {
                        if (ColoredPass[j]) {
                            locColoredCells++;
                            CellList.Add(j);
                            Debug.Assert(Colored[j] == false);
                            Colored[j] = true;
                        }
                    }
                    int LocColoredPass = CellList.Count;

                    int GlobColoredPass = LocColoredPass.MPISum();
                    //Console.WriteLine("Colored in pass {0}: {1}", ColorListsTmp.Count, GlobColoredPass);
                    if (GlobColoredPass <= 0) {
                        DeadlockWatch++;
                        if (DeadlockWatch >= 1000)
                            throw new ApplicationException("Deadlock in parallel coloring.");
                        continue;
                        //Debugger.Launch();
                    }


                    ColorListsTmp.Add(CellList.ToArray());

                    // communicate external lists
                    // ==========================

                    if (gDat.MpiSize > 1) {

                        //Debugger.Launch();

                        var ExchData = new Dictionary<int, List<Tuple<int, int>>>();

                        foreach (int j in CellList) {
                            int[] Neighs_j = Neighs[j];
                            foreach (int jN in Neighs_j) {
                                if (jN >= J) {

                                    int Gl_jN = (int)GlidxExt[jN - J];
                                    int iProc = CellPart.FindProcess(Gl_jN);
                                    int Gl_j = j + CellPart.i0;

                                    if (!ExchData.TryGetValue(iProc, out var ExchData_iProc)) {
                                        ExchData_iProc = new List<Tuple<int, int>>();
                                        ExchData.Add(iProc, ExchData_iProc);
                                    }

                                    ExchData_iProc.Add(new Tuple<int, int>(Gl_j, Gl_jN));
                                }
                            }
                        }

                        var RcvData = SerialisationMessenger.ExchangeData(ExchData);

                        var ExtColor = new Dictionary<int, List<int>>();

                        foreach (var kv in RcvData) {
                            int iProc = kv.Key;
                            var list = kv.Value;

                            foreach (var t in list) {
                                int Gl_j = t.Item1;
                                int Gl_jN = t.Item2;
                                Debug.Assert(CellPart.FindProcess(Gl_j) == iProc);
                                Debug.Assert(CellPart.IsInLocalRange(Gl_jN));

                                int Loc_jN = Gl_jN - CellPart.i0;
                                Debug.Assert(Loc_jN >= 0 && Loc_jN < J);
                                int Loc_j = Gl2LocExt[Gl_j];
                                Debug.Assert(Loc_j >= J && Loc_j < JE);

                                List<int> Neighs_Loc_jN;
                                if (!ExtColor.TryGetValue(Loc_j, out Neighs_Loc_jN)) {
                                    Neighs_Loc_jN = new List<int>();
                                    ExtColor.Add(Loc_j, Neighs_Loc_jN);
                                }

                                Neighs_Loc_jN.Add(Loc_jN);
                            }
                        }


                        int[] T2 = new int[ExtColor.Count];
                        int[][] T3 = new int[ExtColor.Count][];
                        int cnt = 0;
                        foreach (var kv in ExtColor) {
                            T2[cnt] = kv.Key;
                            T3[cnt] = kv.Value.ToArray();
                            cnt++;
                        }

                        ExternalColorListsTmp.Add(T2);
                        ExternalColorListsNeighborsTmp.Add(T3);

                    } else {
                        ExternalColorListsTmp.Add(new int[0]);
                        ExternalColorListsNeighborsTmp.Add(new int[0][]);
                    }


                    // check for loop termination
                    // ==========================
                    Debug.Assert(locColoredCells <= J);
                    int bRunLoc = 0xFFFFFF;
                    if (locColoredCells >= J)
                        bRunLoc = 0;
                    bRun = bRunLoc.MPIMax();
                }

                // store
                // =====
                this.ColorLists = ColorListsTmp.ToArray();
                this.ExternalColorLists = ExternalColorListsTmp.ToArray();
                this.ExternalColorListsNeighbors = ExternalColorListsNeighborsTmp.ToArray();

                // checks
                // ======
#if DEBUG
                {
                    int[] TouchLoc = new int[JE];

                    foreach (int[] _CellList in this.ColorLists) {

                        int[] NeighTouchLoc = new int[JE];

                        foreach (int j in _CellList) {
                            Debug.Assert(j < JE);
                            Debug.Assert(TouchLoc[j] == 0);
                            TouchLoc[j]++;

                            Debug.Assert(NeighTouchLoc[j] == 0);
                            NeighTouchLoc[j]++;

                            foreach (int jn in Neighs[j]) {
                                Debug.Assert(NeighTouchLoc[jn] == 0);
                                NeighTouchLoc[jn]++;
                            }

                        }

                        int[] NeighTouchGlob = NeighTouchLoc.CloneAs();
                        NeighTouchGlob.MPIExchange(gDat);
                        for (int j = J; j < JE; j++) {
                            NeighTouchGlob[j] += NeighTouchLoc[j];
                        }

                        for (int j = 0; j < JE; j++) {
                            Debug.Assert(NeighTouchGlob[j] <= 1);
                        }
                    }

                    int[] TouchGlob = TouchLoc.CloneAs();
                    TouchGlob.MPIExchange(gDat);
                    for (int j = J; j < JE; j++) {
                        TouchGlob[j] += TouchLoc[j];
                    }

                    for (int j = 0; j < JE; j++) {
                        Debug.Assert(TouchGlob[j] == 1);
                    }

                }
#endif

            }


            /// <summary>
            /// computes a approximate linearization of the operator in the form 
            /// \f[
            ///    \mathcal{M} U + \mathcal{B}.
            /// \f]
            /// </summary>
            /// <param name="Matrix">
            /// Output, the approximate Jacobian matrix of the operator is accumulated here
            /// </param>
            /// <param name="AffineOffset">
            /// Output, the operator value in the linearization point
            /// </param>
            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {

                // init locals
                // ===========
                var codMap = XEval.CodomainMapping;
                var domMap = XEval.DomainMapping;
                DGField[] domFields = XEval.DomainFields.Fields.ToArray();
                var U0 = new CoordinateVector(XEval.DomainFields);

                int j0 = XEval.GridData.CellPartitioning.i0;
                int J = XEval.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = XEval.GridData.iLogicalCells.Count;
                int NoOfDomFields = domMap.BasisS.Count;
                int NoOfCodFields = codMap.BasisS.Count;

                int Lout = XEval.CodomainMapping.LocalLength;
                int Lin = domMap.LocalLength;

                int[][] Neighs = XEval.GridData.iLogicalCells.CellNeighbours;

                var lastCodB = codMap.BasisS.Last();
                var lastDomB = domMap.BasisS.Last();

                int NoOfEvals = 0;

                // check Args
                // ==========

                if (!Matrix.RowPartitioning.EqualsPartition(codMap))
                    throw new ArgumentException("Mismatch in matrix row partition.");
                if (!Matrix.ColPartition.EqualsPartition(domMap))
                    throw new ArgumentException("Mismatch in matrix column partition.");
                if (AffineOffset.Count != codMap.LocalLength)
                    throw new ArgumentException("Mismatch in length of affine offset.");

                // evaluate at linearization point
                // ===============================

                double[] F0 = new double[Lout];
                DelParamUpdate(domFields, XEval.Parameters.ToArray());
                XEval.Evaluate(1.0, 0.0, F0);
                AffineOffset.AccV(1.0, F0);
                NoOfEvals++;

                // compute epsilon's
                // =================

                double[] Epsilons = new double[domMap.Ntotal];
                double relEps = this.Eps;
                //double absEps = 1.0e-15; 
                double absEps = 1.0;
                for (int i = 0; i < Lin; i++) {
                    double EpsBase = Math.Abs(U0[i]);
                    if (EpsBase < absEps)
                        EpsBase = absEps;

                    Epsilons[i] = EpsBase * relEps;
                }
                Epsilons.MPIExchange(XEval.GridData);


                // compute directional derivatives
                // ===============================

                double[] U0backup = new double[Lin];
                double[] EvalBuf = new double[Lout];
                if (!domMap.AllBlockSizesEqual)
                    throw new NotSupportedException();
                MultidimensionalArray Buffer = MultidimensionalArray.Create(Lout, domMap.GetBlockLen(domMap.FirstBlock));

                for (int iCellPass = 0; iCellPass < ColorLists.Length; iCellPass++) { // loop over all cell lists...
                    int[] CellList = this.ColorLists[iCellPass];
                    int[] ExtCellList = this.ExternalColorLists[iCellPass];

                    int[] CoordCounter = new int[JE];
                    int[] FieldCounter = new int[JE];

                    int maxNj = 0;
                    foreach (int j in CellList) {
                        int Nj = domMap.GetTotalNoOfCoordinatesPerCell(j);
                        maxNj = Math.Max(Nj, maxNj);
                    }
                    maxNj = maxNj.MPIMax();

                    Buffer.Clear();

                    for (int n = 0; n < maxNj; n++) { // loop over DG coordinates in cell

                        // backup DG coordinates
                        // ---------------------
                        U0backup.SetV(U0);

                        // apply distortions
                        // -----------------
                        int AnyLoc = 0;
                        foreach (int j in CellList) {
                            int iFld = FieldCounter[j];
                            int nFld = CoordCounter[j];
                            if (iFld > NoOfDomFields)
                                continue; // finished with cell 'j'

                            AnyLoc = -1;

                            Debug.Assert(j >= 0 && j < J);
                            int iLoc = domMap.LocalUniqueCoordinateIndex(iFld, j, nFld);
                            Debug.Assert(iLoc >= 0 && iLoc < domMap.LocalLength);

                            double oldVal = domFields[iFld].Coordinates[j, nFld];
                            domFields[iFld].Coordinates[j, nFld] = oldVal + Epsilons[iLoc];
                            Debug.Assert(domFields[iFld].Coordinates[j, nFld] != oldVal);
                        }
                        int AnyGlob = AnyLoc.MPIMin();
                        if (AnyGlob >= 0)
                            break; // finished with entire cell list on all processors

                        // evaluate operator
                        // -------------------
                        EvalBuf.ClearEntries();
                        DelParamUpdate(domFields, XEval.Parameters.ToArray());
                        XEval.Evaluate(1.0, 0.0, EvalBuf);
                        NoOfEvals++;

                        // ------------------------------

                        for (int IntExt = 0; IntExt < 2; IntExt++) {
                            int[] __CellList;
                            switch (IntExt) {
                                case 0: __CellList = CellList; break;
                                case 1: __CellList = ExtCellList; break;
                                default: throw new ApplicationException();
                            }


                            // save results
                            // -------------------------------
                            int cnt = 0;
                            foreach (int _j in __CellList) {
                                int[] Neighs_j; // = Neighs[_j];
                                switch (IntExt) {
                                    case 0: Neighs_j = Neighs[_j]; break;
                                    case 1: Neighs_j = this.ExternalColorListsNeighbors[iCellPass][cnt]; break;
                                    default: throw new ApplicationException();
                                }
                                cnt++;

                                int jCol = _j;

                                int iFldCol = FieldCounter[jCol];
                                int nFldCol = CoordCounter[jCol];
                                if (iFldCol > NoOfDomFields)
                                    continue; // finished with cell

                                int iCol = domMap.LocalUniqueCoordinateIndex(iFldCol, jCol, nFldCol);
                                int i0Col = domMap.LocalUniqueCoordinateIndex(0, jCol, 0);
                                int iRelCol = iCol - i0Col;

                                for (int k = 0; k <= Neighs_j.Length; k++) { // loop over neighbors which are influenced by the distortion
                                    int jRow;
                                    if (k == 0) {
                                        jRow = _j;
                                    } else {
                                        jRow = Neighs_j[k - 1];
                                    }

                                    if (jRow >= J) {
                                        continue; // external cell; should be treated on other proc.
                                    }
                                    int i0Row = codMap.LocalUniqueCoordinateIndex(0, jRow, 0);
                                    int NoOfRows = codMap.GetBlockLen(jRow);

                                    for (int iRelRow = 0; iRelRow < NoOfRows; iRelRow++) {
                                        int iRow = i0Row + iRelRow;

                                        double u1 = EvalBuf[iRow];
                                        double u0 = F0[iRow];
                                        double h = Epsilons[iCol];

                                        double diff = (u1 - u0) / h;
                                        Buffer[iRow, iRelCol] = diff;

                                    }
                                }
                            }

                            // increase counters
                            // ------------------
                            foreach (int j in __CellList) {
                                int iFld = FieldCounter[j];
                                if (iFld > NoOfDomFields)
                                    continue; // finished with cell 'j'

                                int Nj = domMap.BasisS[iFld].GetLength(j);
                                CoordCounter[j]++;
                                if (CoordCounter[j] >= Nj) {
                                    CoordCounter[j] = 0;
                                    FieldCounter[j]++;
                                }

                            }
                        }

                        // restore original DG coordinates
                        // -------------------------------
                        U0.SetV(U0backup);
                    }

                    // save to matrix
                    // --------------

                    for (int IntExt = 0; IntExt < 2; IntExt++) {
                        int[] __CellList;
                        switch (IntExt) {
                            case 0: __CellList = CellList; break;
                            case 1: __CellList = ExtCellList; break;
                            default: throw new ApplicationException();
                        }

                        int cnt = 0;
                        foreach (int _j in __CellList) {
                            int[] Neighs_j; // = Neighs[_j];
                            switch (IntExt) {
                                case 0: Neighs_j = Neighs[_j]; break;
                                case 1: Neighs_j = this.ExternalColorListsNeighbors[iCellPass][cnt]; break;
                                default: throw new ApplicationException();
                            }
                            cnt++;

                            int jCol = _j;
                            int i0Col = domMap.LocalUniqueCoordinateIndex(0, jCol, 0);
                            int iECol = domMap.LocalUniqueCoordinateIndex(NoOfDomFields - 1, jCol, lastDomB.GetLength(jCol) - 1);

                            for (int k = 0; k <= Neighs_j.Length; k++) { // loop over neighbors which are influenced by the distortion
                                int jRow;
                                if (k == 0) {
                                    jRow = _j;
                                } else {
                                    jRow = Neighs_j[k - 1];
                                }

                                if (jRow >= J)
                                    continue; // external cell; should be treated on other proc.


                                int i0Row = domMap.LocalUniqueCoordinateIndex(0, jRow, 0);
                                int iERow = domMap.LocalUniqueCoordinateIndex(NoOfCodFields - 1, jRow, lastCodB.GetLength(jRow) - 1);

                                var Block = Buffer.ExtractSubArrayShallow(new int[] { i0Row, 0 }, new int[] { iERow, iECol - i0Col });

                                Matrix.AccBlock(i0Row + codMap.i0,
                                    //i0Col + domMap.i0, 
                                    domMap.GlobalUniqueCoordinateIndex(0, jCol, 0),
                                    1.0, Block);
                            }
                        }

                    }
                }
                // restore original state before return
                // ====================================
                U0.SetV(U0backup);
                DelParamUpdate(domFields, XEval.Parameters.ToArray());
                //Console.WriteLine("Total number of evaluations: " + NoOfEvals);

                // correct the affine offset
                // =========================

                // actually, the Jacobian is the approx M*(U-U0) + b
                // but we want                          M*U + (b - M*U0)
                // (in this fashion, this approximation also works with BDF schemes *in the same fashion* as other linearizations
                Matrix.SpMV(-1.0, U0, 1.0, AffineOffset);
            }

            /// <summary>
            /// Evaluation at the linearization point
            /// </summary>
            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                //int Lout = Eval.CodomainMapping.LocalLength;

                //double[] F0 = new double[Lout];
                //DelParamUpdate(Eval.DomainFields.Fields.ToArray(), Eval.Parameters.ToArray());
                //Eval.Evaluate(1.0, 0.0, F0);
                //AffineOffset.AccV(1.0, F0);

                BlockMsrMatrix dummy = new BlockMsrMatrix(this.CodomainMapping, this.DomainMapping);
                this.ComputeMatrix(dummy, AffineOffset);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="sgrd"></param>
            /// <param name="subGridBoundaryTreatment"></param>
            public void ActivateSubgridBoundary(CellMask sgrd, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                if (sgrd != null && sgrd.MaskType != MaskType.Logical)
                    throw new ArgumentException("expecting logical mask");
                XEval.ActivateSubgridBoundary(sgrd, subGridBoundaryTreatment);
            }


        }


    }
}
