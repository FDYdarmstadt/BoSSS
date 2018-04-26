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
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// An operator which is specialized in XDG Fields, i.e.
    /// it can have components which couple the phases.
    /// </summary>
    public class XSpatialOperator : SpatialOperator {
        
        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(IList{string},IList{string},Func{int[],int[],int[],int})"/>
        /// </summary>
        public XSpatialOperator(IList<string> domVar, IList<string> codVar, Func<int[], int[], int[], int> QuadOrderFunc)
            : base(domVar, codVar, QuadOrderFunc) {
            ConstructorCommon();
        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(IList{string},IList{string},IList{string},Func{int[],int[],int[],int})"/>
        /// </summary>
        public XSpatialOperator(IList<string> domVar, IList<string> paramVar, IList<string> codVar, Func<int[], int[], int[], int> QuadOrderFunc)
            : base(domVar, paramVar, codVar, QuadOrderFunc) {
            ConstructorCommon();
        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(int,int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XSpatialOperator(int NoOfDomFields, int NoOfParameters, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
            : base(NoOfDomFields, NoOfParameters, NoOfCodomFields, QuadOrderFunc, __varnames) {
            ConstructorCommon();
        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XSpatialOperator(int NoOfDomFields, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
            : base(NoOfDomFields, NoOfCodomFields, QuadOrderFunc, __varnames) {
            ConstructorCommon();
        }

        void ConstructorCommon() {
            GhostEdgesOperator = new FixedOrder_SpatialOperator(base.DomainVar, base.ParameterVar, base.CodomainVar);
            SurfaceElementOperator = new FixedOrder_SpatialOperator(base.DomainVar, base.ParameterVar, base.CodomainVar);
        }

        /// <summary>
        /// Operator working on ghost edges, see e.g.
        /// @article{burman_ghost_2010,                                                             
        ///     title = {Ghost penalty},                                                        
        ///     volume = {348},                                                                 
        ///     issn = {1631073X},                                                              
        ///     url = {http://linkinghub.elsevier.com/retrieve/pii/S1631073X10002827},          
        ///     doi = {10.1016/j.crma.2010.10.006},                                             
        ///     language = {en},                                                                
        ///     number = {21-22},                                                               
        ///     urldate = {2015-09-11},                                                         
        ///     journal = {Comptes Rendus Mathematique},                                        
        ///     author = {Burman, Erik},                                                        
        ///     month = nov,                                                                    
        ///     year = {2010},                                                                  
        ///     pages = {1217--1220}
        /// }                                                
        /// </summary>
        public SpatialOperator GhostEdgesOperator {
            get;
            private set;
        }

        /// <summary>
        /// Non-coupling surface terms; originally intended to implement the flux-form of the surface tension.
        /// </summary>
        public SpatialOperator SurfaceElementOperator {
            get;
            private set;
        }

        /// <summary>
        /// not supported for <see cref="XSpatialOperator"/>, use 
        /// </summary>
        public override IEvaluatorNonLin GetEvaluatorEx(IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap, EdgeQuadratureScheme edgeQrCtx = null, CellQuadratureScheme volQrCtx = null) {
            throw new NotSupportedException("Use specific implementation for XSpatialOperator.");
        }

        /// <summary>
        /// not supported for <see cref="XSpatialOperator"/>, use specific implementation
        /// <see cref="GetMatrixBuilder(LevelSetTracker, UnsetteledCoordinateMapping, IList{DGField}, UnsetteledCoordinateMapping, IDictionary{SpeciesId, QrSchemPair})"/>.
        /// </summary>
        public override IEvaluatorLinear GetMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap, EdgeQuadratureScheme edgeQrCtx = null, CellQuadratureScheme volQrCtx = null) {
            throw new NotSupportedException("Use specific implementation for XSpatialOperator.");
        }
         
        /// <summary>
        /// nix für kleine Kinder
        /// </summary>
        IEvaluatorLinear GetMatrixBuilderBase(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap, EdgeQuadratureScheme edgeQrCtx = null, CellQuadratureScheme volQrCtx = null) {
            return base.GetMatrixBuilder(DomainVarMap, ParameterMap, CodomainVarMap, edgeQrCtx, volQrCtx);
        }


        
        /// <summary>
        /// edge and cell scheme for a certain species
        /// </summary>
        public struct QrSchemPair {

            /// <summary>
            /// if null, a default scheme is chosen for the integration
            /// </summary>
            public EdgeQuadratureScheme EdgeScheme;


            /// <summary>
            /// if null, a default scheme is chosen for the integration
            /// </summary>
            public CellQuadratureScheme CellScheme;
        }

        
        /// <summary>
        /// create a matrix from this operator
        /// </summary>
        public XEvaluatorLinear GetMatrixBuilder(
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap, 
            IDictionary<SpeciesId, QrSchemPair> SpeciesSchemes = null
            ) {

            return new XEvaluatorLinear(this, lsTrk, DomainVarMap, ParameterMap, CodomainVarMap, 
                1, // based on actual level-set tracker state
                SpeciesSchemes);
        }

        /// <summary>
        /// explicit evaluation of the operator
        /// </summary>
        public XEvaluatorNonlin GetEvaluatorEx(
            LevelSetTracker lsTrk,
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap, 
             IDictionary<SpeciesId, QrSchemPair> SpeciesSchemes = null
            ) {
            return new XEvaluatorNonlin(this, lsTrk,
                new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap,
                1,
                SpeciesSchemes);
        }


        /// <summary>
        /// Assembly of matrices for XDG operators
        /// </summary>
        public class XEvaluatorLinear : XEvaluatorBase, IEvaluatorLinear {

            internal XEvaluatorLinear(XSpatialOperator ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory,
                IDictionary<SpeciesId, QrSchemPair> __SpeciesSchemes) :
                base(ownr, lsTrk, DomainVarMap, ParameterMap, CodomainVarMap, TrackerHistory, __SpeciesSchemes) //
            {
                using(var tr = new FuncTrace()) {
                    base.MPITtransceive = true;
                    
                    
                   
                }
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cellScheme, EdgeQuadratureScheme edgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params) {
                var BulkMtxBuilder = base.m_Xowner.GetMatrixBuilderBase(DomainFrame.FrameMap, Params, CodomFrame.FrameMap,
                                edgeScheme, cellScheme);
                Debug.Assert(((EvaluatorBase)BulkMtxBuilder).order == base.order);
                BulkMtxBuilder.MPITtransceive = false;
                SpeciesBulkMtxBuilder.Add(SpeciesId, BulkMtxBuilder);
            }

            /// <summary>
            /// creates a matrix builder
            /// </summary>
            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme nullvolumeScheme, EdgeQuadratureScheme ghostEdgeScheme, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params) {
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
            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme SurfaceElement_volume, EdgeQuadratureScheme SurfaceElement_Edge, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params) {
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
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                using(var tr = new FuncTrace()) {
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

                    if(OnlyAffine == false) {
                        if(!Matrix.RowPartitioning.Equals(base.CodomainMapping))
                            throw new ArgumentException("wrong number of columns in matrix.", "Matrix");
                        if(!Matrix.ColPartition.Equals(base.DomainMapping))
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
                    if(trx != null) {
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
                        foreach(var SpeciesId in ReqSpecies) {
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

                            foreach(var SpeciesBuilder in new[] { SpeciesBulkMtxBuilder, SpeciesGhostEdgeBuilder, SpeciesSurfElmBuilder }) {
                                
                                if(SpeciesBuilder.ContainsKey(SpeciesId)) {
                                    
                                    var builder = SpeciesBulkMtxBuilder[SpeciesId];
                                    builder.OperatorCoefficients = this.SpeciesOperatorCoefficients[SpeciesId];
                                    NotifySpecies(builder.Owner, this.m_lsTrk, SpeciesId);
                                        
                                    if(trx != null) {
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

                    using(new BlockTrace("surface_integration", tr)) {
                        foreach(var tt in this.CouplingRules) {
                            int iLevSet = tt.Item1;
                            var SpeciesA = tt.Item2;
                            var SpeciesB = tt.Item3;
                            var rule = tt.Item4;


                            var MtxBuilder = new LECQuadratureLevelSet<M, V>(GridDat,
                                                             m_Xowner,
                                                             OnlyAffine ? default(M) : Matrix, AffineOffset,
                                                             CodomainMapping, Parameters, DomainMapping,
                                                             lsTrk, iLevSet, new Tuple<SpeciesId, SpeciesId>(SpeciesA, SpeciesB),
                                                             rule);
                            MtxBuilder.time = time;
                            UpdateLevelSetCoefficients(this.SpeciesOperatorCoefficients[SpeciesA], this.SpeciesOperatorCoefficients[SpeciesB]);
                            MtxBuilder.Execute();

#if DEBUG
                            if(Matrix != null && OnlyAffine == false)
                                Matrix.CheckForNanOrInfM();
                            if(AffineOffset != null)
                                GenericBlas.CheckForNanOrInfV(AffineOffset);
#endif

                        }


                        /*
                        if(m_Xowner.ContainesComponentType(typeof(ILevelSetForm))) {


                            var AllSpc = lsTrk.SpeciesIdS;

                            // loop over all possible pairs of species
                            for(int iSpcA = 0; iSpcA < AllSpc.Count; iSpcA++) {
                                var SpeciesA = AllSpc[iSpcA];
                                var SpeciesADom = lsTrk.Regions.GetSpeciesMask(SpeciesA);
                                if(SpeciesADom.NoOfItemsLocally <= 0)
                                    continue;

                                int _iSpcA = Array.IndexOf(ReqSpecies, SpeciesA);

                                for(int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                                    var SpeciesB = AllSpc[iSpcB];

                                    int _iSpcB = Array.IndexOf(ReqSpecies, SpeciesB);
                                    if(_iSpcA < 0 && _iSpcB < 0)
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
                                    for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {

                                        var LsDom = lsTrk.Regions.GetCutCellMask4LevSet(iLevSet);
                                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                        // Check removed since it can cause parallel problems
                                        //if (IntegrationDom.NoOfItemsLocally > 0) {

                                        ICompositeQuadRule<QuadRule> rule;
                                        using(new BlockTrace("QuadRule-compilation", tr)) {
                                            CellQuadratureScheme SurfIntegration = SchemeHelper.GetLevelSetquadScheme(iLevSet, IntegrationDom);
                                            rule = SurfIntegration.Compile(GridDat, order);

                                            if(ruleDiagnosis) {
                                                rule.SumOfWeightsToTextFileVolume(GridDat, string.Format("Levset_{0}.csv", iLevSet));
                                            }
                                        }


                                        //}
                                    }
                                }
                            }
                        }*/
                    }

                    // allow all processes to catch up
                    // -------------------------------
                    if(trx != null) {
                        trx.TransceiveFinish();
                        trx = null;
                    }
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                }
            }

            /// <summary>
            /// nix
            /// </summary>
            protected override void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule) {
            }
        }


        public class XEvaluatorNonlin : XEvaluatorBase, IEvaluatorNonLin {

            internal XEvaluatorNonlin(XSpatialOperator ownr,
                LevelSetTracker lsTrk,
                CoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory,
                IDictionary<SpeciesId, QrSchemPair> __SpeciesSchemes) :
                base(ownr, lsTrk, DomainVarMap, ParameterMap, CodomainVarMap, TrackerHistory, __SpeciesSchemes) {
                this.DomainFields = DomainVarMap;
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
                throw new NotImplementedException();
            }

            protected override void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params) {
                throw new NotImplementedException();
            }

            protected override void ctorLevSetFormIntegrator(int iLevSet, SpeciesId SpeciesA, SpeciesId SpeciesB, ICompositeQuadRule<QuadRule> rule) {
                throw new NotImplementedException();
            }

            protected override void ctorSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params) {
                throw new NotImplementedException();
            }

            protected override void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params) {
                throw new NotImplementedException();
            }

            
        }


        /// <summary>
        /// Common nase-class for linear/nonlinear evaluation of <see cref="XSpatialOperator"/>
        /// </summary>
        abstract public class XEvaluatorBase : SpatialOperator.EvaluatorBase {

            internal XEvaluatorBase(XSpatialOperator ownr,
                LevelSetTracker lsTrk,
                UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
                int TrackerHistory,
                IDictionary<SpeciesId, QrSchemPair> __SpeciesSchemes) :
                base(ownr, DomainVarMap, ParameterMap, CodomainVarMap) //
            {
                using(var tr = new FuncTrace()) {
                    if(!object.ReferenceEquals(base.GridData, lsTrk.GridDat))
                        throw new ArgumentException("grid data mismatch");
                    m_lsTrk = lsTrk;
                    m_Xowner = ownr;
                    SpeciesSchemes = __SpeciesSchemes;
                    ReqSpecies = SpeciesSchemes.Keys.ToArray();
                    base.MPITtransceive = true;

                    var SchemeHelper = lsTrk.GetXDGSpaceMetrics(ReqSpecies, order, TrackerHistory).XQuadSchemeHelper;
                    var TrackerRegions = lsTrk.RegionsHistory[TrackerHistory];

                    ((FixedOrder_SpatialOperator)(m_Xowner.GhostEdgesOperator)).m_Order = base.order;
                    ((FixedOrder_SpatialOperator)(m_Xowner.SurfaceElementOperator)).m_Order = base.order;
                    tr.Info("XSpatialOperator.ComputeMatrixEx quad order: " + order);


                    // compile quadrature rules & create matrix builders for each species 
                    // ------------------------------------------------------------------
                    foreach(var SpeciesId in ReqSpecies) {
                        int iSpecies = Array.IndexOf(ReqSpecies, SpeciesId);

                        // parameters for species
                        // ----------------------
                        DGField[] Params = (from f in (Parameters ?? new DGField[0])
                                            select ((f is XDGField) ? ((XDGField)f).GetSpeciesShadowField(SpeciesId) : f)).ToArray<DGField>();
                        SpeciesParams.Add(SpeciesId, Params);

                        // species frames
                        // --------------

                        var CodomFrame = new FrameBase(TrackerRegions, SpeciesId, base.CodomainMapping, false);
                        var DomainFrame = new FrameBase(TrackerRegions, SpeciesId, base.DomainMapping, true);
                        SpeciesCodomFrame.Add(SpeciesId, CodomFrame);
                        SpeciesDomainFrame.Add(SpeciesId, DomainFrame);

                        // quadrature rules
                        // ----------------
                        EdgeQuadratureScheme edgeScheme;
                        CellQuadratureScheme cellScheme;
                        using(new BlockTrace("QuadRule-compilation", tr)) {

                            var qrSchemes = SpeciesSchemes[SpeciesId];

                            //bool AssembleOnFullGrid = (SubGrid == null);
                            if(qrSchemes.EdgeScheme == null) {
                                edgeScheme = SchemeHelper.GetEdgeQuadScheme(SpeciesId);// AssembleOnFullGrid, SubGridEdgeMask);
                            } else {
                                edgeScheme = qrSchemes.EdgeScheme;
                                //throw new NotSupportedException();
                            }
                            if(qrSchemes.CellScheme == null) {
                                cellScheme = SchemeHelper.GetVolumeQuadScheme(SpeciesId);//, AssembleOnFullGrid, SubGridCellMask);
                            } else {
                                cellScheme = qrSchemes.CellScheme;
                                //throw new NotSupportedException();
                            }
                        }

                        if(ruleDiagnosis) {
                            var edgeRule = edgeScheme.Compile(base.GridData, base.order);
                            var volRule = cellScheme.Compile(base.GridData, base.order);
                            edgeRule.SumOfWeightsToTextFileEdge(base.GridData, string.Format("PhysEdge_{0}.csv", lsTrk.GetSpeciesName(SpeciesId)));
                            volRule.SumOfWeightsToTextFileVolume(base.GridData, string.Format("Volume_{0}.csv", lsTrk.GetSpeciesName(SpeciesId)));
                        }


                        if(m_Xowner.TotalNoOfComponents > 0) {
                            ctorSpeciesIntegrator(SpeciesId, cellScheme, edgeScheme, DomainFrame, CodomFrame, Params);
                        }

                        if(m_Xowner.GhostEdgesOperator.TotalNoOfComponents > 0) {
                            CellQuadratureScheme nullvolumeScheme = new CellQuadratureScheme(false, CellMask.GetEmptyMask(GridData));
                            EdgeQuadratureScheme ghostEdgeScheme = null;
                            ghostEdgeScheme = SchemeHelper.GetEdgeGhostScheme(SpeciesId);//, SubGridEdgeMask);
                            ctorGhostSpeciesIntegrator(SpeciesId, nullvolumeScheme, ghostEdgeScheme, DomainFrame, CodomFrame, Params);
                        }

                        if(m_Xowner.SurfaceElementOperator.TotalNoOfComponents > 0) {
                            EdgeQuadratureScheme SurfaceElement_Edge;
                            CellQuadratureScheme SurfaceElement_volume;
                            SurfaceElement_Edge = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(SpeciesId);
                            SurfaceElement_volume = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(SpeciesId);
                            ctorSurfaceElementSpeciesIntegrator(SpeciesId, SurfaceElement_volume, SurfaceElement_Edge, DomainFrame, CodomFrame, Params);
                        }
                    }

                    // coupling terms
                    // --------------

                    using(new BlockTrace("surface_integration", tr)) {
                        if(m_Xowner.ContainesComponentType(typeof(ILevelSetForm))) {


                            var AllSpc = lsTrk.SpeciesIdS;

                            // loop over all possible pairs of species
                            for(int iSpcA = 0; iSpcA < AllSpc.Count; iSpcA++) {
                                var SpeciesA = AllSpc[iSpcA];
                                var SpeciesADom = lsTrk.Regions.GetSpeciesMask(SpeciesA);
                                if(SpeciesADom.NoOfItemsLocally <= 0)
                                    continue;

                                int _iSpcA = Array.IndexOf(ReqSpecies, SpeciesA);

                                for(int iSpcB = iSpcA + 1; iSpcB < AllSpc.Count; iSpcB++) {
                                    var SpeciesB = AllSpc[iSpcB];

                                    int _iSpcB = Array.IndexOf(ReqSpecies, SpeciesB);
                                    if(_iSpcA < 0 && _iSpcB < 0)
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
                                    for(int iLevSet = 0; iLevSet < NoOfLs; iLevSet++) {

                                        var LsDom = lsTrk.Regions.GetCutCellMask4LevSet(iLevSet);
                                        var IntegrationDom = LsDom.Intersect(SpeciesCommonDom);

                                        // Check removed since it can cause parallel problems
                                        //if (IntegrationDom.NoOfItemsLocally > 0) {

                                        ICompositeQuadRule<QuadRule> rule;
                                        using(new BlockTrace("QuadRule-compilation", tr)) {
                                            CellQuadratureScheme SurfIntegration = SchemeHelper.GetLevelSetquadScheme(iLevSet, IntegrationDom);
                                            rule = SurfIntegration.Compile(GridData, order);

                                            if(ruleDiagnosis) {
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
                    foreach(var SpeciesId in ReqSpecies) {
                        this.SpeciesOperatorCoefficients.Add(SpeciesId,
                            new CoefficientSet() {
                                CellLengthScales = null, // ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.cj,
                                EdgeLengthScales = null, //((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge,
                                UserDefinedValues = new Dictionary<string, object>()
                            });
                    }
                    // m_OperatorCoefficients = new CoefficientSet() {
                    //    CellLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.cj,
                    //    EdgeLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge,
                    //    UserDefinedValues = new Dictionary<string, object>()
                    //};
                }
            }

            /// <summary>
            /// create integrator for bulk phase
            /// </summary>
            abstract protected void ctorSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params);

            /// <summary>
            /// Create integrator for <see cref="XSpatialOperator.GhostEdgesOperator"/>
            /// </summary>
            abstract protected void ctorGhostSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params);
 
            /// <summary>
            /// Create integrator for <see cref="XSpatialOperator.SurfaceElementOperator"/>
            /// </summary>
            abstract protected void ctorSurfaceElementSpeciesIntegrator(SpeciesId SpeciesId, CellQuadratureScheme cqs, EdgeQuadratureScheme eqs, FrameBase DomainFrame, FrameBase CodomFrame, DGField[] Params);
            
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
            protected XSpatialOperator m_Xowner;

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
                for(int iCod = 0; iCod < CodDGdeg.Length; iCod++) {
                    var comps = Owner.EquationComponents[CodNames[iCod]];
                    foreach(var c in comps) {
                        if(c is ILevelSetEquationComponentCoefficient) {
                            var ce = c as ILevelSetEquationComponentCoefficient;

                            int[] DomDGdeg_cd = new int[ce.ArgumentOrdering.Count];
                            for(int i = 0; i < DomDGdeg_cd.Length; i++) {
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
                    foreach(var c in comps) {
                        if(c is IEquationComponentSpeciesNotification) {
                            var ce = c as IEquationComponentSpeciesNotification;
                            ce.SetParameter(sNmn, id);
                        }
                    }
                }
            }


            /// <summary>
            /// Not supported, use <see cref="SpeciesOperatorCoefficients"/>.
            /// </summary>
            public override CoefficientSet OperatorCoefficients {
                get {
                    throw new NotSupportedException("Use per-species implementation.");
                }
                set {
                    throw new NotSupportedException("Use per-species implementation.");
                }
            }

            /// <summary>
            /// Operator coefficients for each species
            /// </summary>
            public Dictionary<SpeciesId, CoefficientSet> SpeciesOperatorCoefficients {
                get;
            }

        }



        internal class SpeciesFrameVector<V> : IList<double>
            where V : IList<double> {

            /// <summary>
            /// ctor.
            /// </summary>
            /// <param name="lsTrk_Regions">see <see cref="LevelSetTracker.Regions"/></param>
            /// <param name="spcId">species which should be framed</param>
            /// <param name="Full">the vector that should be framed</param>
            /// <param name="FullMap"></param>
            public SpeciesFrameVector(LevelSetTracker.LevelSetRegions lsTrk_Regions, SpeciesId spcId, V Full, UnsetteledCoordinateMapping FullMap)
                : this(Full, new FrameBase(lsTrk_Regions, spcId, FullMap, false)) {
            }


            /// <summary>
            /// ctor.
            /// </summary>
            public SpeciesFrameVector(V FullVec, FrameBase __fr) {
                fr = __fr;
                m_Full = FullVec;
            }

            FrameBase fr;

            V m_Full;


            /// <summary>
            /// Mapping for the framed vector.
            /// </summary>
            public UnsetteledCoordinateMapping Mapping {
                get {
                    return fr.FrameMap;
                }
            }

            /// <summary>
            ///  not supported.
            /// </summary>
            public int IndexOf(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public void Insert(int index, double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public void RemoveAt(int index) {
                throw new NotSupportedException();
            }


            /// <summary>
            /// set/get one element
            /// </summary>
            public double this[int i] {
                get {
                    int idx = fr.Frame2Full_Loc(i);
                    return m_Full[idx];
                }
                set {
                    int idx = fr.Frame2Full_Loc(i);
                    m_Full[idx] = value;
                }
            }

            /// <summary>
            /// not supported
            /// </summary>
            public void Add(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// Sets all entries to 0.0.
            /// </summary>
            public void Clear() {
                int L = this.Count;
                for (int i = 0; i < L; i++)
                    this[i] = 0;
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public bool Contains(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// copy to array.
            /// </summary>
            public void CopyTo(double[] array, int arrayIndex) {
                int L = this.Count;
                for (int i = 0; i < L; i++)
                    array[i + arrayIndex] = this[i];
            }

            /// <summary>
            /// number of elements
            /// </summary>
            public int Count {
                get {
                    return Mapping.LocalLength;
                }
            }

            /// <summary>
            /// depends on framed object
            /// </summary>
            public bool IsReadOnly {
                get {
                    return m_Full.IsReadOnly;
                }
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public bool Remove(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// %
            /// </summary>
            public IEnumerator<double> GetEnumerator() {
                throw new NotImplementedException("leck mi - Wolfgang Amadeus Mozart, Köchelverzeichnis KV233 bzw. 382c und d");

            }

            /// <summary>
            /// %
            /// </summary>
            IEnumerator IEnumerable.GetEnumerator() {
                throw new NotImplementedException("leck mi.");
            }
        }

        /// <summary>
        /// This class acts as a frame fore some other matrix, and presents only those entries which are associated
        /// with a given species.
        /// </summary>
        public class SpeciesFrameMatrix<M> : IMutableMatrixEx
            where M : IMutableMatrixEx {

            /// <summary>
            /// Pseudo-implementation.
            /// </summary>
            public object Clone() {
                throw new NotImplementedException();
            }

            /// <summary>
            /// Not Implemented.
            /// </summary>
            public void Clear() {
                throw new NotImplementedException(); 
            }

            /// <summary>
            /// ctor.
            /// </summary>
            /// <param name="full">the full operator matrix, from which the species <paramref name="spcId"/> should e framed (extracted)</param>
            /// <param name="lsTrk_regions">see <see cref="LevelSetTracker.Regions"/></param>
            /// <param name="spcId">the species of interest</param>
            /// <param name="fullMapRow">row mapping for the operator matrix <paramref name="full"/></param>
            /// <param name="fullMapCol">column mapping for the operator matrix <paramref name="full"/></param>
            public SpeciesFrameMatrix(M full, LevelSetTracker.LevelSetRegions lsTrk_regions, SpeciesId spcId, UnsetteledCoordinateMapping fullMapRow, UnsetteledCoordinateMapping fullMapCol)
                : this(full, new FrameBase(lsTrk_regions, spcId, fullMapRow, false), new FrameBase(lsTrk_regions, spcId, fullMapCol, true)) {
            }

            LevelSetTracker.LevelSetRegions m_LsTrk_regions;

            /// <summary>
            /// ctor.
            /// </summary>
            public SpeciesFrameMatrix(M full, FrameBase __RowFrame, FrameBase __ColFrame) {
                m_full = full;
                Debug.Assert(object.ReferenceEquals(__RowFrame.Regions, __ColFrame.Regions));
                Debug.Assert(__RowFrame.Species == __ColFrame.Species);
                RowFrame = __RowFrame;
                ColFrame = __ColFrame;
                m_LsTrk_regions = __RowFrame.Regions;

#if DEBUG
                var grdDat = RowFrame.FullMap.BasisS.First().GridDat;
                int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = grdDat.iLogicalCells.NoOfCells;
                var spc = RowFrame.Species;
                Debug.Assert(RowFrame.Species.Equals(ColFrame.Species));
                var lsTrk = m_LsTrk_regions;
                Basis[] RowBase = RowMapping.BasisS.ToArray();
                Basis[] ColBase = ColMapping.BasisS.ToArray();

                var _AvailableRowIdx = new List<int>();
                for (int j = 0; j < J; j++) {
                    int iSpc = m_LsTrk_regions.GetSpeciesIndex(spc, j);

                    if (iSpc >= 0) {
                        for (int k = 0; k < RowBase.Length; k++) {
                            int N = RowBase[k].GetLength(j);
                            for (int n = 0; n < N; n++) {
                                int iRow = RowFrame.FrameMap.GlobalUniqueCoordinateIndex(k, j, n);
                                _AvailableRowIdx.Add(iRow);
                                Debug.Assert(RowFrame.Frame2Full(iRow) >= 0);
                            }
                        }
                    }
                }
                this.AvailableRowIdx = _AvailableRowIdx.ToArray();


                var _AvailableColIdx = new List<int>();
                for (int j = 0; j < JE; j++) {
                    int iSpc = m_LsTrk_regions.GetSpeciesIndex(spc, j);

                    if (iSpc >= 0) {
                        for (int k = 0; k < ColBase.Length; k++) {
                            int N = ColBase[k].GetLength(j);
                            for (int n = 0; n < N; n++) {
                                int iCol = ColFrame.FrameMap.GlobalUniqueCoordinateIndex(k, j, n);
                                _AvailableColIdx.Add(iCol);
                                Debug.Assert(ColFrame.Frame2Full(iCol) >= 0);
                            }
                        }
                    }
                }
                this.AvailableColIdx = _AvailableColIdx.ToArray();
#endif
            }

            FrameBase RowFrame;
            FrameBase ColFrame;


#if DEBUG
            int[] AvailableRowIdx;
            int[] AvailableColIdx;
#endif

            /// <summary>
            /// row mapping for the framed matrix
            /// </summary>
            public UnsetteledCoordinateMapping RowMapping {
                get {
                    return RowFrame.FrameMap;
                }
            }

            /// <summary>
            /// mpi comm
            /// </summary>
            public MPI_Comm MPI_Comm {
                get {
                    return RowPartitioning.MPI_Comm;
                }
            }

            /// <summary>
            /// column mapping for the framed matrix
            /// </summary>
            public UnsetteledCoordinateMapping ColMapping {
                get {
                    return ColFrame.FrameMap;
                }
            }

            M m_full;


            /// <summary>
            /// get a whole bunch of elements at once
            /// </summary>
            public double[] GetValues(int RowIndex, int[] ColumnIndices) {
                int L = ColumnIndices.Length;
                double[] ret = new double[L];
                for (int j = 0; j < L; j++) {
                    ret[j] = this[RowIndex, ColumnIndices[j]];
                }
                return ret;
            }

            /// <summary>
            /// set a whole bunch of elements at once
            /// </summary>
            public void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues) {
                if (ColumnIndices.Length != newValues.Length)
                    throw new ArgumentException();

                int L = ColumnIndices.Length;
                for (int j = 0; j < L; j++) {
                    this[RowIndex, ColumnIndices[j]] = newValues[j];
                }
            }


            internal int iRowFrame2Full(int iFrame) {
                int iFull = RowFrame.Frame2Full(iFrame);
                return iFull;
            }

            //int iRowFull2Frame(int iFull) {
            //    return RowFrame.Full2Frame(iFull);
            //}

            internal int iColFrame2Full(int iFrame) {
                int iFull = ColFrame.Frame2Full(iFrame);
                return iFull;
            }

            //int iColFull2Frame(int iFull) {
            //    return ColFrame.Full2Frame(iFull);
            //}



            /// <summary>
            /// set/get a specific entry
            /// </summary>
            /// <param name="i">row index</param>
            /// <param name="j">column index</param>
            public double this[int i, int j] {
                get {
                    return m_full[iRowFrame2Full(i), iColFrame2Full(j)];
                }
                set {
                    m_full[iRowFrame2Full(i), iColFrame2Full(j)] = value;
                }
            }

            /// <summary>
            /// Accumulates <paramref name="Block"/>*<paramref name="alpha"/> to this matrix,
            /// at the row/column offset <paramref name="i0"/> resp. <paramref name="j0"/>.
            /// </summary>
            /// <param name="i0">Row offset.</param>
            /// <param name="j0">Column offset.</param>
            /// <param name="alpha">Scaling factor for the accumulation operation.</param>
            /// <param name="Block">Block to accumulate.</param>
            public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block) {
                this.AccBlock(i0, j0, alpha, Block, 1.0);
            }


            /// <summary>
            /// Accumulates a block of entries to this matrix.
            /// </summary>
            /// <param name="i0">Row index offset.</param>
            /// <param name="j0">Column index offset.</param>
            /// <param name="alpha">Scaling factor for the accumulation.</param>
            /// <param name="Block">Block to add.</param>
            /// <param param name="beta">pre-scaling</param>
            public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block, double beta) {
                if (Block.Dimension != 2)
                    throw new ArgumentException();
                int I = Block.NoOfRows;
                int J = Block.NoOfCols;

                //for (int i = 0; i < I; i++)
                //    for (int j = 0; j < J; j++)
                //        this[i0 + i, j0 + j] += alpha * Block[i, j];


                if (I <= 0 || J <= 0)
                    return;

                int NoIBlk = 1;
                int[] i0S = new int[I];
                int[] i0T = new int[I];
                int[] iLT = new int[I];
                i0T[0] = iRowFrame2Full(i0);
                iLT[0] = 1;

                int NoJBlk = 1;
                int[] j0S = new int[J];
                int[] j0T = new int[J];
                int[] jLT = new int[J];
                j0T[0] = iColFrame2Full(j0);
                jLT[0] = 1;

                for (int i = 1; i < I; i++) {
                    int iT = iRowFrame2Full(i0 + i);
                    if (iT == i0T[NoIBlk - 1] + iLT[NoIBlk - 1]) {
                        iLT[NoIBlk - 1]++;
                    } else {
                        i0T[NoIBlk] = iT;
                        iLT[NoIBlk] = 1;
                        i0S[NoIBlk] = i;
                        NoIBlk++;
                    }
                }

                for (int j = 1; j < J; j++) {
                    int jT = iColFrame2Full(j0 + j);
                    if (jT == j0T[NoJBlk - 1] + jLT[NoJBlk - 1]) {
                        jLT[NoJBlk - 1]++;
                    } else {
                        j0T[NoJBlk] = jT;
                        jLT[NoJBlk] = 1;
                        j0S[NoJBlk] = j;
                        NoJBlk++;
                    }
                }


                for (int iBlk = 0; iBlk < NoIBlk; iBlk++) {
                    for (int jBlk = 0; jBlk < NoJBlk; jBlk++) {
                        var SubBlock = Block.ExtractSubArrayShallow(new int[] { i0S[iBlk], j0S[jBlk] }, new int[] { i0S[iBlk] + iLT[iBlk] - 1, j0S[jBlk] + jLT[jBlk] - 1 });
                        //double SubLinf = SubBlock.AbsSum();
                        //if(SubLinf > 0)
                        m_full.AccBlock(i0T[iBlk], j0T[jBlk], alpha, SubBlock, beta);
                    }
                }

            }

            /// <summary>
            /// depends on the framed matrix
            /// </summary>
            public bool OccupationMutable {
                get {
                    return m_full.OccupationMutable;
                }
            }

            /// <summary>
            /// read the value of the diagonal element.
            /// </summary>
            public double GetDiagonalElement(int row) {
                return this[row, row];
            }

            /// <summary>
            /// setting of diagonal element.
            /// </summary>
            public void SetDiagonalElement(int row, double val) {
                this[row, row] = val;
            }

            /// <summary>
            /// partitioning of rows over all MPI processes
            /// </summary>
            public IPartitioning RowPartitioning {
                get {
                    return RowMapping;
                }
            }

            /// <summary>
            /// partitioning of columns over all MPI processes
            /// </summary>
            public IPartitioning ColPartition {
                get {
                    return ColMapping;
                }
            }

            /// <summary>
            /// total number of rows over all MPI processes
            /// </summary>
            public int NoOfRows {
                get {
                    return (int)(RowPartitioning.TotalLength);
                }
            }

            /// <summary>
            /// total number of rows over all MPI processes
            /// </summary>
            public int NoOfCols {
                get {
                    return (int)(ColMapping.TotalLength);
                }
            }

            /// <summary>
            /// not supported
            /// </summary>
            public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
                where VectorType1 : IList<double>
                where VectorType2 : IList<double> {
                throw new NotImplementedException();
            }

            public int GetOccupiedColumnIndices(int RowIndex, ref int[] R) {
                throw new NotImplementedException();
                ////this.Ma

                //int[] FullMtxOccupied = m_full.GetOccupiedColumnIndices(iRowFrame2Full(RowIndex));

                //List<int> ret = new List<int>();
                //for (int i = 0; i < FullMtxOccupied.Length; i++) {
                //    //iColFull2Frame(FullMtxOccupied)
                //    throw new NotImplementedException();
                //}

                //return ret.ToArray();
            }

            public int GetRow(int RowIndex, ref int[] ColumnIndices, ref double[] Values) {
                //public MsrMatrix.MatrixEntry[] GetRow(int RowIndex) {
                int LR = GetOccupiedColumnIndices(RowIndex, ref ColumnIndices);
                //var row = new MsrMatrix.MatrixEntry[Occ.Length];
                if (Values == null || Values.Length < LR)
                    Values = new double[LR];

                for (int i = 0; i < LR; i++) {
                    Values[i] = this[RowIndex, ColumnIndices[i]];
                }
                return LR;
            }
        }


        /// <summary>
        /// see <see cref="SpatialOperator.Commit"/>
        /// </summary>
        public override void Commit() {
            base.Commit();
            GhostEdgesOperator.Commit();
            SurfaceElementOperator.Commit();
        }

        class FixedOrder_SpatialOperator : SpatialOperator {
            public FixedOrder_SpatialOperator(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar)
                : base(__DomainVar, __ParameterVar, __CoDomainVar, null) //
            {
                base.QuadOrderFunction = this.QOF;
            }

            int QOF(int[] A, int[] B, int[] C) {
                return m_Order;
            }

            internal int m_Order;
        }

    }
}
