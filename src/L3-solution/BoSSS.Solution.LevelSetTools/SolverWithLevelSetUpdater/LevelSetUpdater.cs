﻿
using BoSSS.Foundation;
using BoSSS.Foundation.Caching;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Tecplot;
using ilPSP;
using ilPSP.LinSolvers.monkey.CUDA;
using ilPSP.Tracing;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {



    /// <summary>
    /// A pair of two DG-fields for level-set handling
    /// - a (potentially) discontinuous (at cell boundaries) representation, <see cref="DGLevelSet"/>
    /// - a continuous representation, <see cref="C0LevelSet"/>
    /// </summary>
    /// <remarks>
    /// For a well-defined XDG method, continuity of the level-set at cell boundaries is mandatory, otherwise the level-set (-surface) is not closed.
    /// However, typical evolution methods (e.g. scalar convection) do not guarantee continuity. 
    /// (Alternatively, one might use a Continuous Galerkin (CG) method for the evolution, e.g. Streamwise-Upwind-Petrov-Galerkin (SUPG).
    /// The CG discretization in BoSSS, however does not support 3D in general and hanging nodes in 2D and 3D and there are no plans for this in future.)
    /// 
    /// For these reasons, a level-set duality is used:
    /// - a DG representation for the level-set evolution
    /// - a continuous projection for the construction of the XDG space.
    /// 
    /// Numerical experiments showed that it is a bad practice to use the continuity projection for the evolution:
    /// The continuity projection **increases** oscillations in higher modes (i.e. jumps in derivatives grow).
    /// If this is involved into the evolution, these oscillations can be coupled back into the interface over multiple (typically 100) timesteps.
    /// Therefore, also **curvature in surface tension terms should be computed from the discontinuous representation**.
    /// In order to reduce the induced oscillations, the DG polynomial degree is typically higher than the degree of the DG representation.
    /// </remarks>
    public class DualLevelSet {

        /// <summary>
        /// 
        /// </summary>
        public string Identification;

        /// <summary>
        /// index within the level-set tracker 
        /// </summary>
        public int LevelSetIndex;

        /// <summary>
        /// The continuous projection of <see cref="DGLevelSet"/>;
        /// Typically, this is of higher polynomial degree in order to fulfill the additional continuity constraints,
        /// without to much L2-difference to the original <see cref="DGLevelSet"/>.
        /// Continuity projection is performed by the <see cref="LevelSetUpdater"/>.
        /// </summary>
        public LevelSet C0LevelSet;

        /// <summary>
        /// The DG representation which should be used for level-set evolution.
        /// </summary>
        public LevelSet DGLevelSet;

        /// <summary>
        /// the associated tracker
        /// </summary>
        public LevelSetTracker Tracker;
    }

    /// <summary>
    /// Update functionality for parameters fields 
    /// - which are either required for the level-set evolution (cf. passed to <see cref="ILevelSetEvolver.MovePhaseInterface"/>)
    /// - or required by the solver, e.g. interface normals and curvature.
    /// </summary>
    public interface ILevelSetParameter {

        /// <summary>
        /// Global/application-wide naming of the respective parameters, cf. <see cref="VariableNames"/>
        /// </summary>
        IList<string> ParameterNames { get; }

        /// <summary>
        /// Allocation; note that this methods signature **is designed to match** <see cref="DelParameterFactory"/>
        /// </summary>
        (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields);

        /// <summary>
        /// update; 
        /// </summary>
        void LevelSetParameterUpdate(
            DualLevelSet levelSet,
            double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }

    /// <summary>
    /// Driver class for level-set evolution; 
    /// It can manage the evolution of multiple level-sets, were each of them is treated with 
    /// its own <see cref="ILevelSetEvolver"/>
    /// - can update more than one level-set at once
    /// - handles the continuity projection
    /// - initializes a level-set-tracker, <see cref="Tracker"/>
    /// </summary>
    public class LevelSetUpdater : XdgTimestepping.ISlaveTimeIntegrator {

        /// <summary>
        /// internal implementation
        /// </summary>
        class SingleLevelSetUpdater {
            internal ILevelSetEvolver lsMover;

            ContinuityProjection enforcer;

            public DualLevelSet phaseInterface;

            public override string ToString() {
                return phaseInterface.Identification;
            }


            ICollection<ILevelSetParameter> lsParameters;

            public ICollection<ILevelSetParameter> LevelSetParameters {
                get { return lsParameters; }
            }

            readonly Func<int> QuadOrderFunc;

            public SingleLevelSetUpdater(DualLevelSet phaseInterface, ContinuityProjection enforcer, Func<int> quadOrderFunc) {
                this.enforcer = enforcer;
                this.phaseInterface = phaseInterface;
                lsParameters = new List<ILevelSetParameter>();
                this.QuadOrderFunc = quadOrderFunc;
            }

            public void SetLevelSetEvolver(ILevelSetEvolver evolver) {
                if (lsMover != null) {
                    throw new Exception("Only one evolver allowed for each levelSet.");
                }
                lsMover = evolver;
            }

            public void AddLevelSetParameter(ILevelSetParameter parameter) {
                // Check if already registered
                foreach (ILevelSetParameter registeredParameter in lsParameters) {
                    foreach (string registeredName in registeredParameter.ParameterNames) {
                        foreach (string name in parameter.ParameterNames) {
                            if (registeredName == name) {
                                throw new Exception("Parameter can only be registered once.");
                            }
                        }
                    }
                }
                lsParameters.Add(parameter);
            }



            internal double UpdateLevelSet(
                IReadOnlyDictionary<string, DGField> DomainVarFields,
                IReadOnlyDictionary<string, DGField> ParameterVarFields,
                double time,
                double dt,
                double underRelax,
                bool incremental) {
                LevelSet ls = phaseInterface.C0LevelSet;
                LevelSet lsBkUp = ls.CloneAs();

                // Move LevelSet and update Params
                // ===============================
                if (dt > 0 && lsMover != null) {
                    MoveLevelSet(
                        phaseInterface,
                        DomainVarFields,
                        ParameterVarFields,
                        time,
                        dt,
                        underRelax,
                        incremental);
                }

                // Make Continuous
                EnforceContinuity();

                if (dt > 0 && lsMover != null) {

                    // perform reintialization, etc. (depends on the lsMover)
                    // ======================================================
                    bool changed = AfterMoveLevelSet(
                        phaseInterface,
                        DomainVarFields,
                        ParameterVarFields,
                        time,
                        dt,
                        underRelax,
                        incremental);

                    if(changed) 
                        EnforceContinuity();
                }


                // Calculate Residual
                // ==================
                CellMask oldCC = phaseInterface.Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex);
                var newCC = phaseInterface.Tracker.Regions.GetCutCellMask();
                lsBkUp.Acc(-1.0, ls);
                double levSetResidual = lsBkUp.L2Norm(newCC.Union(oldCC));
                return levSetResidual;
            }

            /// <summary>
            /// 1. calls the <see cref="ILevelSetEvolver.MovePhaseInterface"/> method
            /// 2. applies under-relaxation
            /// </summary>
            void MoveLevelSet(
                DualLevelSet phaseInterface,
                IReadOnlyDictionary<string, DGField> DomainVarFields,
                IReadOnlyDictionary<string, DGField> ParameterVarFields,
                double time,
                double dt,
                double underRelax,
                bool incremental) {
                LevelSet dglsBkUp = null;
                if (underRelax < 1.0) {
                    dglsBkUp = phaseInterface.DGLevelSet.CloneAs();
                }
                AssertAllRequiredFieldsArePresent(lsMover.ParameterNames, DomainVarFields, ParameterVarFields);
                AssertAllRequiredFieldsArePresent(lsMover.VariableNames, DomainVarFields, ParameterVarFields);

                lsMover.MovePhaseInterface(
                    phaseInterface,
                    time,
                    dt,
                    incremental,
                    DomainVarFields,
                    ParameterVarFields);

                //UnderRelax
                if (underRelax < 1.0) {
                    LevelSet dgLs = phaseInterface.DGLevelSet;
                    dgLs.Scale(underRelax);
                    dgLs.Acc((1.0 - underRelax), dglsBkUp);
                }
            }

            /// <summary>
            /// Additional routines, the evolver might call after Level-Set Movement, <see cref="ILevelSetEvolver.AfterMovePhaseInterface"/>
            /// E.g. Reinitialization
            /// </summary>
            bool AfterMoveLevelSet(
                DualLevelSet phaseInterface,
                IReadOnlyDictionary<string, DGField> DomainVarFields,
                IReadOnlyDictionary<string, DGField> ParameterVarFields,
                double time,
                double dt,
                double underRelax,
                bool incremental) {


                return lsMover.AfterMovePhaseInterface(
                        phaseInterface,
                        time,
                        dt,
                        incremental,
                        DomainVarFields,
                        ParameterVarFields);
            }


            /// <summary>
            /// 
            /// </summary>
            internal void EnforceContinuity() {

                

                EnforceContinuityWithPreEnforcer(); // performs continuity projection on union of old and new cut-cells
                bool isClosed = IsInterfaceClosed(out _); // in some cases, the continuity projection might modify the level-set so, that it leaves the cut-cell domain

                if(enforcer.myOption != ContinuityProjectionOption.None && !isClosed) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // continuity projection on near-band, because continuity projection on cut-cells left some hole in the level-set
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    
                    Console.WriteLine("Enforce continuity on nearband");
                    //EnforceContinuityWithPreEnforcer(ContinuityProjectionOption.ConstrainedDG);

                    LevelSetTracker Tracker = phaseInterface.Tracker;
                    CellMask Near1 = Tracker.Regions.GetSpeciesRestrictedNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
                    CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

                    enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.C0LevelSet, Near1, PosFF);

                    CellMask RemainingNear1 = Tracker.Regions.GetNearMask4LevSet(phaseInterface.LevelSetIndex, 1).Except(Near1);
                    phaseInterface.C0LevelSet.Clear(RemainingNear1);
                    phaseInterface.C0LevelSet.AccLaidBack(1.0, phaseInterface.DGLevelSet, RemainingNear1);
                }
            }

            /// <summary>
            /// The Pre-Enforcer ensures that the projection is performed on new cut-cells in case of a moving interface
            /// </summary>
            internal void EnforceContinuityWithPreEnforcer(ContinuityProjectionOption ProjOpt = ContinuityProjectionOption.None) {

                LevelSetTracker Tracker = phaseInterface.Tracker;
                CellMask Near1 = Tracker.Regions.GetSpeciesRestrictedNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
                CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

                ContinuityProjection preEnforcer = new ContinuityProjection(
                    phaseInterface.C0LevelSet.Basis,
                    phaseInterface.DGLevelSet.Basis,
                    Tracker.GridDat,
                    ProjOpt);

                LevelSet preCGLevelSet = phaseInterface.C0LevelSet.CloneAs();

                preEnforcer.MakeContinuous(phaseInterface.DGLevelSet, preCGLevelSet, Near1, PosFF);

                using(LevelSetTracker preTracker = new LevelSetTracker(Tracker.GridDat, Tracker.CutCellQuadratureType, 1, [ "A", "B" ], preCGLevelSet)) {
                    preTracker.UpdateTracker(0.0);

                    CellMask CC = preTracker.Regions.GetCutCellMask4LevSet(0);
                    CellMask CCplus = CC.Union(Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex));
                    PosFF = preTracker.Regions.GetLevelSetWing(0, +1).VolumeMask;

                    enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.C0LevelSet, CCplus, PosFF);
                }
            }

            //static int counter = 1;

            /// <summary>
            /// Checks for inner contact points/lines, on the boundary of the cut-cell domain.
            /// This might occur when the continuity projection, or the re-initialization, modifies the level-set in a way 
            /// so that it leaves its original cut-cell domain.
            /// 
            /// Sometimes one need to do the projection on the near-band in order to remove holes in the interface
            /// </summary>
            /// <returns></returns>
            internal bool IsInterfaceClosed(out CellMask AdditionalCellsToIncludeInContinuityProjection) {

                //if(this.phaseInterface.LevelSetIndex == 0) {
                //    AdditionalCellsToIncludeInContinuityProjection = null;
                //    return true;
                //}

                using(var tr = new FuncTrace()) {
                    LevelSet preCGLevelSet = phaseInterface.C0LevelSet.CloneAs();
                    using(LevelSetTracker LsTrk = new LevelSetTracker(phaseInterface.Tracker.GridDat, phaseInterface.Tracker.CutCellQuadratureType, 1, ["A", "B"], preCGLevelSet)) { // to be removed...
                        LsTrk.UpdateTracker(0.0);

                        int iLevSet = 0;
                        //int iLevSet = this.phaseInterface.LevelSetIndex;
                        //LevelSetTracker LsTrk = phaseInterface.Tracker;
                        var gdat = LsTrk.GridDat;
                        int JupL = gdat.iLogicalCells.NoOfLocalUpdatedCells;

                        var ccmask = LsTrk.Regions.GetCutCellMask4LevSet(iLevSet);
                        var ccbitmask = ccmask.GetBitMaskWithExternal();

                        // ======================================================================================================================================================================================
                        // Implementation note:
                        // -----------------------
                        // The goal is to detect whether the level-set leaves the cut-cell-domain `ccmask`
                        // Therefore, we search for intersections of the level-set with the boundary of `ccmask`; 
                        // Note: one MUST NOT use the edge-scheme for boundary elements here, since this scheme is implemented so that the edge integrals for the boundary of `ccmask` are zero!!!
                        // the `_GetSurfaceElement_BoundaryRuleFactory(...)`, on the other hand, provides quadrature rules for each cell and does not care about interior or boundary edges of `ccmask`.
                        // ======================================================================================================================================================================================

                        ICompositeQuadRule<CellBoundaryQuadRule> ccBoundaryRule;
                        {
                            int quadOrder = this.QuadOrderFunc();
                            var factoryHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, QuadOrderFunc(), HistoryIndex: 1).XQuadFactoryHelper;
                            var geom_ccmask = ccmask.ToGeometicalMask();

                            var scheme = new CellBoundaryQuadratureScheme(scaling: null, UseDefaultFactories: false, domain: geom_ccmask); // we can omit scaling, since all we are interested is whether the integral is 0.0 or not
                            foreach(var KrefVol in gdat.iGeomCells.RefElements) {
                                var domain_KrefVol = gdat.GetCells4Refelement(KrefVol).Intersect(geom_ccmask);
                                var factory = factoryHelper._GetSurfaceElement_BoundaryRuleFactory(iLevSet, KrefVol);
                                scheme.AddFactoryDomainPair(factory, domain_KrefVol);
                            }

                            ccBoundaryRule = scheme.Compile(gdat, quadOrder);
                        }


                        bool isClosed = true;
                        BitArray addCells = null;


                        // loop over all cut cells
                        // =======================

                        foreach(var crp in ccBoundaryRule) {
                            var rule = crp.Rule;
                            foreach(int jCell in crp.Chunk.Elements) {
                                Debug.Assert(rule.RefElement == gdat.iGeomCells.GetRefElement(jCell));
                                int NoofFaces = gdat.iGeomCells.GetRefElement(jCell).NoOfFaces;

                                int k0 = 0;
                                for(int iFace = 0; iFace < NoofFaces; iFace++) {
                                    int KF = rule.NumbersOfNodesPerFace[iFace];
                                    double WeightSum = KF == 0 ? 0.0 : rule.Weights.ExtractSubArrayShallow([k0], [k0 + KF - 1]).Sum();
                                    k0 += KF;

                                    if(WeightSum == 0.0)
                                        // emty face - no need to check anything
                                        continue;

                                    // for a non-empty face, we need to check whether all neighbors are **also** cut-cells!
                                    var neighs = gdat.GetGeometricalNeighborCells(jCell, iFace);
                                    foreach(int jNeigh in neighs) {
                                        int jNeighLog = gdat.GetLogicalCellIndex(jNeigh);
                                        if(ccbitmask[jNeighLog] == false) {
                                            // we have an non-empty boundary line (`WeightSum != 0`) at face `iFace` and the neighbor is not part of the cut-cell-domain:
                                            // => at face `iFace` of cell `jCell`, we are leaving the cut-cell domain:
                                            // => interface has a hole, is not closed; the neighbor cell should have been included in the continuity projection

                                            isClosed = false;

                                            if(addCells == null)
                                                addCells = new BitArray(JupL);

                                            //int jLog = gdat.GetLogicalCellIndex(jNeighLog);
                                            if(jNeighLog < JupL)
                                                addCells[jNeighLog] = true;
                                        }
                                    }

                                }
                            }
                        }

                        // return
                        // ======

                        bool isClosedGlob = isClosed.MPIOr();

                        if(addCells == null)
                            addCells = new BitArray(JupL);


                        if(!isClosedGlob) {
                            AdditionalCellsToIncludeInContinuityProjection = new CellMask(gdat, addCells, MaskType.Logical);
                        } else {
                            AdditionalCellsToIncludeInContinuityProjection = null;
                        }

                        //if(isClosedGlob == false) {

                        //    var ccMarker = new SinglePhaseField(new Basis(gdat, 0), "ccMarker");
                        //    foreach(int jCell in ccmask.ItemEnum) {
                        //        ccMarker.SetMeanValue(jCell, 1);
                        //    }
                        //    var outMarker = new SinglePhaseField(new Basis(gdat, 0), "outMarker");
                        //    foreach(int jCell in AdditionalCellsToIncludeInContinuityProjection.ItemEnum) {
                        //        outMarker.SetMeanValue(jCell, 1);
                        //    }

                        //    ccBoundaryRule.SaveToTextFileCellBoundary(gdat, "bdnyrule-" + counter + ".csv", writeHeader: false);
                        //    Tecplot.Tecplot.PlotFields([preCGLevelSet, ccMarker, outMarker], "outshit-" + counter, 0.0, 3);

                        //    counter++;

                        //}

                        return isClosedGlob;

                        /*

                            var gdat = LsTrk.GridDat;

                            if(gdat.Grid.RefElements.Length > 1)
                                throw new NotImplementedException("todo");


                            int order = phaseInterface.C0LevelSet.Basis.Degree;
                            var testMetrics = testTracker.GetXDGSpaceMetrics(testTracker.SpeciesIdS.ToArray(), order);
                            var testFactory = testMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, testTracker.GridDat.iGeomEdges.EdgeRefElements.Single());


                            // the boundary of the cut-cell domain
                            EdgeMask CutCellBoundaryEdgeMask = testTracker.Regions.GetCutCellMask().AllEdges().Except(
                                testTracker.Regions.GetCutCellMask().GetAllInnerEdgesMask()).Except(testTracker.GridDat.BoundaryEdges);
                            EdgeQuadratureScheme CutCellInnerBoundary_Scheme = new EdgeQuadratureScheme(
                                new BoSSS.Foundation.XDG.Quadrature.SurfaceElementEdgeIntegrationMetric(testMetrics.LevelSetData[iLevSet]),
                                UseDefaultFactories: false, domain: CutCellBoundaryEdgeMask);
                            //EdgeQuadratureScheme CutCellInnerBoundary_Scheme = new EdgeQuadratureScheme( 
                            //    UseDefaultFactories: false, domain: CutCellBoundaryEdgeMask);

                            CutCellInnerBoundary_Scheme.AddFactoryDomainPair(testFactory);

                            // integrate over the **boundary** of the cut-cell domain.
                            // if the level-set has no holes due to un-detected cut cells, this integral should be zero!
                            // - the quadrature rules on the boundary of the cut cell domain (`GetSurfaceElement_BoundaryRuleFactory`) must be empty!
                            // - only exception: edges at MPI boundaries; 
                            int D = testTracker.GridDat.SpatialDimension;
                            Vector result_InterProcess = new Vector(D);
                            double result_loc = 0;

                            var GlbProblemList = new List<(long jCell0, long jCell1, int iEdge, Vector normal, Vector pt, double weights_sum, double usedmul, int originProc)>();

                            BitArray failEdgesMask = new BitArray(gdat.Edges.Count);
                            var failQRs = new List<(QuadRule qr, int iEdge, int jCell0, int jCell1, bool lCell0isCut, bool jCell1IsCut)>();
                            var ccBitMask = testTracker.Regions.GetCutCellMask().GetBitMask();

                            EdgeQuadrature.GetQuadrature([D + 1], gdat,
                                CutCellInnerBoundary_Scheme.Compile(gdat, order),
                                delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {
                                    //var edgeNormals = testTracker.GridDat.Edges.NormalsCache.GetNormals_Edge(QR.Nodes, i0, length);
                                    // note: on inter-process-edges the normals point into opposite directions,
                                    //       therefore they should cancel out!

                                    for(int i = 0; i < length; i++) { // loop over edges
                                        int iEdge = i0 + i;
                                        EdgeInfo edgInfo = gdat.Edges.Info[iEdge];

                                        if(edgInfo.HasFlag(EdgeInfo.Boundary))
                                            throw new ArithmeticException("some boundary edge should not be part of the integration domain");

                                        if(edgInfo.HasFlag(EdgeInfo.Interprocess)) {
                                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                            // An inter-process-edge: measure is not zero
                                            // If both cells are correctly identified as cut-cells,
                                            // the normals from both processors cancel out.
                                            // when the sum over all processors is taken.
                                            // We also multiply every edge with some unique number, to avoid two errors to cancel each other.
                                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                            double multiplyer; // 
                                            long jCell0_global, jCell1_global;
                                            {
                                                int jCell0 = gdat.iGeomEdges.LogicalCellIndices[iEdge, 0];
                                                int jCell1 = gdat.iGeomEdges.LogicalCellIndices[iEdge, 1];
                                                if(jCell0 < 0 || jCell0 >= gdat.CellPartitioning.LocalLength)
                                                    throw new ApplicationException("expecting Cell0 to be in local range");
                                                if(jCell1 < gdat.CellPartitioning.LocalLength)
                                                    throw new ApplicationException("expecting Cell1 to be external/ghost");

                                                jCell0_global = jCell0 + gdat.CellPartitioning.i0;
                                                jCell1_global = gdat.iParallel.GlobalIndicesExternalCells[jCell1 - gdat.CellPartitioning.LocalLength];
                                                multiplyer = Math.Max(jCell0_global, jCell1_global);

                                                if(jCell0_global >  jCell1_global) {
                                                    long t = jCell0_global;
                                                    jCell0_global = jCell1_global;
                                                    jCell1_global = t;
                                                }
                                            }

                                            var edgeNormals = gdat.Edges.NormalsCache.GetNormals_Edge(QR.Nodes, iEdge, 1);
                                            var globCoords = QR.Nodes.TransformLocal2Global(gdat, iEdge);
                                            var edgNormal = edgeNormals.GetRowPt(0, 0);
                                            var globCoord = globCoords.GetRowPt(0);

                                            GlbProblemList.Add((jCell0_global, jCell1_global, iEdge, edgNormal, globCoord, QR.Weights.Sum(), multiplyer, gdat.MpiRank));

                                            EvalResult.ExtractSubArrayShallow([i, 0, 0], [i, QR.NoOfNodes - 1, D - 1]).Acc(multiplyer, edgeNormals);


                                        } else {
                                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                            // interior edge within local MPI process
                                            // SurfaceElement_BoundaryRule measure must be 0.0,
                                            // unless the Level-Set has an open end!!!
                                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                            for(int qn = 0; qn < QR.NoOfNodes; qn++) {
                                                // on a closed interface, this integrand should be zero
                                                // (either the quadrature rule is empty or the weights are zero)

                                                EvalResult[i, qn, D] = 1.0;


                                                int jCell0 = gdat.iGeomEdges.CellIndices[iEdge, 0];
                                                int jCell1 = gdat.iGeomEdges.CellIndices[iEdge, 1];

                                                failEdgesMask[iEdge] = true;
                                                failQRs.Add((QR, iEdge, jCell0, jCell1, ccBitMask[jCell0], jCell1 >= 0 ? ccBitMask[jCell1] : false));
                                            }
                                        }
                                    }
                                },
                                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                                    for(int i = 0; i < length; i++) {
                                        //if (resOfIntg != 0) {
                                        //    int[,] cellIdx = LsTrk.GridDat.Edges.CellIndices;
                                        //    long jIn = (rank == 0) ? cellIdx[i0 + i, 0] : cellIdx[i0 + i, 0] + cellPart.i0;
                                        //    long jOut = (rank == 0) ? cellIdx[i0 + i, 1] : cellIdx[i0 + i, 1] - LsTrk.GridDat.Cells.NoOfExternalCells;
                                        //    Console.WriteLine($"proc {rank} - edge {i0 + i} ({jIn}, {jOut}): integration result = {ResultsOfIntegration[i, 0]}");
                                        //}
                                        result_loc += ResultsOfIntegration[i, D];
                                        result_InterProcess += ResultsOfIntegration.GetRowPart(i, 0, D);
                                    }
                                }
                            ).Execute();

                            //ilPSP.Environment.StdoutOnlyOnRank0 = true;

                            result_InterProcess = result_InterProcess.ToArray().MPISum();
                            bool isClosed_Interprocess = result_InterProcess.L2Norm() < 1e-10 * gdat.CellPartitioning.TotalLength;
                            bool isClosed_Locally = result_loc == 0.0;

                            bool isClosed = (isClosed_Interprocess && isClosed_Locally).MPIAnd();


                            if(!isClosed_Interprocess) {
                                tr.Info($"Interface not closed between MPI processes (interprocess result = {result_InterProcess.L2Norm()})");
                                //throw new ArithmeticException($"Interface not closed between MPI processes: local result = {result_loc}, interprocess result = {result_InterProcess.L2Norm()}");
                            }

                            if(!isClosed_Locally) {
                                tr.Info($"Interface not closed: result = {isClosed_Locally}");
                                //(new EdgeMask(gdat, failEdgesMask, MaskType.Geometrical)).SaveToTextFile($"unbalancedEdges.csv");
                                //testFactory.GetQuadRuleSet(CutCellInnerBoundary_Scheme.Domain, order);
                                //throw new ArithmeticException($"Interface not closed: local result = {result_loc}, interprocess result = {result_InterProcess.L2Norm()}");
                            }

                            return isClosed;
                        */
                    }
                }

            }

            public void UpdateParameters(
                IReadOnlyDictionary<string, DGField> DomainVarFields,
                IReadOnlyDictionary<string, DGField> ParameterVarFields,
                double time) {
                int i = 0;
                foreach (ILevelSetParameter parameter in lsParameters) {

                    if (!phaseInterface.C0LevelSet.GridDat.IsAlive())
                        throw new ApplicationException("CG level set on invalidated mesh -- something went wrong during mesh adaptation/load balancing.");
                    if (!phaseInterface.DGLevelSet.GridDat.IsAlive())
                        throw new ApplicationException("CG level set on invalidated mesh -- something went wrong during mesh adaptation/load balancing.");

                    parameter.LevelSetParameterUpdate(
                        phaseInterface,
                        time,
                        DomainVarFields,
                        ParameterVarFields);
                    i++;
                }
            }

            void AssertAllRequiredFieldsArePresent(IList<string> names, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
                if (names != null) {
                    foreach (string name in names) {
                        if (!(DomainVarFields.ContainsKey(name) || ParameterVarFields.ContainsKey(name))) {
                            throw new Exception($"LevelSetUpdater is missing field with name {name}");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Internal by-products, which can be used e.g. for the sake of logging and plotting;
        /// cf. <see cref="ILevelSetEvolver.InternalFields"/>
        /// </summary>
        public IDictionary<string, DGField> InternalFields {
            get {
                var Ret = new Dictionary<string, DGField>();

                foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) {
                    var IP = updater.lsMover?.InternalFields;
                    if (IP != null) {
                        foreach (var kv in IP) {
                            Ret.Add(kv.Key, kv.Value);
                        }
                    }
                }

                return Ret;
            }
        }


        /// <summary>
        /// initialized by the constructor
        /// </summary>
        public LevelSetTracker Tracker {
            get;
            private set;
        }

        Dictionary<string, SingleLevelSetUpdater> lsUpdaters;

        Dictionary<string, DGField> lsParameterFields;


        //Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> GetNamedInputFields;


        /// <summary>
        /// constructor for one level-set, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[], ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, 
            CutCellQuadratureMethod cutCellquadType, Func<int> quadOrderFunc,
            int __NearRegionWidth, string[] _SpeciesTable,
            Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet, string interfaceName, ContinuityProjectionOption continuityMode, int CGLevelSetDegree = -1) {

            LevelSet cgLevelSet = ContinuityProjection.CreateField(
                    dgLevelSet, backgroundGrid, continuityMode, CGLevelSetDegree);
            cgLevelSet.Identification = interfaceName;
            //cgLevelSet.AccLaidBack(1.0, dgLevelSet);

            //this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSet);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(1);

            DualLevelSet levelSet0 = new DualLevelSet {
                Identification = interfaceName,
                LevelSetIndex = 0,
                C0LevelSet = cgLevelSet,
                DGLevelSet = dgLevelSet,
                Tracker = Tracker,
            };
            SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(levelSet0, backgroundGrid, continuityMode, quadOrderFunc);
            lsUpdaters.Add(levelSet0.Identification, singleUpdater);
        }

        /// <summary>
        /// Constructor for two level-sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[,], ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, 
            CutCellQuadratureMethod cutCellquadType, Func<int> quadOrderFunc,
            int __NearRegionWidth, string[,] _SpeciesTable,
            //Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet0, string interfaceName0, ContinuityProjectionOption continuityMode0,
            LevelSet dgLevelSet1, string interfaceName1, ContinuityProjectionOption continuityMode1) {

            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1 };
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1 };
            ContinuityProjectionOption[] options = new ContinuityProjectionOption[] { continuityMode0, continuityMode1 };


            LevelSet[] cgLevelSets = new LevelSet[2];

            for (int i = 0; i < 2; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, options[i]);
                cgLevelSets[i].Identification = interfaceNames[i];
                //cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            //this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(2);
            for (int i = 0; i < 2; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    C0LevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, options[i], quadOrderFunc);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        
        /// <summary>
        /// Constructor for three level-sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[,,], ILevelSet, ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, 
            CutCellQuadratureMethod cutCellquadType, Func<int> quadOrderFunc,
            int __NearRegionWidth, string[,,] _SpeciesTable,
            Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet0, string interfaceName0, ContinuityProjectionOption continuityMode0,
            LevelSet dgLevelSet1, string interfaceName1, ContinuityProjectionOption continuityMode1,
            LevelSet dgLevelSet2, string interfaceName2, ContinuityProjectionOption continuityMode2) {

            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1, dgLevelSet2 };
            LevelSet[] cgLevelSets = new LevelSet[dgLevelSets.Length];
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1, interfaceName2 };
            ContinuityProjectionOption[] options = new ContinuityProjectionOption[] { continuityMode0, continuityMode1, continuityMode2 };

            for (int i = 0; i < dgLevelSets.Length; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, options[i]);
                cgLevelSets[i].Identification = interfaceNames[i];
                //cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            //this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1], cgLevelSets[2]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(3);
            for (int i = 0; i < dgLevelSets.Length; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    C0LevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, options[i], quadOrderFunc);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        /// <summary>
        /// constructor for four level sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[,,,], ILevelSet, ILevelSet, ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, 
            CutCellQuadratureMethod cutCellquadType, Func<int> quadOrderFunc,
            int __NearRegionWidth, string[,,,] _SpeciesTable,
            //Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet0, string interfaceName0, ContinuityProjectionOption continuityMode0,
            LevelSet dgLevelSet1, string interfaceName1, ContinuityProjectionOption continuityMode1,
            LevelSet dgLevelSet2, string interfaceName2, ContinuityProjectionOption continuityMode2,
            LevelSet dgLevelSet3, string interfaceName3, ContinuityProjectionOption continuityMode3) {

            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1, dgLevelSet2, dgLevelSet3 };
            LevelSet[] cgLevelSets = new LevelSet[dgLevelSets.Length];
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1, interfaceName2, interfaceName3 };
            ContinuityProjectionOption[] options = new ContinuityProjectionOption[] { continuityMode0, continuityMode1, continuityMode2, continuityMode3 };


            for (int i = 0; i < dgLevelSets.Length; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, options[i]);
                cgLevelSets[i].Identification = interfaceNames[i];
                //cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            //this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1], cgLevelSets[2], cgLevelSets[3]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(4);
            for (int i = 0; i < dgLevelSets.Length; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    C0LevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, options[i], quadOrderFunc);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        static SingleLevelSetUpdater CreateSingleLevelSetUpdater(DualLevelSet levelSet, GridData grid, ContinuityProjectionOption continuityMode, Func<int> quadOrderFunc) {
            ContinuityProjection enforcer1 = new ContinuityProjection(
                levelSet.C0LevelSet.Basis,
                levelSet.DGLevelSet.Basis,
                grid,
                continuityMode);
            return new SingleLevelSetUpdater(levelSet, enforcer1, quadOrderFunc);
        }


        /// <summary>
        /// current (i.e. after latest update) values of level-set related parameter fields
        /// </summary>
        public IReadOnlyDictionary<string, DGField> Parameters {
            get { return lsParameterFields; }
        }


        /// <summary>
        /// All fluid interfaces handled by this updater
        /// - key: the name of the interface
        /// - value: A pair of a discontinuous (DG) and continuous (CD) field; the zero set of the latter describes the interface. the first one is used to compute the evolution.
        /// </summary>
        public IReadOnlyDictionary<string, DualLevelSet> LevelSets {
            get {
                var levelSets = new Dictionary<string, DualLevelSet>(lsUpdaters.Count);
                foreach (var lsUpater in lsUpdaters) {
                    string lsName = lsUpater.Key;
                    DualLevelSet ls = lsUpater.Value.phaseInterface;
                    levelSets.Add(lsName, ls);
                }
                return levelSets;
            }
        }

        public void AddLevelSetParameter(string levelSetName, ILevelSetParameter levelSetParameter) {
            if (lsUpdaters.TryGetValue(levelSetName, out SingleLevelSetUpdater mover)) {
                mover.AddLevelSetParameter(levelSetParameter);
            }
            else {
                throw new Exception("LevelSet not registered");
            }
        }

        public void AddEvolver(string levelSetName, ILevelSetEvolver evolver) {
            if (lsUpdaters.TryGetValue(levelSetName, out SingleLevelSetUpdater mover)) {
                mover.SetLevelSetEvolver(evolver);
            }
            else {
                throw new Exception("LevelSet not registered");
            }
        }

        /// <summary>
        /// Evolution of all level-sets, fits <see cref="XdgTimestepping.DelUpdateLevelset"/>.
        /// </summary>
        double UpdateLevelSets(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time,
            double dt,
            double underRelax,
            bool incremental) {
            var InnerParameterFields = Combine(ParameterVarFields, this.lsParameterFields);
            double residual = 0;

            //if(dt > 0) {
            //    Console.WriteLine("Testcode jkkjxhvckjhvkjhvbkjxykjvcjxyvnxkcjy nvb kjyvhxckj ykjxyhvndfjfjadhgfvjdsaxcdvfhjdcovk njldlkvfd. vb");
            //}

            UpdateParameters(DomainVarFields, InnerParameterFields, time);
            //if(dt > 0) {
            //    DGField[] LevelSetFields = lsUpdaters.Values.Select(lsu => lsu.phaseInterface.C0LevelSet).Concat(lsUpdaters.Values.Select(lsu => lsu.phaseInterface.DGLevelSet)).ToArray();
            //    Tecplot.Tecplot.PlotFields(LevelSetFields, "LevsetBeforeUpdate", time, 3);
            //}
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) {
                var resi_x = updater.UpdateLevelSet(
                    DomainVarFields,
                    InnerParameterFields,
                    time,
                    dt,
                    underRelax,
                    incremental);

                residual += resi_x.Abs();
            }
            //if(dt > 0) {
            //    DGField[] LevelSetFields = lsUpdaters.Values.Select(lsu => lsu.phaseInterface.C0LevelSet).Concat(lsUpdaters.Values.Select(lsu => lsu.phaseInterface.DGLevelSet)).ToArray();
            //    LevelSetFields = LevelSetFields.Concat(InnerParameterFields.Select(ipf => ipf.Value)).ToArray();

            //    Tecplot.Tecplot.PlotFields(LevelSetFields, "LevsetAfterUpdate", time, 3);
            //    Console.WriteLine("so what");
            //}
            Tracker.UpdateTracker(time + dt, -1, incremental: true);
            UpdateParameters(DomainVarFields, InnerParameterFields, time + dt); // update parameters after change of level-set.
            //Tecplot.Tecplot.PlotFields(new DGField[] { InnerParameterFields["VelocityX@Phi"], InnerParameterFields["VelocityY@Phi"] }, "Velocity@Phi_AfterSecondUpdate", time, 3);
            //Console.WriteLine("... done");

            return residual;
        }


        /// <summary>
        /// 
        /// </summary>
        public void EnforceContinuity() {
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) {
                updater.EnforceContinuity();
            }
        }


        /// <summary>
        /// Allocation
        /// </summary>
        public void InitializeParameters(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time = 0.0) {
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) {
                InitializeParameters(updater.LevelSetParameters, DomainVarFields, ParameterVarFields);
            }
            var InnerParameterFields = Combine(ParameterVarFields, this.lsParameterFields);
            //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(DomainVarFields.Values, InnerParameterFields.Values), "beforeUpdateParametersInitial", 0.0, 2);
            UpdateParameters(DomainVarFields, InnerParameterFields, 0.0);
            //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(DomainVarFields.Values, InnerParameterFields.Values), "afterUpdateParametersInitial", 0.0, 2);


        }

        void InitializeParameters(
            ICollection<ILevelSetParameter> parameters,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            if (lsParameterFields == null)
                lsParameterFields = new Dictionary<string, DGField>();
            foreach (ILevelSetParameter parameter in parameters) {
                LinkedList<string> notFound = new LinkedList<string>();
                foreach (string pName in parameter.ParameterNames) {
                    if (!ParameterVarFields.ContainsKey(pName) || ParameterVarFields[pName] == null) {
                        notFound.AddLast(pName);
                    }
                }
                if (notFound.Count > 0) {
                    (string name, DGField field)[] parameterFields = parameter.ParameterFactory(DomainVarFields);
                    while (notFound.Count > 0) {
                        string current = notFound.First.Value;
                        notFound.RemoveFirst();

                        for (int i = 0; i < parameterFields.Length; ++i) {
                            if (parameterFields[i].name == current) {
                                lsParameterFields.Add(parameterFields[i].name, parameterFields[i].field);
                            }
                        }
                    }
                }
            }
        }

        IReadOnlyDictionary<string, DGField> Combine(IReadOnlyDictionary<string, DGField> a, IReadOnlyDictionary<string, DGField> b) {
            Dictionary<string, DGField> combination = new Dictionary<string, DGField>(a.Count + b.Count + 5);
            foreach (var entry in a) {
                combination.Add(entry.Key, entry.Value);
            }
            foreach (var entry in b) {
                combination.Add(entry.Key, entry.Value);
            }
            return combination;
        }

        /// <summary>
        /// 
        /// </summary>
        void UpdateParameters(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time) {
            int cnt = 0;
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) {
                updater.UpdateParameters(DomainVarFields, ParameterVarFields, time);
                cnt++;
            }
        }


        /// <summary>
        /// <see cref="XdgTimestepping.ISlaveTimeIntegrator.Update"/>
        /// </summary>
        public double Update(
            IList<string> domainVar, DGField[] CurrentVariables,
            IList<string> parameterVar, DGField[] CurrentParameters,
            double time, double dt, double UnderRelax, bool incremental) {

            //(IReadOnlyDictionary<string, DGField> DomainVarFields,
            //IReadOnlyDictionary<string, DGField> ParameterVarFields) = GetNamedInputFields(CurrentState);

            var DomainVarFields = new Dictionary<string, DGField>(CurrentVariables.Length);
            for (int iVar = 0; iVar < CurrentVariables.Length; iVar++) {
                if (!CurrentVariables[iVar].GridDat.IsAlive())
                    throw new ApplicationException("Trying to work on field with invalidated grid object.");
                DomainVarFields.Add(domainVar[iVar], CurrentVariables[iVar]);
            }
            var ParameterVarFields = new Dictionary<string, DGField>(CurrentParameters.Count());
            for (int iVar = 0; iVar < CurrentParameters.Count(); iVar++) {
                if (!CurrentParameters[iVar].GridDat.IsAlive())
                    throw new ApplicationException("Trying to work on field with invalidated grid object.");
                ParameterVarFields.Add(parameterVar[iVar], CurrentParameters[iVar]);
            }


            return this.UpdateLevelSets(DomainVarFields, ParameterVarFields, time, dt, UnderRelax, incremental);
        }

        public void Push() {

        }

        public void Pop() {

        }
    }
}
