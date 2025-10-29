﻿
using BoSSS.Foundation;
using BoSSS.Foundation.Caching;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.LinSolvers.monkey.CUDA;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
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

            public SingleLevelSetUpdater(DualLevelSet phaseInterface, ContinuityProjection enforcer) {
                this.enforcer = enforcer;
                this.phaseInterface = phaseInterface;
                lsParameters = new List<ILevelSetParameter>(10);
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

                //Move LevelSet and update Params
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

                //Make Continuous
                EnforceContinuity();

                if (dt > 0 && lsMover != null) {
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


                //Calculate Residual
                CellMask oldCC = phaseInterface.Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex);
                var newCC = phaseInterface.Tracker.Regions.GetCutCellMask();
                lsBkUp.Acc(-1.0, ls);
                double levSetResidual = lsBkUp.L2Norm(newCC.Union(oldCC));
                return levSetResidual;
            }

            /// <summary>
            /// 1. calls the <see cref="ILevelSetEvolver.MovePhaseInterface"/> method
            /// 2. calls the continuity projection
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
            /// Additional routines, the Evolver might call after LS Movement
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

                if(lsMover.AfterMovePhaseInterface != null)
                    return lsMover.AfterMovePhaseInterface.Invoke(
                        phaseInterface,
                        time,
                        dt,
                        incremental,
                        DomainVarFields,
                        ParameterVarFields);

                return false;

            }


            //internal void EnforceContinuity(bool enforceOnNearband = false) {
            //    LevelSetTracker Tracker = phaseInterface.Tracker;
            //    CellMask Near1 = Tracker.Regions.GetSpeciesRestrictedNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
            //    CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

            //    //enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);

            //    ContinuityProjection preEnforcer = new ContinuityProjection(
            //        phaseInterface.CGLevelSet.Basis,
            //        phaseInterface.DGLevelSet.Basis,
            //        Tracker.GridDat,
            //        ContinuityProjectionOption.None);
            //    LevelSet preCGLevelSet = phaseInterface.CGLevelSet.CloneAs();
            //    preEnforcer.MakeContinuous(phaseInterface.DGLevelSet, preCGLevelSet, Near1, PosFF);
            //    LevelSetTracker preTracker = new LevelSetTracker(Tracker.GridDat, Tracker.CutCellQuadratureType, 1, new string[] { "A", "B" }, preCGLevelSet);
            //    preTracker.UpdateTracker(0.0);

            //    CellMask CC = preTracker.Regions.GetCutCellMask4LevSet(0);
            //    CellMask CCplus = CC.Union(Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex));
            //    PosFF = preTracker.Regions.GetLevelSetWing(0, +1).VolumeMask;
            //    preTracker.Dispose();

            //    bool FalseContactline = false;

            //    if(!enforceOnNearband) {
            //        enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, CCplus, PosFF);

            //        EdgeMask CCplusBnd = CCplus.AllEdges().Intersect(CellMask.GetFullMask(Tracker.GridDat, MaskType.Logical).GetAllInnerEdgesMask().Except(CCplus.GetAllInnerEdgesMask()));
            //        foreach(int iEdge in CCplusBnd.ItemEnum) {
            //            Tracker.GridDat.Edges.GetRefElement(iEdge).GetNodeSet(5, out NodeSet TestNodes, out _, out _);
            //            MultidimensionalArray phiIn = MultidimensionalArray.Create(1, TestNodes.NoOfNodes);
            //            MultidimensionalArray phiOut = MultidimensionalArray.Create(1, TestNodes.NoOfNodes);
            //            phaseInterface.CGLevelSet.EvaluateEdge(iEdge, 1, TestNodes, phiIn, phiOut);

            //            if(phiIn.Max() * phiIn.Min() < 0 || phiOut.Max() * phiOut.Min() < 0) {
            //                FalseContactline = true;
            //            }
            //        }

            //        FalseContactline = FalseContactline.MPIOr();

            //        if(FalseContactline) {
            //            Console.WriteLine("Error in continuity projection, extending computation to nearband!");
            //            EnforceContinuity(true);
            //        }
            //    } else {
            //        enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);

            //        EdgeMask NearBnd = Near1.AllEdges().Intersect(CellMask.GetFullMask(Tracker.GridDat, MaskType.Logical).GetAllInnerEdgesMask().Except(Near1.GetAllInnerEdgesMask()));
            //        foreach(int iEdge in NearBnd.ItemEnum) {
            //            Tracker.GridDat.Edges.GetRefElement(iEdge).GetNodeSet(5, out NodeSet TestNodes, out _, out _);
            //            MultidimensionalArray phiIn = MultidimensionalArray.Create(1, TestNodes.NoOfNodes);
            //            MultidimensionalArray phiOut = MultidimensionalArray.Create(1, TestNodes.NoOfNodes);
            //            phaseInterface.CGLevelSet.EvaluateEdge(iEdge, 1, TestNodes, phiIn, phiOut);

            //            if(phiIn.Max() * phiIn.Min() < 0) {
            //                FalseContactline = true;
            //            }
            //        }

            //        FalseContactline = FalseContactline.MPIOr();

            //        if(FalseContactline) {
            //            throw new ApplicationException("Continuity projection failed, cannot recover!");
            //        }
            //    }
            //}


            /// <summary>
            /// 
            /// </summary>
            internal void EnforceContinuity() {

                EnforceContinuityWithPreEnforcer();

                if(enforcer.myOption != ContinuityProjectionOption.None && !IsInterfaceClosed()) {
                    Console.WriteLine("Enforce continuity on nearband");
                    //EnforceContinuityWithPreEnforcer(ContinuityProjectionOption.ConstrainedDG);

                    LevelSetTracker Tracker = phaseInterface.Tracker;
                    CellMask Near1 = Tracker.Regions.GetSpeciesRestrictedNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
                    CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

                    enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.C0LevelSet, Near1, PosFF);
                }
            }

            /// <summary>
            /// The preEnforcer ensures that the projection is performed on new cut-cells in case of a moving interface
            /// </summary>
            /// <param name="ProjOpt"></param>
            internal void EnforceContinuityWithPreEnforcer(ContinuityProjectionOption ProjOpt = ContinuityProjectionOption.None) {

                LevelSetTracker Tracker = phaseInterface.Tracker;
                CellMask Near1 = Tracker.Regions.GetSpeciesRestrictedNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
                CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

                //enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);

                ContinuityProjection preEnforcer = new ContinuityProjection(
                    phaseInterface.C0LevelSet.Basis,
                    phaseInterface.DGLevelSet.Basis,
                    Tracker.GridDat,
                    ProjOpt);

                LevelSet preCGLevelSet = phaseInterface.C0LevelSet.CloneAs();

                preEnforcer.MakeContinuous(phaseInterface.DGLevelSet, preCGLevelSet, Near1, PosFF);
                LevelSetTracker preTracker = new LevelSetTracker(Tracker.GridDat, Tracker.CutCellQuadratureType, 1, new string[] { "A", "B" }, preCGLevelSet);
                preTracker.UpdateTracker(0.0);

                CellMask CC = preTracker.Regions.GetCutCellMask4LevSet(0);
                CellMask CCplus = CC.Union(Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex));
                PosFF = preTracker.Regions.GetLevelSetWing(0, +1).VolumeMask;


                enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.C0LevelSet, CCplus, PosFF);
                preTracker.Dispose();
            }

            /// <summary>
            /// Checks for inner contact points/lines,
            /// which may occur when an actually cut cell is **not** detected as cut!
            /// 
            /// Sometimes (so far only 2D) one need to do the projection on the nearband in order to remove holes in the interface
            /// </summary>
            /// <returns></returns>
            internal bool IsInterfaceClosed() {
                const int iLevSet = 0;
                LevelSetTracker LsTrk = phaseInterface.Tracker;
                LevelSet preCGLevelSet = phaseInterface.C0LevelSet.CloneAs();
                LevelSetTracker testTracker = new LevelSetTracker(LsTrk.GridDat, LsTrk.CutCellQuadratureType, 1, new string[] { "A", "B" }, preCGLevelSet);
                testTracker.UpdateTracker(0.0);

                int order = phaseInterface.C0LevelSet.Basis.Degree;
                var testMetrics = testTracker.GetXDGSpaceMetrics(testTracker.SpeciesIdS.ToArray(), order);
                var testFactory = testMetrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(iLevSet, testTracker.GridDat.Grid.RefElements[0]);
                

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
                // if the level-set has no holes due to un-detected cut cells, this integral should be zero
                double result = 0.0;
                int D = testTracker.GridDat.SpatialDimension;
                EdgeQuadrature.GetQuadrature([ 1 ], testTracker.GridDat,
                    CutCellInnerBoundary_Scheme.Compile(testTracker.GridDat, order),
                    delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                        for (int i = 0; i < length; i++) { // 
                            EdgeInfo edgInfo = testTracker.GridDat.Edges.Info[i0 + i];
                            double edgSign = 1.0;
                            if (edgInfo.HasFlag(EdgeInfo.Interprocess)) {
                                // on an inter-process-edge, if both cells are correctly identified as cut-cells
                                // the integral contributions from both processors cancel out 
                                // when the sum over all processors is taken.

                                double[] edgNormal = testTracker.GridDat.Edges.NormalsForAffine.ExtractSubArrayShallow(i0 + i, -1).To1DArray();
                                edgSign = edgNormal.Sum();
                            }
                            for (int qn = 0; qn < QR.NoOfNodes; qn++) {
                                // on a closed interface, this integrand should be zero
                                // (either the quadrature rule is empty or the weights are zero)
                                EvalResult[i, qn, 0] = 1.0 * edgSign;
                            }
                        }
                        //EvalResult.SetAll(1.0);

                    },
                    delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < length; i++) {
                            //if (resOfIntg != 0) {
                            //    int[,] cellIdx = LsTrk.GridDat.Edges.CellIndices;
                            //    long jIn = (rank == 0) ? cellIdx[i0 + i, 0] : cellIdx[i0 + i, 0] + cellPart.i0;
                            //    long jOut = (rank == 0) ? cellIdx[i0 + i, 1] : cellIdx[i0 + i, 1] - LsTrk.GridDat.Cells.NoOfExternalCells;
                            //    Console.WriteLine($"proc {rank} - edge {i0 + i} ({jIn}, {jOut}): integration result = {ResultsOfIntegration[i, 0]}");
                            //}
                            result += ResultsOfIntegration[i, 0];
                        }
                    }
                ).Execute();
                testTracker.Dispose();

                //ilPSP.Environment.StdoutOnlyOnRank0 = true;

                result = result.MPISum();
                bool isClosed = Math.Abs(result) < 1e-10;

                if(!isClosed) {
                    Console.Error.WriteLine($"Interface not closed: result = {result}");
                    //throw new ArithmeticException($"Interface not closed: result = {result}");
                }

                return isClosed;
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
        public LevelSetUpdater(GridData backgroundGrid, CutCellQuadratureMethod cutCellquadType,
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
            SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(levelSet0, backgroundGrid, continuityMode);
            lsUpdaters.Add(levelSet0.Identification, singleUpdater);

            //Tracker.UpdateTracker(0.0);
        }

        /// <summary>
        /// Constructor for two level-sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[,], ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, CutCellQuadratureMethod cutCellquadType,
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
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, options[i]);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }

            //Tracker.UpdateTracker(0.0);
        }

        /// <summary>
        /// Constructor for three level-sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[,,], ILevelSet, ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, CutCellQuadratureMethod cutCellquadType,
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
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, options[i]);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        /// <summary>
        /// constructor for four level sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, CutCellQuadratureMethod, int, string[,,,], ILevelSet, ILevelSet, ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, CutCellQuadratureMethod cutCellquadType,
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
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, options[i]);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        static SingleLevelSetUpdater CreateSingleLevelSetUpdater(DualLevelSet levelSet, GridData grid, ContinuityProjectionOption continuityMode) {
            ContinuityProjection enforcer1 = new ContinuityProjection(
                levelSet.C0LevelSet.Basis,
                levelSet.DGLevelSet.Basis,
                grid,
                continuityMode);
            return new SingleLevelSetUpdater(levelSet, enforcer1);
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

            UpdateParameters(DomainVarFields, InnerParameterFields, time);
            //Tecplot.Tecplot.PlotFields( new DGField[] {lsUpdaters["Phi"].phaseInterface.DGLevelSet, lsUpdaters["Phi"].phaseInterface.CGLevelSet, }, "LevsetBeforeUpdate", time, 3);
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
            //Tecplot.Tecplot.PlotFields( new DGField[] { lsUpdaters["Phi"].phaseInterface.DGLevelSet, lsUpdaters["Phi"].phaseInterface.CGLevelSet, }, "LevsetAfterUpdate", time, 3);
            //Tecplot.Tecplot.PlotFields(new DGField[] { InnerParameterFields["VelocityX@Phi"], InnerParameterFields["VelocityY@Phi"] }, "Velocity@Phi_AfterFirstUpdate", time, 3);
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
