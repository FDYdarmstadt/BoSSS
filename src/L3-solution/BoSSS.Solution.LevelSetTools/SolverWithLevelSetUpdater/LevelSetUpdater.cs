﻿
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {


    /// <summary>
    /// A pair of two DG-fields for level-set handling
    /// - a (potentially) discontinuous (at cell boundaries) representation, <see cref="DGLevelSet"/>
    /// - a continuous representation, <see cref="CGLevelSet"/>
    /// </summary>
    /// <remarks>
    /// For a well-defined XDG method, continuity of the level-set at cell boundaries is mandatory.
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
        public LevelSet CGLevelSet;

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
                if(lsMover != null) {
                    throw new Exception("Only one evolver allowed for each levelSet.");
                }
                lsMover = evolver;
            }

            public void AddLevelSetParameter(ILevelSetParameter parameter) {
                // Check if already registered
                foreach(ILevelSetParameter registeredParameter in lsParameters) {
                    foreach(string registeredName in registeredParameter.ParameterNames) {
                        foreach(string name in parameter.ParameterNames) {
                            if(registeredName == name) {
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
                LevelSet ls = phaseInterface.CGLevelSet;
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
                    AfterMoveLevelSet(
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
                if(underRelax < 1.0) {
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
                if(underRelax < 1.0) {
                    LevelSet dgLs = phaseInterface.DGLevelSet;
                    dgLs.Scale(underRelax);
                    dgLs.Acc((1.0 - underRelax), dglsBkUp);
                }
            }

            /// <summary>
            /// Additional routines, the Evolver might call after LS Movement
            /// E.g. Reinitialization
            /// </summary>
            void AfterMoveLevelSet(
                DualLevelSet phaseInterface,
                IReadOnlyDictionary<string, DGField> DomainVarFields,
                IReadOnlyDictionary<string, DGField> ParameterVarFields,
                double time,
                double dt,
                double underRelax,
                bool incremental) {

                if(lsMover.AfterMovePhaseInterface != null)
                    lsMover.AfterMovePhaseInterface.Invoke(
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
                LevelSetTracker Tracker = phaseInterface.Tracker;
                CellMask Near1 = Tracker.Regions.GetSpeciesRestrictedNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
                CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

                //enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);

                ContinuityProjection preEnforcer = new ContinuityProjection(
                    phaseInterface.CGLevelSet.Basis,
                    phaseInterface.DGLevelSet.Basis,
                    Tracker.GridDat,
                    ContinuityProjectionOption.None);
                LevelSet preCGLevelSet = phaseInterface.CGLevelSet.CloneAs();
                preEnforcer.MakeContinuous(phaseInterface.DGLevelSet, preCGLevelSet, Near1, PosFF);
                LevelSetTracker preTracker = new LevelSetTracker(Tracker.GridDat, Tracker.CutCellQuadratureType, 1, new string[] { "A", "B" }, preCGLevelSet);
                preTracker.UpdateTracker(0.0);

                CellMask CC = preTracker.Regions.GetCutCellMask4LevSet(0);
                CellMask CCplus = CC.Union(Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex));
                PosFF = preTracker.Regions.GetLevelSetWing(0, +1).VolumeMask;

                enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, CCplus, PosFF);
                preTracker.Dispose();
            }


            public void UpdateParameters(
                IReadOnlyDictionary<string, DGField> DomainVarFields,
                IReadOnlyDictionary<string, DGField> ParameterVarFields,
                double time) {
                int i = 0;
                foreach(ILevelSetParameter parameter in lsParameters) {

                    if(!phaseInterface.CGLevelSet.GridDat.IsAlive())
                        throw new ApplicationException("CG level set on invalidated mesh -- something went wrong during mesh adaptation/load balancing.");
                    if(!phaseInterface.DGLevelSet.GridDat.IsAlive())
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
                if(names != null) {
                    foreach(string name in names) {
                        if(!(DomainVarFields.ContainsKey(name) || ParameterVarFields.ContainsKey(name))) {
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

                foreach(SingleLevelSetUpdater updater in lsUpdaters.Values) {
                    var IP = updater.lsMover?.InternalFields;
                    if(IP != null) {
                        foreach(var kv in IP) {
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


        Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> GetNamedInputFields;


        /// <summary>
        /// constructor for one level-set, <see cref="LevelSetTracker.LevelSetTracker(GridData, XQuadFactoryHelper.MomentFittingVariants, int, string[], ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[] _SpeciesTable,
            Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet, string interfaceName, ContinuityProjectionOption continuityMode) {

            LevelSet cgLevelSet = ContinuityProjection.CreateField(
                    dgLevelSet, backgroundGrid, continuityMode);
            cgLevelSet.Identification = interfaceName;
            //cgLevelSet.AccLaidBack(1.0, dgLevelSet);

            this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSet);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(1);

            DualLevelSet levelSet0 = new DualLevelSet {
                Identification = interfaceName,
                LevelSetIndex = 0,
                CGLevelSet = cgLevelSet,
                DGLevelSet = dgLevelSet,
                Tracker = Tracker,
            };
            SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(levelSet0, backgroundGrid, continuityMode);
            lsUpdaters.Add(levelSet0.Identification, singleUpdater);

            //Tracker.UpdateTracker(0.0);
        }

        /// <summary>
        /// Constructor for two level-sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, XQuadFactoryHelper.MomentFittingVariants, int, string[,], ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[,] _SpeciesTable,
            Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet0, string interfaceName0, LevelSet dgLevelSet1, string interfaceName1,
            ContinuityProjectionOption continuityMode) {

            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1 };
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1 };

            LevelSet[] cgLevelSets = new LevelSet[2];

            for (int i = 0; i < 2; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, continuityMode);
                cgLevelSets[i].Identification = interfaceNames[i];
                //cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(2);
            for (int i = 0; i < 2; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    CGLevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, continuityMode);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }

            //Tracker.UpdateTracker(0.0);
        }

        /// <summary>
        /// Constructor for three level-sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, XQuadFactoryHelper.MomentFittingVariants, int, string[,,], ILevelSet, ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[,,] _SpeciesTable, 
            Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet0, string interfaceName0, LevelSet dgLevelSet1, string interfaceName1, LevelSet dgLevelSet2, string interfaceName2,
            ContinuityProjectionOption continuityMode) {
            
            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1, dgLevelSet2 };
            LevelSet[] cgLevelSets = new LevelSet[dgLevelSets.Length];
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1, interfaceName2 };

            for(int i = 0; i < dgLevelSets.Length; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, continuityMode);
                cgLevelSets[i].Identification = interfaceNames[i];
                //cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1], cgLevelSets[2]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(3);
            for(int i = 0; i < dgLevelSets.Length; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    CGLevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, continuityMode);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        /// <summary>
        /// constructor for four level sets, <see cref="LevelSetTracker.LevelSetTracker(GridData, XQuadFactoryHelper.MomentFittingVariants, int, string[,,,], ILevelSet, ILevelSet, ILevelSet, ILevelSet)"/>
        /// </summary>
        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[,,,] _SpeciesTable,
            Func<DGField[], (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields)> __GetNamedInputFields,
            LevelSet dgLevelSet0, string interfaceName0, LevelSet dgLevelSet1, string interfaceName1, LevelSet dgLevelSet2, string interfaceName2, LevelSet dgLevelSet3, string interfaceName3,
            ContinuityProjectionOption continuityMode) {
            
            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1, dgLevelSet2, dgLevelSet3 };
            LevelSet[] cgLevelSets = new LevelSet[dgLevelSets.Length];
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1, interfaceName2, interfaceName3 };

            for(int i = 0; i < dgLevelSets.Length; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, continuityMode);
                cgLevelSets[i].Identification = interfaceNames[i];
                //cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            this.GetNamedInputFields = __GetNamedInputFields;
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1], cgLevelSets[2], cgLevelSets[3]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(4);
            for(int i = 0; i < dgLevelSets.Length; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    CGLevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, continuityMode);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }
        }

        static SingleLevelSetUpdater CreateSingleLevelSetUpdater(DualLevelSet levelSet, GridData grid, ContinuityProjectionOption continuityMode) {
            ContinuityProjection enforcer1 = new ContinuityProjection(
                levelSet.CGLevelSet.Basis,
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
                foreach(var lsUpater in lsUpdaters) {
                    string lsName = lsUpater.Key;
                    DualLevelSet ls = lsUpater.Value.phaseInterface;
                    levelSets.Add(lsName, ls);
                }
                return levelSets;
            }
        }

        public void AddLevelSetParameter(string levelSetName, ILevelSetParameter levelSetParameter) {
            if(lsUpdaters.TryGetValue(levelSetName, out SingleLevelSetUpdater mover)) {
                mover.AddLevelSetParameter(levelSetParameter);
            } else {
                throw new Exception("LevelSet not registered");
            }
        }

        public void AddEvolver(string levelSetName, ILevelSetEvolver evolver) {
            if(lsUpdaters.TryGetValue(levelSetName, out SingleLevelSetUpdater mover)) {
                mover.SetLevelSetEvolver(evolver);
            } else {
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

            //Console.WriteLine("Plotting ...");
            //Tecplot.Tecplot.PlotFields( new DGField[] {InnerParameterFields["VelocityX@Phi"], InnerParameterFields["VelocityY@Phi"] }, "Velocity@Phi_BeforeFirstUpdate", time, 3);
            UpdateParameters(DomainVarFields, InnerParameterFields, time);
            //Tecplot.Tecplot.PlotFields( new DGField[] {lsUpdaters["Phi"].phaseInterface.DGLevelSet, lsUpdaters["Phi"].phaseInterface.CGLevelSet, }, "LevsetBeforeUpdate", time, 3);
            foreach(SingleLevelSetUpdater updater in lsUpdaters.Values) {
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
            foreach(SingleLevelSetUpdater updater in lsUpdaters.Values) {
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
            
            if(lsParameterFields == null)
                lsParameterFields = new Dictionary<string, DGField>();
            foreach(ILevelSetParameter parameter in parameters) {
                LinkedList<string> notFound = new LinkedList<string>();
                foreach(string pName in parameter.ParameterNames) {
                    if(!ParameterVarFields.ContainsKey(pName) || ParameterVarFields[pName] == null) {
                        notFound.AddLast(pName);
                    }
                }
                if(notFound.Count > 0) {
                    (string name, DGField field)[] parameterFields = parameter.ParameterFactory(DomainVarFields);
                    while(notFound.Count > 0) {
                        string current = notFound.First.Value;
                        notFound.RemoveFirst();

                        for(int i = 0; i < parameterFields.Length; ++i) {
                            if(parameterFields[i].name == current) {
                                lsParameterFields.Add(parameterFields[i].name, parameterFields[i].field);
                            }
                        }
                    }
                }
            }
        }

        IReadOnlyDictionary<string, DGField> Combine(IReadOnlyDictionary<string, DGField> a, IReadOnlyDictionary<string, DGField> b) {
            Dictionary<string, DGField> combination = new Dictionary<string, DGField>(a.Count + b.Count + 5);
            foreach(var entry in a) {
                combination.Add(entry.Key, entry.Value);
            }
            foreach(var entry in b) {
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
            foreach(SingleLevelSetUpdater updater in lsUpdaters.Values) {
                updater.UpdateParameters(DomainVarFields, ParameterVarFields, time);
                cnt++;
            }
        }


        /// <summary>
        /// <see cref="XdgTimestepping.ISlaveTimeIntegrator.Update"/>
        /// </summary>
        public double Update(DGField[] CurrentState, double time, double dt, double UnderRelax, bool incremental) {

            (IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) = GetNamedInputFields(CurrentState);


            return this.UpdateLevelSets(DomainVarFields, ParameterVarFields, time, dt, UnderRelax, incremental);
        }

        public void Push() {
            
        }

        public void Pop() {
            
        }
    }
}
