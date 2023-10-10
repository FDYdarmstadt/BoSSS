using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.Tecplot;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using ilPSP.Tracing;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    
    /// <summary>
    /// The major contribution of this class in addition to the base class is that 
    /// it formalizes the evolution of the Level-Sets by using a <see cref="LevelSetUpdater"/>, cf. <see cref="LsUpdater"/>
    /// </summary>
    /// <remarks>
    /// - created by Lauritz Beck, dec 2020
    /// - further modified and reworked by F Kummer, jan 2021
    /// </remarks>
    public abstract class SolverWithLevelSetUpdater<T> : XdgApplicationWithSolver<T> where T : SolverWithLevelSetUpdaterControl, new() {
        
        /// <summary>
        /// 
        /// </summary>
        public LevelSetUpdater LsUpdater;

        /// <summary>
        /// Quadrature order for everything;
        /// This central computation of the quadrature order should ensure that the cut-cell quadrature rules are only 
        /// constructed for a single order.
        /// </summary>
        abstract public int QuadOrder();

        protected override MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration entry will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[15][];
                for(int iLevel = 0; iLevel < Math.Min(15, configs.Length); iLevel++) { // after level 14, the same config is used over an over;
                    var configsLevel = new List<MultigridOperator.ChangeOfBasisConfig>();

                    AddMultigridConfigLevel(configsLevel, iLevel);

                    configs[iLevel] = configsLevel.ToArray();
                }
                return configs;
            }
        }

        /// <summary>
        /// Configuration of operator pre-pre-conditioning (not a typo), cf. <see cref="MultigridOperatorConfig"/>.
        /// </summary>
        protected abstract void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel);



        /// <summary>
        /// Instantiation of the spatial operator;
        /// Can only be called once per gird lifetime (until <see cref=""/>
        /// </summary>
        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {
            // fails on a second call, if the LsUpdater is already configured.
            // access `base.XOperator`
            XSpatialOperatorMk2 xOperator = GetOperatorInstance(D, LsUpdater);
            return xOperator;
        }

        /// <summary>
        /// Instantiation of the spatial operator; 
        /// </summary>
        protected abstract XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater);

        /// <summary>
        /// Setup of the level-set system and the tracker
        /// </summary>
        protected override LevelSetTracker InstantiateTracker() {
            LsUpdater = InstantiateLevelSetUpdater();

            // register all managed LevelSets
            foreach (DualLevelSet LevSet in LsUpdater.LevelSets.Values) {
                base.RegisterField(LevSet.CGLevelSet);
                base.RegisterField(LevSet.DGLevelSet);
            }

            // register internal fields, e.g. extension velocity etc.
            foreach (var field in LsUpdater.InternalFields.Values) {
                base.RegisterField(field);
            }

            return LsUpdater.Tracker;
        }


        /// <summary>
        /// Number of different interfaces 
        /// </summary>
        protected abstract int NoOfLevelSets {
            get;
        }

        /// <summary>
        /// Species table for the initialization of the <see cref="LevelSetTracker"/>,
        /// see <see cref="LevelSetTracker.SpeciesTable"/>.
        /// Dimension of array must be equal to <see cref="NoOfLevelSets"/>.
        /// </summary>
        protected abstract Array SpeciesTable {
            get;
        }

        /// <summary>
        /// Predefined level-set names; this can be overridden, but is not recommended.
        /// The recommended practice for an app is to 
        /// override <see cref="NoOfLevelSets"/>; then, default names for the level-set-fields are chosen.
        /// </summary>
        protected virtual (string ContLs, string DgLs)[] LevelSetNames {
            get {
                var ret = new ValueTuple<string, string>[NoOfLevelSets];
                for(int i = 0; i < NoOfLevelSets; i++) {
                    ret[i] = (VariableNames.LevelSetCGidx(i), VariableNames.LevelSetDGidx(i));
                }
                return ret;
            }
        }

        /// <summary>
        /// The evolution velocity for the <paramref name="iLevSet"/>-th level-set;
        /// - can be null, if the level-set should not be moved
        /// - The <see cref="ILevelSetParameter.ParameterNames"/> must comply with the following convention:
        ///   <see cref="VariableNames.AsLevelSetVariable"/>( s , <see cref="VariableNames.VelocityVector"/> ),
        ///   where s is the first item of <see cref="LevelSetNames"/>
        /// </summary>
        abstract protected ILevelSetParameter GetLevelSetVelocity(int iLevSet);

        /// <summary>
        /// boundary condition mapping, mainly required for the Stokes extension, where a velocity boundary condition is required.
        /// </summary>
        protected abstract IncompressibleBoundaryCondMap GetBcMap();
        
        /// <summary>
        /// Instantiate the level-set-system (fields for storing, evolution operators, ...) 
        /// Before creating XDG-fields one need to
        /// initialize the level-set fields <see cref="InitializeLevelSets"/>
        /// </summary>
        protected virtual LevelSetUpdater InstantiateLevelSetUpdater() {
            if(!this.GridData.IsAlive())
                throw new ApplicationException("invalid grid -- most likely something went wrong during mesh adaptation/redistribution");
            int D = this.Grid.SpatialDimension;
            var lsNames = this.LevelSetNames;
            //ISpatialOperator test = this.Operator; hier ist noch kein OP
            int NoOfLevelSets = lsNames.Length;
            if(NoOfLevelSets != this.NoOfLevelSets)
                throw new ApplicationException();
            //bool isRestart = Control.RestartInfo != null;


            // phase 1: create DG level-sets
            // ======================================
            LevelSet[] DGlevelSets = new LevelSet[NoOfLevelSets];
            for (int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                var LevelSetCG = lsNames[iLevSet].ContLs;
                var LevelSetDG = lsNames[iLevSet].DgLs;

                int levelSetDegree = Control.FieldOptions[LevelSetCG].Degree;    // need to change naming convention of old XNSE_Solver

                switch (Control.Get_Option_LevelSetEvolution(iLevSet)) {
                    case LevelSetEvolution.Fourier:
                        FourierLevelSet fourierLevelSetDG = new FourierLevelSet(this.Control.FourierLevSetControl, new Basis(GridData, levelSetDegree), LevelSetDG);
                        DGlevelSets[iLevSet] = fourierLevelSetDG;
                        break;
                    case LevelSetEvolution.Prescribed:
                    case LevelSetEvolution.StokesExtension:
                    case LevelSetEvolution.FastMarching:
                    case LevelSetEvolution.Phasefield:
                    case LevelSetEvolution.None: 
                    case LevelSetEvolution.SplineLS: {
                        LevelSet levelSetDG = new LevelSet(new Basis(GridData, levelSetDegree), LevelSetDG);
                        DGlevelSets[iLevSet] = levelSetDG;
                        break;
                    }
                    case LevelSetEvolution.RigidObject: {
                        DGlevelSets[iLevSet] = SetRigidLevelSet(new Basis(GridData, levelSetDegree), LevelSetDG);
                        break;
                    }
                    default:
                        throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

            }


            // phase 2: create updater
            // =======================
            LevelSetUpdater lsUpdater;
            switch (NoOfLevelSets) {
                case 1:
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, Control.LS_TrackerWidth,
                    (string[])this.SpeciesTable,
                    this.GetLsUpdaterInputFields,
                    DGlevelSets[0], lsNames[0].ContLs, Control.LSContiProjectionMethod);
                break;

                case 2:
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, Control.LS_TrackerWidth,
                    (string[,])this.SpeciesTable,
                    this.GetLsUpdaterInputFields,
                    DGlevelSets[0], lsNames[0].ContLs, DGlevelSets[1], lsNames[1].ContLs, Control.LSContiProjectionMethod);
                break;

                default:
                throw new NotImplementedException("Unsupported number of level-sets: " + NoOfLevelSets);
            }


            // phase 3: instantiate evolvers
            // ============================
            for (int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                var LevelSetCG = lsNames[iLevSet].ContLs;
                var LevelSetDG = lsNames[iLevSet].DgLs;

                // create evolver:
                switch(Control.Get_Option_LevelSetEvolution(iLevSet)) {
                    case LevelSetEvolution.Fourier: {
                        FourierLevelSet ls = (FourierLevelSet)DGlevelSets[iLevSet];
                        var fourier = new FourierEvolver(
                            VariableNames.LevelSetCG,
                            ls,
                            Control.FourierLevSetControl,
                            Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Curvature].Degree);

                        lsUpdater.AddEvolver(LevelSetCG, fourier);             
                        break;
                    }
                    case LevelSetEvolution.FastMarching: {
                        var fastMarcher = new FastMarchingEvolver(LevelSetCG, QuadOrder(), D, Control.FastMarchingReInitPeriod);
                        lsUpdater.AddEvolver(LevelSetCG, fastMarcher);
                        break;
                    }
                    case LevelSetEvolution.StokesExtension: {
                        ILevelSetEvolver stokesExtEvo;
                        if (LevelSetHandling == LevelSetHandling.Coupled_Iterative) {
                            stokesExtEvo = new ImplicitStokesExtensionEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(),
                            this.Control.AgglomerationThreshold, this.GridData);
                        } else {
                            stokesExtEvo = new StokesExtensionEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(),
                            this.Control.AgglomerationThreshold, this.GridData,
                            ReInitPeriod: Control.ReInitPeriod);
                        }
                        lsUpdater.AddEvolver(LevelSetCG, stokesExtEvo);
                        break;
                    }
                    case LevelSetEvolution.Phasefield: {
                        var PhasefieldEvolver = new PhasefieldEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(), this.Control,
                            this.Control.AgglomerationThreshold, this.GridData);

                        lsUpdater.AddEvolver(LevelSetCG, PhasefieldEvolver);
                        break;
                    }
                    case LevelSetEvolution.SplineLS: {
                        int nodeCount = 30;
                        Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                        var SplineEvolver = new SplineLevelSetEvolver(LevelSetCG, (GridData)(this.GridData));
                        lsUpdater.AddEvolver(LevelSetCG, SplineEvolver);
                        break;
                    }
                    case LevelSetEvolution.Prescribed: {
                        var prescrEvo = new PrescribedEvolver(this.Control.InitialValues_EvaluatorsVec[LevelSetCG]);
                        lsUpdater.AddEvolver(LevelSetCG, prescrEvo);
                        break;
                    }
                    case LevelSetEvolution.RigidObject: {
                        var rigidEvolver = EvolveRigidLevelSet();
                        lsUpdater.AddEvolver(LevelSetCG, rigidEvolver);
                        break;
                    }
                    case LevelSetEvolution.None: {
                        break;
                    }
                    default:
                    throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

                // add velocity parameter:
                var levelSetVelocity = GetLevelSetVelocity(iLevSet);
                if(levelSetVelocity != null) {
                    if(!ArrayTools.ListEquals(levelSetVelocity.ParameterNames,
                        BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(LevelSetCG, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)))) {
                        throw new ApplicationException($"Parameter names for the level-set velocity provider for level-set #{iLevSet} ({LevelSetCG}) does not comply with convention.");
                    }
                    lsUpdater.AddLevelSetParameter(LevelSetCG, levelSetVelocity);
                }

            }

            // return
            // ======
            return lsUpdater;
        }

        /*
        /// <summary>
        /// Cell-performance classes:
        /// cell performance class equals number of species present in that cell
        /// </summary>
        protected override void GetCellPerformanceClasses(out int NoOfClasses, out int[] CellPerfomanceClasses, int TimeStepNo, double physTime) {
            //throw new NotImplementedException("Dynamic Load Balancing not yet fully implemented - Look at Application MpiRedistributeAndMeshAdapt!");
            NoOfClasses = this.LsTrk.TotalNoOfSpecies;
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            CellPerfomanceClasses = new int[J];
            for(int j = 0; j<J; j++) {
                CellPerfomanceClasses[j] = this.LsTrk.Regions.GetNoOfSpecies(j) - 1;
                //Console.WriteLine("No of Species in cell {0} : {1}", j, CellPerfomanceClasses[j]);
            }
        }
        */

        /// <summary>
        /// Corresponding to <see cref="LevelSetEvolution"/> initialization of LevelSetDG
        /// and projection on continuous LevelSetCG
        /// calls <see cref="LevelSetTracker.UpdateTracker(double, int, bool, int[])">
        /// </summary>
        protected virtual void InitializeLevelSets(LevelSetUpdater lsUpdater, double time) {

            var lsNames = this.LevelSetNames;
            int NoOfLevelSets = lsNames.Length;
            if (NoOfLevelSets != lsUpdater.LevelSets.Count)
                throw new ApplicationException();


            for (int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                string LevelSetDG = lsNames[iLevSet].DgLs;
                string LevelSetCG = lsNames[iLevSet].ContLs;
                DualLevelSet pair = LsUpdater.LevelSets[LevelSetCG];

                int levelSetDegree = Control.FieldOptions[LevelSetCG].Degree;    // need to change naming convention of old XNSE_Solver
                if (levelSetDegree != pair.DGLevelSet.Basis.Degree)
                    throw new ApplicationException();


                ScalarFunction Phi_InitialValue = null;
                if (Control.InitialValues_EvaluatorsVec.TryGetValue(LevelSetCG, out var scalarFunctionTimeDep)) {
                    Phi_InitialValue = scalarFunctionTimeDep.SetTime(0.0);
                }

                switch (Control.Get_Option_LevelSetEvolution(iLevSet)) {
                    case LevelSetEvolution.Fourier: {
                        pair.DGLevelSet.Clear();
                        if (Phi_InitialValue != null)
                            pair.DGLevelSet.ProjectField(Control.InitialValues_EvaluatorsVec[LevelSetCG].SetTime(time));
                        break;
                    }
                    case LevelSetEvolution.Prescribed:
                    case LevelSetEvolution.StokesExtension:
                    case LevelSetEvolution.FastMarching:
                    case LevelSetEvolution.None: {
                        pair.DGLevelSet.Clear();
                        if (Phi_InitialValue != null)
                            pair.DGLevelSet.ProjectField(Control.InitialValues_EvaluatorsVec[LevelSetCG].SetTime(time));
                        break;
                    }
                    case LevelSetEvolution.Phasefield: {
                        pair.DGLevelSet.Clear();
                        if (Phi_InitialValue != null)
                            pair.DGLevelSet.ProjectField(Control.InitialValues_EvaluatorsVec[LevelSetCG].SetTime(time));
                        
                        break;
                    }
                    case LevelSetEvolution.SplineLS: {
                        int nodeCount = 30;
                        Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                        SplineLevelSet SplineLevelSet = new SplineLevelSet(Control.Phi0Initial, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG, nodeCount);
                        if (time != 0.0)
                            Console.WriteLine("Warning: no time dependent initial value");
                        pair.DGLevelSet = SplineLevelSet;
                        break;
                    }
                    case LevelSetEvolution.RigidObject:
                        break;
                    default:
                        throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

                if (pair.DGLevelSet.L2Norm() == 0.0) {
                    Console.WriteLine($"Level-Set field {LevelSetCG} is **exactly** zero: setting entire field to -1.");
                    pair.DGLevelSet.AccConstant(-1.0);
                }

                pair.CGLevelSet.Clear();
                pair.CGLevelSet.AccLaidBack(1.0, pair.DGLevelSet);

            }

            LsUpdater.Tracker.UpdateTracker(time); // update the tracker **before** pushing

            LsUpdater.Tracker.PushStacks();

            LsUpdater.Tracker.UpdateTracker(time);

            //LsUpdater.Tracker.PushStacks();

        }


        /// <summary>
        /// <see cref="XdgTimestepping.LevelSetHandling"/>
        /// </summary>
        protected override LevelSetHandling LevelSetHandling {
            get {
                return this.Control.Timestepper_LevelSetHandling;
            }
        }


        /// <summary>
        /// Used to provide information about rigid objects to the levelSetUpdater in case of <see cref="LevelSetEvolution.RigidObject"/>
        /// Overwrite this to enter information.
        /// </summary>
        protected virtual RigidObjectLevelSet SetRigidLevelSet(Basis Basis, string Name) {
            throw new NotImplementedException();
        }


        /// <summary>
        /// Used to provide information about rigid objects to the levelSetEvolver in case of <see cref="LevelSetEvolution.RigidObject"/>
        /// Overwrite this to enter information.
        /// </summary>
        protected virtual RigidObjectLevelSetEvolver EvolveRigidLevelSet() {
            throw new NotImplementedException();
        }

        /// <summary>
        /// The base implementation <see cref="Solution.Application{T}.SetInitial"/>
        /// must be overridden, since it does not preform the continuity projection, see <see cref="DualLevelSet"/>,
        /// but it may overwrite the continuous level set.
        ///
        /// This implementation, however, ensures continuity of the level-set at the cell boundaries.
        /// </summary>
        protected override void SetInitial(double t) {
            base.SetInitial(t); // base implementation does not considers the DG/CG pair.
            this.InitializeLevelSets(LsUpdater, t);
        }

        /*
        /// <summary>
        /// - Matches <see cref="DelUpdateLevelset"/>, used by the <see cref="ApplicationWithSolver{T}.Timestepping"/> to advance the interfaces
        /// - Uses the <see cref="LsUpdater"/>
        /// </summary>
        public double UpdateLevelset(DGField[] domainFields, double time, double dt, double UnderRelax, bool incremental) {
            using(var tr = new FuncTrace()) {
                var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Length);
                for(int iVar = 0; iVar < domainFields.Length; iVar++) {
                    DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
                }

                var parameterFields = Timestepping.Parameters;

                var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
                for(int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                    ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
                }
                double residual = LsUpdater.UpdateLevelSets(DomainVarsDict, ParameterVarsDict, time, dt, UnderRelax, incremental);
                tr.Info("Residual of level-set update: " + residual);
                return 0.0;
            }
        }
        */


        (IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) GetLsUpdaterInputFields(DGField[] domainFields) {
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Length);
            for(int iVar = 0; iVar < domainFields.Length; iVar++) {
                if(!domainFields[iVar].GridDat.IsAlive())
                    throw new ApplicationException("Trying to work on field with invalidated grid object.");
                if(!object.ReferenceEquals(domainFields[iVar].GridDat, this.GridData))
                    throw new ApplicationException("Grid data object mismatch");
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;

            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for(int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                if(!parameterFields[iVar].GridDat.IsAlive())
                    throw new ApplicationException("Trying to work on field with invalidated grid object.");
                if(!object.ReferenceEquals(parameterFields[iVar].GridDat, this.GridData))
                    throw new ApplicationException("Grid data object mismatch");

                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            return (DomainVarsDict, ParameterVarsDict);
        }


        public override ISlaveTimeIntegrator GetLevelSetUpdater() {
            if(this.LsUpdater == null)
                throw new ApplicationException();
            return this.LsUpdater;
        }


        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            base.CreateEquationsAndSolvers(L);

            // Level Set Parameters
            var domainFields = CurrentState.Fields;
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Count);
            for (int iVar = 0; iVar < domainFields.Count; iVar++) {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = base.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }

            // check if all objects have been updated correctly after mesh adaptation or redistribution
            if(!this.GridData.IsAlive())
                throw new ApplicationException("running on old grid.");
            if(!object.ReferenceEquals(this.GridData, LsTrk.GridDat))
                throw new ApplicationException("LevelSetTracker linked to old grid.");
            if(!object.ReferenceEquals(this.LsTrk, LsUpdater.Tracker))
                throw new ApplicationException("LevelSetUpdater linked to old LevelSetTracker.");
            foreach(var kv in LsUpdater.LevelSets) {
                if(!object.ReferenceEquals(this.GridData, kv.Value.CGLevelSet.GridDat)) {
                    throw new ApplicationException($"CG Level set {kv.Key} linked to old GridData");
                }
                if(!object.ReferenceEquals(this.GridData, kv.Value.DGLevelSet.GridDat)) {
                    throw new ApplicationException($"DG Level set {kv.Key} linked to old GridData");
                }
            }


            LsUpdater.InitializeParameters(DomainVarsDict, ParameterVarsDict);

            foreach (var f in LsUpdater.Parameters.Values) {
                base.RegisterField(f);
            }

            // enforce continuity
            // ------------------

            if (L == null) {
                var pair1 = LsUpdater.LevelSets.First().Value;
                var oldCoords1 = pair1.DGLevelSet.CoordinateVector.ToArray();
                this.LsUpdater.Update(this.CurrentState.Fields.ToArray(), restartTime, 0.0, 1.0, false); // enforces the continuity projection upon the initial level set
                double dist1 = pair1.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
                if (dist1 != 0)
                    throw new Exception("illegal modification of DG level-set when evolving for dt = 0.");
                this.LsUpdater.Update(this.CurrentState.Fields.ToArray(), restartTime, 0.0, 1.0, false); // und doppelt hält besser ;)
                double dist2 = pair1.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
                if (dist2 != 0)
                    throw new Exception("illegal modification of DG level-set when evolving for dt = 0.");
            }
#if TEST
            var MPIrankField = new SinglePhaseField(new Basis(this.GridData, 0), "MPIRank");
            MPIrankField.AccConstant(this.MPIRank);
            base.RegisterField(MPIrankField, IOListOption.Always);

            var CostClusterField = new SinglePhaseField(new Basis(this.GridData, 0), "CostCluster");
            var MaskSpcA = LsTrk.Regions.GetSpeciesMask("A");
            var VoidMask = CellMask.Complement(MaskSpcA);

            CostClusterField.AccNoOfSpecies(1.0, LsTrk, 2);
            CostClusterField.AccConstant(2, VoidMask);
            base.RegisterField(CostClusterField, IOListOption.Always);
#endif
        }


        // Hack to set the correct time for the levelset tracker on restart
        double restartTime = 0.0;


        public override void PostRestart(double time, TimestepNumber timestep) {
            base.PostRestart(time, timestep);

            // Set DG LevelSet by CG LevelSet, if for some reason only the CG is loaded
            if (this.LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet.L2Norm() == 0.0 && this.LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet.L2Norm() != 0.0)
                this.LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet.AccLaidBack(1.0, this.LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet);

            // set restart time, used later in the intial tracker updates
            restartTime = time;

            // push stacks, otherwise we get a problem when updating the tracker, parts of the xdg fields are cleared or something
            this.LsUpdater.Tracker.PushStacks();

        }
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="control"></param>
        public override void Init(AppControl control) {

            void KatastrophenPlot(DGField[] dGFields,string Tag = "") {

                List<DGField> allfields = new();
                allfields.AddRange(dGFields);

                if (string.IsNullOrEmpty(Tag)) {
                    Tag = "AgglomerationKatastrophe";
                }
                
                foreach (var f in this.RegisteredFields) {
                    if (!allfields.Contains(f, (a, b) => object.ReferenceEquals(a, b)))
                        allfields.Add(f);
                }

                Tecplot.Tecplot.PlotFields(allfields, Tag, 0.0, 0);

                //if (Tag.Length > 0) 
                //    Tecplot.Tecplot.PlotFields(allfields, Tag + "AgglomerationKatastrophe_HighRes", 0.0, 2);

            }
            AgglomerationAlgorithm.Katastrophenplot = KatastrophenPlot;
            AgglomerationAlgorithm.PlotAgglomeration = control.PlotAgglomeration;

            base.Init(control);
        }

    }
}
