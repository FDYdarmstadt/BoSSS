using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;


namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    
    /// <summary>
    /// The major contribution of this class in addition to the base class is that 
    /// it formalizes the evolution of the Level-Sets by using a <see cref="LevelSetUpdater"/>, cf. <see cref="LevelSetUpdater"/>
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
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for(int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    var configsLevel = new List<MultigridOperator.ChangeOfBasisConfig>();

                    AddMultigridConfigLevel(configsLevel);

                    configs[iLevel] = configsLevel.ToArray();
                }
                return configs;
            }
        }

        /// <summary>
        /// Configuration of operator pre-pre-conditioning (not a typo), cf. <see cref="MultigridOperatorConfig"/>.
        /// </summary>
        protected abstract void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel);


        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {
            XSpatialOperatorMk2 xOperator = GetOperatorInstance(D, LsUpdater);
            return xOperator;
        }

        protected abstract XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater);

        protected override LevelSetTracker InstantiateTracker() {
            LsUpdater = InstantiateLevelSetUpdater();
            return LsUpdater.Tracker;
        }


        /// <summary>
        /// Number of different interfaces 
        /// </summary>
        protected abstract int NoOfLevelSets {
            get;
        }

        /// <summary>
        /// Predefined level-set names; this can be overridden, but is not recommended.
        /// The recommended practice for an app is to 
        /// override <see cref="NoOfLevelSets"/>; then, default names for the level-set-fields are chosen.
        /// </summary>
        protected virtual (string ContLs, string DgLs)[] LevelSetNames {
            get {
                if(NoOfLevelSets == 1) {
                    return new[] { (VariableNames.LevelSetCG, VariableNames.LevelSetDG) };
                } else {
                    throw new NotImplementedException();
                }
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
        /// Sets up the level-set-system (fields for storing, evolution operators, ...) before
        /// XDG-fields can be created.
        /// </summary>
        protected virtual LevelSetUpdater InstantiateLevelSetUpdater() {
            int D = this.Grid.SpatialDimension;
            var lsNames = this.LevelSetNames;
            int NoOfLevelSets = lsNames.Length;
            if(NoOfLevelSets != this.NoOfLevelSets)
                throw new ApplicationException();

            // phase 1: initialization of level-sets
            // ======================================
            LevelSet[] DGlevelSets = new LevelSet[NoOfLevelSets];
            for(int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                var LevelSetCG = lsNames[iLevSet].ContLs;
                var LevelSetDG = lsNames[iLevSet].DgLs;

                int levelSetDegree = Control.FieldOptions[VariableNames.LevelSetCG].Degree;    // need to change naming convention of old XNSE_Solver
                

                switch(Control.Option_LevelSetEvolution) {
                    case LevelSetEvolution.Fourier: {
                        //if(Control.EnforceLevelSetConservation) {
                        //    throw new NotSupportedException("mass conservation correction currently not supported");
                        //}
                        FourierLevelSet fourierLevelSet = new FourierLevelSet(Control.FourierLevSetControl, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                        fourierLevelSet.ProjectField(Control.InitialValues_Evaluators[LevelSetCG]);
                        DGlevelSets[iLevSet] = fourierLevelSet;
                        break;
                    }
                    case LevelSetEvolution.FastMarching: {
                        LevelSet levelSetDG = new LevelSet(new Basis(GridData, levelSetDegree), LevelSetDG);
                        levelSetDG.ProjectField(Control.InitialValues_Evaluators[LevelSetCG]);
                        DGlevelSets[iLevSet] = levelSetDG;
                        break;
                    }
                    case LevelSetEvolution.StokesExtension: {
                        LevelSet levelSetDG = new LevelSet(new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG);
                        levelSetDG.ProjectField(Control.InitialValues_Evaluators[VariableNames.LevelSetCG]);
                        DGlevelSets[iLevSet] = levelSetDG;
                        break;
                    }
                    case LevelSetEvolution.SplineLS: {
                        int nodeCount = 30;
                        Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                        SplineLevelSet SplineLevelSet = new SplineLevelSet(Control.Phi0Initial, new Basis(GridData, levelSetDegree), VariableNames.LevelSetDG, nodeCount);
                        DGlevelSets[iLevSet] = SplineLevelSet;
                        break;
                    }
                    case LevelSetEvolution.None: {
                        LevelSet levelSet1 = new LevelSet(new Basis(GridData, levelSetDegree), LevelSetDG);
                        levelSet1.ProjectField(Control.InitialValues_Evaluators[LevelSetCG]);
                        DGlevelSets[iLevSet] = levelSet1;
                        break;
                    }
                    default:
                    throw new NotImplementedException($"Unknown option for level-set evolution: {Control.Option_LevelSetEvolution}");
                }

                if(DGlevelSets[iLevSet].L2Norm() == 0.0) {
                    Console.WriteLine($"Level-Set field {LevelSetCG} is **exactly** zero: setting entire field to -1.");
                    DGlevelSets[iLevSet].AccConstant(-1.0);
                }
            }

            // phase 2: initialization of updater
            // ==================================
            LevelSetUpdater lsUpdater;
            switch(NoOfLevelSets) {
                case 1:
                lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, 
                    new string[] { "A", "B" }, DGlevelSets[0], lsNames[0].ContLs);
                break;

                default:
                throw new NotImplementedException("Unsupported number of level-sets: " + NoOfLevelSets);
            }


            // phase 3: init of evolvers
            // =========================
            for(int iLevSet = 0; iLevSet < this.LevelSetNames.Length; iLevSet++) {
                var LevelSetCG = lsNames[iLevSet].ContLs;
                var LevelSetDG = lsNames[iLevSet].DgLs;

                // create evolver:
                switch(Control.Option_LevelSetEvolution) {
                    case LevelSetEvolution.Fourier: {
                        break;
                    }
                    case LevelSetEvolution.FastMarching: {
                        var fastMarcher = new FastMarchingEvolver(LevelSetCG, QuadOrder(), D);
                        lsUpdater.AddEvolver(LevelSetCG, fastMarcher);
                        break;
                    }
                    case LevelSetEvolution.StokesExtension: {
                        var stokesExtEvo = new StokesExtensionEvolver(LevelSetCG, QuadOrder(), D,
                            GetBcMap(),
                            this.Control.AgglomerationThreshold, this.GridData);
                        lsUpdater.AddEvolver(LevelSetCG, stokesExtEvo);
                        break;
                    }
                    case LevelSetEvolution.SplineLS: {
                        int nodeCount = 30;
                        Console.WriteLine("Achtung, Spline node count ist hart gesetzt. Was soll hier hin?");
                        var SplineEvolver = new SplineLevelSetEvolver(LevelSetCG, (GridData)(this.GridData));
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
                    lsUpdater.AddLevelSetParameter(VariableNames.LevelSetCG, levelSetVelocity);
                }

            }

            // return
            // ======
            return lsUpdater;
        }


        protected override LevelSetHandling LevelSetHandling {
            get {
                return this.Control.Timestepper_LevelSetHandling;
            }
        }


        

        /// <summary>
        /// The base implementation <see cref="Solution.Application{T}.SetInitial"/>
        /// must be overridden, since it does not preform the continuity projection, see <see cref="DualLevelSet"/>,
        /// but it may overwrite the continuous level set.
        ///
        /// This implementation, however, ensures continuity of the level-set at the cell boundaries.
        /// </summary>
        protected override void SetInitial() {
            base.SetInitial(); // base implementation does not considers the DG/CG pair.

            foreach(var NamePair in this.LevelSetNames) {
                string LevelSetCG = NamePair.ContLs;
                
                // we just overwrite the DG-level-set, continuity projection is set later when the operator is fully set-up
                var pair1 = LsUpdater.LevelSets[LevelSetCG];
                pair1.DGLevelSet.ProjectField(Control.InitialValues_Evaluators[LevelSetCG]);
            }
        }



        /// <summary>
        /// - Matches <see cref="DelUpdateLevelset"/>, used by the <see cref="ApplicationWithSolver{T}.Timestepping"/> to advance the interfaces
        /// - Uses the <see cref="LsUpdater"/>
        /// </summary>
        public override double UpdateLevelset(DGField[] domainFields, double time, double dt, double UnderRelax, bool incremental) {
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Length);
            for (int iVar = 0; iVar < domainFields.Length; iVar++) {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;

            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            double residual = LsUpdater.UpdateLevelSets(DomainVarsDict, ParameterVarsDict, time, dt, UnderRelax, incremental);
            Console.WriteLine("Residual of level-set update: " + residual);
            return 0.0;
        }

        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            base.CreateEquationsAndSolvers(L);

            var domainFields = CurrentState.Fields;
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Count);
            for (int iVar = 0; iVar < domainFields.Count; iVar++) {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++) {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            LsUpdater.InitializeParameters(DomainVarsDict, ParameterVarsDict);
            
            // enforce continuity
            // ------------------
            
            var pair1 = LsUpdater.LevelSets.First().Value;
            var oldCoords1 = pair1.DGLevelSet.CoordinateVector.ToArray();
            UpdateLevelset(this.CurrentState.Fields.ToArray(), 0.0, 0.0, 1.0, false); // enforces the continuity projection upon the initial level set
            double dist1 = pair1.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
            if(dist1 != 0)
                throw new Exception("illegal modification of DG level-set when evolving for dt = 0.");
            UpdateLevelset(this.CurrentState.Fields.ToArray(), 0.0, 0.0, 1.0, false); // und doppelt hält besser ;)
            double dist2 = pair1.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
            if(dist2 != 0)
                throw new Exception("illegal modification of DG level-set when evolving for dt = 0.");
        }
    }
}
