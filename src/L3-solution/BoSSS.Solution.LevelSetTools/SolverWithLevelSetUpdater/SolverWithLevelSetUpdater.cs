using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;


namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public abstract class SolverWithLevelSetUpdater<T> : XdgApplicationWithSolver<T> where T : AppControlSolver, new() {
        
        public LevelSetUpdater LsUpdater;

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

        protected abstract LevelSetUpdater InstantiateLevelSetUpdater();

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {
            XSpatialOperatorMk2 xOperator = GetOperatorInstance(D, LsUpdater);
            return xOperator;
        }

        protected abstract XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater);

        protected override LevelSetTracker InstantiateTracker() {
            LsUpdater = InstantiateLevelSetUpdater();

            // register all managed LevelSets
            foreach (DualLevelSet LevSet in LsUpdater.LevelSets.Values) {
                base.RegisterField(LevSet.CGLevelSet);
                base.RegisterField(LevSet.DGLevelSet);
            }

            return LsUpdater.Tracker;
        }

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
            UpdateLevelset(this.CurrentState.Fields.ToArray(), restartTime, 0.0, 1.0, false); // enforces the continuity projection upon the initial level set
            double dist1 = pair1.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
            if(dist1 != 0)
                throw new Exception("illegal modification of DG level-set when evolving for dt = 0.");
            UpdateLevelset(this.CurrentState.Fields.ToArray(), restartTime, 0.0, 1.0, false); // und doppelt hält besser ;)
            double dist2 = pair1.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
            if(dist2 != 0)
                throw new Exception("illegal modification of DG level-set when evolving for dt = 0.");
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

    }
}
