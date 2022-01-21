using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Linq;


namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public abstract class SimpleSolverWithLevelSetUpdater<T> : XdgApplicationWithSolver<T> where T : AppControlSolver, new() {

        public LevelSetUpdater LsUpdater;

        protected override MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig {
            get {
                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++) {
                    var configsLevel = new List<MultigridOperator.ChangeOfBasisConfig>();

                    AddMultigridConfigLevel(configsLevel);

                    configs[iLevel] = configsLevel.ToArray();
                }
                return configs;
            }
        }

        protected abstract void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel);

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {
            XSpatialOperatorMk2 xOperator = GetOperatorInstance(D, LsUpdater);
            return xOperator;
        }

        protected abstract XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater);

        protected abstract LevelSetUpdater InstantiateLevelSetUpdater();
        
        protected override LevelSetTracker InstantiateTracker() {
            LsUpdater = InstantiateLevelSetUpdater();
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

        public override ISlaveTimeIntegrator GetLevelSetUpdater() {
            if (this.LsUpdater == null)
                throw new ApplicationException();
            return this.LsUpdater;
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

            foreach (var f in LsUpdater.Parameters.Values) {
                base.RegisterField(f);
            }

            this.LsUpdater.Update(Operator.DomainVar, this.CurrentState.Fields.ToArray(), Operator.ParameterVar, this.Timestepping.Parameters.ToArray(),
                0.0, 0.0, 1.0, false); // enforces the continuity projection upon the initial level set
        }
    }
}

