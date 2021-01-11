﻿using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;


namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public abstract class SolverWithLevelSetUpdater<T> : XdgApplicationWithSolver<T> where T : AppControlSolver, new() {
        
        public LevelSetUpdater LsUpdater;

        protected override MultigridOperator.ChangeOfBasisConfig[][] MultigridOperatorConfig
        {
            get
            {
                // set the MultigridOperator configuration for each level:
                // it is not necessary to have exactly as many configurations as actual multigrid levels:
                // the last configuration enty will be used for all higher level
                MultigridOperator.ChangeOfBasisConfig[][] configs = new MultigridOperator.ChangeOfBasisConfig[3][];
                for (int iLevel = 0; iLevel < configs.Length; iLevel++)
                {
                    var configsLevel = new List<MultigridOperator.ChangeOfBasisConfig>();

                    AddMultigridConfigLevel(configsLevel);

                    configs[iLevel] = configsLevel.ToArray();
                }
                return configs;
            }
        }

        protected abstract void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel);

        protected abstract LevelSetUpdater InstantiateLevelSetUpdater();

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {
            XSpatialOperatorMk2 xOperator = GetOperatorInstance(D, LsUpdater);
            return xOperator;
        }

        protected abstract XSpatialOperatorMk2 GetOperatorInstance(int D, LevelSetUpdater levelSetUpdater);

        protected override LevelSetTracker InstantiateTracker() {
            LsUpdater = InstantiateLevelSetUpdater();
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
            
            var pair = LsUpdater.LevelSets.First().Value;
            var oldCoords1 = pair.DGLevelSet.CoordinateVector.ToArray();
            UpdateLevelset(this.CurrentState.Fields.ToArray(), 0.0, 0.0, 1.0, false); // enforces the continuity projection upon the initial level set
            double dist1 = pair.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
            if(dist1 != 0)
                throw new Exception("illegal modification of DG level-set");

            UpdateLevelset(this.CurrentState.Fields.ToArray(), 0.0, 0.0, 1.0, false); // und doppelt hält besser ;)
            double dist2 = pair.DGLevelSet.CoordinateVector.L2Distance(oldCoords1);
            if(dist2 != 0)
                throw new Exception("illegal modification of DG level-set");
        }                    
    }
}
