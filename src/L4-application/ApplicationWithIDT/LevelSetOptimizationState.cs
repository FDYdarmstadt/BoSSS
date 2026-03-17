using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ApplicationWithIDT.OptiLevelSets;
using System.Collections.Generic;
using ilPSP.Utils;

namespace ApplicationWithIDT {
    /// <summary>
    /// Container class for level-set optimization state variables.
    /// Encapsulates all state related to level-set optimization to enable
    /// supporting multiple level-sets in future extensions.
    /// </summary>
    public class LevelSetOptimizationState {
        /// <summary>
        /// The level-set(s) to be optimized - currently singular but extensible to multiple
        /// </summary>
        public LevelSet LsTBO { get; set; }

        /// <summary>
        /// The optimizer object for the level-set
        /// </summary>
        public IOptiLevelSet LevelSetOpti { get; set; }

        /// <summary>
        /// Backup of the optimizer object
        /// </summary>
        public IOptiLevelSet LevelSetOptiBackup { get; set; }

        /// <summary>
        /// Grid used for level-set optimization
        /// </summary>
        public IGrid LevelSetOptiGrid { get; set; }

        /// <summary>
        /// History of level-set optimization parameters
        /// </summary>
        public List<double[]> LevelSetOptiParams { get; set; }

        /// <summary>
        /// Constructor initializes collections
        /// </summary>
        public LevelSetOptimizationState() {
            LevelSetOptiParams = new List<double[]>();
        }

        /// <summary>
        /// Get the current parameter array from the optimizer
        /// </summary>
        public double[] GetCurrentParams() {
            return LevelSetOpti?.GetParamsAsArray();
        }

        /// <summary>
        /// Set parameters on the optimizer from an array
        /// </summary>
        public void SetParams(double[] parameters) {
            if (LevelSetOpti != null && parameters != null) {
                for (int i = 0; i < parameters.Length && i < LevelSetOpti.GetLength(); i++) {
                    LevelSetOpti.SetParam(i, parameters[i]);
                }
            }
        }

        /// <summary>
        /// Add current parameters to history
        /// </summary>
        public void AddCurrentParamsToHistory() {
            var currentParams = GetCurrentParams();
            if (currentParams != null) {
                LevelSetOptiParams.Add(currentParams);
            }
        }

        /// <summary>
        /// Project the optimizer state onto the target DG level-set
        /// </summary>
        public void ProjectOptimizerOntoLevelSet() {
            LevelSetOpti?.ProjectOntoLevelSet(LsTBO);
        }

        /// <summary>
        /// Get total DOF count for optimization
        /// </summary>
        public int GetOptimizationDOFCount() {
            return LevelSetOpti?.GetLength() ?? 0;
        }

        /// <summary>
        /// Create backup of current state
        /// </summary>
        public void CreateBackup() {
            if (LevelSetOpti != null) {
                LevelSetOptiBackup = (IOptiLevelSet)LevelSetOpti.Clone();
            }
        }

        /// <summary>
        /// Restore from backup
        /// </summary>
        public void RestoreFromBackup() {
            if (LevelSetOptiBackup != null) {
                LevelSetOpti = (IOptiLevelSet)LevelSetOptiBackup.Clone();
            }
        }
    }
}