using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ApplicationWithIDT.OptiLevelSets;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP.Utils;

namespace ApplicationWithIDT {
    /// <summary>
    /// Aggregation layer for optimizer DOFs to support multiple level set optimization.
    /// Provides unified helpers for SQP code when multiple IOptiLevelSet objects are involved.
    /// </summary>
    public class LevelSetOptimizationAggregator {
        private List<LevelSetOptimizationState> optimizationStates;

        /// <summary>
        /// Constructor 
        /// </summary>
        public LevelSetOptimizationAggregator() {
            optimizationStates = new List<LevelSetOptimizationState>();
        }

        /// <summary>
        /// Add an optimization state to be managed by this aggregator
        /// </summary>
        public void AddOptimizationState(LevelSetOptimizationState state) {
            if (state != null) {
                optimizationStates.Add(state);
            }
        }

        /// <summary>
        /// Remove all optimization states
        /// </summary>
        public void ClearOptimizationStates() {
            optimizationStates.Clear();
        }

        /// <summary>
        /// Get total level-set DOF count across all optimizers
        /// </summary>
        public int GetTotalDOFCount() {
            return optimizationStates.Sum(state => state.GetOptimizationDOFCount());
        }

        /// <summary>
        /// Map global phi-index to (optimized level-set block, local param index)
        /// </summary>
        /// <param name="globalIndex">Global parameter index</param>
        /// <returns>Tuple of (block index, local index) or null if not found</returns>
        public (int blockIndex, int localIndex)? MapGlobalToLocal(int globalIndex) {
            int currentOffset = 0;
            for (int blockIndex = 0; blockIndex < optimizationStates.Count; blockIndex++) {
                int blockSize = optimizationStates[blockIndex].GetOptimizationDOFCount();
                if (globalIndex >= currentOffset && globalIndex < currentOffset + blockSize) {
                    return (blockIndex, globalIndex - currentOffset);
                }
                currentOffset += blockSize;
            }
            return null; // Invalid global index
        }

        /// <summary>
        /// Project all optimizer objects onto their DG level-sets
        /// </summary>
        public void ProjectAllOptimizerOntoLevelSets() {
            foreach (var state in optimizationStates) {
                state.ProjectOptimizerOntoLevelSet();
            }
        }

        /// <summary>
        /// Clone/backup all optimizer parameter blocks
        /// </summary>
        public void CreateBackupOfAllOptimizers() {
            foreach (var state in optimizationStates) {
                state.CreateBackup();
            }
        }

        /// <summary>
        /// Restore all optimizer parameter blocks from backup
        /// </summary>
        public void RestoreAllOptimizersFromBackup() {
            foreach (var state in optimizationStates) {
                state.RestoreFromBackup();
            }
        }

        /// <summary>
        /// Get all phi-DOFs as a global step vector
        /// </summary>
        public double[] GetGlobalPhiDOFs() {
            var result = new double[GetTotalDOFCount()];
            int offset = 0;
            
            foreach (var state in optimizationStates) {
                var localParams = state.GetCurrentParams();
                if (localParams != null) {
                    Array.Copy(localParams, 0, result, offset, localParams.Length);
                    offset += localParams.Length;
                }
            }
            
            return result;
        }

        /// <summary>
        /// Set all phi-DOFs from a global step vector
        /// </summary>
        public void SetGlobalPhiDOFs(double[] globalStep) {
            if (globalStep == null) return;
            
            int offset = 0;
            foreach (var state in optimizationStates) {
                int blockSize = state.GetOptimizationDOFCount();
                if (offset + blockSize <= globalStep.Length) {
                    var localParams = new double[blockSize];
                    Array.Copy(globalStep, offset, localParams, 0, blockSize);
                    state.SetParams(localParams);
                    offset += blockSize;
                }
            }
        }

        /// <summary>
        /// Accumulate step to all optimizer parameters
        /// </summary>
        /// <param name="globalStep">Global step vector to accumulate</param>
        public void AccumulateGlobalPhiStep(double[] globalStep) {
            if (globalStep == null) return;
            
            int offset = 0;
            foreach (var state in optimizationStates) {
                int blockSize = state.GetOptimizationDOFCount();
                if (offset + blockSize <= globalStep.Length && state.LevelSetOpti != null) {
                    for (int i = 0; i < blockSize; i++) {
                        state.LevelSetOpti.AccToParam(i, globalStep[offset + i]);
                    }
                    offset += blockSize;
                }
            }
        }

        /// <summary>
        /// Get parameter at global index
        /// </summary>
        public double GetParameterAt(int globalIndex) {
            var mapping = MapGlobalToLocal(globalIndex);
            if (mapping.HasValue) {
                var (blockIndex, localIndex) = mapping.Value;
                return optimizationStates[blockIndex].LevelSetOpti?.GetParam(localIndex) ?? 0.0;
            }
            return 0.0;
        }

        /// <summary>
        /// Set parameter at global index
        /// </summary>
        public void SetParameterAt(int globalIndex, double value) {
            var mapping = MapGlobalToLocal(globalIndex);
            if (mapping.HasValue) {
                var (blockIndex, localIndex) = mapping.Value;
                optimizationStates[blockIndex].LevelSetOpti?.SetParam(localIndex, value);
            }
        }

        /// <summary>
        /// Accumulate to parameter at global index
        /// </summary>
        public void AccumulateToParameterAt(int globalIndex, double value) {
            var mapping = MapGlobalToLocal(globalIndex);
            if (mapping.HasValue) {
                var (blockIndex, localIndex) = mapping.Value;
                optimizationStates[blockIndex].LevelSetOpti?.AccToParam(localIndex, value);
            }
        }

        /// <summary>
        /// Add current parameters to history for all optimization states
        /// </summary>
        public void AddAllCurrentParamsToHistory() {
            foreach (var state in optimizationStates) {
                state.AddCurrentParamsToHistory();
            }
        }

        /// <summary>
        /// Get count of optimization states
        /// </summary>
        public int GetOptimizationStateCount() {
            return optimizationStates.Count;
        }

        /// <summary>
        /// Get optimization state at index (for compatibility)
        /// </summary>
        public LevelSetOptimizationState GetOptimizationState(int index) {
            if (index >= 0 && index < optimizationStates.Count) {
                return optimizationStates[index];
            }
            return null;
        }
    }
}