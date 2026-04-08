using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ApplicationWithIDT.OptiLevelSets;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP.Utils;
using ilPSP.LinSolvers;
using ilPSP;

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
        internal void AddOptimizationState(LevelSetOptimizationState state) {
            if (state != null) {
                optimizationStates.Add(state);
            }
        }

        /// <summary>
        /// Remove all optimization states
        /// </summary>
        public void ClearOptimizationStates() {
            foreach ( var state in optimizationStates ) {
                state.LsTBO.Clear();
            }
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
        /// <returns>Tuple of (state index, local index)
        /// - 1st entry (state index): index into <see cref="optimizationStates"/>
        /// - 2nd entry (parameter index into level set) : local parameter index, c.f. <see cref="IOptiLevelSet.GetParam(int)"/>, <see cref="IOptiLevelSet.SetParam(int, double)"/>
        /// </returns>
        public (int stateIndex, int localIndex) MapGlobalToLocal(int globalIndex) {
            int currentOffset = 0;
            for ( int blockIndex = 0; blockIndex < optimizationStates.Count; blockIndex++ ) {
                int blockSize = optimizationStates[blockIndex].GetOptimizationDOFCount();
                if ( globalIndex >= currentOffset && globalIndex < currentOffset + blockSize ) {
                    return (blockIndex, globalIndex - currentOffset);
                }
                currentOffset += blockSize;
            }
            throw new IndexOutOfRangeException(); // Invalid global index
        }

        /// <summary>
        /// Map local phi-index to (optimized level-set block, local param index) global one
        /// </summary>
        /// <param name="globalIndex">Global parameter index</param>
        /// <returns>Tuple of (block index, local index) or null if not found</returns>
        public int MapLocalToGlobal(int blockIndex, int localIndex) {
            int globalIndex = 0;

            for ( int bi = 0; bi < blockIndex; bi++ )
                globalIndex += optimizationStates[bi].GetOptimizationDOFCount();

            globalIndex += localIndex;
            
            return globalIndex;
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
        /// Project all optimizer objects onto their DG level-sets
        /// </summary>
        public void ProjectAllOptimizerOntoLevelSets(LevelSet[] projectionTargets) {
            if ( projectionTargets.Length != optimizationStates.Count )
                throw new ArgumentException();
            int i = 0;
            foreach ( var state in optimizationStates ) {
                state.ProjectOptimizerOntoLevelSet(projectionTargets[i]);
                i++;
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
        /// Backup DG level-set fields for all optimization states.
        /// </summary>
        public void CreateBackupOfAllLevelSets() {
            foreach(var state in optimizationStates) {
                state.CreateLevelSetBackup();
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
        /// Restore DG level-set fields for all optimization states from backup.
        /// </summary>
        public void RestoreAllLevelSetsFromBackup() {
            foreach(var state in optimizationStates) {
                state.RestoreLevelSetFromBackup();
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
            var (blockIndex, localIndex) = MapGlobalToLocal(globalIndex);
            return optimizationStates[blockIndex].LevelSetOpti.GetParam(localIndex);
        }

        /// <summary>
        /// Set parameter at global index
        /// </summary>
        public void SetParameterAt(int globalIndex, double value) {
            var (blockIndex, localIndex) = MapGlobalToLocal(globalIndex);
            optimizationStates[blockIndex].LevelSetOpti.SetParam(localIndex, value);
        }

        /// <summary>
        /// Accumulate to parameter at global index
        /// </summary>
        public void AccumulateToParameterAt(int globalIndex, double value) {
            var (blockIndex, localIndex) = MapGlobalToLocal(globalIndex);
            optimizationStates[blockIndex].LevelSetOpti.AccToParam(localIndex, value);
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
        /// Get parameter name at global index
        /// </summary>
        public string GetParameterNameAt(int globalIndex) {
            var (blockIndex, localIndex) = MapGlobalToLocal(globalIndex);
            return optimizationStates[blockIndex].LevelSetOpti.GetParamName(localIndex);
        }

        /// <summary>
        /// Compute norm of parameter vector
        /// </summary>
        public double Norm(double[] paramVector) {
            if (paramVector == null) return 0.0;
            
            double norm = 0.0;
            int offset = 0;
            
            foreach (var state in optimizationStates) {
                int blockSize = state.GetOptimizationDOFCount();
                if (offset + blockSize <= paramVector.Length && state.LevelSetOpti != null) {
                    var subVector = new double[blockSize];
                    Array.Copy(paramVector, offset, subVector, 0, blockSize);
                    double blockNorm = state.LevelSetOpti.Norm(subVector);
                    norm += blockNorm * blockNorm;
                }
                offset += blockSize;
            }
            
            return Math.Sqrt(norm);
        }

        /// <summary>
        /// Get count of optimization states
        /// </summary>
        public int GetOptimizationStateCount() {
            return optimizationStates.Count;
        }


        public MsrMatrix GetRegMatrix() {
            MsrMatrix RegMatrix_0 = optimizationStates[0].LevelSetOpti.GetRegMatrix();

            for(int i = 1; i < optimizationStates.Count; i++) {
                int I0 = RegMatrix_0.RowPartitioning.LocalLength;
                int J0 = RegMatrix_0.ColPartition.LocalLength;
                long i0 = RegMatrix_0.RowPartitioning.i0;
                long j0 = RegMatrix_0.ColPartition.i0;

                MsrMatrix RegMatrix_i = optimizationStates[i].LevelSetOpti.GetRegMatrix();
                int I1 = RegMatrix_i.RowPartitioning.LocalLength;
                int J1 = RegMatrix_i.ColPartition.LocalLength;
                

                MsrMatrix RegMatrixCat = new MsrMatrix(I0 + I1, J0 + J1);
                RegMatrix_0.AccSubMatrixTo(1.0, RegMatrixCat, default(long[]), I0.ForLoop(i => i + i0), default(long[]), J0.ForLoop(j => j + j0));
                RegMatrix_i.AccSubMatrixTo(1.0, RegMatrixCat, default(long[]), I1.ForLoop(i => i + I0 + i0), default(long[]), J1.ForLoop(j => j + J0 + j0));
                RegMatrix_0 = RegMatrixCat;
            }

            return RegMatrix_0;
        }
    }
}