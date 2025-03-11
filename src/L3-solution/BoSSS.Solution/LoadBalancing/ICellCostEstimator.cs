/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution.LoadBalancing {

    /// <summary>
    /// Provides a measure to estimate the runtime cost (aka. weight) of all cells in the computational grid.
    /// These weights (<see cref="GetEstimatedCellCosts"/>) 
    /// are then uses by <see cref="LoadBalancer.GetNewPartitioning"/>
    /// to compute a cell partitioning (i.e., which cell is assigned to which MPI process),
    /// which assigns approximately the same weight to each MPI process.
    /// </summary>
    public interface ICellCostEstimator : ICloneable {

        void Init(IApplication app);

        /// <summary>
        /// Updates the internal cell costs
        /// </summary>
        void UpdateEstimates(IApplication app);


        /// <summary>
        /// The estimated cost of each individual cell on this process, for each cluster.
        /// 
        /// Why multiple weights, i.e., clusters?
        /// Very often, one has multiple kinds of cells (aka. classes or clusters, like, e.g., cut and un-cut cells)
        /// where the computational costs vary by one or two magnitudes.
        /// Furthermore, the factor between this costs
        /// - is not constant and depends on many settings, e.g. DG polynomial degree and
        /// - varies for different parts of the algorithm (e.g., for matrix assembly, a cut cell might be 10 times more expensive 
        ///   than an un-cut cell; for the linear solver, this factor might be only 2).
        /// For these reasons, one would like to balance not only the total weight, but also the weight within each cluster.
        ///
        /// If multiple weights per cell are given, the load balancer in
        /// <see cref="LoadBalancer.GetNewPartitioning"/> tries not only to balance the total weight 
        /// (sum over all weights) across the MPI processors, 
        /// but it also tries to balance the weight in each cluster
        /// (i.e., the sum of all weights over all cells, **for each cluster**, is roughly the same for each MPI process.)
        /// This is also referred to as multi-constraint optimization. 
        /// </summary>
        /// <returns>
        /// cell weights for multi-constraint partitioning, where multiple weights are assigned to each cells.
        /// - 1st index: cell cluster/constraint index (only 0 for a single constraint)
        /// - 2nd index: correlates with local cell index.
        /// </returns>
        int[][] GetEstimatedCellCosts();
    }

    
    /// <summary>
    /// Extension methods for <see cref="ICellCostEstimator"/>
    /// </summary>
    public static class ICellCostEstimatorExtensions {

        public static double[] EstimatedLocalCost(this ICellCostEstimator estimator) {
            int[][] cellCost = estimator.GetEstimatedCellCosts();
            double[] ret = new double[cellCost.Length];
            for(int i = 0; i < ret.Length; i++) {
                ret[i] = cellCost[i]?.Sum() ?? 0.0;
            }
            return ret; 
        }

        /// <summary>
        /// Estimates the cost imbalance between all processes using the given
        /// <paramref name="estimator"/>
        /// </summary>
        /// <param name="estimator"></param>
        /// <returns>
        /// A number between 0 and 1 which represents an estimate of the load
        /// imbalance in percent
        /// </returns>
        public static double[] ImbalanceEstimate(this ICellCostEstimator estimator) {
            MPICollectiveWatchDog.Watch();

            double[] localCost = estimator.EstimatedLocalCost();
            int N = localCost.Length;

            double[] globalMin, globalMax;
            {
                // compute MPI min and max in one pass, by using max = -min( -x_i)
                // this saves a lot of expensive MPI communication
                double[] globalMinMax = localCost.CloneAs();
                localCost.ScaleV(-1.0);
                globalMinMax = globalMinMax.Cat(localCost);
                globalMinMax = globalMinMax.MPIMin();

                globalMin = globalMinMax.GetSubVector(0, N);
                globalMax = globalMinMax.GetSubVector(N, N);
                globalMax.ScaleV(-1.0);
            }


            double[] imbalance = new double[N];
            for(int i = 0; i < N; i++) {
                double maxCost = globalMax[i];
                double minCost = globalMin[i];
                imbalance[i] = (maxCost - minCost) / Math.Max(double.Epsilon, maxCost);
            }

            return imbalance;
        }
    }
    
}