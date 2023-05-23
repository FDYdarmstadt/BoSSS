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
    /// Provides a measure to estimate the runtime cost of a given cell in
    /// the computational grid
    /// </summary>
    public interface ICellCostEstimator : ICloneable {

        void Init(IApplication app);

        /// <summary>
        /// Updates the internal cell costs
        /// </summary>
        void UpdateEstimates(IApplication app);

      /*
        /// <summary>
        /// The estimated total cost of all cells on this process
        /// </summary>
        double EstimatedLocalCost {
            get;
        }
      */
        /// <summary>
        /// The estimated cost of each individual cell on this process;
        /// 
        /// 
        /// </summary>
        /// <returns>
        /// cell weights for multi-constraint partitioning, where multiple weights are assigned to each cells.
        /// - 1st index: constraint index (only 0 for a single constraint)
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