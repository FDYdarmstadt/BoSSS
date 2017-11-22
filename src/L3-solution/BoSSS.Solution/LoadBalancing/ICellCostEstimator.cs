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
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution {

    /// <summary>
    /// Provides a measure to estimate the runtime cost of a given cell in
    /// the computational grid
    /// </summary>
    public interface ICellCostEstimator {

        /// <summary>
        /// Updates <see cref="EstimatedLocalCost"/> and
        /// <see cref="GetEstimatedCellCosts"/> for the given mapping
        /// cells and performance classes
        /// </summary>
        /// <param name="performanceClassCount"></param>
        /// <param name="cellToPerformanceClassMap"></param>
        void UpdateEstimates(int performanceClassCount, int[] cellToPerformanceClassMap);

        /// <summary>
        /// The total number of performance classes
        /// </summary>
        int CurrentPerformanceClassCount {
            get;
        }

        /// <summary>
        /// The estimated total cost of all cells on this process
        /// </summary>
        double EstimatedLocalCost {
            get;
        }

        /// <summary>
        /// The estimated cost of each individual cell on this process
        /// </summary>
        /// <returns></returns>
        int[] GetEstimatedCellCosts();
    }

    /// <summary>
    /// Extension methods for <see cref="ICellCostEstimator"/>
    /// </summary>
    public static class ICellCostEstimatorExtensions {

        /// <summary>
        /// Estimates the cost imbalance between all processes using the given
        /// <paramref name="estimator"/>
        /// </summary>
        /// <param name="estimator"></param>
        /// <returns>
        /// A number betweeen 0 and 1 which represents an estimate of the load
        /// imbalance in percent
        /// </returns>
        public static double ImbalanceEstimate(this ICellCostEstimator estimator) {
            MPICollectiveWatchDog.Watch();

            double localCost = estimator.EstimatedLocalCost;
            double[] allCosts = localCost.MPIAllGather();

            double minCost = allCosts.Min();
            double maxCost = allCosts.Max();
            double imbalance = (maxCost - minCost) / Math.Max(double.Epsilon, maxCost);

            return imbalance;
        }
    }
}