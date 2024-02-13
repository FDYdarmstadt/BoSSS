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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution.LoadBalancing {

    /// <summary>
    /// In an MPI-parallel run, this class can be used to create an MPI
    /// partition of the grid according a given cell cost estimation, so that
    /// the compute load is balanced evenly among the MPI processors.
    /// </summary>
    public class LoadBalancer {

        ///// <summary>
        ///// A factory used to update
        ///// <see cref="CurrentCellCostEstimators"/> if required
        ///// </summary>
        //private List<Func<IApplication, int, ICellCostEstimator>> cellCostEstimatorFactories;

        /// <summary>
        /// A set of models that estimates the costs of different cells
        /// </summary>
        public ICellCostEstimator[] CellCostEstimators {
            get;
            private set;
        }

        /// <summary>
        /// Indicates whether load balancer has already suggest a new partitioning
        /// before. If this value is true, the load balancer is allowed to
        /// rearrange cells more aggressively
        /// </summary>
        private bool isFirstRepartitioning = true;

        /// <summary>
        /// Constructor.
        /// </summary>
        public LoadBalancer(IEnumerable<ICellCostEstimator> cellCostEstimators, IApplication app) {
            this.CellCostEstimators = cellCostEstimators.ToArray();
            foreach (var e in cellCostEstimators)
                e.Init(app);
        }

        /// <summary>
        /// Returns a new grid partition based on the performance model.
        /// </summary>
        /// <param name="app"></param>
        /// <param name="TimestepNo"></param>
        /// <param name="gridPartType">Grid partitioning type.</param>
        /// <param name="PartOptions"></param>
        /// <param name="imbalanceThreshold">
        /// See <see cref="Control.AppControl.DynamicLoadBalancing_ImbalanceThreshold"/>.
        /// </param>
        /// <param name="Period">
        /// See <see cref="Control.AppControl.DynamicLoadBalancing_Period"/>.
        /// </param>
        /// <param name="redistributeAtStartup"></param>
        /// <param name="TimestepNoRestart"></param>
        /// <returns></returns>
        public int[] GetNewPartitioning(IApplication app, int TimestepNo, GridPartType gridPartType, string PartOptions, double imbalanceThreshold, int Period, bool redistributeAtStartup, TimestepNumber TimestepNoRestart) {
            // Create new model if number of cell classes has changed
            using (var tr = new FuncTrace()) {
                //tr.InfoToConsole = true;
                tr.Info($"Computing Partition of using {gridPartType}...");
               
                bool performPertationing;

                if (TimestepNo == 0 || (TimestepNoRestart != null && TimestepNo == TimestepNoRestart.MajorNumber)) {
                    tr.Info("redistributeAtStartup = " + redistributeAtStartup);
                    tr.Info("TimestepNo = " + TimestepNo);
                    tr.Info("TimestepNoRestart != null = " + (TimestepNoRestart != null));
                    tr.Info("TimestepNoRestart.MajorNumber = " + TimestepNoRestart.MajorNumber);

                    performPertationing = redistributeAtStartup;
                } else {
                    tr.Info("Period = " + Period);
                    tr.Info("TimestepNo % Period = " + (TimestepNo % Period));

                    performPertationing = (Period > 0 && TimestepNo % Period == 0);
                }

                if (!performPertationing) {
                    tr.Info("No new partition will be computed.");
                    return null;
                }

                bool imbalanceTooLarge = CheckImbalance(app, imbalanceThreshold);
                if (!imbalanceTooLarge) {
                    tr.Info("Imbalance is not sufficiently high for load balancing!");
                    return null;
                }

                IList<int[]> allCellCosts = new List<int[]>();
                foreach (var estimator in CellCostEstimators) {
                    int J = app.Grid.CellPartitioning.LocalLength;

                    var costs = estimator.GetEstimatedCellCosts();
                    foreach (int[] cc in costs) {
                        if (cc.Length != J) {
                            throw new ApplicationException($"Illegal cell cost list returned by cell cost estimator {estimator}; expecting a length of {J}, but got {cc?.Length}");
                        }
                    }
                    allCellCosts.AddRange(costs);
                }
                if (allCellCosts == null || allCellCosts.All(c => c == null)) {
                    return null;
                }


                int[] result;
                switch (gridPartType) {
                    case GridPartType.METIS:
                        int.TryParse(PartOptions, out int noOfPartitioningsToChooseFrom);
                        noOfPartitioningsToChooseFrom = Math.Max(1, noOfPartitioningsToChooseFrom);
                        result = ((GridCommons)(app.Grid)).ComputePartitionMETIS(allCellCosts);
                        isFirstRepartitioning = false;
                        break;

                    case GridPartType.ParMETIS:
                        // Do full ParMETIS run on first repartitioning since
                        // initial partitioning may be _really_ bad
                        if (isFirstRepartitioning) {
                            result = ((GridCommons)(app.Grid)).ComputePartitionParMETIS(allCellCosts);
                            isFirstRepartitioning = false;
                        } else {
                            // Refinement currently deactivate because it behaves
                            // strangely when large numbers of cells should be
                            // repartitioned
                            //result = Grid.ComputePartitionParMETIS(cellCosts, refineCurrentPartitioning: true);
                            result = ((GridCommons)(app.Grid)).ComputePartitionParMETIS(allCellCosts);
                        }
                        break;

                    case GridPartType.clusterHilbert:
                        return ((GridCommons)(app.Grid)).ComputePartitionHilbert(localcellCosts: allCellCosts, Functype: 0);

                    case GridPartType.Hilbert:
                        return ((GridCommons)(app.Grid)).ComputePartitionHilbert(localcellCosts: allCellCosts, Functype: 1);

                    case GridPartType.none:
                        result = IndexBasedPartition(allCellCosts.Single());
                        break;

                    case GridPartType.OtherSession:
                    case GridPartType.Predefined:
                        return null;

                    default:
                        throw new NotImplementedException();
                }


                if (result.Length != app.Grid.CellPartitioning.LocalLength) {
                    throw new ApplicationException($"Load balancing '{gridPartType}' computed invalid partitioning on rank {app.Grid.MyRank}; length mismatch ({result.Length} vs {app.Grid.CellPartitioning.LocalLength}).");
                }

                {
                    int myRank = app.Grid.MyRank;
                    long ToOther = 0;
                    for (int j = 0; j < result.Length; j++) {
                        if (result[j] != myRank)
                            ToOther++;
                    }

                    long ToOtherGlobal = ToOther.MPISum();
                    tr.Info("Number of cells which go to other processors: " + ToOther);
                }

                //int[] CurrentPart = new int[result.Length];
                //CurrentPart.SetAll(app.MPIRank);
                //CurrentPart.SaveToTextFile($"OldPart{counter}.txt");
                //result.SaveToTextFile($"NewPart{counter}.txt");


                return result;
            }
        }

        /// <summary>
        /// Checks the imbalance w.r.t. pre-defined estimators
        /// </summary>
        /// <param name="app"></param>
        /// <param name="imbalanceThreshold"></param>
        /// <returns> False: Imbalance is tolerable. True: Imbalance is too big.</returns>
        public bool CheckImbalance(IApplication app, double imbalanceThreshold) {
            using (var tr = new FuncTrace()) {
#if DEBUG
                tr.InfoToConsole = true;
#endif
                var imbalanceEstimates = new List<double>();
                foreach (var estimator in CellCostEstimators) {
                    estimator.UpdateEstimates(app);
                    imbalanceEstimates.AddRange(estimator.ImbalanceEstimate());
                }
                bool imbalanceTooLarge = false;


                for (int i = 0; i < imbalanceEstimates.Count; i++) {
                    if (imbalanceEstimates[i].IsNaNorInf())
                        throw new ArithmeticException("NaN/Inf in imbalance estimate: " + imbalanceEstimates[i]);
                    imbalanceTooLarge |= (imbalanceEstimates[i] > imbalanceThreshold);
                }

                if (imbalanceTooLarge) {
                    tr.Info( $"At least one runtime imbalance estimate ({String.Join(", ", imbalanceEstimates.Select(e => String.Format("{0:P1}", e)))}) was above configured threshold ({imbalanceThreshold:P1}); attempting repartitioning");
                } else {
                    tr.Info($"All imbalance estimate ({String.Join(", ", imbalanceEstimates.Select(e => String.Format("{0:P1}", e)))}) are within threshold ({imbalanceThreshold:P1}); no repartitioning necessary");
                }
                return imbalanceTooLarge;
            }
        }

        static int[] IndexBasedPartition(int[] Cost) {
            int J = Cost.Length;
            int[] AccCost = Cost.CloneAs();

            for (int j = 0; j < J - 1; j++) {
                AccCost[j + 1] += AccCost[j];
            }

            int MpiRank, MpiSize;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MpiRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MpiSize);

            int[] GlobAccCost = new int[MpiSize];
            unsafe {
                fixed (int* pGlobAccCost = GlobAccCost) {
                    int locAccCost = AccCost[J - 1];

                    csMPI.Raw.Allgather(
                        (IntPtr)(&locAccCost), 1, csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pGlobAccCost, 1, csMPI.Raw._DATATYPE.INT,
                        csMPI.Raw._COMM.WORLD);
                }
            }

            for (int proc = 0; proc < MpiSize - 1; proc++) {
                GlobAccCost[proc + 1] += GlobAccCost[proc];
            }

            int glbOffset = MpiRank > 0 ? GlobAccCost[MpiRank - 1] : 0;
            for (int j = 0; j < J; j++) {
                AccCost[j] += glbOffset;
            }


            double TotalCost = GlobAccCost[MpiSize - 1];
            double[] Ranges = new double[MpiSize + 1];
            for (int proc = 1; proc <= MpiSize; proc++) {
                Ranges[proc] = (TotalCost * proc) / MpiSize;
            }
            Ranges[MpiSize] = Math.Max(Ranges[MpiSize], TotalCost);

            int[] R = new int[J];

            int rank = 0;
            for (int j = 0; j < J; j++) {
                Debug.Assert(j == 0 || AccCost[j - 1] < AccCost[j]);
                while (!(Ranges[rank] < AccCost[j] && AccCost[j] <= Ranges[rank + 1]) && (rank < MpiSize - 1)) {
                    rank++;
                }
                Debug.Assert(Ranges[rank] < AccCost[j] && AccCost[j] <= Ranges[rank + 1]);
                R[j] = rank;
            }

            return R;
        }
    }
}
