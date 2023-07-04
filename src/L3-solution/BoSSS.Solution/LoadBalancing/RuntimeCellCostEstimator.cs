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
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution.LoadBalancing {

    /// <summary>
    /// In an MPI-parallel run, this class can be used to obtain a model of how
    /// load is distributed to MPI processes.
    /// </summary>
    [Serializable]
    public class RuntimeCellCostEstimator : CellTypeBasedEstimator {

        /// <summary>
        /// <see cref="MaxNoOfTimesteps"/>
        /// </summary>
        private int m_MaxNoOfTimesteps = 20;

        /// <summary>
        /// Number of time-steps which are used to obtain the performance model.
        /// </summary>
        public int MaxNoOfTimesteps {
            get {
                return m_MaxNoOfTimesteps;
            }
        }

        /// <summary>
        /// Call paths into the instrumentation tree (see
        /// <see cref="Tracer.Root"/>), from which the runtime model will be
        /// derived.
        /// </summary>
        private string[][] InstrumentationPaths;

        /// <summary>
        /// Runtime measurements for the paths in
        /// <see cref="InstrumentationPaths"/>, for the most recent call to
        /// <see cref="UpdateLocalTimes"/>.
        /// </summary>
        [NonSerialized]
        private double[] ActualInstrumentationTimeStamps;

        /// <summary>
        /// Old values of <see cref="ActualInstrumentationTimeStamps"/>.
        /// </summary>
        [NonSerialized]
        private double[] LastInstrumentationTimeStamps;


        int CurrentPerformanceClassCount;

        /// <summary>
        /// Constructor.
        /// </summary>
        public RuntimeCellCostEstimator(string[][] instrumentationPaths) {
            this.InstrumentationPaths = instrumentationPaths;
        }

        /// <summary>
        /// 
        /// </summary>
        override public void UpdateEstimates(IApplication app) {

            var cellToPerformanceClassMap = base.CellClassifier.ClassifyCells(app);
            int performanceClassCount = cellToPerformanceClassMap.Max().MPIMax();

            CurrentPerformanceClassCount = performanceClassCount + 1;
            currentCellToPerformanceClassMap = cellToPerformanceClassMap;
            int J = m_app.GridData.CellPartitioning.LocalLength;

            UpdateLocalTimes();
            CallCount++;

            // don't measure warm-up effects
            if (CallCount >= 3) {

                // Locally count cells for each performance class
                double[] cellCountPerClass = new double[CurrentPerformanceClassCount];
                for (int j = 0; j < J; j++) {
                    int performanceClass = cellToPerformanceClassMap[j];
                    cellCountPerClass[performanceClass]++;
                }

                PushData(cellCountPerClass, EstimatedLocalCost);
            }
        }

        [NonSerialized]
        double EstimatedLocalCost;

        /// <summary>
        /// Solves a linear regression cost model
        /// </summary>
        /// <returns></returns>
        override public int[][] GetEstimatedCellCosts() {
            MPICollectiveWatchDog.Watch();

            int noOfEstimates = localCellCounts.Count;
            Debug.Assert(noOfEstimates == localRunTimeEstimates.Count);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSize);

            MultidimensionalArray SendBuf = MultidimensionalArray.Create(noOfEstimates, CurrentPerformanceClassCount + 1);
            for (int i = 0; i < noOfEstimates; i++) {
                for (int j = 0; j < CurrentPerformanceClassCount; j++) {
                    SendBuf[i, j] = localCellCounts[i][j];
                }

                SendBuf[i, CurrentPerformanceClassCount] = localRunTimeEstimates[i];
            }

            MultidimensionalArray RecvBuf = MultidimensionalArray.Create(MpiSize, noOfEstimates, CurrentPerformanceClassCount + 1);

            unsafe {
                fixed (double* pSendBuf = SendBuf.Storage, pRecvBuf = RecvBuf.Storage) {
                    csMPI.Raw.Allgather(
                        (IntPtr)pSendBuf, SendBuf.Length, csMPI.Raw._DATATYPE.DOUBLE,
                        (IntPtr)pRecvBuf, SendBuf.Length, csMPI.Raw._DATATYPE.DOUBLE,
                        csMPI.Raw._COMM.WORLD);
                }
            }

#if DEBUG
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            for (int i = 0; i < noOfEstimates; i++) {
                for (int j = 0; j < CurrentPerformanceClassCount; j++) {
                    Debug.Assert(RecvBuf[rank, i, j] == SendBuf[i, j]);
                    Debug.Assert(RecvBuf[rank, i, j] == localCellCounts[i][j]);
                }

                Debug.Assert(RecvBuf[rank, i, CurrentPerformanceClassCount] == SendBuf[i, CurrentPerformanceClassCount]);
                Debug.Assert(RecvBuf[rank, i, CurrentPerformanceClassCount] == localRunTimeEstimates[i]);
            }
#endif
            
            // TO DO: MAKE EFFICIENT
            int noOfRows = MpiSize * noOfEstimates;
            int noOfColumns = CurrentPerformanceClassCount;
            MultidimensionalArray matrix = MultidimensionalArray.Create(noOfRows, noOfColumns);
            MultidimensionalArray rhs = MultidimensionalArray.Create(noOfRows, 1);
            for (int i = 0; i < MpiSize; i++) {
                for (int j = 0; j < noOfEstimates; j++) {
                    for (int k = 0; k < CurrentPerformanceClassCount; k++) {
                        matrix[i * noOfEstimates + j, k] = RecvBuf[i, j, k];
                    }

                    rhs[i * noOfEstimates + j, 0] = RecvBuf[i, j, CurrentPerformanceClassCount];
                }
            }

            

            if (MpiSize > 1 && rhs.Length > MpiSize) {
                MultidimensionalArray costs = MultidimensionalArray.Create(CurrentPerformanceClassCount, 1);
                matrix.LeastSquareSolve(costs, rhs);
                //reducedMatrix.LeastSquareSolve(costs, reducedRHS);

                double minCost = costs.Min();
                double maxCost = costs.Max();
                int[] intCost = costs.Storage.Select(
                    cost => Math.Max(1, (int)Math.Round(1e3 * (cost / maxCost)))).ToArray();

                for (int i = 0; i < CurrentPerformanceClassCount; i++) {
                    Console.WriteLine(
                        "Cost of cell type {0}: \tabs:{1:0.##E-00} \trel1: {2:0.0%}\trel2: {3:0.0%}\tint: {4}",
                        i,
                        costs[i, 0],
                        costs[i, 0] / minCost,
                        costs[i, 0] / maxCost,
                        intCost[i]);
                }

                int[] cellCosts = new int[currentCellToPerformanceClassMap.Length];
                for (int j = 0; j < cellCosts.Length; j++) {
                    cellCosts[j] = intCost[currentCellToPerformanceClassMap[j]];
                }
                return new int[][] { cellCosts };
            } else {
                return new int[0][];
            }
        }

        /// <summary>
        /// Extracts the total time of some call path (in the application stack)#
        /// from the tracing infrastructure, <see cref="Tracer.Root"/>.
        /// </summary>
        private double GetTime(string[] Path) {
            return GetTimeRecursive(new MethodCallRecord[] { Tracer.Root }, Path);
        }

        private double GetTimeRecursive(IEnumerable<MethodCallRecord> mcrS, string[] Path) {
            double Time = 0.0;
            foreach (var mcr in mcrS) {
                if (Path.Length <= 0) {
                    Time += mcr.TimeSpentInMethod.TotalSeconds;
                } else {
                    var childs = mcr.FindChildren(Path[0]);
                    Time += GetTimeRecursive(childs, Path.Skip(1).ToArray());
                }
            }

            return Time;
        }

        private void UpdateLocalTimes() {
            int I = InstrumentationPaths.Length;

            LastInstrumentationTimeStamps = ActualInstrumentationTimeStamps;
            ActualInstrumentationTimeStamps = new double[I];

            for (int i = 0; i < I; i++) {
                ActualInstrumentationTimeStamps[i] = GetTime(InstrumentationPaths[i]);
            }

            EstimatedLocalCost = ActualInstrumentationTimeStamps.Sum();
            if (LastInstrumentationTimeStamps != null) {
                EstimatedLocalCost -= LastInstrumentationTimeStamps.Sum();
            }
        }

        [NonSerialized]
        private List<double[]> localCellCounts = new List<double[]>();

        [NonSerialized]
        private List<double> localRunTimeEstimates = new List<double>();
        
        /// <summary>
        /// Pushes the latest timing measurements
        /// </summary>
        private void PushData(double[] cellCountPerClass, double cost) {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSize);

            localCellCounts.Add(cellCountPerClass);
            localRunTimeEstimates.Add(cost);

            while (localCellCounts.Count > MaxNoOfTimesteps) {
                Debug.Assert(localCellCounts.Count == localRunTimeEstimates.Count);
                localCellCounts.RemoveAt(0);
                localRunTimeEstimates.RemoveAt(0);
            }
        }

        [NonSerialized]
        private int CallCount = 0;

        [NonSerialized]
        private int[] currentCellToPerformanceClassMap;

        /// <summary>
        /// Serialization constructor
        /// </summary>
        private RuntimeCellCostEstimator() {

        }


        public override object Clone() {
            return new RuntimeCellCostEstimator() {
                CellClassifier = this.CellClassifier.CloneAs(),
                m_MaxNoOfTimesteps = this.m_MaxNoOfTimesteps,
                InstrumentationPaths = this.InstrumentationPaths.Select(sa => sa.Select(s => s.CloneAs()).ToArray()).ToArray()
            };
        }
    }
}