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
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution {

    /// <summary>
    /// In an MPI-parallel run, this class can be used to obtain a model of how
    /// load is distributed to MPI processes.
    /// </summary>
    public class RuntimeCellCostEstimator : ICellCostEstimator {

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
        private double[] ActualInstrumentationTimeStamps;

        /// <summary>
        /// Old values of <see cref="ActualInstrumentationTimeStamps"/>.
        /// </summary>
        private double[] LastInstrumentationTimeStamps;

        /// <summary>
        /// Difference between the sum of
        /// <see cref="ActualInstrumentationTimeStamps"/> and
        /// <see cref="LastInstrumentationTimeStamps"/>.
        /// </summary>
        public double EstimatedLocalCost {
            get;
            private set;
        }

        public int CurrentPerformanceClassCount {
            get;
            private set;
        }

        /// <summary>
        /// Constructor.
        /// </summary>
        public RuntimeCellCostEstimator(string[][] instrumentationPaths) {
            this.InstrumentationPaths = instrumentationPaths;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="performanceClassCount"></param>
        /// <param name="cellToPerformanceClassMap"></param>
        public void UpdateEstimates(int performanceClassCount, int[] cellToPerformanceClassMap) {
            CurrentPerformanceClassCount = performanceClassCount;
            currentCellToPerformanceClassMap = cellToPerformanceClassMap;
            int J = cellToPerformanceClassMap.Length;

            UpdateLocalTimes();
            CallCount++;

            // don't measure warm-up effects
            if (CallCount >= 3) {

                // Locally count cells for each performance class
                double[] cellCountPerClass = new double[performanceClassCount];
                for (int j = 0; j < J; j++) {
                    int performanceClass = cellToPerformanceClassMap[j];
                    cellCountPerClass[performanceClass]++;
                }

                PushData(cellCountPerClass, EstimatedLocalCost);
            }
        }

        /// <summary>
        /// Solves a linear regression cost model
        /// </summary>
        /// <returns></returns>
        public int[] GetEstimatedCellCosts() {
            MPICollectiveWatchDog.Watch();

            int noOfEstimates = localCellCounts.Count;
            Debug.Assert(noOfEstimates == localRunTimeEstimates.Count);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSize);

            MultidimensionalArray SendBuf = MultidimensionalArray.Create(
                noOfEstimates, CurrentPerformanceClassCount + 1);
            for (int i = 0; i < noOfEstimates; i++) {
                for (int j = 0; j < CurrentPerformanceClassCount; j++) {
                    SendBuf[i, j] = localCellCounts[i][j];
                }

                SendBuf[i, CurrentPerformanceClassCount] = localRunTimeEstimates[i];
            }

            MultidimensionalArray RecvBuf = MultidimensionalArray.Create(
                MpiSize, noOfEstimates, CurrentPerformanceClassCount + 1);

            unsafe {
                fixed (double* pSendBuf = &SendBuf.Storage[0], pRecvBuf = &RecvBuf.Storage[0]) {
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


            //// REMOVE IDENTICAL ROWS
            //List<int> obsoleteRows = new List<int>();
            //List<double> averages = new List<double>();
            //for (int row0 = 0; row0 < matrix.NoOfRows; row0++) {
            //    if (obsoleteRows.Contains(row0)) {
            //        continue;
            //    }

            //    List<int> identicalRows = new List<int>();
            //    for (int row1 = row0 + 1; row1 < matrix.NoOfRows; row1++) {
            //        bool rowIsIdentical = true;
            //        for (int col = 0; col < matrix.NoOfCols; col++) {
            //            rowIsIdentical &= (matrix[row0, col] == matrix[row1, col]);
            //        }

            //        if (rowIsIdentical) {
            //            identicalRows.Add(row1);
            //        }
            //    }

            //    obsoleteRows.AddRange(identicalRows);
            //    double average = (rhs[row0, 0] + identicalRows.Sum(i => rhs[i, 0])) / (1 + identicalRows.Count);
            //    averages.Add(average);
            //}

            //MultidimensionalArray reducedMatrix = MultidimensionalArray.Create(
            //    matrix.NoOfRows - obsoleteRows.Count(), matrix.NoOfCols);
            //MultidimensionalArray reducedRHS = MultidimensionalArray.Create(
            //    reducedMatrix.NoOfRows, 1);
            //int reducedRow = 0;
            //foreach (int row in Enumerable.Range(0, matrix.NoOfRows).Except(obsoleteRows)) {
            //    for (int col = 0; col < reducedMatrix.NoOfCols; col++) {
            //        reducedMatrix[reducedRow, col] = matrix[row, col];
            //    }
            //    reducedRHS[reducedRow, 0] = averages[reducedRow];

            //    reducedRow++;
            //}

            

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
                return cellCosts;
            } else {
                return null;
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

        private List<double[]> localCellCounts = new List<double[]>();

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

        private int CallCount = 0;

        private int[] currentCellToPerformanceClassMap;
    }
}