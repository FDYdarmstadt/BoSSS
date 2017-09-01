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

using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System.Linq;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Helper Class for LTS: Performs the local Adams-Bashforth time integration on a specific sub-grid
    /// </summary>
    public class ABevolve : AdamsBashforth {

        /// <summary>
        /// Stores History of DGCoordinates
        /// </summary>
        public Queue<double[]> HistoryDGCoordinate {
            get;
            private set;
        }

        internal protected SubGrid sgrd;

        internal protected Queue<double[]> historyTimePerCell;

        internal double[][] ABCoefficientsPerCell;

        private int numOfEdges;

        public Queue<double[]> HistoryBndFluxes {
            get;
            private set;
        }

        public double[] completeBndFluxes;

        private double[] currentBndFluxes;

        /// <summary>
        /// Stores the local cell indices of each boundary sub-grid.
        /// Avoids MPI communication within each call
        /// </summary>
        protected int[] jSub2jCell;

        internal bool adaptive;

        public Queue<double[]> HistoryABCoefficients {
            get;
            private set;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters">
        /// optional parameter fields, can be null if <paramref name="spatialOp"/> contains no parameters;
        /// must match the parameter field list of <paramref name="spatialOp"/>, see <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="order">Same order than for LTS</param>
        /// <param name="sgrd">Sub-grid in which the local Adams-Bashforth step is evaluated</param>
        /// <remarks>The result of the local sub-step is saved in historyDGC, not directly in m_DGCoordinates</remarks>
        public ABevolve(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, bool adaptive = false, SubGrid sgrd = null)
                : base(spatialOp, Fieldsmap, Parameters, order, null, sgrd) {
            this.sgrd = sgrd;
            HistoryDGCoordinate = new Queue<double[]>(order);
            RungeKuttaScheme = null; // Instance of RungeKutta not needed 
            jSub2jCell = sgrd.SubgridIndex2LocalCellIndex;

            numOfEdges = Fieldsmap.Fields.First().GridDat.iGeomEdges.Count;
            HistoryBndFluxes = new Queue<double[]>(order - 1);
            completeBndFluxes = new double[numOfEdges * Fieldsmap.Count];

            historyTimePerCell = new Queue<double[]>(order - 1);

            this.adaptive = adaptive;
        }


        /// <summary>
        /// Performs one local sub-step
        /// </summary>
        /// <param name="dt">size of time step</param>
        public override double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {
                // Checking History
                if (HistoryDGCoordinate.Count >= order)
                    HistoryDGCoordinate.Dequeue();
                if (HistoryChangeRate.Count >= order) {
                    HistoryChangeRate.Dequeue();
                    HistoryBndFluxes.Dequeue();
                }

                // Compute AB Coefficents                
                if (adaptive) {
                    if (historyTimePerCell.Count > order - 1)
                        historyTimePerCell.Dequeue();

                    ABCoefficientsPerCell = new double[sgrd.LocalNoOfCells][];
                    for (int cell = 0; cell < sgrd.LocalNoOfCells; cell++) {
                        double[] historyTimeArray = new double[order];
                        int i = 0;
                        foreach (double[] historyPerCell in historyTimePerCell) {
                            historyTimeArray[i] = historyPerCell[cell];
                            i++;
                        }
                        historyTimeArray[i] = m_Time;
                        ABCoefficientsPerCell[cell] = ComputeCoefficients(dt, historyTimeArray);
                    }
                } else {
                    double[] historyTimeArray = HistoryTime.ToArray();
                    ABCoefficients = ComputeCoefficients(dt, historyTimeArray);
                }

                CurrentChangeRate = new double[Mapping.LocalLength];
                currentBndFluxes = new double[numOfEdges * Mapping.Fields.Count];
                ComputeChangeRate(CurrentChangeRate, m_Time, 0, currentBndFluxes);

                ////////////// Calculate complete change rate with older steps
                MakeABStep();

                // Update DGCoordinates for History
                // (Hint: calls extra function to gives the ability to modify the update procedure, i.e in IBM case)
                double[] upDGC = ComputesUpdatedDGCoordinates(CompleteChangeRate);

                // Keeps track of histories
                HistoryDGCoordinate.Enqueue(upDGC);
                HistoryChangeRate.Enqueue(CurrentChangeRate);
                HistoryBndFluxes.Enqueue(currentBndFluxes);

                UpdateTimeHistory(dt);

                //HistoryDGCoordinate.Dequeue();
                //HistoryChangeRate.Dequeue();
                //HistoryBndFluxes.Dequeue();
                //if (adaptive)
                //    historyTimePerCell.Dequeue();

                // Update local sub-grid time
                m_Time = m_Time + dt;
            }
            return dt;
        }

        protected virtual void MakeABStep() {
            CompleteChangeRate = new double[Mapping.LocalLength];

            if (adaptive) {
                CompleteChangeRate = new double[Mapping.LocalLength];
                for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                    int cell = jSub2jCell[j];
                    // cell in sgrd
                    // f== each field
                    // n== basis polynomial
                    foreach (DGField f in Mapping.Fields) {
                        for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                            int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                            CompleteChangeRate[index] = CompleteChangeRate[index] + ABCoefficientsPerCell[j][0] * CurrentChangeRate[index];
                            int i = 1;
                            foreach (double[] oldRate in HistoryChangeRate) {
                                CompleteChangeRate[index] = CompleteChangeRate[index] + ABCoefficientsPerCell[j][i] * oldRate[index];
                                i++;
                            }
                        }
                    }
                }
            } else {
                //y  <--  alpha*x + y
                BLAS.daxpy(CompleteChangeRate.Length, ABCoefficients[0], CurrentChangeRate, 1, CompleteChangeRate, 1);

                BLAS.daxpy(completeBndFluxes.Length, ABCoefficients[0], currentBndFluxes, 1, completeBndFluxes, 1);

                // calculate completeChangeRate
                int i = 1;
                foreach (double[] oldRate in HistoryChangeRate) {
                    BLAS.daxpy(CompleteChangeRate.Length, ABCoefficients[i], oldRate, 1, CompleteChangeRate, 1);
                    i++;
                }

                i = 1;
                foreach (double[] oldRate in HistoryBndFluxes) {
                    BLAS.daxpy(completeBndFluxes.Length, ABCoefficients[i], oldRate, 1, completeBndFluxes, 1);
                    i++;
                }
            }
        }

        protected virtual void UpdateTimeHistory(double dt) {
            if (adaptive) {
                double[] currentTime = new double[sgrd.LocalNoOfCells];
                for (int i = 0; i < currentTime.Length; i++) {
                    currentTime[i] = m_Time;
                }
                historyTimePerCell.Enqueue(currentTime);
                if (historyTimePerCell.Count > order - 1) { //eventuell nicht nötig, später mal überprüfen
                    historyTimePerCell.Dequeue();
                }
            } else {
                HistoryTime.Enqueue(m_Time + dt); // --> später noch schön machen und ähnlich wie timeHistoryPerCell, einen eintrag weglassen, weil das eh m_time ist!
                if (HistoryTime.Count >= order) //eventuell nicht nötig, später mal überprüfen
                    HistoryTime.Dequeue();
            }
        }

        /// <summary>
        /// Computes the new intermediate DGCoordinates. 
        /// Important: It doesn't changes the DGCoordinates <see cref="ExplicitEuler.DGCoordinates"/>.
        /// </summary>
        /// <param name="completeChangeRate">complete ChangeRate of a sub-grid for one sub-step</param>
        /// <returns>intermediate DGCoordinates as array</returns>
        protected virtual double[] ComputesUpdatedDGCoordinates(double[] completeChangeRate) {
            // Standard case: Just add completeChangeRate to DGCoordinates as array
            double[] upDGC = new double[Mapping.LocalLength];
            DGCoordinates.CopyTo(upDGC, 0);
            upDGC = OrderValuesBySgrd(upDGC);
            BLAS.daxpy(upDGC.Length, -1, OrderValuesBySgrd(completeChangeRate), 1, upDGC, 1);
            return upDGC;
        }

        /// <summary>
        /// Takes a double[] with results for the global grid and gives an
        /// array with only entries for the specific sub-grid
        /// </summary>
        /// <param name="results">Result for the complete Grid</param>
        /// <returns>Array with entries only for the sgrd-cells</returns>
        private double[] OrderValuesBySgrd(double[] results) {
            double[] ordered = new double[Mapping.LocalLength];

            for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                int cell = jSub2jCell[j];
                // cell in sgrd
                // f== each field
                // n== basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        ordered[index] = results[index];
                    }
                }
            }
            return ordered;
        }


        private double[] ZerosValuesOfSgrd(double[] results) {
            double[] ordered = new double[Mapping.LocalLength];
            Array.Copy(results, ordered, results.Length);

            for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                int cell = jSub2jCell[j];
                // cell in sgrd
                // f== each field
                // n== basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        ordered[index] = 0.0;
                    }
                }
            }
            return ordered;
        }
    }
}
