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
    /// Helper class for LTS: Performs Adams-Bashforth time integration on a specific sub-grid
    /// </summary>
    public class ABevolve : AdamsBashforth {

        /// <summary>
        /// History of the DG Coordiantes
        /// </summary>
        public Queue<double[]> HistoryDGCoordinate {
            get;
            private set;
        }

        /// <summary>
        /// Subgrid of all cells belonging to the cluster
        /// </summary>
        internal protected SubGrid ABSubGrid;

        /// <summary>
        /// Time history of each cell
        /// </summary>
        internal protected Queue<double[]> historyTimePerCell;

        /// <summary>
        /// Coefficients for Adams-Bashforth time integration of each cell
        /// </summary>
        internal double[][] ABCoefficientsPerCell;

        /// <summary>
        /// Number of edges
        /// </summary>
        private int numOfEdges;

        /// <summary>
        /// History of the fluxes across cell boundaries
        /// </summary>
        public Queue<double[]> HistoryBoundaryFluxes {
            get;
            private set;
        }

        /// <summary>
        /// Complete boundary fluxes across cell boundaries
        /// </summary>
        public double[] CompleteBoundaryFluxes {
            get;
            private set;
        }

        /// <summary>
        /// Current boundary fluxes
        /// </summary>
        private double[] currentBndFluxes;

        /// <summary>
        /// Stores the local cell indices of each boundary sub-grid.
        /// Avoids MPI communication within each call.
        /// </summary>
        protected int[] jSub2jCell;

        /// <summary>
        /// Adaptive mode is used if <see cref="AdamsBashforthLTS.reclusteringInterval"/> != 0
        /// </summary>
        internal bool adaptive;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap"><see cref="CoordinateMapping"/> of the fields</param>
        /// <param name="Parameters">Optional parameter fields, can be null if <paramref name="spatialOp"/> contains no parameters. Must match the parameter field list of <paramref name="spatialOp"/>, see <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="order">Order of the LTS algorithm</param>
        /// <param name="adaptive"><see cref="adaptive"/></param>
        /// <param name="sgrd">Sub-grid in which the local Adams-Bashforth step is evaluated</param>
        /// <remarks>Result of the local sub-step is saved in historyDGC, not directly in m_DGCoordinates</remarks>
        public ABevolve(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, bool adaptive = false, SubGrid sgrd = null)
                : base(spatialOp, Fieldsmap, Parameters, order, null, sgrd) {
            this.ABSubGrid = sgrd;
            HistoryDGCoordinate = new Queue<double[]>(order);
            RungeKuttaScheme = null; // Instance of RungeKutta not needed 
            jSub2jCell = sgrd.SubgridIndex2LocalCellIndex;

            numOfEdges = Fieldsmap.Fields.First().GridDat.iGeomEdges.Count;
            HistoryBoundaryFluxes = new Queue<double[]>(order - 1);
            CompleteBoundaryFluxes = new double[numOfEdges * Fieldsmap.Count];

            historyTimePerCell = new Queue<double[]>(order - 1);

            this.adaptive = adaptive;
        }


        /// <summary>
        /// Performs a sub-step
        /// </summary>
        /// <param name="dt">size of time step</param>
        public override double Perform(double dt) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {
                // Checking History
                if (HistoryDGCoordinate.Count >= order)
                    HistoryDGCoordinate.Dequeue();
                if (HistoryChangeRate.Count >= order) {
                    HistoryChangeRate.Dequeue();
                    HistoryBoundaryFluxes.Dequeue();
                }

                // Compute AB Coefficents                
                if (adaptive) {
                    if (historyTimePerCell.Count > order - 1)
                        historyTimePerCell.Dequeue();

                    ABCoefficientsPerCell = new double[ABSubGrid.LocalNoOfCells][];
                    for (int cell = 0; cell < ABSubGrid.LocalNoOfCells; cell++) {
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

                double[] upDGC;
                CurrentChangeRate = new double[Mapping.LocalLength];
                currentBndFluxes = new double[numOfEdges * Mapping.Fields.Count];
                ComputeChangeRate(CurrentChangeRate, m_Time, 0, currentBndFluxes);

                ////////////// Calculate complete change rate with older steps
                MakeABStep();

                // Update DGCoordinates for History
                // (Hint: calls extra function to gives the ability to modify the update procedure, i.e in IBM case)
                upDGC = ComputesUpdatedDGCoordinates(CompleteChangeRate);

                // Keeps track of histories
                HistoryDGCoordinate.Enqueue(upDGC);
                HistoryChangeRate.Enqueue(CurrentChangeRate);
                HistoryBoundaryFluxes.Enqueue(currentBndFluxes);
                UpdateTimeHistory(dt);

                // Update local sub-grid time
                m_Time = m_Time + dt;
            }
            return dt;
        }

        /// <summary>
        /// Perfoms the actual Adams-Bashforth time integration sub-step in a cluster
        /// </summary>
        protected virtual void MakeABStep() {
            CompleteChangeRate = new double[Mapping.LocalLength];

            if (adaptive) {
                CompleteChangeRate = new double[Mapping.LocalLength];
                for (int j = 0; j < ABSubGrid.LocalNoOfCells; j++) {
                    int cell = jSub2jCell[j];
                    // cell = global cell index
                    // f = each field
                    // n = basis polynomial
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

                BLAS.daxpy(CompleteBoundaryFluxes.Length, ABCoefficients[0], currentBndFluxes, 1, CompleteBoundaryFluxes, 1);

                // calculate completeChangeRate
                int i = 1;
                foreach (double[] oldRate in HistoryChangeRate) {
                    BLAS.daxpy(CompleteChangeRate.Length, ABCoefficients[i], oldRate, 1, CompleteChangeRate, 1);
                    i++;
                }

                i = 1;
                foreach (double[] oldRate in HistoryBoundaryFluxes) {
                    BLAS.daxpy(CompleteBoundaryFluxes.Length, ABCoefficients[i], oldRate, 1, CompleteBoundaryFluxes, 1);
                    i++;
                }
            }
        }

        /// <summary>
        /// Updates the time history of the entire cluster (if <see cref="adaptive"/>== false), or
        /// of each cell (if <see cref="adaptive"/>== true)
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void UpdateTimeHistory(double dt) {
            if (adaptive) {
                double[] currentTime = new double[ABSubGrid.LocalNoOfCells];
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
        /// Important: It does not change the DGCoordinates <see cref="ExplicitEuler.DGCoordinates"/>.
        /// </summary>
        /// <param name="completeChangeRate">Complete ChangeRate of a cluster for one sub-step</param>
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
        /// <param name="results">Result for the complete grid</param>
        /// <returns>Array with entries only for the sub-grid cells</returns>
        private double[] OrderValuesBySgrd(double[] results) {
            double[] ordered = new double[Mapping.LocalLength];

            for (int j = 0; j < ABSubGrid.LocalNoOfCells; j++) {
                int cell = jSub2jCell[j];
                // cell in the sub-grid
                // f = each field
                // n = basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        ordered[index] = results[index];
                    }
                }
            }
            return ordered;
        }
    }
}
