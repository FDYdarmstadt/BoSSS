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
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using MPI.Wrappers;

namespace BoSSS.Solution {
    /// <summary>
    /// Abstract bass class for the implementation of various time step constraints, 
    /// e.g. CFL condition, capillary wave constraint
    /// </summary>
    public abstract class TimeStepConstraint {

        /// <summary>
        /// Information about the grid
        /// </summary>
        protected readonly IGridData gridData;

        /// <summary>
        /// Nodes for the evaluation of the time step constraint per reference element
        /// </summary>
        private NodeSet[] evaluationPoints;

        public double dtMin, dtMax, dtFraction, Endtime;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="gridData"></param>
        public TimeStepConstraint(IGridData gridData, double dtMin, double dtMax, double dtFraction, double Endtime) {
            this.gridData = gridData;
            this.dtMin = dtMin;
            this.dtMax = dtMax;
            this.dtFraction = dtFraction;
            this.Endtime = Endtime;
        }

        /// <summary>
        /// Standard evaluation points for the CFL criterion. Here, we choose
        /// the vertices of the reference element plus the cell center.
        /// </summary>
        public NodeSet[] EvaluationPoints {
            get {
                if (evaluationPoints == null) {
                    var KrefS = gridData.iGeomCells.RefElements;
                    int D = gridData.SpatialDimension;

                    evaluationPoints = new NodeSet[KrefS.Length];
                    for (int i = 0; i < KrefS.Length; i++) {
                        var Kref = KrefS[i];
                        int N = Kref.NoOfVertices + 1;

                        NodeSet vertices = new NodeSet(Kref, N, D);
                        vertices.SetSubArray(
                            Kref.Vertices,
                            new int[] { 0, 0 },
                            new int[] { N - 2, D - 1 });
                        vertices.LockForever();

                        evaluationPoints[i] = vertices;
                    }
                }

                return evaluationPoints;
            }
        }


        /// <summary>
        /// Determines the admissible step-size within
        /// <paramref name="subGrid"/> by calling
        /// <see cref="GetLocalStepSize"/> for all contiguous chunks of cells.
        /// </summary>
        /// <param name="subGrid">
        /// The sub-grid for which the step size shall be determined. If null
        /// is given, the full domain will be considered.
        /// </param>
        /// <returns>
        /// The maximum step size dictated by the CFL restriction in the given
        /// <paramref name="subGrid"/>.
        /// </returns>
        /// <remarks>
        /// This is a collective call which requires all processes to
        /// synchronize.
        /// </remarks>
        public double GetGloballyAdmissibleStepSize(SubGrid subGrid = null) {
            MPICollectiveWatchDog.Watch();
            subGrid = subGrid ?? new SubGrid(CellMask.GetFullMask(gridData));

            double maxTimeStep = double.MaxValue;
            double globalMaxTimeStep;
            System.Exception e = null;
            try {
                foreach (Chunk chunk in subGrid.VolumeMask) {
                    maxTimeStep = Math.Min(
                        maxTimeStep,
                        GetLocalStepSize(chunk.i0, chunk.Len));
                }
            } catch (System.Exception ee) {
                e = ee;
            }
            e.ExceptionBcast();

            unsafe
            {
                csMPI.Raw.Allreduce(
                    (IntPtr)(&maxTimeStep),
                    (IntPtr)(&globalMaxTimeStep),
                    1,
                    csMPI.Raw._DATATYPE.DOUBLE,
                    csMPI.Raw._OP.MIN,
                    csMPI.Raw._COMM.WORLD);
            }

            return dtFraction * globalMaxTimeStep;
        }

        /// <summary>
        /// Implement this method by determining the maximum admissible
        /// step-size in the local chunk of cells from <paramref name="i0"/> to
        /// <paramref name="i0"/> + <paramref name="Length"/>.
        /// </summary>
        /// <param name="i0">First cell index</param>
        /// <param name="Length">Number of cells</param>
        /// <returns>
        /// The maximum admissible step-size for explicit time integrators in
        /// the local chunk of cells from <paramref name="i0"/> to
        /// <paramref name="i0"/> + <paramref name="Length"/>.
        /// </returns>
        public abstract double GetLocalStepSize(int i0, int Length);
    }
}
