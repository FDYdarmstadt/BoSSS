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
using System.Text;
using BoSSS.Foundation;
using ilPSP.Tracing;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Utilities common to several incompressible solvers
    /// </summary>
    static public class SolverUtils {

        /// <summary>
        /// computes the index of the 
        /// 0-th DG coordinate (see <see cref="SetRefPtPressure_Matrix"/>, <see cref="SetRefPtPressure_Rhs"/>)
        /// from a position <paramref name="Point"/> vector.
        /// </summary>
        /// <param name="map">the mapping which defines how cell indices translate.</param>
        /// <param name="iVar">the index of the pressure variable in the mapping <paramref name="map"/>.</param>
        /// <param name="Point">position vector</param>
        static public int GetIndexOfPressureReferencePoint(double[] Point, UnsetteledCoordinateMapping map, int iVar) {
            using(new FuncTrace()) {
                var GridDat = map.GridDat;
                
                Basis PressureBasis = map.BasisS[iVar];
                int D = GridDat.SpatialDimension;

                long GlobalID, GlobalIndex;
                bool IsInside, onthisProc;
                GridDat.LocatePoint(Point, out GlobalID, out GlobalIndex, out IsInside, out onthisProc);
                
                int iRowGl = -111;
                if(onthisProc) {
                    int jCell = (int)GlobalIndex - GridDat.CellPartitioning.i0;
                    iRowGl = (int)map.GlobalUniqueCoordinateIndex_FromGlobal(iVar, GlobalIndex, 0);
                }
                iRowGl = iRowGl.MPIMax();

                return iRowGl;
            }
        }

        /// <summary>
        /// Set reference point for pressure by manipulating <paramref name="Matrix"/>, i.e.
        /// mean value for pressure correction (i.e. zeroth order polynomial coefficient)
        /// is set to zero in cell with GlobalID of zero, which is given by <paramref name="IndexRefPtPressure"/>.
        /// </summary>
        /// <param name="Matrix"></param>
        /// <param name="IndexRefPtPressure"></param>
        public static void SetRefPtPressure_Matrix(IMutableMatrixEx Matrix, int IndexRefPtPressure) {
            int i0 = Matrix.RowPartitioning.i0;
            int LocalLength = Matrix.RowPartitioning.LocalLength;

            if((IndexRefPtPressure >= i0) && (IndexRefPtPressure < i0 + LocalLength)) {
                Matrix[IndexRefPtPressure, IndexRefPtPressure] = 1.0;
                
                int[] colIdx = Matrix.GetOccupiedColumnIndices(IndexRefPtPressure);
                foreach (int col in colIdx) {
                    if (col != IndexRefPtPressure)
                        Matrix[IndexRefPtPressure, col] = 0.0;
                }

                for(int row = 0; row < LocalLength; row++)
                    if((row + i0) != IndexRefPtPressure)
                        Matrix[row + i0, IndexRefPtPressure] = 0.0;
            } else {
                for(int row = 0; row < LocalLength; row++)
                    Matrix[row + i0, IndexRefPtPressure] = 0.0;
            }
        }

        /// <summary>
        /// Set reference point for pressure by manipulating <paramref name="Rhs"/>.
        /// </summary>
        /// <param name="Rhs"></param>
        /// <param name="IndexRefPtPressure"></param>
        /// <param name="i0"></param>        
        public static void SetRefPtPressure_Rhs(double[] Rhs, int IndexRefPtPressure, int i0) {
            if((IndexRefPtPressure >= i0) && (IndexRefPtPressure < i0 + Rhs.Length)) {
                int LocalIndexRefPtPressure = IndexRefPtPressure - i0;
                Rhs[LocalIndexRefPtPressure] = 0.0;
            }
        }
    }
}
