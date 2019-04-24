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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CNS.IBM {

    /// <summary>
    /// Static utility class for common IBM methods
    /// </summary>
    public static class IBMUtility {

        /// <summary>
        /// Specialized version of
        /// <see cref="BlockDiagonalMatrix.SpMV{VectorType1, VectorType2}(double, VectorType1, double, VectorType2)"/>
        /// where each block corresponds to an individual cell and where the
        /// multiplication shall only be performed over the cells in a given
        /// <paramref name="mask"/>
        /// </summary>
        /// <typeparam name="VectorType1"></typeparam>
        /// <typeparam name="VectorType2"></typeparam>
        /// <param name="M"></param>
        /// <param name="alpha"></param>
        /// <param name="a"></param>
        /// <param name="beta"></param>
        /// <param name="acc"></param>
        /// <param name="mask"></param>
        internal static void SubMatrixSpMV<VectorType1, VectorType2>(BlockMsrMatrix M, double alpha, VectorType1 a, double beta, VectorType2 acc, CellMask mask)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> {

            if (acc.Count != M.RowPartitioning.LocalLength)
                throw new ArgumentException("mismatch between accumulator and number of rows, i.e. row partition.", "acc");
            if (a.Count != M.ColPartition.LocalLength)
                throw new ArgumentException("mismatch between input vector 'a' and number of columns, i.e. column partition.", "a");

            //int M1 = M.NoOfColsPerBlock;
            //int N1 = M.NoOfRowsPerBlock;
            double[] aBuf = null;// new double[M1];
            MultidimensionalArray Mbuf = null;

            int Glob_cell0 = mask.GridData.CellPartitioning.i0;

            foreach (int cell in mask.ItemEnum) {
                //int i0 = cell * N1; // row
                //int j0 = cell * M1; // col
                int i0Glob = M._RowPartitioning.GetBlockI0(cell + Glob_cell0);
                int j0Glob = M._ColPartitioning.GetBlockI0(cell + Glob_cell0);
                int i0 = i0Glob - M._RowPartitioning.i0;
                int j0 = j0Glob - M._ColPartitioning.i0;
                int M1 = M._RowPartitioning.GetBlockLen(cell + Glob_cell0);
                int N1 = M._ColPartitioning.GetBlockLen(cell + Glob_cell0);

                if (aBuf == null || aBuf.Length < M1)
                    aBuf = new double[M1];

                for (int j = 0; j < M1; j++) {
                    aBuf[j] = a[j + j0];
                }

                if (Mbuf == null || Mbuf.GetLength(0) != M1 || Mbuf.GetLength(1) != N1)
                    Mbuf = MultidimensionalArray.Create(M1, N1);

                M.ReadBlock(i0Glob, j0Glob, Mbuf);

                for (int i = 0; i < N1; i++) {
                    double _acc = 0;

                    for (int j = 0; j < M1; j++) {
                        _acc += Mbuf[i,j] * aBuf[j];
                    }

                    acc[i + i0] = acc[i + i0] * beta + _acc * alpha;
                }
            }

        }

        internal static MultidimensionalArray GetBlock(BlockMsrMatrix Mtx, int jCell) {
            Debug.Assert(Mtx._RowPartitioning.FirstBlock + jCell == Mtx._ColPartitioning.FirstBlock + jCell);
            int jCellGlb = Mtx._RowPartitioning.FirstBlock + jCell;
            int i0Glb = Mtx._RowPartitioning.GetBlockI0(jCellGlb);
            int j0Glb = Mtx._ColPartitioning.GetBlockI0(jCellGlb);
            int Mblk = Mtx._RowPartitioning.GetBlockLen(jCellGlb);
            int Nblk = Mtx._ColPartitioning.GetBlockLen(jCellGlb);
            MultidimensionalArray ret = MultidimensionalArray.Create(Mblk, Nblk);
            Mtx.ReadBlock(i0Glb, j0Glb, ret);
            return ret;
        }

        internal static double GetMassMatrixDeterminant(ImmersedSpeciesMap speciesMap, CoordinateMapping mapping, CellMask cellMask) {
            BlockMsrMatrix massMatrix = speciesMap.GetMassMatrixFactory(mapping).MassMatrix;
            Basis maxBasis = mapping.BasisS.ElementAtMax(b => b.Degree);

            MultidimensionalArray subMatrix = MultidimensionalArray.Create(
                cellMask.NoOfItemsLocally + 1, maxBasis.Length, cellMask.NoOfItemsLocally + 1, maxBasis.Length);
            int cellIndex = 0;
            foreach (Chunk chunk in cellMask) {
                for (int i = 0; i < chunk.Len; i++) {
                    //IMatrix block = massMatrix.GetBlock(i + chunk.i0);
                    MultidimensionalArray block = GetBlock(massMatrix, i + chunk.i0);
                    for (int j = 0; j < maxBasis.Length; j++) {
                        for (int k = 0; k < maxBasis.Length; k++) {
                            subMatrix[cellIndex, j, cellIndex, k] = block[j, k];
                        }
                    }

                    cellIndex++;
                }
            }

            for (int j = 0; j < maxBasis.Length; j++) {
                subMatrix[cellIndex, j, cellIndex, j] = 1.0;
            }

            subMatrix = subMatrix.ResizeShallow(
                (cellMask.NoOfItemsLocally + 1) * maxBasis.Length,
                (cellMask.NoOfItemsLocally + 1) * maxBasis.Length);

            return subMatrix.cond();
        }

        internal static double GetMaxLocalMassMatrixDeterminant(ImmersedSpeciesMap speciesMap, CoordinateMapping mapping, CellMask cellMask, out int maxCondCell) {
            BlockMsrMatrix massMatrix = speciesMap.GetMassMatrixFactory(mapping).MassMatrix;
            Basis maxBasis = mapping.BasisS.ElementAtMax(b => b.Degree);

            MultidimensionalArray subMatrix = MultidimensionalArray.Create(
                maxBasis.Length, maxBasis.Length);

            double maxCond = 0.0;
            maxCondCell = -1;
            foreach (Chunk chunk in cellMask) {
                for (int i = 0; i < chunk.Len; i++) {
                    //IMatrix block = massMatrix.GetBlock(i + chunk.i0);
                    MultidimensionalArray block = GetBlock(massMatrix, i + chunk.i0);
                    for (int j = 0; j < maxBasis.Length; j++) {
                        for (int k = 0; k < maxBasis.Length; k++) {
                            subMatrix[j, k] = block[j, k];
                        }
                    }

                    double cond = subMatrix.cond();
                    if (cond > maxCond) {
                        maxCond = cond;
                        maxCondCell = i + chunk.i0;
                    }
                }
            }

            return maxCond;
        }
    }
}
