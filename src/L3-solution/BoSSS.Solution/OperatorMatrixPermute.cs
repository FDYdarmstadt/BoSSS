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
using BoSSS.Foundation.Grid;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Solution.Utils {
    
    
    /// <summary>
    /// Extension functions for operator matrices
    /// </summary>
    public static class OperatorMatrixPermute {


        /// <summary>
        /// A low-Performance method to permute some operator Matrix according to the GlobalID-permutation;
        /// </summary>
        /// <param name="M"></param>
        /// <param name="grd"></param>
        /// <param name="RowMap"></param>
        /// <param name="ColMap"></param>
        /// <returns></returns>
        /// <remarks>
        /// Low-Performance implementation, only for debugging purposes.
        /// </remarks>
        public static MsrMatrix ResortMatrix(this MsrMatrix M, GridData grd, UnsetteledCoordinateMapping RowMap, UnsetteledCoordinateMapping ColMap) {
            if (!M.RowPartitioning.Equals(RowMap)) {
                throw new ArgumentException();
            }
            if (!M.ColPartition.Equals(ColMap)) {
                throw new ArgumentException();
            }


            int[] RowPerm = ComputePermutation(RowMap, grd);
            int[] ColPerm = ComputePermutation(ColMap, grd);

            // permute columns
            MsrMatrix M2 = ColumnPermute(M, ColPerm);

            // permute rows
            MsrMatrix M3 = ColumnPermute(M2.Transpose(), RowPerm).Transpose();

            // return
            return M3;
        }

        /// <summary>
        /// A low-Performance method to permute the affine part of some operator according to the GlobalID-permutation;
        /// </summary>
        /// <param name="Affine"></param>
        /// <param name="grd"></param>
        /// <param name="RowMap"></param>
        /// <returns></returns>
        ///  <remarks>
        /// Low-Performance implementation, only for debugging purposes.
        /// </remarks>
        public static IList<double> ResortAffine(this IList<double> Affine, GridData grd, UnsetteledCoordinateMapping RowMap) {
            if (RowMap.LocalLength != Affine.Count) {
                throw new ArgumentException();
            }            

            // Copy values of Affine to diagonal matrix
            MsrMatrix Matrix = new MsrMatrix(RowMap, RowMap);
            int i0 = RowMap.i0;

            for (int row = 0; row < RowMap.LocalLength; row++)
                Matrix[i0 + row, i0 + row] = Affine[row];

            // Permute matrix
            MsrMatrix MatrixPermuted = Matrix.ResortMatrix(grd, RowMap, RowMap);

            // Copy values from permuted matrix to permuted affine
            IList<double> AffinePermuted = new double[Affine.Count];             

            for (int row = 0; row < RowMap.LocalLength; row++)
                AffinePermuted[row] = MatrixPermuted[i0 + row, i0 + row];

            return AffinePermuted;
        }

        private static MsrMatrix ColumnPermute(MsrMatrix M, int[] ColPerm) {
            MsrMatrix M2;
            M2 = new MsrMatrix(M.RowPartitioning, M.ColPartition);

            int[] ColIdx = null;
            double[] MtxVals = null;
            int LR;

            for (int i = M2.RowPartitioning.i0; i < M2.RowPartitioning.iE; i++) {
                //var row = M.GetRow(i);
                LR = M.GetRow(i, ref ColIdx, ref MtxVals);

                for(int lr = 0; lr < LR; lr++) {
                    int iCol = ColIdx[lr];
                    int icol_Targ = ColPerm[iCol];
                    M2[i, icol_Targ] = MtxVals[lr];
                }
            }
            return M2;
        }

        /// <summary>
        /// Computes the permutation for the given <paramref name="RowMap"/>
        /// </summary>
        /// <param name="RowMap"></param>
        /// <param name="grd"></param>
        /// <returns></returns>
        public static int[] ComputePermutation(UnsetteledCoordinateMapping RowMap, GridData grd) {
            int J = grd.Grid.NoOfUpdateCells;
            int[] N = RowMap.BasisS.Select(x => x.MaximalLength).ToArray();
            int NT = N.Sum();
            int Gamma = RowMap.BasisS.Count;
            var GlobalIDs = grd.CurrentGlobalIdPermutation.Values;

            int[] PermutationTable;
            PermutationTable = new int[RowMap.LocalLength];
            Debug.Assert(PermutationTable.Length == J*NT);
            for (int j = 0; j < J; j++) {

                int gid = (int)GlobalIDs[j];
                int iTarg0 = gid*NT;

                int i0 = (int)RowMap.LocalUniqueCoordinateIndex(0, j, 0);
                for (int f = 0; f < Gamma; f++) {
                    for (int n = 0; n < N[f]; n++) {
                        int i = (int)RowMap.LocalUniqueCoordinateIndex(f, j, n);
                        int iTarg = iTarg0 + (i - i0);
                        PermutationTable[i] = iTarg;
                    }
                }
            }


            int[] globalPermutationTable = new int[grd.Grid.NumberOfCells*NT];
            unsafe {
                int[] displ = RowMap.GetI0s().Select(x => (int)x).ToArray();
                int[] rcvCnt = new int[RowMap.MpiSize];
                for (int i = 0; i < RowMap.MpiSize; i++)
                    rcvCnt[i] = displ[i + 1] - displ[i];

                fixed (int* pSnd = PermutationTable, pRcv = globalPermutationTable, pDispl = displ, pRcvCnt = rcvCnt) {
                    csMPI.Raw.Allgatherv((IntPtr)pSnd, PermutationTable.Length, csMPI.Raw._DATATYPE.INT,
                        (IntPtr)pRcv, (IntPtr)pRcvCnt, (IntPtr)pDispl, csMPI.Raw._DATATYPE.INT,
                        csMPI.Raw._COMM.WORLD);
                }
            }

            return globalPermutationTable;
        }
    }
}
