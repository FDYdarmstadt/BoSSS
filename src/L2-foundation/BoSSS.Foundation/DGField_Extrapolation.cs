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
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation {

    partial class DGField {

        /// <summary>
        /// Polynomial extrapolation between cells.
        /// </summary>
        /// <param name="CellPairs">
        /// A list of which cells should be extrapolated to which cells. <br/>
        /// 1st index: enumeration;<br/>
        /// 2nd index, in the range of 0 to 1: 0: cell to extrapolate from (source), i.e. values in this cell remain unchanged.
        /// 1: cell to extrapolate to (target), the values in this cell will be the polynomial continuation of the "cell to extrapolate from".
        /// </param>
        /// <param name="scl">
        /// optional scaling for the target cells
        /// </param>
        /// <param name="beta">
        /// optional input scaling scaling for the target cells
        /// </param>        
        virtual public void CellExtrapolation<T, U>(int[,] CellPairs, T scl = default(T), U beta = default(U))
            where T : IList<double>
            where U : IList<double> //
        {
            //MPICollectiveWatchDog.Watch();

            if (CellPairs.GetLength(1) != 2)
                throw new ArgumentOutOfRangeException("second dimension is expected to be 2!");

            //if (!SurpressMPIUpdate)
            //    this.MPIExchange();

            int Esub = CellPairs.GetLength(0);
            Basis B = this.Basis;
            int N = B.Length;
            int JE = this.GridDat.iLogicalCells.NoOfCells;
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if (scl != null) {
                if (scl.Count != Esub)
                    throw new ArgumentException();
            }

            var M = MultidimensionalArray.Create(Esub, N, N);
            B.GetExtrapolationMatrices(CellPairs, M);
            double[] V = null, W = null;
            for (int esub = 0; esub < Esub; esub++) {
                var M_tmp = M.ExtractSubArrayShallow(esub, -1, -1);

                int jCell0 = CellPairs[esub, 0];
                int jCell1 = CellPairs[esub, 1];
                if (jCell0 < 0 || jCell0 >= JE)
                    throw new ArgumentOutOfRangeException();
                if (jCell1 < 0 || jCell1 >= J)
                    throw new ArgumentOutOfRangeException();

                if (!this.GridDat.iGeomCells.IsCellAffineLinear(jCell0))
                    throw new NotSupportedException("Currently not supported for curved cells.");
                if (!this.GridDat.iGeomCells.IsCellAffineLinear(jCell1))
                    throw new NotSupportedException("Currently not supported for curved cells.");

                double _scl = 1.0;
                if (scl != null)
                    _scl = scl[esub];

                double _beta = 0.0;
                if (beta != null)
                    _beta = beta[esub];

                if(esub == 0) {
                    V = new double[this.Coordinates.NoOfCols];
                    W = new double[this.Coordinates.NoOfCols];
                }

                {
                    this.Coordinates.GetRow(jCell0, V);
                    this.Coordinates.GetRow(jCell1, W);
                    
                    for (int n = 0; n < N; n++) {

                        double acc0 = 0;
                        for (int m = 0; m < N; m++) {
                            acc0 += V[m] * M_tmp[n, m];
                        }

                        W[n] = W[n] * _beta + acc0 * _scl;
                    }

                    this.Coordinates.SetRow(jCell1, W);
                }
            }
        }


        /// <summary>
        /// Polynomial extrapolation between cells.
        /// </summary>
        /// <param name="ExtrapolateTo">contains all cells for which extrapolated values should be computed</param>
        /// <param name="ExtrapolateFrom">contains all cells with "known values"</param>
        /// <returns>
        /// The number of cells (locally) for which the algorithm was not able to extrapolate a value
        /// </returns>
        virtual public int CellExtrapolation(CellMask ExtrapolateTo, CellMask ExtrapolateFrom) {
            MPICollectiveWatchDog.Watch();
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int NoOfNeigh;
            int[] NeighIdx;
            int[][] CN = this.GridDat.iLogicalCells.CellNeighbours;


            // mark all cells in which species 'Id' is known
            // ---------------------------------------------
            BitArray ValueIsKnown = ExtrapolateFrom.GetBitMaskWithExternal().CloneAs();
            CellMask _ExtrapolateTo = ExtrapolateTo.Except(ExtrapolateFrom);
            BitArray NeedToKnowValue = _ExtrapolateTo.GetBitMask();
            this.Clear(_ExtrapolateTo);

            // repeat until (for species 'Id' the DOF' of) all cells
            // that contain (at least a fraction of) species 'Id' are known...
            // ------------------------------------------------------------------
            int NoOfCells_ToExtrapolate = _ExtrapolateTo.NoOfItemsLocally;
            int NoOfCells_ExtrapolatedInSweep = 1;
            int sweepcnt = 0;

            List<double> scaling = new List<double>();
            List<int[]> CellPairs = new List<int[]>();
            BitArray cells_mod_in_sweep = new BitArray(J);
            while (true) {


                // MPI-parallel evaluation of termination criterion
                // ------------------------------------------------

                int bool_LocalRun = ((NoOfCells_ToExtrapolate > 0) && (NoOfCells_ExtrapolatedInSweep > 0)) ? 1 : 0;
                int bool_GlobalRun = 0;
                unsafe
                {
                    csMPI.Raw.Allreduce((IntPtr)(&bool_LocalRun), (IntPtr)(&bool_GlobalRun), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                }

                if (bool_GlobalRun <= 0)
                    // finished on all MPI processors
                    break;

                // MPI-update before local sweep: necessary for consistent result of alg.
                this.MPIExchange();
                if (sweepcnt > 0)
                    ValueIsKnown.MPIExchange(this.GridDat);



                // local work
                // ----------
                if (bool_LocalRun > 0) {


                    sweepcnt++;
                    NoOfCells_ExtrapolatedInSweep = 0;

                    scaling.Clear();
                    CellPairs.Clear();

                    cells_mod_in_sweep.SetAll(false);

                    for (int j = 0; j < J; j++) {

                        // determine whether there is something to do for cell 'j' or not ...
                        bool needToKnowSpecies = NeedToKnowValue[j];
                        bool _mustbeExtrapolated = needToKnowSpecies && !ValueIsKnown[j];

                        if (_mustbeExtrapolated) {
                            // continuation for this cell is needed
                            // ++++++++++++++++++++++++++++++++++++


                            // try to find a neighbour
                            NeighIdx = CN[j].CloneAs();
                            NoOfNeigh = NeighIdx.Length;
                            int FoundNeighs = 0;
                            for (int nn = 0; nn < NoOfNeigh; nn++) {
                                if (ValueIsKnown[NeighIdx[nn]]) {
                                    // bingo
                                    FoundNeighs++;
                                } else {
                                    NeighIdx[nn] = -1; // not usable
                                }
                            }

                            if (FoundNeighs <= 0)
                                // hope for better luck in next sweep
                                continue;

                            //Array.Clear(u2, 0, N);

                            for (int nn = 0; nn < NoOfNeigh; nn++) {
                                if (NeighIdx[nn] < 0)
                                    continue;

                                int _2 = j;            // cell to extrapolate TO
                                int _1 = NeighIdx[nn]; // cell to extrapolate FROM

                                CellPairs.Add(new int[] { _1, _2 });
                                double ooFoundNeighs = 1.0 / (double)FoundNeighs;
                                scaling.Add(ooFoundNeighs);
                            }

                            cells_mod_in_sweep[j] = true;
                            NoOfCells_ExtrapolatedInSweep++;
                        }
                    }

                    int E = CellPairs.Count;
                    int[,] _CellPairs = new int[E, 2];
                    for (int e = 0; e < E; e++)
                        _CellPairs.SetRow(e, CellPairs[e]);

                    double[] preScale = new double[scaling.Count];
                    preScale.SetAll(1.0);

                    this.CellExtrapolation(_CellPairs, scaling, preScale);

                    for (int j = 0; j < J; j++) {
                        if (cells_mod_in_sweep[j] == true) {
                            ValueIsKnown[j] = true;
                        }
                    }

                    NoOfCells_ToExtrapolate -= NoOfCells_ExtrapolatedInSweep;
                }
            }

            // return
            // ------

            return NoOfCells_ToExtrapolate;
        }

    }
}
