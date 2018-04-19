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
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.XDG {
    public partial class XDGField {

        delegate void EvaluateInternalSignature(int j0, int L, NodeSet NS, Basis basis, MultidimensionalArray Coördinates, MultidimensionalArray ResultAcc, double ResultPreScale);
        delegate void Picker(MultidimensionalArray R, int offset, int m, int SpcInd, MultidimensionalArray[] SR);

        private void GenericEval(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale,
            EvaluateInternalSignature EvalFunc,
            ref MultidimensionalArray[] Buffer, Func<int, int[]> BufferDim,
            Picker p) {
            LevelSetTracker trk = m_CCBasis.Tracker;

            if (Len > result.GetLength(0) + ResultCellindexOffset)
                throw new ArgumentOutOfRangeException("mismatch between Len and 0-th length of result");

            int j = 0;
            while (j < Len) {
                int L = Math.Min(trk.Regions.m_LenToNextChange[j0 + j], Len - j);

                ReducedRegionCode ReducedRegionCode;
                int NoOfSpecies;
                NoOfSpecies = trk.Regions.GetNoOfSpecies(j0 + j, out ReducedRegionCode);
                if (NoOfSpecies == 1) {
                    // next "L" cells are single-phase -> vectorized evaluation
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    EvaluateStd(j0 + j, L, NodeSet, result, ResultCellindexOffset + j, ResultPreScale, EvalFunc);
                } else {
                    // multiphase cell
                    // +++++++++++++++
                    L = 1; // 'EvaluateMultiphase' is not vectorized
                    EvaluateMultiphase(j0 + j, NodeSet, ReducedRegionCode, NoOfSpecies, result, ResultCellindexOffset + j, ResultPreScale, EvalFunc,
                        ref Buffer, BufferDim, p);
                }
                j += L;
            }
        }


        /// <summary>
        /// used by <see cref="EvaluateMultiphase"/>
        /// </summary>
        MultidimensionalArray[] _levSetVals = new MultidimensionalArray[4];

        /// <summary>
        /// used by <see cref="EvaluateMultiphase"/>
        /// </summary>
        double[] _levSetSign = new double[4];

        /// <summary>
        /// evaluation of field in cut- or near - cells;
        /// </summary>
        private void EvaluateMultiphase(int jCell, NodeSet NS, ReducedRegionCode ReducedRegionCode, int NoOfSpecies, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale,
            EvaluateInternalSignature EvalFunc,
            ref MultidimensionalArray[] Buffer, Func<int, int[]> BufferDim,
            Picker p) {
            LevelSetTracker trk = this.Basis.Tracker;
            int NoOfLevSets = trk.LevelSets.Count;
            int M = NS.NoOfNodes;

            {
                //var resultAcc = result.ExtractSubArrayShallow(new int[] { ResultCellindexOffset + j, 0 }, new int[] { ResultCellindexOffset + j, M - 1 });

                ushort RegionCode = trk.Regions.m_LevSetRegions[jCell];

                bool OnlyOneSpecies = true;
                for (int i = 0; i < NoOfLevSets; i++) {
                    int dst = LevelSetTracker.DecodeLevelSetDist(RegionCode, i);
                    if (dst == 0) {
                        // cut-cell with respect or levelset #i
                        // +++++++++++++++++++++++++++++++++++++

                        OnlyOneSpecies = false;
                        _levSetVals[i] = trk.DataHistories[i].Current.GetLevSetValues(NS, jCell, 1);
                    } else {
                        // near-cell with respect or levelset #i
                        // +++++++++++++++++++++++++++++++++++++

                        // find the species index...
                        _levSetVals[i] = null;
                        _levSetSign[i] = dst;

                    }
                    //_levSetSign2[i] = dst;
                }


                if (OnlyOneSpecies) {
                    // all level sets are 'NEAR'
                    // +++++++++++++++++++++++++

                    LevelSetSignCode levset_bytecode = LevelSetSignCode.ComputeLevelSetBytecode(_levSetSign);
                    int SpecInd = trk.GetSpeciesIndex(ReducedRegionCode, levset_bytecode);

                    // evaluate species ...
                    Evaluate_ithSpecies(jCell, NS, result, ResultCellindexOffset, ResultPreScale, SpecInd, EvalFunc);


                } else {
                    // at least one cell is 'CUT'
                    // ++++++++++++++++++++++++++

                    // allocate buffers/evaluate all species ...
                    if (Buffer.Length < NoOfSpecies) {
                        int oldL = Buffer.Length;
                        Array.Resize(ref Buffer, NoOfSpecies);
                        for (int i = oldL; i < NoOfSpecies; i++) {
                            Buffer[i] = MultidimensionalArray.Create(BufferDim(M));
                        }
                    }
                    for (int i = 0; i < NoOfSpecies; i++) {
                        if (Buffer[i].Dimension > 1 && Buffer[i].GetLength(1) != M)
                            Buffer[i].Allocate(BufferDim(M));

                        this.Evaluate_ithSpecies(jCell, NS, Buffer[i], 0, 0.0, i, EvalFunc);
                    }

                    // ... and pick the right species for the Node
                    for (int m = 0; m < M; m++) { // loop over nodes


                        for (int i = 0; i < NoOfLevSets; i++) {
                            if (_levSetVals[i] != null)
                                _levSetSign[i] = _levSetVals[i][0, m];
                        }

                        LevelSetSignCode levset_bytecode = LevelSetSignCode.ComputeLevelSetBytecode(_levSetSign);
                        int SpecInd = trk.GetSpeciesIndex(ReducedRegionCode, levset_bytecode);

                        p(result, ResultCellindexOffset, m, SpecInd, Buffer);
                    }
                }
            }
        }



        MultidimensionalArray m_CoordinateBuffer;

        private void Evaluate_ithSpecies(int j0, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale, int iSpecies, EvaluateInternalSignature EvalFunc) {
            int M = this.Basis.NonX_Basis.Length;
            Debug.Assert(M ==  this.Basis.DOFperSpeciesPerCell);
            int K = NodeSet.NoOfNodes;        // number of nodes
            Debug.Assert(result.GetLength(1) == K, "rank 1 is assumed to correlate with node set");

            if (m_CoordinateBuffer == null)
                m_CoordinateBuffer = MultidimensionalArray.Create(1, M);

            int N0 = iSpecies*M;
            for (int n = 0; n < M; n++)
                m_CoordinateBuffer[0, n] = this.m_Coordinates[j0, n + N0];

            int LL = result.Dimension;
            int[] I0 = new int[LL], IE = new int[LL];
            I0[0] = ResultCellindexOffset;
            I0[1] = 0;
            IE[0] = ResultCellindexOffset;
            IE[1] = K - 1;
            for (int ll = 2; ll < LL; ll++) {
                I0[ll] = 0;
                IE[ll] = result.GetLength(ll) - 1;
            }

            var resAcc = result.ExtractSubArrayShallow(I0, IE);

            EvalFunc(j0, 1, NodeSet,
                this.Basis.NonX_Basis,
                m_CoordinateBuffer, resAcc, ResultPreScale);
        }

        /// <summary>
        /// evaluation in standard cells (single phase cells)
        /// </summary>
        private void EvaluateStd(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale, EvaluateInternalSignature EvalFunc) {

            int M = this.Basis.NonX_Basis.Length;
            int K = NodeSet.NoOfNodes;        // number of nodes
            Debug.Assert(result.GetLength(1) == K, "rank 1 is assumed to correlate with node set");

            double[] CoördBase = m_Coordinates.m_BaseStorage;
            var Coörd = MultidimensionalArray.CreateWrapper(CoördBase, this.GridDat.iLogicalCells.NoOfCells, M);

            var C = Coörd.ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, M - 1 });

            int LL = result.Dimension;
            int[] I0 = new int[LL], IE = new int[LL];
            I0[0] = ResultCellindexOffset;
            I0[1] = 0;
            IE[0] = ResultCellindexOffset + Len - 1;
            IE[1] = K - 1;
            for (int ll = 2; ll < LL; ll++) {
                I0[ll] = 0;
                IE[ll] = result.GetLength(ll) - 1;
            }

            var resAcc = result.ExtractSubArrayShallow(I0, IE);

            EvalFunc(j0, Len, NodeSet,
                this.Basis.NonX_Basis,
                C, resAcc, ResultPreScale);
        }
    }
}
