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
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Tools for updating operator matrices under level-set movement
    /// </summary>
    public class TimeSteppingUtils {


        /// <summary>
        /// Tools for updating operator matrices under level-set movement, when ordering of DOFs (the mapping) changes.
        /// </summary>
        public static void OperatorLevelSetUpdate(LevelSetTracker LsTrk,
            BlockMsrMatrix[] OpMtx, UnsetteledCoordinateMapping RowMap, UnsetteledCoordinateMapping ColMap) {

            using (new FuncTrace()) {
                BoSSS.Foundation.Grid.IGridData GridData = LsTrk.GridDat;
                int JE = GridData.iLogicalCells.NoOfCells;
                int Jup = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int cell_j0 = GridData.CellPartitioning.i0;

                // build index mappings
                // ====================
                int[,] _old2NewRows, _old2NewCols;
                {
                    List<Tuple<int, int>> old2NewRows = new List<Tuple<int, int>>();
                    List<Tuple<int, int>> old2NewCols = new List<Tuple<int, int>>();
                    List<Tuple<int, int>> old2NewColsExt = new List<Tuple<int, int>>();

                    // loop over cells
                    for (int j = 0; j < JE; j++) {
                        int[] RowIdxUpdate;
                        if (j < Jup)
                            RowIdxUpdate = MappingUpdate(LsTrk, j, RowMap, false);
                        else
                            RowIdxUpdate = null;

                        int[] ColIdxUpdate = MappingUpdate(LsTrk, j, ColMap, false);

                        Debug.Assert((j >= Jup) || ((RowIdxUpdate != null) == (ColIdxUpdate != null)));
                        if (RowIdxUpdate != null) {
                            //Debug.Assert(RowMap.GetBlockLen(cell_j0 + j) == RowMap.MaxTotalNoOfCoordinatesPerCell);
                            //Debug.Assert(RowIdxUpdate.Length == RowMap.GetBlockLen(cell_j0 + j));

                            int i0 = RowMap.GetBlockI0(cell_j0 + j);

                            int II = RowIdxUpdate.Length;
                            for (int ii = 0; ii < II; ii++) {
                                int iOld = ii + i0;
                                int iNew = RowIdxUpdate[ii] + i0;
                                if (RowIdxUpdate[ii] >= 0) {
                                    Debug.Assert(old2NewRows.Count == 0 || old2NewRows[old2NewRows.Count - 1].Item1 < iOld);
                                    old2NewRows.Add(new Tuple<int, int>(iOld, iNew));
                                }
                            }
                        }

                        if (ColIdxUpdate != null) {
                            //Debug.Assert(ColMap.MinTotalNoOfCoordinatesPerCell == ColMap.MaxTotalNoOfCoordinatesPerCell);
                            //Debug.Assert(ColMap.GetBlockLen(cell_j0 + j) == ColMap.MaxTotalNoOfCoordinatesPerCell);
                            //Debug.Assert(ColIdxUpdate.Length == ColMap.GetBlockLen(cell_j0 + j));

                            // mapping is old-->new, i.e. new[IdxUpdate[k]] = old[k] for all k

                            //int i0 = ColMap.GetBlockI0(cell_j0 + j);
                            //int jLoc;
                            //if (j < Jup)
                            //    jLoc = j;
                            //else 
                            //    jLoc = LsTrk.GridDat.iParallel.Global2LocalIdx[cell_j0 + j];
                            int i0 = ColMap.GlobalUniqueCoordinateIndex(0, j, 0);
                            //Debug.Assert()
                            
                            int II = ColIdxUpdate.Length;
                            for (int ii = 0; ii < II; ii++) {
                                int iOld = ii + i0;
                                int iNew = ColIdxUpdate[ii] + i0;
                                if (ColIdxUpdate[ii] >= 0) {
                                    if (j < Jup) {
                                        Debug.Assert(old2NewCols.Count == 0 || old2NewCols[old2NewCols.Count - 1].Item1 < iOld);
                                        old2NewCols.Add(new Tuple<int, int>(iOld, iNew));
                                    } else {
                                        old2NewColsExt.Add(new Tuple<int, int>(iOld, iNew));
                                    }
                                }
                            }
                        }
                    }

                    // data conversion
                    int I = old2NewRows.Count;
                    _old2NewRows = new int[I, 2];
                    for (int i = 0; i < I; i++) {
                        var t = old2NewRows[i];
                        _old2NewRows[i, 0] = t.Item1;
                        _old2NewRows[i, 1] = t.Item2;
                    }
                    int J = old2NewCols.Count + old2NewColsExt.Count;
                    old2NewColsExt.Sort((a, b) => a.Item1 - b.Item1);
                    _old2NewCols = new int[J, 2];
                    int ja = 0, jb = 0;
                    for (int j = 0; j < J; j++) {
                        Tuple<int, int> t;
                        if (ja < old2NewCols.Count && jb < old2NewColsExt.Count) {
                            Tuple<int, int> t1, t2;
                            t1 = old2NewCols[ja];
                            t2 = old2NewColsExt[jb];
                            Debug.Assert(t1.Item1 != t2.Item1);
                            if (t1.Item1 < t2.Item1) {
                                t = t1;
                                ja++;
                            } else {
                                t = t2;
                                jb++;
                            }
                        } else if (ja < old2NewCols.Count) {
                            t = old2NewCols[ja];
                            ja++;
                        } else {
                            Debug.Assert(jb < old2NewColsExt.Count);
                            Debug.Assert(ja >= old2NewCols.Count);
                            t = old2NewColsExt[jb];
                            jb++;
                        }
                        

                        _old2NewCols[j, 0] = t.Item1;
                        _old2NewCols[j, 1] = t.Item2;
                        Debug.Assert(j == 0 || _old2NewCols[j - 1, 0] < _old2NewCols[j, 0]);
                    }
                    Debug.Assert(ja == old2NewCols.Count);
                    Debug.Assert(jb == old2NewColsExt.Count);
                }

                // update matrices
                // ===============
                for (int iMtx = 0; iMtx < OpMtx.Length; iMtx++) {
                    if(OpMtx[iMtx] != null)
                        OpMtx[iMtx] = OpMtx[iMtx].RecyclePermute(RowMap, ColMap, _old2NewRows, _old2NewCols);
                }
            }
        }

        /// <summary>
        /// Tools for updating RHS under level-set movement, when ordering of DOFs (the mapping) changes.
        /// </summary>
        public static void OperatorLevelSetUpdate(LevelSetTracker LsTrk, double[] Affine, UnsetteledCoordinateMapping RowMap) {
            using (new FuncTrace()) {
                if (RowMap.LocalLength != Affine.Length)
                    throw new ArgumentException();
                if (!object.ReferenceEquals(LsTrk.GridDat, RowMap.GridDat))
                    throw new ArgumentException();

                BoSSS.Foundation.Grid.IGridData GridData = LsTrk.GridDat;


                int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int GAMMA = RowMap.NoOfVariables;
                //int DELTA = ColMap.NoOfVariables;

                //int ColBlockSize = ColMap.MaxTotalNoOfCoordinatesPerCell;
                int RowBlockSize = RowMap.MaxTotalNoOfCoordinatesPerCell;

                int Gj0 = GridData.CellPartitioning.i0;

                int NoOfSpc = LsTrk.TotalNoOfSpecies;

                double[] OldAffine = new double[RowBlockSize];
                Dictionary<int, int[]> jNeighs = new Dictionary<int, int[]>(); // Key:     global cell index
                                                                               //                                                                values:  re-sorting for species in this cell
                Tuple<int, int[], double[]>[] oldRows = new Tuple<int, int[], double[]>[RowBlockSize];

                // loop over cells
                for (int j = 0; j < J; j++) {
                    int[] RowIdxUpdate = MappingUpdate(LsTrk, j, RowMap, false);

                    // affine vector update
                    // --------------------
                    if (Affine != null && RowIdxUpdate != null) {
                        // some update occurred 
                        // +++++++++++++++++++++++

                        int i0 = RowMap.LocalUniqueCoordinateIndex(0, j, 0);
                        Debug.Assert(i0 + RowBlockSize <= RowMap.LocalLength);
                        Debug.Assert(i0 >= 0);

                        for (int i = 0; i < RowBlockSize; i++) {
                            OldAffine[i] = Affine[i + i0];
                            Affine[i + i0] = 0.0;
                        }

                        for (int i = 0; i < RowBlockSize; i++) {
                            int k = RowIdxUpdate[i];

                            if (k < 0) {
                                if (OldAffine[i] != 0.0)
                                    throw new ApplicationException("Cannot kill non-null entry.");
                            } else {
                                Affine[i0 + k] = OldAffine[i];
                            }
                        }
                    }
                }
            }
        }

        static double RowNorm(Tuple<int, int[], double[]> row) {
            int L = row.Item1;
            int[] idxs = row.Item2;
            double[] vals = row.Item3;

            double acc = 0;
            for(int l = 0; l < L; l++) {
                Debug.Assert(idxs[l] >= 0);
                acc += Math.Abs(vals[l]);
            }
#if DEBUG
            for (int l = L; l < Math.Max(idxs.Length, vals.Length); l++) {
                Debug.Assert(idxs[l] < 0);
            }
#endif
            return acc;
        }
        
        /// <summary>
        /// Computes a mapping from the new, cell local DG coordinate index to the previous/old DG coordinate index, or vice-versa.
        /// </summary>
        /// <param name="jCell">
        /// Local cell index.
        /// </param>
        /// <param name="New2Old">
        /// If true, a mapping from new species index to previous/old species index is returned,
        /// if false the other way around.
        /// </param>
        /// <param name="LsTrk">
        /// </param>
        /// <param name="Map">
        /// DG coordinate mapping.
        /// </param>
        /// <returns>
        /// Null, if there is no change in species ordering in cell <paramref name="jCell"/>
        /// from the previous level-set tracker state to the actual. 
        /// 
        /// Otherwise, if <paramref name="New2Old"/> is true, an array of the same length 
        /// as the number of degrees-of-freedom currently in cell <paramref name="jCell"/>.
        ///  - index: current DOF index in cell <paramref name="jCell"/>.
        ///  - content: index of this DOF with respect to the previous level-set-tracker state.
        /// If <paramref name="New2Old"/> false, the other way around.
        /// </returns>
        public static int[] MappingUpdate(LevelSetTracker LsTrk, int jCell, UnsetteledCoordinateMapping Map, bool New2Old) {
            int[] SpeciesMap = SpeciesUpdate(LsTrk, jCell, New2Old);
            if (SpeciesMap == null)
                return null;

            int NoVar = Map.NoOfVariables;

            XDGBasis[] XBasiseS = new XDGBasis[NoVar];
            Basis[] BasiseS = new Basis[NoVar];
            int[] Ns = new int[NoVar];

            for (int iVar = 0; iVar < NoVar; iVar++) {
                Basis b = Map.BasisS[iVar];
                XBasiseS[iVar] = b as XDGBasis;
                if (XBasiseS[iVar] != null)
                    BasiseS[iVar] = XBasiseS[iVar].NonX_Basis;
                else
                    BasiseS[iVar] = b;
                Ns[iVar] = BasiseS[iVar].Length;
            }

            //List<int> oldIdx = new List<int>();
            //List<int> newIdx = new List<int>();
            int[] IdxMap = new int[Map.MaxTotalNoOfCoordinatesPerCell];
            ArrayTools.SetAll(IdxMap, int.MinValue);

            bool AnyRelocation = false;

            int i0 = Map.LocalUniqueCoordinateIndex(0, jCell, 0);

            for (int iVar = 0; iVar < NoVar; iVar++) {

                int i0Var = Map.LocalUniqueCoordinateIndex(iVar, jCell, 0);
                int N = Ns[iVar];
                int offset = i0Var - i0;

                if (XBasiseS[iVar] != null) {
                    // XDG-field

                    for (int i = 0; i < SpeciesMap.Length; i++) {
                        int ii = SpeciesMap[i];
                        for (int n = 0; n < N; n++) {
                            int idxS = offset + n + i * N;
                            int idxT = ii >= 0 ? offset + n + ii * N : int.MinValue;
                            Debug.Assert(IdxMap[idxS] < 0);
                            IdxMap[idxS] = idxT;

                            AnyRelocation |= (idxS != idxT);
                                
                        }
                    }
                } else {
                    // single-phase-field
                    for (int n = 0; n < N; n++) {
                        int idx = offset + n;
                        //oldIdx.Add(idx);
                        //newIdx.Add(idx);
                        Debug.Assert(IdxMap[idx] < 0);
                        IdxMap[idx] = idx;
                    }
                }
            }

            return IdxMap;
            //if (AnyRelocation)
            //    return IdxMap;
            //else
            //    return null;

        }


        /// <summary>
        /// Computes a mapping from the new species index to the old species index, or vice-versa.
        /// </summary>
        /// <param name="jCell">
        /// Local cell index.
        /// </param>
        /// <param name="New2Old">
        /// If true, a mapping from new species index to previous/old species index is returned,
        /// if false the other way around.
        /// </param>
        /// <param name="LsTrk">
        /// </param>
        /// <returns>
        /// Null, if there is no change in species ordering in cell <paramref name="jCell"/>
        /// from the previous level-set tracker state to the actual. 
        /// 
        /// Otherwise, if <paramref name="New2Old"/> is true, an array of the same length 
        /// as the number of species currently in cell <paramref name="jCell"/>.
        ///  - index: current species index in cell <paramref name="jCell"/>.
        ///  - content: index of this species with respect to the previous level-set-tracker state.
        /// If <paramref name="New2Old"/> false, the other way around.
        /// </returns>
        public static int[] SpeciesUpdate(LevelSetTracker LsTrk, int jCell, bool New2Old) {
            ushort[] OldCode = LsTrk.RegionsHistory[0].RegionsCode;
            ushort[] NewCode = LsTrk.RegionsHistory[1].RegionsCode;

            // Only if the return value of this function is positive, a mapping from current species index to previous species index
            // - index: current species index in cell <paramref name="jCell"/>.
            // - content: index of this species with respect to the previous level-set-tracker state.
            int[] UpdateMap = null;
            bool AnyUpdate = false;

            ushort oc = OldCode[jCell];
            ushort nc = NewCode[jCell];

            if (!New2Old) {
                // swap old and new 
                ushort b = oc;
                oc = nc;
                nc = b;
            }

            ReducedRegionCode oldRRC, newRRC;
            int OldNoOfSpc = LsTrk.GetNoOfSpeciesByRegionCode(oc, out oldRRC);
            int NewNoOfSpc = LsTrk.GetNoOfSpeciesByRegionCode(nc, out newRRC);


            for (int new_iSpc = 0; new_iSpc < NewNoOfSpc; new_iSpc++) {
                SpeciesId spId = LsTrk.GetSpeciesIdFromIndex(newRRC, new_iSpc);
                int old_iSpc = LsTrk.GetSpeciesIndex(oldRRC, spId); // if the species was not present in the previous state, this should be negative

                if (old_iSpc != new_iSpc) {
                    if (AnyUpdate == false) {
                        // init the update map: no change up to index 'new_iSpc'
                        AnyUpdate = true;
                        UpdateMap = new int[NewNoOfSpc];
                        for (int i = 0; i < new_iSpc; i++)
                            UpdateMap[i] = i;
                    }
                    // index 'new_iSpc' maps to old species index
                    UpdateMap[new_iSpc] = old_iSpc;
                }
            }

            if (NewNoOfSpc < OldNoOfSpc) {
                // the number of species is reduced - this also counts as change

                if (AnyUpdate == false) {
                    // init the update map: no change up to index 'new_iSpc'
                    AnyUpdate = true;
                    UpdateMap = new int[NewNoOfSpc];
                    for (int i = 0; i < NewNoOfSpc; i++)
                        UpdateMap[i] = i;
                }
            }

            // return
            return UpdateMap;
        }


        /// <summary>
        /// Diagnostic output.
        /// </summary>
        /// <param name="cout"></param>
        /// <param name="time"></param>
        public static void PrintCellSpecisTable(LevelSetTracker LsTrk, TextWriter cout, double time) {
            int J = LsTrk.GridDat.Cells.NoOfLocalUpdatedCells;

            cout.WriteLine("Species at time {0}: ==============", time);
            Console.WriteLine("j\tDist\t#Spc\t");

            for (int j = 0; j < J; j++) {
                // cell index
                cout.Write(j);
                cout.Write("\t");

                // level-set distance
                int dist = LevelSetTracker.DecodeLevelSetDist(LsTrk.Regions.RegionsCode[j], 0);
                cout.Write(dist);
                cout.Write("\t");

                // number of species in cell
                ReducedRegionCode rrc;
                int NoOfSpc = LsTrk.Regions.GetNoOfSpecies(j, out rrc);
                cout.Write(NoOfSpc);
                cout.Write("\t");

                // species sequence
                for (int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                    var SpId = LsTrk.GetSpeciesIdFromIndex(rrc, iSpc);
                    var SpNm = LsTrk.GetSpeciesName(SpId);

                    cout.Write(SpNm);
                    if (iSpc < (NoOfSpc - 1))
                        cout.Write(",");
                }
                cout.Write("\t");

                // new 2 old
                int[] N2O = TimeSteppingUtils.SpeciesUpdate(LsTrk, j, true);
                cout.Write("new2old: ");
                if (N2O != null) {
                    for (int i = 0; i < N2O.Length; i++) {
                        cout.Write(N2O[i]);
                        if (i < (N2O.Length - 1))
                            cout.Write(",");
                    }
                } else {
                    cout.Write("-");
                }
                cout.Write("\t");

                // old 2 new
                int[] O2N = TimeSteppingUtils.SpeciesUpdate(LsTrk, j, false);
                cout.Write("old2new: ");
                if (O2N != null) {
                    for (int i = 0; i < O2N.Length; i++) {
                        cout.Write(O2N[i]);
                        if (i < (O2N.Length - 1))
                            cout.Write(",");
                    }
                } else {
                    cout.Write("-");
                }
                cout.Write("\t");

                // end-of-line
                cout.WriteLine();
            }

            cout.WriteLine("---------------------");
        }

        /// <summary>
        /// Diagnostic output.
        /// </summary>
        static public void PrintMappingUpdate(LevelSetTracker LsTrk, UnsetteledCoordinateMapping Mapping) {
            int JE = LsTrk.GridDat.Cells.NoOfCells;

            for (int j = 0; j < JE; j++) {
                Console.Write(j);
                Console.Write(":\t");

                int[] mu = TimeSteppingUtils.MappingUpdate(LsTrk, j, Mapping, false);

                if (mu == null) {
                    Console.Write("-");
                } else {
                    for (int k = 0; k < mu.Length; k++) {
                        if (mu[k] >= 0)
                            Console.Write(mu[k]);
                        else
                            Console.Write("x");
                        if (k < mu.Length - 1)
                            Console.Write(",");
                    }
                }
                Console.WriteLine();
            }
        }



    }
}
