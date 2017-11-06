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
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using System.Collections;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    partial class LevelSetTracker {

        /// <summary>
        /// A summary of information about a level set
        /// </summary>
        public class LevelSetRegionsInfo {

            LevelSetTracker m_owner;

            /// <summary>
            /// Constructor
            /// </summary>
            internal LevelSetRegionsInfo(LevelSetTracker owner) {
                m_owner = owner;
            }

            /// <summary>
            /// The 'version' of the data, i.e. a number that is always increased
            /// when the data changes (see <see cref="Update"/>).
            /// </summary>
            public int Version {
                get;
                internal set;
            }

            /// <summary>
            /// Level set region codes for locally updated and external cells
            /// </summary>
            public ushort[] LevelSetRegions {
                get {
                    return m_LevSetRegions;
                }
            }

            /// <summary>
            /// For each species: The sub-grid of cells populated with this species
            /// </summary>
            Dictionary<SpeciesId, SubGrid> m_SpeciesSubGrids;

            /// <summary>
            /// Level set region code,
            /// for locally updated and external cells
            /// </summary>
            public ushort[] m_LevSetRegions;

            /// <summary>
            /// For each cell <em>j</em>, the number of cells that follow with the same 
            /// region code as cell <em>j</em>.<br/>
            /// index: cell index <em>j</em>;
            /// </summary>
            public int[] m_LenToNextChange;

            /// <summary>
            /// Update of <see cref="m_LenToNextChange"/>.
            /// </summary>
            public void Recalc_LenToNextchange() {
                int JA = m_owner.GridDat.iLogicalCells.NoOfCells;
                m_LenToNextChange = new int[JA];
                m_LenToNextChange[JA - 1] = 1;
                ushort regionCd = m_LevSetRegions[JA - 1];
                for (int j = JA - 2; j >= 0; j--) {
                    if (m_LevSetRegions[j] != regionCd) {
                        m_LenToNextChange[j] = 1;
                        regionCd = m_LevSetRegions[j];
                    } else {
                        m_LenToNextChange[j] = m_LenToNextChange[j + 1] + 1;
                    }

                }

            }

            /// <summary>
            /// Returns the subgrid of those cells which are cut by level set No. <paramref name="LevSetIdx"/>.
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public SubGrid GetCutCellSubgrid4LevSet(int LevSetIdx) {
                return GetNearFieldSubgrid4LevSet(LevSetIdx, 0);
            }

            /// <summary>
            /// Returns the subgrid that contains all cells that are cut by any of the zero -- level sets
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public SubGrid GetCutCellSubGrid() {
                return GetNearFieldSubgrid(0);
            }

            /// <summary>
            /// returns the subgrid of those cells which are cut by level set No. <paramref name="LevSetIdx"/>.
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public CellMask GetCutCellMask4LevSet(int LevSetIdx) {
                return GetNearMask4LevSet(LevSetIdx, 0);
            }

            /// <summary>
            /// gets the subgrid that contains all cells that are cut by any of the zero -- level sets
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public CellMask GetCutCellMask() {
                return GetNearFieldMask(0);
            }

            /// <summary>
            /// Caches the values returned by <see cref="GetNearFieldSubgrid4LevSet"/>;
            /// index: width of subgrid around the cut cells
            /// </summary>
            CellMask[] m_NearMask;

            /// <summary>
            /// gets a mask of width <paramref name="FieldWidth"/>, around any level set
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public CellMask GetNearFieldMask(int FieldWidth) {

                if (FieldWidth > m_owner.NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.NearRegionWidth + ".", "FieldWidth");

                if (this.m_owner.LevelSets.Count == 1)
                    // if there is only one Level Set, no need to separate between
                    // cells-cut-by-any-level-set and cells-cut-by-a-specific-level-set
                    return GetNearMask4LevSet(0, FieldWidth);

                if (m_NearMask == null) {
                    m_NearMask = new CellMask[m_owner.m_NearRegionWidth + 1];
                }

                if (m_NearMask[FieldWidth] == null) {
                    int J = m_owner.m_gDat.Cells.NoOfLocalUpdatedCells;
                    BitArray ba = new BitArray(J, false);
                    //int NoL = m_LevelSets.Length;

                    for (int j = 0; j < J; j++) {
                        ushort code = m_LevSetRegions[j];

                        for (int levSetIdx = 0; levSetIdx < m_owner.m_LevelSets.Length; levSetIdx++) {
                            int dist = LevelSetTracker.DecodeLevelSetDist(code, levSetIdx);
                            if (Math.Abs(dist) <= FieldWidth) {
                                ba[j] = true;
                                continue;
                            }
                        }
                    }

                    m_NearMask[FieldWidth] = new CellMask(m_owner.m_gDat, ba);
                }

                return m_NearMask[FieldWidth];
            }

            /// <summary>
            /// caches the values returned by <see cref="GetNearFieldSubgrid4LevSet"/>;
            /// <br/>
            /// index: width of subgrid around the cut cells
            /// </summary>
            SubGrid[] m_NearField;

            /// <summary>
            /// gets a subgrid of width <paramref name="FieldWidth"/>, around any level set
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public SubGrid GetNearFieldSubgrid(int FieldWidth) {
                MPICollectiveWatchDog.Watch();

                if (FieldWidth > m_owner.m_NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.m_NearRegionWidth + ".", "FieldWidth");

                if (m_owner.m_LevelSets.Length == 1)
                    // if there is only one Level Set, no need to separate between
                    // cells-cut-by-any-level-set and cells-cut-by-a-specific-level-set
                    return GetNearFieldSubgrid4LevSet(0, FieldWidth);

                if (m_NearField == null) {
                    m_NearField = new SubGrid[m_owner.m_NearRegionWidth + 1];
                }

                if (m_NearField[FieldWidth] == null) {
                    m_NearField[FieldWidth] = new SubGrid(GetNearFieldMask(FieldWidth));
                }

                return m_NearField[FieldWidth];
            }

            /// <summary>
            /// caches the values returned by <see cref="GetNearFieldSubgrid4LevSet"/>;
            /// <br/>
            /// 1st index: Level Set - index <br/>
            /// 2nd index: width of subgrid around the cut cells
            /// </summary>
            SubGrid[,] m_NearField4LevelSet;

            /// <summary>
            /// gets a subgrid of width <paramref name="FieldWidth"/>, around level set No. <paramref name="levSetIdx"/>;
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public SubGrid GetNearFieldSubgrid4LevSet(int levSetIdx, int FieldWidth) {
                MPICollectiveWatchDog.Watch();

                if (FieldWidth > m_owner.m_NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.m_NearRegionWidth + ".", "FieldWidth");
                if (levSetIdx < 0 || levSetIdx >= this.m_owner.m_LevelSets.Length)
                    throw new IndexOutOfRangeException();


                if (m_NearField4LevelSet == null || m_NearField4LevelSet.GetLength(1) != this.m_owner.m_NearRegionWidth) {
                    m_NearField4LevelSet = new SubGrid[this.m_owner.m_LevelSets.Length, this.m_owner.m_NearRegionWidth + 1];
                }

                if (m_NearField4LevelSet[levSetIdx, FieldWidth] == null) {
                    // create subgrid
                    m_NearField4LevelSet[levSetIdx, FieldWidth] = new SubGrid(GetNearMask4LevSet(levSetIdx, FieldWidth));
                }

                return m_NearField4LevelSet[levSetIdx, FieldWidth];
            }

            /// <summary>
            /// caches the values returned by <see cref="GetNearFieldSubgrid4LevSet"/>;
            /// <br/>
            /// 1st index: Level Set - index <br/>
            /// 2nd index: width of subgrid around the cut cells
            /// </summary>
            CellMask[,] m_NearMask4LevelSet;

            /// <summary>
            /// gets a cell-mask of width <paramref name="FieldWidth"/>, around level set No. <paramref name="levSetIdx"/>;
            /// </summary>
            /// <remarks>
            /// In contrast to <see cref="SubGrid"/>'s, <see cref="CellMask"/>'s are not MPI-collective, and should therefore be preferred.
            /// </remarks>
            public CellMask GetNearMask4LevSet(int levSetIdx, int FieldWidth) {


                if (FieldWidth > m_owner.m_NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.m_NearRegionWidth + ".", "FieldWidth");
                if (levSetIdx < 0 || levSetIdx >= this.m_owner.m_LevelSets.Length)
                    throw new IndexOutOfRangeException();


                if (m_NearMask4LevelSet == null || m_NearMask4LevelSet.GetLength(1) != m_owner.m_NearRegionWidth) {
                    m_NearMask4LevelSet = new CellMask[m_owner.m_LevelSets.Length, m_owner.m_NearRegionWidth + 1];
                }

                if (m_NearMask4LevelSet[levSetIdx, FieldWidth] == null) {
                    // create subgrid

                    int J = m_owner.m_gDat.Cells.NoOfLocalUpdatedCells;
                    BitArray ba = new BitArray(J, false);
                    //int NoL = m_LevelSets.Length;

                    for (int j = 0; j < J; j++) {
                        ushort code = m_LevSetRegions[j];

                        int dist = LevelSetTracker.DecodeLevelSetDist(code, levSetIdx);
                        if (Math.Abs(dist) <= FieldWidth) {
                            ba[j] = true;
                            continue;
                        }
                    }

                    m_NearMask4LevelSet[levSetIdx, FieldWidth] = new CellMask(m_owner.GridDat, ba);
                }

                return m_NearMask4LevelSet[levSetIdx, FieldWidth];
            }

            /// <summary>
            /// Equivalent to <see cref="GetSpeciesSubGrid(string)"/>
            /// </summary>
            public SubGrid GetSpeciesSubGrid(SpeciesId specId) {
                MPICollectiveWatchDog.Watch();

                if (m_SpeciesSubGrids == null)
                    m_SpeciesSubGrids = new Dictionary<SpeciesId, SubGrid>();

                if (!m_SpeciesSubGrids.ContainsKey(specId)) {
                    CellMask cm = GetSpeciesMask(specId);
                    m_SpeciesSubGrids.Add(specId, new SubGrid(cm));
                } 
                return this.m_SpeciesSubGrids[specId];
            }

            SortedDictionary<int, SubGrid> m_LevelSetWings;

            /// <summary>
            /// returns the subgrid containing all cells in which the sign of the level set function #<paramref name="LevelSetIndex"/>
            /// is at least in one point equal to <paramref name="sign"/>.
            /// </summary>
            public SubGrid GetLevelSetWing(int LevelSetIndex, double sign) {
                MPICollectiveWatchDog.Watch();
                if (sign == 0.0)
                    throw new ArgumentException("must be either positive or negative");
                if (LevelSetIndex < 0 || LevelSetIndex >= m_owner.m_LevelSets.Length)
                    throw new IndexOutOfRangeException("invalid level set index");

                int _sign = Math.Sign(sign);
                int iwing = _sign * (LevelSetIndex + 1);

                if (m_LevelSetWings == null)
                    m_LevelSetWings = new SortedDictionary<int, SubGrid>();

                if (!m_LevelSetWings.ContainsKey(iwing)) {

                    int J = m_owner.GridDat.Cells.NoOfLocalUpdatedCells;
                    BitArray mask = new BitArray(J);
                    for (int j = 0; j < J; j++) {
                        int dist = DecodeLevelSetDist(m_LevSetRegions[j], LevelSetIndex);
                        mask[j] = (dist * _sign >= 0);
                    }

                    SubGrid LevelSetWing = new SubGrid(new CellMask(m_owner.GridDat, mask));
                    m_LevelSetWings.Add(iwing, LevelSetWing);
                }

                return m_LevelSetWings[iwing];
            }

            /// <summary>
            /// Creates an subgrid which only contains cells
            /// that are at least partly occupied by species
            /// <paramref name="speciesName"/>;
            /// </summary>
            /// <param name="speciesName">
            /// The name of the species
            /// </param>
            /// <returns>
            /// A subgrid containing only cells with at least a portion of
            /// species <paramref name="speciesName"/>
            /// </returns>
            public SubGrid GetSpeciesSubGrid(string speciesName) {
                SpeciesId id = m_owner.GetSpeciesId(speciesName);
                return GetSpeciesSubGrid(id);// this.SubGrids[id];
            }

            /// <summary>
            /// see <see cref="GetSpeciesMask(SpeciesId)"/>
            /// </summary>
            /// <param name="speciesName"></param>
            /// <returns></returns>
            public CellMask GetSpeciesMask(string speciesName) {
                SpeciesId id = m_owner.GetSpeciesId(speciesName);
                return GetSpeciesMask(id);
            }

            Dictionary<SpeciesId, CellMask> m_SpeciesMask;

            /// <summary>
            /// Creates an execution mask (volume mask) which only contains cells
            /// that are at least partly occupied by species
            /// <paramref name="speciesId"/>;
            /// </summary>
            /// <remarks>
            /// Note: The method is so complex because we do <b>not</b> want to
            /// consider the whole narrow-band because it may contain additional
            /// cells on the "other" side of the level set!
            /// </remarks>
            /// <param name="speciesId">
            /// The name of the species
            /// </param>
            /// <returns>
            /// A volume mask containing only cells with at least a portion of
            /// species <paramref name="speciesId"/>
            /// </returns>
            public CellMask GetSpeciesMask(SpeciesId speciesId) {
                if (m_SpeciesMask == null)
                    m_SpeciesMask = new Dictionary<SpeciesId, CellMask>();
                if (!m_SpeciesMask.ContainsKey(speciesId)) {

                    int J = m_owner.m_gDat.Grid.NoOfUpdateCells;
                    BitArray mask = new BitArray(J);
                    LevelSetSignCode[] signCodes = m_owner.GetLevelSetSignCodes(speciesId);

                    if (signCodes.Length == 0) {
                        throw new ArgumentException("Unknown species " + speciesId, "speciesName");
                    }

                    for (int jCell = 0; jCell < J; jCell++) {
                        bool matchesOne = IsSpeciesPresentInCell(speciesId, jCell);

                        // mask[i] = matches would also do it but a set operation on
                        // a BitMask is heavier than this ugly additional if statement
                        if (matchesOne) {
                            mask[jCell] = true;
                        }
                    }

                    m_SpeciesMask.Add(speciesId, new CellMask(m_owner.m_gDat, mask));
                }

                return m_SpeciesMask[speciesId];
            }

            /// <summary>
            /// Tests if a species  \f$\mathfrak{s}\f$ is actually \f$ \emph{present} \f$ in some cell \f$K_{\text{\tt jCell}}\f$,
            /// i.e. if 
            /// the measure of the species in the cell is positive, i.e. 
            /// \f[
            ///   \int_{K_{\text{\tt jCell}} \cap \mathfrak{s} } 1 \dV > 0 
            /// \f]
            /// </summary>
            /// <param name="speciesId">The id of species \mathfrak{s}.</param>
            /// <param name="jCell">a cell index</param>
            /// <returns></returns>
            /// <remarks>
            /// Note that 
            /// a species may not be 'present' in some cell 
            /// although 
            /// the index of a species (e.g. obtained by <see cref="GetSpeciesIndex(SpeciesId,int)"/>)
            /// may be positive.
            /// In the near-field however, some species index may be allocated (on the 'other' side of the level-set),
            /// but the species may not be present. 
            /// </remarks>
            public bool IsSpeciesPresentInCell(SpeciesId speciesId, int jCell) {
                bool matchesOne = false;

                LevelSetSignCode[] signCodes = m_owner.GetLevelSetSignCodes(speciesId);

                // Check if at least one of the sign codes matches the
                // situation in the current cell
                for (int k = 0; k < signCodes.Length; k++) {
                    bool matches = true;
                    for (int j = 0; j < m_owner.m_LevelSets.Length; j++) {
                        int sign = Math.Sign(DecodeLevelSetDist(m_LevSetRegions[jCell], j));

                        // Cell is cut, both signs exist and thus the mask
                        // matches for sure
                        if (sign == 0) {
                            continue;
                        }

                        int signCodeEntry = (signCodes[k].val & 0x1 << j) >> j;

                        // Code translation: 2 * {0, 1} - 1 = {-1, 1}
                        if (sign != 2 * signCodeEntry - 1) {
                            matches = false;
                            break;
                        }
                    }

                    matchesOne = matchesOne || matches;
                }
                return matchesOne;
            }

            /// <summary>
            /// Outputs a piecewise constant field which is 0 outside the narrow band
            /// 1+width on cut cells and decreasing on the near-cells
            /// </summary>
            public SinglePhaseField ToDGField(){
                SinglePhaseField TrackerField = new SinglePhaseField(new Basis(m_owner.GridDat, 0));
                // decrement loop: 
                for (int width = 0; width <= m_owner.NearRegionWidth ; width++) {
                    TrackerField.AccConstant(1, this.GetNearFieldMask(width));
                }                    
                return TrackerField;
            }
        }
    }
}