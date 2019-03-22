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
using BoSSS.Foundation.Grid.Classic;
using System.Linq;
using BoSSS.Foundation.Comm;

namespace BoSSS.Foundation.XDG {

    partial class LevelSetTracker {

        /// <summary>
        /// A summary of information about a level set
        /// </summary>
        public class LevelSetRegions : ICloneable {

            LevelSetTracker m_owner;
         

            /// <summary>
            /// Constructor
            /// </summary>
            internal LevelSetRegions(LevelSetTracker owner) {
                m_owner = owner;
                m_ColorMap4Spc = new Dict_ColorMap4Spc(this);
            }

            /// <summary>
            /// The 'version' of the data, i.e. a number that is always increased
            /// when the data changes (see <see cref="LevelSetTracker.UpdateTracker"/>).
            /// </summary>
            public int Version {
                get;
                internal set;
            }

            /// <summary>
            /// Level set region codes for locally updated and external cells
            /// </summary>
            public ushort[] RegionsCode {
                get {
                    return m_LevSetRegions;
                }
            }

            /// <summary>
            /// A color map for each species; a color is a positive integer. Each topologically, simply connected part
            /// of a species is painted in a unique color. 
            /// - key: species ID
            /// - index into value: local cell index
            /// - each entry value: the unique color of the respective cell; 0 if the species is not present in the respective cell
            /// </summary>
            public IReadOnlyDictionary<SpeciesId,int[]> ColorMap4Spc {
                get {
                    return m_ColorMap4Spc;
                }
            }

            Dict_ColorMap4Spc m_ColorMap4Spc;

            class Dict_ColorMap4Spc : IReadOnlyDictionary<SpeciesId, int[]> {
                public Dict_ColorMap4Spc(LevelSetRegions __owner) {
                    m_owner = __owner;
                }

                LevelSetRegions m_owner;

                public int[] this[SpeciesId key] {
                    get {
                        if (!ContainsKey(key))
                            throw new KeyNotFoundException("Unknown Species");

                        if(!m_internal.TryGetValue(key, out int[] R)) {
                            R = m_owner.UpdateColoring(key);
                            m_internal.Add(key, R);
                        }

                        return R;
                    }

                }

                Dictionary<SpeciesId, int[]> m_internal = new Dictionary<SpeciesId, int[]>();

                internal Dict_ColorMap4Spc CloneNonShallow(LevelSetRegions __owner) {
                    var R = new Dict_ColorMap4Spc(__owner);
                    foreach(var kv in m_internal) {
                        R.m_internal.Add(kv.Key, kv.Value.CloneAs());
                    }
                    return R;
                }


                public IEnumerable<SpeciesId> Keys {
                    get {
                        return m_owner.m_owner.SpeciesIdS;
                    }
                }

                public IEnumerable<int[]> Values {
                    get {
                        var R = new List<int[]>();
                        foreach(var s in this.Keys) {
                            R.Add(this[s]);
                        }
                        return R;
                    }
                }

                public int Count {
                    get {
                        return this.Keys.Count();
                    }
                }

                public bool ContainsKey(SpeciesId key) {
                    return m_owner.SpeciesIdS.Contains(key);
                }

                IEnumerator<KeyValuePair<SpeciesId, int[]>> _GetEnumerator() {
                    var R = new List<KeyValuePair<SpeciesId, int[]>>();//[NoSpc];
                    foreach(var key in this.Keys) {
                        R.Add( new KeyValuePair<SpeciesId, int[]>(key, this[key]));
                    }
                    return R.GetEnumerator();
                }

                public IEnumerator<KeyValuePair<SpeciesId, int[]>> GetEnumerator() {
                    return this._GetEnumerator();
                }

                public bool TryGetValue(SpeciesId key, out int[] value) {
                    if (ContainsKey(key)) {
                        value = this[key];
                        return true;
                    } else {
                        value = null;
                        return false;
                    }

                }

                IEnumerator IEnumerable.GetEnumerator() {
                    return this._GetEnumerator();
                }
            }


            LevelSetRegions GetPreviousRegion() {
                int L = 1 - m_owner.RegionsHistory.GetPopulatedLength();
                
                for(int i = 1; i >= L; i--) {
                    if(object.ReferenceEquals(this, m_owner.RegionsHistory[i])) {
                        if (i > L)
                            return m_owner.RegionsHistory[i - 1];
                        else
                            return null;
                    } 
                }
                throw new ApplicationException();
            }

            /// <summary>
            /// - verifies the consistency of a color map among MPI processors
            /// - in future versions, it may be DEBUG-only
            /// </summary>
            private static void VerifyColoring(IGridData gdat, int[] ColorMap) {
                int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                int Je = gdat.iLogicalCells.NoOfExternalCells + J;

                if (ColorMap.Length != Je)
                    throw new ArgumentException();
                
                //ColorMap.MPIExchange(gdat);

                yx<y

            }

            private int[] UpdateColoring(SpeciesId SpId) {
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                int Je = this.GridDat.iLogicalCells.NoOfExternalCells + J;

                // paint on local processor
                // ========================
                int[] ColorMap = new int[J];
                var UsedColors = new HashSet<int>();
                {
                    int[] oldColorMap = GetPreviousRegion()?.ColorMap4Spc[SpId];
                    bool Incremental = oldColorMap != null;
                    if(Incremental) {
                        VerifyColoring(this.GridDat, oldColorMap);
                    }

                    CellMask SpMask = this.GetSpeciesMask(SpId);
                    BitArray SpBitMask = SpMask.GetBitMask();

                    int ColorCounter = 1;
                    var Part = new List<int>(); // all cells which form one part.
                    for (int j = 0; j < J; j++) { // sweep over cells...
                        if (SpBitMask[j] && ColorMap[j] == 0) {
                            Part.Clear();
                            int CurrentColor = ColorCounter;
                            bool ColorNegogiable = true;
                            ColorCounter = RecursiveColoring(this.GridDat, SpBitMask, j, ref CurrentColor, ColorMap, oldColorMap, ref ColorNegogiable, Part, UsedColors);
                            UsedColors.Add(CurrentColor);
                            Debug.Assert(ColorCounter > CurrentColor);
                        }
                    }
                }

                // parallelization
                // ===============
                int LocColors = UsedColors.Max();
                var ColorPart = new Partitioning(LocColors);
                int ColorOffset = ColorPart.i0;

                ColorMap.MPIExchange(GridDat);
                
                int[,] Edge2Cell = GridDat.iLogicalEdges.CellIndices;
                int NoEdg = Edge2Cell.GetLength(0);
                for(int iEdge = 0; iEdge < NoEdg; iEdge++) {
                    Debug.Assert(Edge2Cell[iEdge, 0] < J, "The external/ghost cell is expected to be the OUT-cell.");
                    int Cell1 = Edge2Cell[iEdge, 1];
                    if (Cell1 >= J) {
                        // reached an MPI boundary
                        int Cell0 = Edge2Cell[iEdge, 0];

                        int Color0 = ColorMap[Cell0];
                        int Color1 = ColorMap[Cell1];

                        if(Color0 != 0 && Color1 != 0 && Color0 != Color1) {
                            // need to do something...
                            
                            if(Color0 < 0 && Color1 > 0) {
                                // external color is fixed, internal color can be re-painted


                            } else if(Color0 > 0 && Color1 < 0) {
                                // internal color is fixed, but external color can be re-painted -> do nothing

                            } else if(Color0 < 0 && Color1 < 0) {
                                // both colors can be re-painted -> pick the minimum

                            } else if(Color0 > 0 && Color1 > 0) {
                                // no color can be re-painted -> this is a collision, but repaint in minimum color


                            } else {
                                Debug.Assert(false, "should never reach this point.");
                            }
                            

                        }
                    }
                }


                return ColorMap;
            }

            private static int RecursiveColoring(IGridData g, BitArray Msk, int j, ref int Color, int[] ColorMap, int[] oldColorMap, ref bool ColorNegotiable, List<int> Part, HashSet<int> UsedColors) {
                Debug.Assert(Msk[j] == true, "illegal to call on non-occupied cells");
                int J = g.iLogicalCells.NoOfLocalUpdatedCells;
                bool incremental = oldColorMap != null;

                int NextColor = Color + 1;

                if (incremental) {
                    int oldColor_j = oldColorMap[j];

                    // we care about the old colors
                    if (oldColor_j != 0) {

                        if (oldColor_j != Color) {
                            if (ColorNegotiable) {
                                // colors in previous map should match old colors
                                // for the current part, we are still not fixed in terms of color - we can re-paint the cells painted so far

                                if (!UsedColors.Contains(oldColor_j)) {
                                    // repainting is allowed, the old color (from previous time-step) is not yet used in this time-step
                                    Color = oldColor_j;
                                    foreach (int jk in Part) {
                                        ColorMap[jk] = oldColor_j;
                                    }
                                } else {
                                    // this is a topology change/a split
                                    jsdklsdjakldjk
                                }

                                NextColor = Math.Max(NextColor, Color + 1);

                                ColorNegotiable = false;
                            } else {
                                // this is a topology change/a merge
                                NextColor = Math.Max(NextColor, oldColorMap[j] + 1);
                                nklansxnakjxnjk
                            }
                        }
                    }
                }

                ColorMap[j] = Color;
                Part.Add(j);


                int[] jNeigh = g.iLogicalCells.CellNeighbours[j];


                foreach (int jN in jNeigh) {
                    if (Msk[jN] == false)
                        // Neighbor cell does not contain species -> end of recursion
                        continue;

                    if (j > J)
                        // external cell -> attention
                        continue;

                    if (ColorMap[jN] > 0) {
                        // already colored -> end of recursion
                        if (ColorMap[jN] != Color)
                            throw new ApplicationException("error in Algorithm.");
                        continue;
                    }

                    int recNextColor = RecursiveColoring(g, Msk, jN, ref Color, ColorMap, oldColorMap, ref ColorNegotiable, Part, UsedColors);
                    NextColor = Math.Max(NextColor, recNextColor);
                }

                return NextColor;
            }


            /// <summary>
            /// For each species: The sub-grid of cells populated with this species
            /// </summary>
            Dictionary<SpeciesId, SubGrid> m_SpeciesSubGrids;

            /// <summary>
            /// Level set region code,
            /// for locally updated and external cells
            /// </summary>
            internal ushort[] m_LevSetRegions;

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
                int JA = m_owner.GridDat.iLogicalCells.Count;
                m_LenToNextChange = new int[JA];
                m_LenToNextChange[JA - 1] = 1;
                ushort regionCd = m_LevSetRegions[JA - 1];
                for(int j = JA - 2; j >= 0; j--) {
                    if(m_LevSetRegions[j] != regionCd) {
                        m_LenToNextChange[j] = 1;
                        regionCd = m_LevSetRegions[j];
                    } else {
                        m_LenToNextChange[j] = m_LenToNextChange[j + 1] + 1;
                    }

                }

            }


            internal void InvalidateCaches() {
                this.m_LevelSetWings = null;
                this.m_NearField = null;
                this.m_NearField4LevelSet = null;
                this.m_NearMask = null;
                this.m_NearMask4LevelSet = null;
                this.m_SpeciesMask = null;
                this.m_SpeciesSubGrids = null;
                this.m_ColorMap4Spc = new Dict_ColorMap4Spc(this);
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

                if(FieldWidth > m_owner.NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.NearRegionWidth + ".", "FieldWidth");

                if(this.m_owner.LevelSets.Count == 1)
                    // if there is only one Level Set, no need to separate between
                    // cells-cut-by-any-level-set and cells-cut-by-a-specific-level-set
                    return GetNearMask4LevSet(0, FieldWidth);

                if(m_NearMask == null) {
                    m_NearMask = new CellMask[m_owner.m_NearRegionWidth + 1];
                }

                if(m_NearMask[FieldWidth] == null) {
                    int J = m_owner.m_gDat.Cells.NoOfLocalUpdatedCells;
                    BitArray ba = new BitArray(J, false);
                    //int NoL = m_LevelSets.Length;

                    for(int j = 0; j < J; j++) {
                        ushort code = m_LevSetRegions[j];

                        for(int levSetIdx = 0; levSetIdx < m_owner.NoOfLevelSets; levSetIdx++) {
                            int dist = LevelSetTracker.DecodeLevelSetDist(code, levSetIdx);
                            if(Math.Abs(dist) <= FieldWidth) {
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

                if(FieldWidth > m_owner.m_NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.m_NearRegionWidth + ".", "FieldWidth");

                if(m_owner.NoOfLevelSets == 1)
                    // if there is only one Level Set, no need to separate between
                    // cells-cut-by-any-level-set and cells-cut-by-a-specific-level-set
                    return GetNearFieldSubgrid4LevSet(0, FieldWidth);

                if(m_NearField == null) {
                    m_NearField = new SubGrid[m_owner.m_NearRegionWidth + 1];
                }

                if(m_NearField[FieldWidth] == null) {
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

                if(FieldWidth > m_owner.m_NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.m_NearRegionWidth + ".", "FieldWidth");
                if(levSetIdx < 0 || levSetIdx >= this.m_owner.NoOfLevelSets)
                    throw new IndexOutOfRangeException();


                if(m_NearField4LevelSet == null || m_NearField4LevelSet.GetLength(1) != this.m_owner.m_NearRegionWidth) {
                    m_NearField4LevelSet = new SubGrid[this.m_owner.NoOfLevelSets, this.m_owner.m_NearRegionWidth + 1];
                }

                if(m_NearField4LevelSet[levSetIdx, FieldWidth] == null) {
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


                if(FieldWidth > m_owner.m_NearRegionWidth)
                    throw new ArgumentException("Near-" + FieldWidth + " cannot be acquired, because this tracker is set to detect at most Near-" + m_owner.m_NearRegionWidth + ".", "FieldWidth");
                if(levSetIdx < 0 || levSetIdx >= this.m_owner.NoOfLevelSets)
                    throw new IndexOutOfRangeException();


                if(m_NearMask4LevelSet == null || m_NearMask4LevelSet.GetLength(1) != (m_owner.m_NearRegionWidth + 1)) {
                    m_NearMask4LevelSet = new CellMask[m_owner.NoOfLevelSets, m_owner.m_NearRegionWidth + 1];
                }
                

                if(m_NearMask4LevelSet[levSetIdx, FieldWidth] == null) {
                    // create subgrid

                    int J = m_owner.m_gDat.Cells.NoOfLocalUpdatedCells;
                    BitArray ba = new BitArray(J, false);
                    //int NoL = m_LevelSets.Length;

                    for(int j = 0; j < J; j++) {
                        ushort code = m_LevSetRegions[j];

                        int dist = LevelSetTracker.DecodeLevelSetDist(code, levSetIdx);
                        if(Math.Abs(dist) <= FieldWidth) {
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

                if(m_SpeciesSubGrids == null)
                    m_SpeciesSubGrids = new Dictionary<SpeciesId, SubGrid>();

                if(!m_SpeciesSubGrids.ContainsKey(specId)) {
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
                if(sign == 0.0)
                    throw new ArgumentException("must be either positive or negative");
                if(LevelSetIndex < 0 || LevelSetIndex >= m_owner.m_LevelSetHistories.Count)
                    throw new IndexOutOfRangeException("invalid level set index");

                int _sign = Math.Sign(sign);
                int iwing = _sign * (LevelSetIndex + 1);

                if(m_LevelSetWings == null)
                    m_LevelSetWings = new SortedDictionary<int, SubGrid>();

                if(!m_LevelSetWings.ContainsKey(iwing)) {

                    int J = m_owner.GridDat.Cells.NoOfLocalUpdatedCells;
                    BitArray mask = new BitArray(J);
                    for(int j = 0; j < J; j++) {
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
                if(m_SpeciesMask == null)
                    m_SpeciesMask = new Dictionary<SpeciesId, CellMask>();
                if(!m_SpeciesMask.ContainsKey(speciesId)) {

                    int J = m_owner.m_gDat.Grid.NoOfUpdateCells;
                    BitArray mask = new BitArray(J);
                    LevelSetSignCode[] signCodes = m_owner.GetLevelSetSignCodes(speciesId);

                    if(signCodes.Length == 0) {
                        throw new ArgumentException("Unknown species " + speciesId, "speciesName");
                    }

                    for(int jCell = 0; jCell < J; jCell++) {
                        bool matchesOne = IsSpeciesPresentInCell(speciesId, jCell);

                        // mask[i] = matches would also do it but a set operation on
                        // a BitMask is heavier than this ugly additional if statement
                        if(matchesOne) {
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
                for(int k = 0; k < signCodes.Length; k++) {
                    bool matches = true;
                    for(int j = 0; j < m_owner.NoOfLevelSets; j++) {
                        int sign = Math.Sign(DecodeLevelSetDist(m_LevSetRegions[jCell], j));

                        // Cell is cut, both signs exist and thus the mask
                        // matches for sure
                        if(sign == 0) {
                            continue;
                        }

                        int signCodeEntry = (signCodes[k].val & 0x1 << j) >> j;

                        // Code translation: 2 * {0, 1} - 1 = {-1, 1}
                        if(sign != 2 * signCodeEntry - 1) {
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
            public SinglePhaseField ToDGField() {
                SinglePhaseField TrackerField = new SinglePhaseField(new Basis(m_owner.GridDat, 0));
                // decrement loop: 
                for(int width = 0; width <= m_owner.NearRegionWidth; width++) {
                    TrackerField.AccConstant(1, this.GetNearFieldMask(width));
                }
                return TrackerField;
            }

            /// <summary>
            /// returns the (possible) number of species in cell <paramref name="j"/>;
            /// </summary>
            /// <param name="j">
            /// local cell index;
            /// </param>
            /// <param name="ReducedRegionCode">
            /// on exit, the reduced region code for cell <paramref name="j"/>
            /// in an 3-adic representation; This number 
            /// is later on required as an input for <see cref="GetSpeciesIndex(SpeciesId,int)"/>;
            /// </param>
            /// <returns></returns>
            /// <remarks>
            /// Here, three states for each of the for level sets are considered:
            /// <list type="bullet">
            ///   <item>positive far (FAR+)</item>
            ///   <item>negative far (FAR-)</item>
            ///   <item>positive near, negative near or cut</item>
            /// </list>
            /// This implies, that also for cells in the near region, memory is allocated for more than one
            /// species.
            /// </remarks>
            public int GetNoOfSpecies(int j, out ReducedRegionCode ReducedRegionCode) {
                ushort celJ = m_LevSetRegions[j];
                return this.m_owner.GetNoOfSpeciesByRegionCode(celJ, out ReducedRegionCode);
            }


            /// <summary>
            /// this function is the inverse to <see cref="GetSpeciesIndex(SpeciesId,int)"/>;
            /// </summary>
            /// <param name="SpeciesIndex"></param>
            /// <returns></returns>
            /// <remarks>
            /// this function is the inverse to <see cref="GetSpeciesIndex(SpeciesId, int)"/>
            /// </remarks>
            /// <param name="jCell">
            /// local cell index.
            /// </param>
            public SpeciesId GetSpeciesIdFromIndex(int jCell, int SpeciesIndex) {
                ReducedRegionCode rrc;
                int NoOfSpc = GetNoOfSpecies(jCell, out rrc);
                if(SpeciesIndex >= NoOfSpc) {
                    SpeciesId invalid;
                    invalid.cntnt = int.MinValue;
                    return invalid;
                } else {
                    return m_owner.GetSpeciesIdFromIndex(rrc, SpeciesIndex);
                }
            }


            /// <summary>
            /// Returns the index of species '<paramref name="_SpeciesId"/>' in the cell '<paramref name="jCell"/>'.
            /// </summary>
            /// <param name="jCell">
            /// a local cell index
            /// </param>
            /// <param name="_SpeciesId">
            /// species identification
            /// </param>
            public int GetSpeciesIndex(SpeciesId _SpeciesId, int jCell) {
                ReducedRegionCode rrc;
                int NoOfSpc = this.GetNoOfSpecies(jCell, out rrc);
                return m_owner.GetSpeciesIndex(rrc, _SpeciesId);
            }

            /// <summary>
            /// returns the (possible) number of species in cell <paramref name="j"/>;
            /// </summary>
            /// <param name="j">
            /// local cell index;
            /// </param>
            /// <returns></returns>
            public int GetNoOfSpecies(int j) {
                ReducedRegionCode dummy;
                int NoOf = GetNoOfSpecies(j, out dummy);
                return NoOf;
            }

            /// <summary>
            /// Equal to <see cref="LevelSetTracker.SpeciesIdS"/>.
            /// </summary>
            public IList<SpeciesId> SpeciesIdS {
                get {
                    return m_owner.SpeciesIdS;
                }
            }

            /// <summary>
            /// See <see cref="LevelSetTracker.SpeciesTable"/>.
            /// </summary>
            public Array SpeciesTable {
                get {
                    return m_owner.SpeciesTable;
                }
            }

            /// <summary>
            /// Equal to <see cref="LevelSetTracker.GetSpeciesName"/>.
            /// </summary>
            public String GetSpeciesName(SpeciesId id) {
                return m_owner.GetSpeciesName(id);
            }

            /// <summary>
            /// Equal to <see cref="LevelSetTracker.GetSpeciesName"/>.
            /// </summary>
            public SpeciesId GetSpeciesId(string SpeciesName) {
                return m_owner.GetSpeciesId(SpeciesName);
            }

            /// <summary>
            /// Non-Shallow copy
            /// </summary>
            /// <returns></returns>
            public object Clone() {
                var L = new LevelSetRegions(this.m_owner);

                L.m_LenToNextChange = this.m_LenToNextChange.CloneAs();
                L.m_LevSetRegions = this.m_LevSetRegions.CloneAs();
                L.Version = this.Version;
                L.m_ColorMap4Spc = this.m_ColorMap4Spc.CloneNonShallow(L);
                
                return L;
            }


            /// <summary>
            /// Link to the underlying background grid of the XDG discretization.
            /// </summary>
            public GridData GridDat {
                get {
                    return m_owner.GridDat;
                }
            }
        }
    }
}