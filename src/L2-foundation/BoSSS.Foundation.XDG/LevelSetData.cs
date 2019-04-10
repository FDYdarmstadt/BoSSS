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
using ilPSP.Utils;
using MPI.Wrappers;

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
            public IReadOnlyDictionary<SpeciesId, int[]> ColorMap4Spc {
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
                        MPICollectiveWatchDog.Watch();
                        if (!ContainsKey(key))
                            throw new KeyNotFoundException("Unknown Species");

                        if (!m_internal.TryGetValue(key, out int[] R)) {
                            R = m_owner.UpdateColoring(key);
                            m_internal.Add(key, R);
                        }

                        return R;
                    }

                }

                Dictionary<SpeciesId, int[]> m_internal = new Dictionary<SpeciesId, int[]>();

                internal Dict_ColorMap4Spc CloneNonShallow(LevelSetRegions __owner) {
                    var R = new Dict_ColorMap4Spc(__owner);
                    foreach (var kv in m_internal) {
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
                        foreach (var s in this.Keys) {
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
                    foreach (var key in this.Keys) {
                        R.Add(new KeyValuePair<SpeciesId, int[]>(key, this[key]));
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

                for (int i = 1; i >= L; i--) {
                    if (object.ReferenceEquals(this, m_owner.RegionsHistory[i])) {
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
                return;

                //
                // Remark: this test is not correct int the following case:
                // - assume a part is shared
                // - on one MPI process, it falls apart into two parts
                // - these will appear a separate parts with non-unique color; this raises a false assertion.

                /*
                int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                int Je = gdat.iLogicalCells.NoOfExternalCells + J;

                if (ColorMap.Length != Je)
                    throw new ArgumentException();

                //  check if color map has been correctly MPI exchanged
                //  ====================================================
                {
                    var ColorMapClone = ColorMap.CloneAs();
                    ColorMapClone.MPIExchange(gdat);

                    for (int j = J; j < Je; j++) {
                        if (ColorMapClone[j] != ColorMap[j])
                            throw new ApplicationException("ColorMap has not been synchronized correctly.");
                    }
                }

                //  Local Checks:
                //  * connected domains are painted in the same color
                //  * locally used colors are unique
                //  ====================================================

                var ColorsOfIsolatedPars = new HashSet<int>(); // all colors of parts that are isolated on this processor, i.e. they *do not* extend into the external cells.
                var ColorsOfAllParts = new HashSet<int>();

                BitArray CheckedCells = new BitArray(Je);
                for (int j = 0; j < J; j++) { 
                    // we only check local cells, because for external we don't have all neighborship info,
                    // e.g. we don't know which external cells are neighbors to which other cells
                    int Color = ColorMap[j];
                    if (Color < 0)
                        throw new ApplicationException("Negative color in cell " + j + " (number of cells " + J + ", excluding external, MPI rank " + gdat.MpiRank + ") - color has not been fixed.");
                    if (CheckedCells[j] == false && Color != 0) {
                        bool NonIsolated = CheckColorRecursive(ColorMap, j, Color, CheckedCells, gdat);

                        if (!ColorsOfAllParts.Add(Color)) {
                            
                            //Debugger.Launch();
                            //using (var stw = new System.IO.StreamWriter("verdammte-zelle.csv")) {
                            //    for (int i = 0; i < Je; i++) {
                            //        var cen = gdat.GlobalNodes.GetValue_Cell(Grid.RefElements.Square.Instance.Center, i, 1);
                            //        double x = cen[0, 0, 0];
                            //        double y = cen[0, 0, 1];
                            //        stw.WriteLine("{0}\t{1}\t{2}", x, y, ColorMap[i]);
                            //    }
                            //}
                            
                            throw new ApplicationException("color " + Color + " is non-unique, cell #" + j + " (number of cells " + J + ", excluding external, MPI rank " + gdat.MpiRank + ")");
                        }

                        if (!NonIsolated) {
                            // isolated part
                            ColorsOfIsolatedPars.Add(Color);
                        }
                    }
                }
                Debug.Assert(ColorsOfIsolatedPars.IsSubsetOf(ColorsOfAllParts));

                for (int j = 0; j < J; j++) {
                    // further algorithm check: are all colored cells checked?
                    Debug.Assert((CheckedCells[j] == true) || (ColorMap[j] == 0));
                    Debug.Assert((CheckedCells[j] == true) == (ColorMap[j] != 0));
                }

                // check global uniqueness
                // =======================

                //
                // Assumption: (a) AND (b) equal (c), where:
                //  (a) colors of isolated parts are globally (over all MPI processors) unique
                //  (b) colors of all parts are locally (only on current MPI processor) unique (must be checked also for external cells)
                //  (c) all parts are globally unique
                //

                {
                    var SendData = new Dictionary<int, int[]>();
                    if (gdat.MpiRank != 0)
                        SendData.Add(0, ColorsOfIsolatedPars.ToArray());
                    var CollectedData = SerialisationMessenger.ExchangeData(SendData);

                    if (gdat.MpiRank == 0) {
                        var ColorsOfIsolatedPars_globally = new HashSet<int>();
                        ColorsOfIsolatedPars_globally.AddRange(ColorsOfIsolatedPars);

                        for (int iRnk = 1; iRnk < gdat.MpiSize; iRnk++) {
                            int[] UsedColors = CollectedData[iRnk];
                            foreach (int Color in UsedColors) {
                                if (!ColorsOfIsolatedPars_globally.Add(Color))
                                    throw new ApplicationException("color " + Color + " is globally non-unique; found a second time on processor " + iRnk);
                            }
                        }
                    }
                }
                */
            }

            /*
            private static bool CheckColorRecursive(int[] ColorMap, int j, int Color, BitArray CheckedCells, IGridData gdat) {
                Debug.Assert(ColorMap[j] != 0, "Recursion error.");  // an error in this algorithm -> debug assertion 
                Debug.Assert(CheckedCells[j] == false); // detto
                if (ColorMap[j] != Color) { // error in the data to check -> 
                    throw new ApplicationException("Color mismatch/ color change within one connected part.");
                }
                int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                CheckedCells[j] = true;

                if (j >= J) {
                    // external cell - no further recursion
                    return true;
                }

                bool R = false;

                int[] Neighs_j = gdat.iLogicalCells.CellNeighbours[j];
                foreach (int jN in Neighs_j) {
                    if (ColorMap[jN] == 0)
                        continue;
                    if (CheckedCells[jN] == true)
                        continue;
                    R |= CheckColorRecursive(ColorMap, jN, Color, CheckedCells, gdat);
                }

                return R;
            }
            */


            private int[] UpdateColoring(SpeciesId SpId) {

                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                int Je = this.GridDat.iLogicalCells.NoOfExternalCells + J;
                MPICollectiveWatchDog.Watch();

                // paint on local processor
                // ========================
                int[] ColorMap = new int[Je];
                int ColorCounter;
                {
                    CellMask SpMask = this.GetSpeciesMask(SpId); 
                    BitArray SpBitMask = SpMask.GetBitMaskWithExternal();

                    //int ColorCounter = this.GridDat.CellPartitioning.i0; // (more-than) worst case estimation of used colors on previous processors
                    ColorCounter = 1;
                    //var Part = new List<int>(); // all cells which form one part.
                    for (int j = 0; j < J; j++) { // sweep over local cells...
                        if (SpBitMask[j] && ColorMap[j] == 0) {
                            //Part.Clear();
                            int CurrentColor = ColorCounter;
                            bool IsIsolated = true;
                            RecursiveColoring(this.GridDat, SpBitMask, j, CurrentColor, ColorMap, ref IsIsolated);

                            ColorCounter++;
                            Debug.Assert(ColorCounter > CurrentColor);
                        }
                    }
                }


                // parallelization, pt 1, make new colors *globally* unique
                // ========================================================
                {
                    int LocColors = ColorCounter - 1;
                    var ColorPart = new Partitioning(LocColors);
                    int ColorOffset = ColorPart.i0;

                    for (int j = 0; j < J; j++) {
                        if (ColorMap[j] != 0) {
                            int Color_j = ColorMap[j];
                            Color_j += ColorOffset;
                            ColorMap[j] = Color_j;
                        }

                        Debug.Assert((ColorMap[j] == 0) || (ColorMap[j] > ColorPart.i0));
                        Debug.Assert((ColorMap[j] == 0) || (ColorMap[j] <= ColorPart.iE));
                    }
                }



                // parallelization, pt 2, synchronize shared parts
                // ===============================================



                ColorMap.MPIExchange(GridDat);

                {
                    /*

                    long[] extCells = GridDat.iParallel.GlobalIndicesExternalCells;
                    var CellPart = GridDat.CellPartitioning;

                    int NoOfConflicts;
                    do {
                        NoOfConflicts = 0;

                        int[,] Edge2Cell = GridDat.iLogicalEdges.CellIndices;
                        int NoEdg = Edge2Cell.GetLength(0);
                        for (int iEdge = 0; iEdge < NoEdg; iEdge++) {
                            Debug.Assert(Edge2Cell[iEdge, 0] < J, "The external/ghost cell is expected to be the OUT-cell.");
                            int Cell1 = Edge2Cell[iEdge, 1];
                            if (Cell1 >= J) {
                                // reached an MPI boundary
                                int Cell0 = Edge2Cell[iEdge, 0];

                                int Color0 = ColorMap[Cell0];
                                int Color1 = ColorMap[Cell1];

                                if (Color0 != 0 && Color1 != 0 && Color0 != Color1) {
                                    // need to do something...
                                    NoOfConflicts++;


                                    // Rule: the lower MPI rank determines the color
                                    long GlobIdx = extCells[Cell1 - J];
                                    Debug.Assert(GlobIdx < CellPart.i0 || GlobIdx >= CellPart.iE);
                                    bool ChangeColor = GlobIdx < CellPart.i0;

                                    if (ChangeColor) {
                                        // re-paint my stuff int 'Color1'
                                        RepaintRecursive(Color1, ColorMap, Cell0, this.GridDat);
                                    }

                                   
                                }
                            }
                        }

                        NoOfConflicts = NoOfConflicts.MPISum();

                        if (NoOfConflicts > 0)
                            ColorMap.MPIExchange(GridDat);

                    } while (NoOfConflicts > 0); // if a part is shared by more than 2 processors, multiple iterations might be necessary
                    */

                    // data structure to store color equality
                    // --------------------------------------

                    var locColEq = new Dictionary<int, HashSet<int>>(); // key: some locally used color value; value: all colors that should be equal (includes also the key)
                    void AddEqPairing(int Color0, int Color1) { // add some color pair to 'locColEq'

                        HashSet<int> equalCols;
                        if (!locColEq.TryGetValue(Color0, out equalCols)) {
                            if (!locColEq.TryGetValue(Color1, out equalCols)) {
                                equalCols = new HashSet<int>();
                                locColEq.Add(Color1, equalCols);
                            }
                            locColEq.Add(Color0, equalCols);
                        } else {
                            if(!locColEq.ContainsKey(Color1)) {
                                locColEq.Add(Color1, equalCols);
                            }
                        }
                        equalCols.Add(Color0);
                        equalCols.Add(Color1);

                        Debug.Assert(object.ReferenceEquals(locColEq[Color0], equalCols));
                        Debug.Assert(object.ReferenceEquals(locColEq[Color1], equalCols));

#if DEBUG
                        // check data integrity:
                        foreach (var kv in locColEq) { // for each color in 'locColEq'
                            int col = kv.Key;
                            var _equalCols = kv.Value;

                            Debug.Assert(_equalCols.Contains(col)); // equal colors must contain key value itself
                            foreach (var pairCol in _equalCols) { // all equal colors...
                                if (pairCol != col) {
                                    Debug.Assert(locColEq.ContainsKey(pairCol)); // must also be in the dictionary...
                                    Debug.Assert(object.ReferenceEquals(locColEq[pairCol], _equalCols)); // ...and the equality set must be tha same.
                                }
                            }
                        }
#endif
                    }


                    // collect all colors that should be locally equal
                    // -----------------------------------------------
                    int[,] Edge2Cell = GridDat.iLogicalEdges.CellIndices;
                    int NoEdg = Edge2Cell.GetLength(0);
                    {
                        //long[] extCells = GridDat.iParallel.GlobalIndicesExternalCells;

                        for (int iEdge = 0; iEdge < NoEdg; iEdge++) {
                            Debug.Assert(Edge2Cell[iEdge, 0] < J, "The external/ghost cell is expected to be the OUT-cell.");
                            int Cell1 = Edge2Cell[iEdge, 1];
                            if (Cell1 >= J) {
                                // reached an MPI boundary
                                int Cell0 = Edge2Cell[iEdge, 0];

                                int Color0 = ColorMap[Cell0];
                                int Color1 = ColorMap[Cell1];

                                if (Color0 != 0 && Color1 != 0 && Color0 != Color1) {
                                    // Color0 should be equal to color1

                                    Debug.Assert(Color0 > 0);
                                    Debug.Assert(Color1 > 0);

                                    AddEqPairing(Color0, Color1);
                                }
                            }
                        }
                    }

                    // synchronize equalities
                    // ----------------------
                    var globColEqArr = locColEq.MPIGatherO(0);
                    Dictionary<int, int> Remappings;
                    if (GridDat.MpiRank == 0) {
                        var globColEq = globColEqArr[0];
                        Debug.Assert(object.ReferenceEquals(globColEq, locColEq));

                        for (int rnk = 1; rnk < GridDat.MpiSize; rnk++) {
                            var rnkColEq = globColEqArr[rnk];

                            foreach (var kv in rnkColEq) {
                                var col = kv.Key;
                                var eqCols = kv.Value;
                                foreach (int eqCol in eqCols) {
                                    if (eqCol != col)
                                        AddEqPairing(col, eqCol);
                                }
                            }
                        }

                        Remappings = new Dictionary<int, int>();
                        foreach (var EqSet in globColEq.Values) {
                            int Cnew = EqSet.Min(); // from equal colors, we pick the minimum
                            foreach (int otherCol in EqSet) {
                                if (otherCol != Cnew) {
                                    if (Remappings.ContainsKey(otherCol)) {
                                        Debug.Assert(Remappings[otherCol] == Cnew);
                                    } else {
                                        Remappings.Add(otherCol, Cnew);
                                    }
                                }
                            }
                        }

                    } else {
                        Remappings = null;
                    }
                    Remappings = Remappings.MPIBroadcast(0);

                    // re-paint shared parts
                    // ---------------------

                    if (Remappings.Count > 0) {
                        for (int iEdge = 0; iEdge < NoEdg; iEdge++) {
                            Debug.Assert(Edge2Cell[iEdge, 0] < J, "The external/ghost cell is expected to be the OUT-cell.");
                            int Cell1 = Edge2Cell[iEdge, 1];
                            if (Cell1 >= J) {
                                // reached an MPI boundary
                                int Cell0 = Edge2Cell[iEdge, 0];
                                int Color0 = ColorMap[Cell0];


                                if (Color0 != 0 && Remappings.ContainsKey(Color0)) {
                                    int cSoll = Remappings[Color0];
                                    if (cSoll != Color0) {
                                        // part has to be re-painted

                                        RepaintRecursive(cSoll, ColorMap, Cell0, this.GridDat);
                                    }
                                }
                            }
                        }
                        ColorMap.MPIExchange(GridDat);
                    }

                }
                VerifyColoring(this.GridDat, ColorMap);

                


                // correlate with old colors
                // =========================

                if(GetPreviousRegion() != null) {
                    
                    // build dictionary: (new color --> old color(s)) which can be send across MPI boundaries
                    // --------------------------------------------------------------------------------------

                    var ColorRecord = new List<int>();
                    var oldColors = new HashSet<int>();

                    void AddColorRecord(int NewColor) {
                        ColorRecord.Add(NewColor);
                        ColorRecord.Add(oldColors.Count);
                        ColorRecord.AddRange(oldColors);
                    }

                    int[] oldColorMap = GetPreviousRegion().ColorMap4Spc[SpId];
                    VerifyColoring(this.GridDat, oldColorMap);
                    
                    BitArray marker = new BitArray(Je);
                    for (int j = 0; j < J; j++) {
                        int newColor = ColorMap[j]; // beim dritten durchgang stehen hier nur nullen 
                        if (newColor != 0 && marker[j] == false) {
                            oldColors.Clear();
                            FindColorsRecursive(oldColors, marker, j, newColor, ColorMap, oldColorMap, this.GridDat);
                            Debug.Assert(oldColors.Contains(0) == false);
                            AddColorRecord(newColor);
                        }
                    }

                    // map new colors to final colors 
                    // ------------------------------

                    // this is done "serially", i.e. only on processor 0; a parallel approach is tricky;
                    // the following part *does not* scale, but i hope it will have 
                    // no effect in the foreseeable future (fk, 27mar19)
                    int[] rcvCounts = ColorRecord.Count.MPIGather(0); // something is wrong here!
                    int[] CollectedColorRecord = ColorRecord.ToArray().MPIGatherv(rcvCounts);

                    Dictionary<int, int> new2finallyNewColor; 
                    if (GridDat.MpiRank == 0) {
                        new2finallyNewColor = new Dictionary<int, int>();
                        var finallyNewColors = new HashSet<int>();

                        int MaxOldColor = 0;
                        int cnt = 0;
                        while (cnt < CollectedColorRecord.Length) {
                            //int NewColor = CollectedColorRecord[cnt]; cnt++;
                            cnt++;
                            int NoOfOldColors = CollectedColorRecord[cnt]; cnt++;

                            for (int kkk = 0; kkk < NoOfOldColors; kkk++) {
                                int OldColor = CollectedColorRecord[cnt]; cnt++;
                                MaxOldColor = Math.Max(MaxOldColor, OldColor);
                            }
                        }

                        cnt = 0;
                        int finallyNewCounter = MaxOldColor + 1;
                        while (cnt < CollectedColorRecord.Length) {
                            int NewColor = CollectedColorRecord[cnt]; cnt++;
                            int NoOfOldColors = CollectedColorRecord[cnt]; cnt++;

                            if (!new2finallyNewColor.ContainsKey(NewColor)) {
                                int PickedColor = -1;
                                for (int kkk = 0; kkk < NoOfOldColors; kkk++) {
                                    int OldColor = CollectedColorRecord[cnt + kkk];
                                    if (!finallyNewColors.Contains(OldColor)) {
                                        PickedColor = OldColor;
                                        break;
                                    }
                                }

                                bool newAssigned = false;
                                if (PickedColor < 0) {
                                    PickedColor = finallyNewCounter;
                                    finallyNewCounter++;
                                    newAssigned = true;
                                }


                                if (NoOfOldColors <= 0) {
                                    // this is a new-born part
                                }

                                if (NoOfOldColors > 1) {
                                    // this is a merge                                  
                                }

                                if (NoOfOldColors > 0 && newAssigned) {
                                    // this is a split from an other part
                                }


                                Debug.Assert(PickedColor > 0);
                                bool bAdd = finallyNewColors.Add(PickedColor);
                                if (bAdd == false)
                                    throw new ApplicationException("error in algorithm; should never happen.");
                                new2finallyNewColor.Add(NewColor, PickedColor);
                            }

                            cnt += NoOfOldColors;
                        }

                        new2finallyNewColor.MPIBroadcast(0);
                    } else {
                        new2finallyNewColor = null;
                        new2finallyNewColor = new2finallyNewColor.MPIBroadcast(0); // this could be done with a more efficient bcast, but for the moment, just fuck it!
                    }
#if DEBUG
                    {
                        int[] FinalColors = new2finallyNewColor.Values.ToArray();
                        for(int i = 0; i < FinalColors.Length; i++) {
                            for(int j =  i + 1; j < FinalColors.Length; j++) {
                                if (FinalColors[i] == FinalColors[j])
                                    throw new ApplicationException();
                            }
                        }

                    }
#endif
                    /*
                    if (counter >= 26) {
                        //int[] oldColorMap = GetPreviousRegion()?.ColorMap4Spc[SpId];

                        var allLocalCol = new HashSet<int>();
                        SinglePhaseField farbe1 = new SinglePhaseField(new Basis(this.GridDat, 0), "farbe1");
                        SinglePhaseField farbe0 = new SinglePhaseField(new Basis(this.GridDat, 0), "farbe0");
                        for (int j = 0; j < J; j++) {
                            farbe1.SetMeanValue(j, ColorMap[j]);
                            allLocalCol.Add(ColorMap[j]);

                            if (oldColorMap != null)
                                farbe0.SetMeanValue(j, oldColorMap[j]);
                        }


                        MultiphaseCellAgglomerator.Katastrophenplot(new DGField[] { (SinglePhaseField)(this.m_owner.LevelSets[0]), farbe1, farbe0 });
                        Debugger.Launch();
                    }
                    */

                    // now, re-paint all new colors in the final new color
                    // ---------------------------------------------------

                    int PrevNewColor = -1, PrevFinalColor = -1;
                    for (int j = 0; j < J; j++) {
                        int NewColor = ColorMap[j];
                        if(NewColor != 0) {

                            int FinalColor;
                            // just a very primitive cache to save lookups in the dict
                            if (NewColor != PrevNewColor) {
                                if (!new2finallyNewColor.ContainsKey(NewColor))
                                    continue;

                                FinalColor = new2finallyNewColor[NewColor];
                                PrevFinalColor = FinalColor;
                                PrevNewColor = NewColor;
                            } else {
                                FinalColor = PrevFinalColor;
                            }

                            // re-paint if necessary
                            if(NewColor != FinalColor) {
                                RepaintRecursive(FinalColor, ColorMap, j, this.GridDat);
                            }
                        }
                    }

                    ColorMap.MPIExchange(GridDat);
                }

                // check & return
                // ===============

                VerifyColoring(this.GridDat, ColorMap);
                return ColorMap;
            }

            private static void FindColorsRecursive(HashSet<int> OldColors, BitArray marker, int j, int NewColor, int[] ColorMap, int[] OldColorMap, IGridData gdat) { 
                int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = gdat.iLogicalCells.NoOfExternalCells + J;
                Debug.Assert(ColorMap.Length == JE);

                Debug.Assert(ColorMap[j] != 0);
                Debug.Assert(ColorMap[j] == NewColor);
                Debug.Assert(NewColor != 0);
                Debug.Assert(marker[j] == false);

                marker[j] = true;
                int OldColor = OldColorMap[j];
                Debug.Assert(OldColor >= 0);
                if(OldColor > 0)
                    OldColors.Add(OldColorMap[j]);
                
                if (j >= J)
                    return; // end of recursion

                foreach(int jN in gdat.iLogicalCells.CellNeighbours[j]) {
                    if (marker[jN])
                        continue;
                    if (ColorMap[jN] == 0)
                        continue;
                    FindColorsRecursive(OldColors, marker, jN, NewColor, ColorMap, OldColorMap, gdat);
                }
            }

            private static void RepaintRecursive(int NewColor, int[] ColorMap, int j, IGridData gdat) {
                int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = gdat.iLogicalCells.NoOfExternalCells + J;
                Debug.Assert(ColorMap.Length == JE);

                Debug.Assert(ColorMap[j] != 0);
                Debug.Assert(ColorMap[j] != NewColor);
                Debug.Assert(NewColor != 0);
                ColorMap[j] = NewColor;

                if (j >= J)
                    return; // end of recursion

                foreach(int jN in gdat.iLogicalCells.CellNeighbours[j]) {
                    if (ColorMap[jN] == NewColor)
                        continue;
                    if (ColorMap[jN] == 0)
                        continue;
                    RepaintRecursive(NewColor, ColorMap, jN, gdat);
                }
            }

            private static void RecursiveColoring(IGridData g, BitArray Msk, int j, int Color, int[] ColorMap, ref bool IsIsolated) {
                Debug.Assert(Msk[j] == true, "illegal to call on non-occupied cells");
                int J = g.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = g.iLogicalCells.NoOfExternalCells + J;
                Debug.Assert(Msk.Length == JE);
                Debug.Assert(ColorMap.Length == JE);
                //Debug.Assert(oldColorMap == null || oldColorMap.Length == JE);
                //bool incremental = oldColorMap != null;

                //int NextColor = Color + 1;

                        /*
                if (incremental) {
                    int oldColor_j = oldColorMap[j];

                    // we care about the old colors
                    if (oldColor_j != 0) {
                        OldColors.Add(oldColor_j);

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

                                }

                                NextColor = Math.Max(NextColor, Color + 1);

                                ColorNegotiable = false;
                            } else {
                                // this is a topology change/a merge
                                NextColor = Math.Max(NextColor, oldColorMap[j] + 1);

                            }
                        }
                    }
                }
                        */

                ColorMap[j] = Color;
                //Part.Add(j);

                if (j >= J)
                    // external cell -> no further recursion
                    //return NextColor;
                    return;

                int[] jNeigh = g.iLogicalCells.CellNeighbours[j];


                foreach (int jN in jNeigh) {
                    if (Msk[jN] == false)
                        // Neighbor cell does not contain species -> end of recursion
                        continue;

                    if (jN >= J) {
                        // external cell -> no further recursion
                        IsIsolated = false;
                    }

                    if (ColorMap[jN] > 0 && jN < J) {
                        // note: there may be the case that a part splits in two on the local MPI processor,
                        //       but the part is connected through another processor; 
                        //       thats why we can't test external cells.

                        // already colored -> end of recursion
                        if (ColorMap[jN] != Color) {
                            /*
                            using (var stw = new System.IO.StreamWriter("megafut.csv")) {
                                //foreach (int i in new[] { j, jN }) {
                                for(int i = 0; i < JE; i++) {
                                    var cen = g.GlobalNodes.GetValue_Cell(Grid.RefElements.Square.Instance.Center, i, 1);
                                    double x = cen[0, 0, 0];
                                    double y = cen[0, 0, 1];
                                    stw.WriteLine("{0}\t{1}\t{2}", x, y, ColorMap[i]);
                                }

                            }
                            */

                            throw new ApplicationException("error in Algorithm, cell " + jN); // Debug.Assert would also be fine, *if* our homies would ever run DEBUG
                        }
                        continue;
                    }

                    RecursiveColoring(g, Msk, jN, Color, ColorMap, ref IsIsolated);
                    //int recNextColor = RecursiveColoring(g, Msk, jN, ref Color, ColorMap, oldColorMap, ref ColorNegotiable, Part, UsedColors, ref IsIsolated);
                    //NextColor = Math.Max(NextColor, recNextColor);
                }

                //return NextColor;
                return;
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