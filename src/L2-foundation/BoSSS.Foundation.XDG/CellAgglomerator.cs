﻿/* =======================================================================
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

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Cell agglomeration for a single species;
    /// This class is responsible for applying a given agglomeration graph onto operator matrices and residual vectors,
    /// see <see cref="ManipulateMatrixAndRHS{M, T}"/> and <see cref="ManipulateRHS{T}"/>.
    /// The agglomeration graph itself is computed by <see cref="AgglomerationAlgorithm"/>.
    /// </summary>
    public class CellAgglomerator {


        /// <summary>
        /// defines a pair (source and target) of cells affected by agglomeration
        /// </summary>
        /// <remarks>
        /// A pair is local when the source cell is local (OwnerRank4Source) but still needs to be present on target ranks (OwnerRank4Target)
        /// </remarks>
        [Serializable]
        public struct AgglomerationPair : IEquatable<AgglomerationPair> {

            /// <summary>
            /// local cell index of agglomeration target cell (the cell where volume is added) 
            /// </summary>
            public int jCellTarget;

            /// <summary>
            /// MPI rank of the process which owns cell <see cref="jCellTarget"/>
            /// </summary>
            public int OwnerRank4Target;

            /// <summary>
            /// local cell index of agglomeration source cell (the cell which is eliminated)
            /// </summary>
            public int jCellSource;

            /// <summary>
            /// MPI rank of the process which owns cell <see cref="jCellSource"/>
            /// </summary>
            public int OwnerRank4Source;

            /// <summary>
            /// level-0 agglomerations have to be done first, level=1 agglomerations second, ...
            /// </summary>
            public int AgglomerationLevel;

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public override int GetHashCode() {
                return jCellSource ^ jCellTarget ^ AgglomerationLevel;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="other"></param>
            /// <returns></returns>
            public bool Equals(AgglomerationPair other) {
                if (object.ReferenceEquals(this, other)) {
                    return true;
                } else if (this.GetHashCode() != other.GetHashCode()) {
                    return false;
                }

                bool result = (this.jCellSource == other.jCellSource)
                    && (this.jCellTarget == other.jCellTarget)
                    && (this.OwnerRank4Target == other.OwnerRank4Target)
                    && (this.OwnerRank4Source == other.OwnerRank4Source)
                    && (this.AgglomerationLevel == other.AgglomerationLevel);
                return result;
            }

            /// <summary>
            /// 
            /// </summary>
            public override string ToString() {
                return $"({jCellSource} [rnk {OwnerRank4Source}] -> {jCellTarget} [rnk {OwnerRank4Target}], Level {AgglomerationLevel})";
            }
        }

        /// <summary>
        /// General information about the structure of the agglomeration, see <see cref="AggInfo"/>
        /// </summary>
        public class AgglomerationInfo : IEquatable<AgglomerationInfo> {

            /// <summary>
            /// The MPI-global maximum over all <see cref="AgglomerationPair.AgglomerationLevel"/> 
            /// for all <see cref="AgglomerationPairs"/>.
            /// </summary>
            public int MaxLevel;

            /// <summary>
            /// true if any agglomeration occurs across MPI boundaries
            /// </summary>
            public bool InterProcessAgglomeration;

            /// <summary>
            /// Edges between agglomerated cells, i.e. pair-wise edges between agglomeration source- and target-cells
            /// </summary>
            public EdgeMask AgglomerationEdges;

            /// <summary>
            /// all cells that are eliminated (aka. 'agglomeration sources', see
            /// <see cref="AgglomerationPair.jCellSource"/>) by agglomeration
            /// (does NOT include agglomeration target cells);
            /// </summary>
            public CellMask SourceCells;

            private CellMask m_TargetCells = null;

            /// <summary>
            /// all cells that are receivers (aka. 'agglomeration targets', see
            /// <see cref="AgglomerationPair.jCellTarget"/>) by agglomeration
            /// (does NOT include agglomeration source cells);
            /// </summary>
            public CellMask TargetCells {
                get {
                    if (m_TargetCells == null) {
                        var grdDat = SourceCells.GridData;
                        int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
                        BitArray msk = new BitArray(J);
                        foreach (var pair in this.AgglomerationPairs) {
                            if (pair.jCellTarget < J) {
                                msk[pair.jCellTarget] = true;
                            }
                        }
                        m_TargetCells = new CellMask(grdDat, msk);
                    }
                    return m_TargetCells;
                }
            }

            /// <summary>
            /// Cell pairs for agglomeration
            /// </summary>
            public AgglomerationPair[] AgglomerationPairs;

            private CellMask m_AllAffectedCells = null;

            /// <summary>
            /// all cells affected by agglomeration (source and target cells)
            /// </summary>
            public CellMask AllAffectedCells {
                get {
                    if (m_AllAffectedCells == null) {
                        var grdDat = SourceCells.GridData;
                        int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
                        BitArray msk = new BitArray(J);
                        foreach (var pair in this.AgglomerationPairs) {
                            if (pair.jCellSource < J)
                                msk[pair.jCellSource] = true;
                            if (pair.jCellTarget < J)
                                msk[pair.jCellTarget] = true;
                        }
                        m_AllAffectedCells = new CellMask(grdDat, msk);
                    }
                    return m_AllAffectedCells;
                }
            }

            /// <summary>
            /// all edges which are connected to at least one cell in
            /// <see cref="SourceCells"/>
            /// </summary>
            public EdgeMask SourceCellsEdges;

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public override int GetHashCode() {
                return this.AgglomerationEdges.GetHashCode();
            }

            /// <summary>
            /// 
            /// </summary>
            public bool Equals(AgglomerationInfo other) {
                if (other == null) {
                    return false;
                } else if (object.ReferenceEquals(this, other)) {
                    return true;
                }

                bool result = (this.AgglomerationPairs.SequenceEqual(other.AgglomerationPairs));
                return result;
            }

            /// <summary>
            /// 
            /// </summary>
            public override string ToString() {
                using (var stw = new StringWriter()) {
                    foreach (var p in AgglomerationPairs) {
                        stw.WriteLine(p.ToString());
                    }

                    return stw.ToString();
                }
            }
        }

        /// <summary>
        /// Grid data.
        /// </summary>
        public GridData GridDat {
            get;
            private set;
        }

        /// <summary>
        /// Agglomeration info.
        /// </summary>
        public CellAgglomerator.AgglomerationInfo AggInfo {
            get;
            private set;
        }

        /// <summary>
        /// Number of agglomerations on all MPI processors.
        /// </summary>
        public int TotalNumberOfAgglomerations {
            get;
            private set;
        }



        /// <summary>
        /// The cell agglomerator v2 (without explicit cycle detection)
        /// </summary>
        /// <param name="g">grid</param>
        /// <param name="AgglomerationPairs">
        /// Local pairs
        /// <see cref="AgglomerationPair"/>)
        /// </param>
        /// <remarks>
        /// The cell agglomerator checks possible problems in the agglomeration pairs coming from <see cref="AgglomerationAlgorithm"/> (Mk3) and then exchange local pairs with targetRanks.
        /// </remarks>
        public CellAgglomerator(GridData g, IEnumerable<AgglomerationPair> AgglomerationPairs) {
            using (new FuncTrace()) {

                // check and init
                // ==============

                MPICollectiveWatchDog.Watch();
                this.GridDat = g;
                int J = g.Cells.NoOfLocalUpdatedCells;
                int JE = g.Cells.Count;
                int mpiRank = g.MpiRank;
                int InterProcessAgglomeration = 0; // used as boolean

                // check cell index validity
                {
                    BitArray sourceUnique = new BitArray(J);
                    foreach (var AggPair in AgglomerationPairs) {
                        int jAggTarget = AggPair.jCellTarget;
                        int jAggSource = AggPair.jCellSource;

                        if (jAggSource < 0 || jAggSource >= J)
                            throw new IndexOutOfRangeException($"Agglomeration source cell (cell index {jAggSource}) must be in the range of locally updated cells.");
                        if (jAggTarget < 0 || jAggTarget >= JE)
                            throw new IndexOutOfRangeException($"Agglomeration target (cell index {jAggTarget}) is not a valid cell index.");
                        if (jAggTarget >= J)
                            InterProcessAgglomeration += 1;

                        if (sourceUnique[jAggSource])
                            throw new ArgumentException("Illegal agglomeration graph: Local cell " + jAggSource + " is occurring in at least two agglomeration pairs as source.");
                        sourceUnique[jAggSource] = true;
                    }
                }

                int NoLocalAggPairs = AgglomerationPairs.Count();
                InterProcessAgglomeration = InterProcessAgglomeration.MPISum();
                List<AgglomerationPair> NeighborPairs = new List<AgglomerationPair>();

                // Collect the interprocess agg. pairs directed to this processor
                if (InterProcessAgglomeration > 0) {
                    int NoOfChangedPairs = DoAggPairsMPIexchangeForOwners(g, AgglomerationPairs.ToList(), ref NeighborPairs);
                    AgglomerationPairs = AgglomerationPairs.Union(NeighborPairs).Distinct().ToList();
                    Debug.Assert(NoOfChangedPairs == InterProcessAgglomeration);
                }

                var AggPairs = AgglomerationPairs.ToList();

                int[] Cells2Aggpairs = new int[JE];
                Cells2Aggpairs.SetAll(-1);
                for (int i = 0; i < AggPairs.Count; i++) {
                    int jAggSrc = AggPairs[i].jCellSource;
                    if (Cells2Aggpairs[jAggSrc] >= 0)
                        throw new ArgumentException("Error in agglomeration graph: a source cell '" + jAggSrc + "' appears twice.");

                    Cells2Aggpairs[jAggSrc] = i;
                }

                void check_Cells2Aggpairs(string inf) {
                    for (int iPair = 0; iPair < AggPairs.Count; iPair++) {
                        if (AggPairs[iPair].jCellSource >= 0) { // ... pair is **not** marked for deletion
                            if (Cells2Aggpairs[AggPairs[iPair].jCellSource] != iPair)
                                throw new ApplicationException($"internal data structure corrupted (Cells2Aggpairs, 1, {inf}): iPair = {iPair}, pair = {AggPairs[iPair]}, Cells2Aggpairs[{AggPairs[iPair].jCellSource}] = {Cells2Aggpairs[AggPairs[iPair].jCellSource]}");
                        }
                    }

                    for (int j = 0; j < JE; j++) {
                        if (Cells2Aggpairs[j] >= 0) {
                            int iPair = Cells2Aggpairs[j];
                            if (AggPairs[iPair].jCellSource != j)
                                throw new ApplicationException($"internal data structure corrupted (Cells2Aggpairs, 2, {inf}): j = {j}, Cells2Aggpairs[{j}] = {iPair}, pair = {AggPairs[iPair]}");
                        }
                    }
                }

                var Level = new List<int>(AggPairs.Count);
                for (int i = 0; i < AggPairs.Count; i++)
                    Level.Add(-1);

                int[] LevelAtCells = new int[JE];
                LevelAtCells.SetAll(-1);

                // determine agglomeration levels by traversing the agglomeration graph
                // --------------------------------------------------------------------
                int Updates = 1;
                int Sweep = 0;
                while (Updates != 0) { // due to MPI parallelization, we have to do multiple sweeps
                                       //                    until nothing changes anymore...
                    Updates = 0; // used like a bool

                    for (int iPair = 0; iPair < AggPairs.Count; iPair++) { // loop over pairs...
                        int CurrentPairIndex = iPair;
                        int CurrentLevel = -1;

                        if (AggPairs[iPair].jCellSource < 0 || AggPairs[iPair].jCellTarget < 0)
                            continue; // pair has been deleted


                        while (CurrentPairIndex >= 0) { // traverse the agglomeration chain...
                            var CurrentPair = AggPairs[CurrentPairIndex];
                            int jAggSrc = CurrentPair.jCellSource;
                            int jAggTrg = CurrentPair.jCellTarget;



                            // update level index
                            int NextLevel = Math.Max(CurrentLevel + 1, LevelAtCells[jAggSrc]);
                            if (NextLevel <= Level[CurrentPairIndex]) {
                                // no need to further traverse this branch
                                break;
                            } else {
                                Updates++; // something changed in this sweep
                            }
                            Level[CurrentPairIndex] = NextLevel;
                            LevelAtCells[jAggSrc] = NextLevel;

                            // move to next pair in chain
                            CurrentPairIndex = Cells2Aggpairs[jAggTrg];
                            CurrentLevel = NextLevel;

                        }
                    }

                    VectorTransceiver_Ext.MPIExchange<int[], int>(LevelAtCells, g);
                    check_Cells2Aggpairs("MPI sweep " + Sweep);
                    Updates = MPIExtensions.MPIMax(Updates);
                    Sweep++;
                }

                int MaxLevel = Level.Count() > 0 ? Level.Max() : 0;
                MaxLevel = MPIExtensions.MPIMax(MaxLevel);

                // re-assign the agglomeration levels (notice AgglomerationPairs are struct (value type variable)
                for (int iPair = 0; iPair < AggPairs.Count; iPair++) { // loop over pairs...
                    var currentPair = AggPairs[iPair];
                    currentPair.AgglomerationLevel = Level[iPair];
                    AggPairs[iPair] = currentPair;
                }

                // mark and check all edges which are used for agglomeration
                // =========================================================

                BitArray AgglomerationEdgesBitMask = new BitArray(g.Edges.Count, false);
                int[,] E2C = g.Edges.CellIndices;
                int[][] C2E = g.Cells.Cells2Edges;
                {
                    // Since we are forming agg groups, too, the agg edges can be between a source and a target as well as between two source cells with the same target, which need to be also removed.
                    // First, we need to distinguish individual agg groups
                    int[] targetCells = AggPairs.Select(p => p.jCellTarget).Distinct().ToArray(); //each group can be identified w.r.t. target cells, which should be unique.
                    foreach (int targetCell in targetCells) { //for each group 
                        var linkedCells = AggPairs.Where(p => p.jCellTarget == targetCell).Select(p => p.jCellSource).ToList(); //gather the source cells of that group
                        linkedCells.Add(targetCell); //target + sources

                        foreach (int linkedCell in linkedCells) { //for each cell in the agg group
                            // exclude non-local cells (if they share an edge with a local cell, then the edge would be marked during the loop of local cell)
                            if (linkedCell >= J)
                                continue;

                            var edges = C2E[linkedCell];
                            foreach (var edge in edges) {
                                int iEdge = Math.Abs(edge) - 1;
                                Debug.Assert(linkedCell == E2C[iEdge, 0] || linkedCell == E2C[iEdge, 1]);

                                int jCellNeighbor = edge > 0 ? E2C[iEdge, 1] : E2C[iEdge, 0];
                                if (linkedCells.Contains(jCellNeighbor)) {
                                    if (iEdge >= 0 && iEdge < AgglomerationEdgesBitMask.Count) { //check if edge is local
                                        AgglomerationEdgesBitMask[iEdge] = true;
                                    }
                                }
                            }
                        }
                    }
                }

                // Save the agglomeration info
                // ============================
                {
                    var ai = new CellAgglomerator.AgglomerationInfo();
                    this.AggInfo = ai;
                    ai.MaxLevel = MaxLevel;

                    BitArray SourceCellBitMask = new BitArray(J);
                    foreach (var pair in AgglomerationPairs) {
                        int jAggSrc = pair.jCellSource;
                        if (jAggSrc < J && jAggSrc >= 0)
                            SourceCellBitMask[jAggSrc] = true;


                    }
                    ai.SourceCells = new CellMask(g, SourceCellBitMask);

                    BitArray SourceCellsEdgesBitMask = new BitArray(g.Edges.Count);
                    foreach (int jCell in ai.SourceCells.ItemEnum) {
                        foreach (int i in C2E[jCell]) {
                            int iEdg = Math.Abs(i) - 1;
                            SourceCellsEdgesBitMask[iEdg] = true;
                        }
                    }
                    ai.SourceCellsEdges = new EdgeMask(g, SourceCellsEdgesBitMask);

                    ai.InterProcessAgglomeration = InterProcessAgglomeration != 0;

                    ai.AgglomerationEdges = new EdgeMask(g, AgglomerationEdgesBitMask);
                    ai.AgglomerationPairs = AggPairs.ToArray();

                    // at this point, because jAggSource must be local, the agglomerations are unique over all MPI processors.
                    this.TotalNumberOfAgglomerations = NoLocalAggPairs.MPISum();
                }
            }
        }


        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="AgglomerationPairs">
        /// Each tuple represents one agglomeration operation:
        ///  - 1st entry: source cell index, i.e. cell which will be removed due to agglomeration; must be in the range of locally updated cells.
        ///  - 2nd entry: target cell index, i.e. cell which will be enlarged due to agglomeration
        /// </param>
        /// <remarks>
        /// The cell agglomerator eliminates (resp. tries to eliminate ;) cycles in the agglomeration graph, since 
        /// it cannot be guaranteed that the <see cref="AgglomerationAlgorithm"/> provides cycle-free agglomeration graphs.
        /// </remarks>
        public CellAgglomerator(GridData g, IEnumerable<(int jSource, int jTarget)> AgglomerationPairs) {

            using (new FuncTrace()) {

                // check and init
                // ==============

                MPICollectiveWatchDog.Watch();

                this.GridDat = g;
                int J = g.Cells.NoOfLocalUpdatedCells;
                int JE = g.Cells.Count;
                int mpiRank = g.MpiRank;

                // check cell index validity
                {
                    BitArray sourceUnique = new BitArray(J);
                    foreach(var AggPair in AgglomerationPairs) {
                        int jAggTarget = AggPair.Item2;
                        int jAggSource = AggPair.Item1;

                        if(jAggSource < 0 || jAggSource >= J)
                            throw new IndexOutOfRangeException($"Agglomeration source cell (cell index {jAggSource}) must be in the range of locally updated cells.");
                        if(jAggTarget < 0 || jAggTarget >= JE)
                            throw new IndexOutOfRangeException($"Agglomeration target (cell index {jAggTarget}) is not a valid cell index.");

                        if(sourceUnique[jAggSource])
                            throw new ArgumentException("Illegal agglomeration graph: Local cell " + jAggSource + " is occurring in at least two agglomeration pairs as source.");
                        sourceUnique[jAggSource] = true;

                    }
                }

                List<(int jSource, int jTarget)> AggPairs = new List<(int, int)>(AgglomerationPairs);
                AgglomerationPairs = null;





                // MPI exchange
                // ============

                // we have to **copy** all agglomeration pairs
                // with a target cell rank on a different processor
                // to the respective processor
                // I.e., some pairs are not unique anymore, but may co-exist on more than one processor.

                int InterProcessAgglomeration; // used as boolean
                List<int> TargetCellMpiRank; // index: correlates with `AggPairs`
                List<int> SourceCellMpiRank; // index: correlates with `AggPairs`
                InterProcessAgglomeration = DoAggPairsMPIexchange(g, AggPairs, out TargetCellMpiRank, out SourceCellMpiRank);




                // mark and check all edges which are used for agglomeration
                // =========================================================

                BitArray AgglomerationEdgesBitMask = new BitArray(g.Edges.Count);
                int[,] E2C = g.Edges.CellIndices;
                int[][] C2E = g.Cells.Cells2Edges;
                {
                    foreach(var pair in AggPairs) {
                        int j1 = pair.jSource;
                        int j2 = pair.jTarget;

                        Debug.Assert(j1 < J || j2 < J);
                        if(j1 >= J) {
                            int a = j2;
                            j2 = j1;
                            j1 = a;
                        }

                        try {
                            int i = C2E[j1].Single(delegate (int k) {
                                int iEdge = Math.Abs(k) - 1;
                                //                          => E2C[Math.Abs(k) - 1, 0] == jAggSrc &&
                                if(k > 0)
                                    return (E2C[iEdge, 0] == j1 && E2C[iEdge, 1] == j2);
                                else if(k < 0)
                                    return (E2C[iEdge, 1] == j1 && E2C[iEdge, 0] == j2);
                                else
                                    throw new ApplicationException();

                            });

                            AgglomerationEdgesBitMask[Math.Abs(i) - 1] = true;
                        } catch(InvalidOperationException) {
                            throw new ArgumentException("Found agglomeration pair which are not neighbor cells.");
                        }
                    }
                }


                // find agglomeration levels
                // =========================

                // level-0 agglomerations have to be done first, level=1 -- agglomerations second, ...

                var Level = new List<int>();
                int MaxLevel;
                {

                    // create mapping
                    // --------------

                    // cell index -> agglomeration pair index.
                    int[] Cells2Aggpairs = new int[JE];
                    Cells2Aggpairs.SetAll(-1);
                    for(int i = 0; i < AggPairs.Count; i++) {
                        int jAggSrc = AggPairs[i].Item1;
                        if(Cells2Aggpairs[jAggSrc] >= 0)
                            throw new ArgumentException("Error in agglomeration graph: a source cell '" + jAggSrc + "' appears twice.");

                        Cells2Aggpairs[jAggSrc] = i;
                    }

                    void check_Cells2Aggpairs(string inf) {
                        for(int iPair = 0; iPair < AggPairs.Count; iPair++) {
                            if(AggPairs[iPair].jSource >= 0) { // ... pair is **not** marked for deletion
                                if(Cells2Aggpairs[AggPairs[iPair].jSource] != iPair)
                                    throw new ApplicationException($"internal data structure corrupted (Cells2Aggpairs, 1, {inf}): iPair = {iPair}, pair = {AggPairs[iPair]}, Cells2Aggpairs[{AggPairs[iPair].jSource}] = {Cells2Aggpairs[AggPairs[iPair].jSource]}");
                            }
                        }

                        for(int j = 0; j < JE; j++) {
                            if(Cells2Aggpairs[j] >= 0) {
                                int iPair = Cells2Aggpairs[j];
                                if(AggPairs[iPair].jSource != j) 
                                    throw new ApplicationException($"internal data structure corrupted (Cells2Aggpairs, 2, {inf}): j = {j}, Cells2Aggpairs[{j}] = {iPair}, pair = {AggPairs[iPair]}");
                            }
                        }
                    }

                    check_Cells2Aggpairs("init");

                    BitArray CycleDetection = new BitArray(JE, false);
                    BitArray MustDeletePair = new BitArray(JE, false);
                    var CycleDetection_sets = new List<int>();

                    void clear_CycleDetection() { // allows O(1) clearing
                        foreach(int i in CycleDetection_sets)
                            CycleDetection[i] = false;
                        CycleDetection_sets.Clear();
                    }
                    void set_CycleDetection(int z) {
                        CycleDetection_sets.Add(z);
                        CycleDetection[z] = true;
                    }



                    // reduce agglomeration levels locally
                    // -----------------------------------

                    // if e.g. we have locally the agglomeration pairs (notation: (jCellSource -> jCellTarget))
                    // (123 -> 234), (234 -> 32), which form an agglomeration chain,
                    // we can replace that by (123 -> 32), (234 ->32) and avoid an agglomeration level, which is expensive.
                    // Thus, agglomeration levels are only necessary, if agglomeration chains cross processor boundaries.


                    for(int iPair = 0; iPair < AggPairs.Count; iPair++) { // loop over pairs...
                        int jSource = AggPairs[iPair].Item1;
                        int jTarget = AggPairs[iPair].Item2;
                        if(jSource < 0 || jTarget < 0)
                            continue; // pair has been filtered

                        if(jSource >= J)
                            // do not reduce over MPI boundaries...
                            // only tricky w.r.t. Extrapolate( ... ), the other stuff would work, however...
                            continue;

                        if(jTarget < J && Cells2Aggpairs[jTarget] >= 0)
                            clear_CycleDetection();

                        while(jTarget < J && Cells2Aggpairs[jTarget] >= 0) { // traverse the agglomeration chain to find its end...
                            //                                                  ... but only as long as the target cell is a local one.
                            if(CycleDetection[jTarget] == true) {

                                // the pair `Cells2Aggpairs[jTarget]` should be removed
                                int currentPair = Cells2Aggpairs[jTarget];

                                // mark as invalid, i.e. set members to illegal values; The pair will be removed later;
                                Cells2Aggpairs[jTarget] = -1;
                                AggPairs[currentPair] = (int.MinValue, int.MinValue);
                                SourceCellMpiRank[currentPair] = -1;
                                TargetCellMpiRank[currentPair] = -1;

                                MustDeletePair[jTarget] = true;

                                //int jAggSrcGlob = (int)GridDat.Parallel.GetGlobalCellIndex(jSource);
                                //int jAggTrgGlob = (int)GridDat.Parallel.GetGlobalCellIndex(jTarget);
                                //Console.WriteLine($"Reduced level for CurrentPair={currentPair}, jAggSrc={jSource}, jAggTrg={jTarget}, jAggSrcGlob={jAggSrcGlob}, jAggTrgGlob={jAggTrgGlob}, rank={mpiRank} offset={g.CellPartitioning.GetI0Offest(mpiRank)}");

                                break;
                            } else {
                                set_CycleDetection(jTarget);
                                int nextPair = Cells2Aggpairs[jTarget];

                                jTarget = AggPairs[nextPair].Item2;
                            }
                        }

                        if(jTarget >= J)
                            // do not reduce over MPI boundaries...
                            // only tricky w.r.t. Extrapolate( ... ), the other stuff would work, however...
                            continue;

                        if(jTarget != AggPairs[iPair].jTarget) {
                            // any local agglomeration chain can be skipped:
                            // we do this by redefining the target.
                            AggPairs[iPair] = (jSource, jTarget);
                        }
                    }



                    void ClearExternalPairs() {
                        MustDeletePair.MPIExchange(g);

                        for(int j = 0; j < J; j++) {
                            if(MustDeletePair[j]) {
                                int currentPair = Cells2Aggpairs[j];
                                if (currentPair >=0 )
                                    if(AggPairs[currentPair].jSource >= 0 || AggPairs[currentPair].jTarget >= 0)
                                        throw new ApplicationException($"Error in Algorithm: " + AggPairs[currentPair].ToString());
                            }
                        }


                        for(int j = J; j < JE; j++) {
                            if(MustDeletePair[j]) {
                                // the pair `Cells2Aggpairs[j]` should be removed
                                int currentPair = Cells2Aggpairs[j];

                                if (currentPair >=0) { // is already assigned as illegal (can happen in other processors)
                                                       // mark as invalid, i.e. set members to illegal values; The pair will be removed later;
                                    //Console.WriteLine($"Queried for CurrentPair={currentPair}, j={j},  AggPairs[currentPair]: ({AggPairs[currentPair].jSource} , {AggPairs[currentPair].jTarget}) on rank viz. ({SourceCellMpiRank[currentPair]} , {TargetCellMpiRank[currentPair]} )  rank={mpiRank} offset={g.CellPartitioning.GetI0Offest(mpiRank)}");

                                    Cells2Aggpairs[j] = -1;
                                    AggPairs[currentPair] = (int.MinValue, int.MinValue);
                                    SourceCellMpiRank[currentPair] = -1;
                                    TargetCellMpiRank[currentPair] = -1;
                                } 
                            }
                        }
                        MustDeletePair.SetAll(false);
                    }

                    ClearExternalPairs();
                    check_Cells2Aggpairs("loc");

                    // init
                    // ----
                    Level = new List<int>(AggPairs.Count);
                    for(int i = 0; i < AggPairs.Count; i++)
                        Level.Add(-1);
                    //Level.SetAll(-1);

                    int[] LevelAtCells = new int[JE];
                    LevelAtCells.SetAll(-1);

                    // determine agglomeration levels by traversing the agglomeration graph
                    // --------------------------------------------------------------------

                    int Updates = 1;
                    int Sweep = 0;
                    while(Updates != 0) { // due to MPI parallelization, we have to do multiple sweeps
                                          //                    until nothing changes anymore...
                        Updates = 0; // used like a bool

                        for(int iPair = 0; iPair < AggPairs.Count; iPair++) { // loop over pairs...
                            int CurrentPairIndex = iPair;
                            int CurrentLevel = -1;

                            if(AggPairs[iPair].jSource < 0 || AggPairs[iPair].jTarget < 0)
                                continue; // pair has been deleted

                            clear_CycleDetection();

                            while(CurrentPairIndex >= 0) { // traverse the agglomeration chain...
                                var CurrentPair = AggPairs[CurrentPairIndex];
                                int jAggSrc = CurrentPair.jSource;
                                int jAggTrg = CurrentPair.jTarget;
                                
                                if(CycleDetection[jAggSrc] == true) {
                                    // the pair `Cells2Aggpairs[jAggTrg]` should be removed
                                    int currentPair = Cells2Aggpairs[jAggTrg];

                                    //int jAggSrcGlob = (int)GridDat.Parallel.GetGlobalCellIndex(jAggSrc);
                                    //int jAggTrgGlob = (int)GridDat.Parallel.GetGlobalCellIndex(jAggTrg);
                                    //Console.WriteLine($"Cycle detected for CurrentPair={currentPair}, jAggSrc={jAggSrc}, jAggTrg={jAggTrg}, jAggSrcGlob={jAggSrcGlob}, jAggTrgGlob={jAggTrgGlob} rank={mpiRank} offset={g.CellPartitioning.GetI0Offest(mpiRank)}");

                                    // mark as invalid, i.e. set members to illegal values; The pair will be removed later;
                                    Cells2Aggpairs[jAggTrg] = -1;
                                    SourceCellMpiRank[currentPair] = -1;
                                    TargetCellMpiRank[currentPair] = -1;
                                    AggPairs[currentPair] = (int.MinValue, int.MinValue);
                                    Level[currentPair] = -1;

                                    MustDeletePair[jAggTrg] = true;

                                    break;

                                    //throw new ArgumentException("Cycle in agglomeration graph (2).");
                                } else {
                                    set_CycleDetection(jAggSrc);

                                    // update level index
                                    int NextLevel = Math.Max(CurrentLevel + 1, LevelAtCells[jAggSrc]);
                                    if(NextLevel <= Level[CurrentPairIndex]) {
                                        // no need to further traverse this branch
                                        break;
                                    } else {
                                        Updates++; // something changed in this sweep
                                    }
                                    Level[CurrentPairIndex] = NextLevel;
                                    LevelAtCells[jAggSrc] = NextLevel;

                                    // move to next pair in chain
                                    CurrentPairIndex = Cells2Aggpairs[jAggTrg];
                                    CurrentLevel = NextLevel;
                                }
                            }
                        }

                        VectorTransceiver_Ext.MPIExchange<int[], int>(LevelAtCells, g);
                        ClearExternalPairs();
                        check_Cells2Aggpairs("MPI sweep " + Sweep);
                        Updates = MPIExtensions.MPIMax(Updates);
                        Sweep++;
                    }

                    MaxLevel = Level.Count() > 0 ? Level.Max() : 0;
                    MaxLevel = MPIExtensions.MPIMax(MaxLevel);
                }

                // Save the agglomeration info
                // ============================
                {
                    var ai = new CellAgglomerator.AgglomerationInfo();
                    this.AggInfo = ai;

                    ai.MaxLevel = MaxLevel;

                    BitArray SourceCellBitMask = new BitArray(J);
                    foreach(var pair in AggPairs) {
                        int jAggSrc = pair.Item1;
                        if(jAggSrc < J && jAggSrc >= 0)
                            SourceCellBitMask[jAggSrc] = true;
                    }
                    ai.SourceCells = new CellMask(g, SourceCellBitMask);

                    BitArray SourceCellsEdgesBitMask = new BitArray(g.Edges.Count);
                    foreach(int jCell in ai.SourceCells.ItemEnum) {
                        foreach(int i in C2E[jCell]) {
                            int iEdg = Math.Abs(i) - 1;
                            SourceCellsEdgesBitMask[iEdg] = true;
                        }
                    }
                    ai.SourceCellsEdges = new EdgeMask(g, SourceCellsEdgesBitMask);

                    ai.InterProcessAgglomeration = InterProcessAgglomeration != 0;


                    //ai.SourceCells = new CellMask(g, SourceCellBitMask);
                    ai.AgglomerationEdges = new EdgeMask(g, AgglomerationEdgesBitMask);

                    ai.AgglomerationPairs = new CellAgglomerator.AgglomerationPair[AggPairs.Count];
                    int ii = 0;
                    int NoOfLocalPairs = 0;
                    int myMpiRank = g.MpiRank;
                    for(int i = 0; i < AggPairs.Count; i++) {
                        if(AggPairs[i].jSource >= 0 && AggPairs[i].jTarget >= 0) { // only take non-deleted pairs
                            ai.AgglomerationPairs[ii].jCellSource = AggPairs[i].jSource;
                            ai.AgglomerationPairs[ii].jCellTarget = AggPairs[i].jTarget;
                            ai.AgglomerationPairs[ii].OwnerRank4Target = TargetCellMpiRank[i];
                            ai.AgglomerationPairs[ii].OwnerRank4Source = SourceCellMpiRank[i];
                            ai.AgglomerationPairs[ii].AgglomerationLevel = Level[i];


                            if(SourceCellMpiRank[ii] == myMpiRank)
                                NoOfLocalPairs++;

                            ii++;
                        }
                    }

                    // at this point, because jAggSource must be local, the agglomerations are unique over all MPI processors.
                    this.TotalNumberOfAgglomerations = NoOfLocalPairs.MPISum();

                    if(ii < AggPairs.Count)
                        Array.Resize(ref ai.AgglomerationPairs, ii);
                    //ai.AgglomerationPairs.SaveToTextFileDebugUnsteady("Agglomerator", ".txt");
                }
            }
        }

        /// <summary>
        /// Exchange the agglomeration pairs with their owners (for target cells), a restricted variant of <see cref="AgglomerationAlgorithm.DoAggPairsMPIexchangeForGhostCells"/>
        /// </summary>
        /// <param name="grdDat"></param>
        /// <param name="LocalAggPairs"></param>
        /// <param name="AggPairsDirectedToThisRank"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private static int DoAggPairsMPIexchangeForOwners(GridData grdDat, List<AgglomerationPair> LocalAggPairs, ref List<AgglomerationPair> AggPairsDirectedToThisRank) {
            int InterProcessAgglomeration = 0;
            int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
            Partitioning CellPart = grdDat.CellPartitioning;
            long j0 = CellPart.i0;
            long[] GidxExt = grdDat.Parallel.GlobalIndicesExternalCells;
            int MyRank = grdDat.MpiRank;

            // Send data
            Dictionary<int, List<CellAgglomerator.AgglomerationPair>> _SendData = new Dictionary<int, List<CellAgglomerator.AgglomerationPair>>(Math.Max(1, (int)Math.Round(((double)(LocalAggPairs.Count)) * 1.03)));
            foreach (var AggPair in LocalAggPairs) {
                int jAggSource = AggPair.jCellSource;
                int jAggTarget = AggPair.jCellTarget;

                int targetRank = AggPair.OwnerRank4Target;
                int sourceRank = AggPair.OwnerRank4Source;


                if (jAggSource >= J) {
                    throw new Exception("Source cells should be in local cell index");
                }

                // is target in another process?
                if (jAggTarget >= J || targetRank != MyRank) { //if so we need to send it
                    InterProcessAgglomeration++;

                    //if (AggPair.AgglomerationLevel < 1)
                    //    throw new Exception("Inter-process agglomeration can not be zero-th level"); //seems that our agg. algorithm does not work properly

                    // source cells are always in local processor
                    int jAggSourceGlob = (int)j0 + jAggSource;

                    //find rank of jTarget
                    long jAggTargetGlob = (GidxExt[jAggTarget - J]);

#if DEBUG
                    int jAggTargetRank = CellPart.FindProcess(jAggTargetGlob);
                    Debug.Assert(jAggTargetRank == targetRank);
#endif

                    List<CellAgglomerator.AgglomerationPair> SendDataList;
                    if (!_SendData.TryGetValue(targetRank, out SendDataList)) {
                        SendDataList = new List<CellAgglomerator.AgglomerationPair>();
                        _SendData.Add(targetRank, SendDataList);
                    }

                    SendDataList.Add(new CellAgglomerator.AgglomerationPair() {
                        jCellSource = jAggSourceGlob,
                        jCellTarget = (int)jAggTargetGlob,
                        OwnerRank4Source = AggPair.OwnerRank4Source,
                        OwnerRank4Target = AggPair.OwnerRank4Target,
                        AgglomerationLevel = AggPair.AgglomerationLevel,
                    });
                }
            }

            Dictionary<int, CellAgglomerator.AgglomerationPair[]> SendData = new Dictionary<int, CellAgglomerator.AgglomerationPair[]>(_SendData.Count);
            foreach (var kv in _SendData) {
                SendData.Add(kv.Key, kv.Value.ToArray());
            }
            _SendData = null;

            InterProcessAgglomeration = MPIExtensions.MPIMax(InterProcessAgglomeration);
            if (InterProcessAgglomeration > 0) {
                // Receive
                var RcvData = SerialisationMessenger.ExchangeData(SendData);
                foreach (var kv in RcvData) {
                    int rcvMpiRank = kv.Key;
                    var ReceivedAggPairs = kv.Value;

                    foreach (var rap in ReceivedAggPairs) {
                        // receive pairs and convert back to local coordinates...
                        long jGlbAggTarget = rap.jCellTarget;
                        long jGlbAggSource = rap.jCellSource;

                        Debug.Assert(CellPart.IsInLocalRange(jGlbAggTarget), "Agglomeration target is expected to be in local cell range.");
                        int jAggTarget = (int)(jGlbAggTarget - j0);

                        Debug.Assert(!CellPart.IsInLocalRange(jGlbAggSource), $"Agglomeration source is expected to be outside of the local cell range. proc-{grdDat.MpiRank}");
                        int jAggSource = grdDat.Parallel.Global2LocalIdx[jGlbAggSource]; //checked((int)(jGlbAggSource - j0));

                        AggPairsDirectedToThisRank.Add(new CellAgglomerator.AgglomerationPair() {
                            jCellSource = jAggSource,
                            jCellTarget = jAggTarget,
                            OwnerRank4Source = rap.OwnerRank4Source,
                            OwnerRank4Target = rap.OwnerRank4Target,
                            AgglomerationLevel = rap.AgglomerationLevel
                        });

                    }
                }
            }
            AggPairsDirectedToThisRank = AggPairsDirectedToThisRank.Distinct().ToList();
            return InterProcessAgglomeration;
        }


        private static int DoAggPairsMPIexchange(GridData g, List<(int jSource, int jTarget)> AggPairs, out List<int> TargetCellMpiRank, out List<int> SourceCellMpiRank) {

            int InterProcessAgglomeration = 0;
            int J = g.iLogicalCells.NoOfLocalUpdatedCells;
            Partitioning CellPart = g.CellPartitioning;
            long j0 = CellPart.i0;
            int mpiRank = g.MpiRank;

            long[] GidxExt = g.Parallel.GlobalIndicesExternalCells;

            TargetCellMpiRank = new List<int>(Math.Max(1, (int)Math.Round(((double)(AggPairs.Count)) * 1.03)));
            SourceCellMpiRank = new List<int>(TargetCellMpiRank.Capacity);

            Dictionary<int, List<(long jGlbSrc, long jGlbTarg)>> _SendData = new Dictionary<int, List<(long, long)>>();
            foreach(var AggPair in AggPairs) {
                // determine all agglomeration pairs which have their target cell on an other processor

                int jAggTarget = AggPair.Item2;
                int jAggSource = AggPair.Item1;

                if(jAggTarget >= J) {
                    // target cell on different MPI process
                    InterProcessAgglomeration = 1;

                    long jGlbTarg = (GidxExt[jAggTarget - J]);
                    int targRank = CellPart.FindProcess(jGlbTarg);

                    TargetCellMpiRank.Add(targRank);

                    List<(long, long)> SendTo;
                    if(!_SendData.TryGetValue(targRank, out SendTo)) {
                        SendTo = new List<(long, long)>();
                        _SendData.Add(targRank, SendTo);
                    }

                    SendTo.Add((j0 + jAggSource, jGlbTarg));

                } else {
                    TargetCellMpiRank.Add(mpiRank);
                }

                SourceCellMpiRank.Add(mpiRank);
            }

            Dictionary<int, (long jGlbSrc, long jGlbTarg)[]> SendData = new Dictionary<int, (long, long)[]>();
            foreach(var kv in _SendData) {
                SendData.Add(kv.Key, kv.Value.ToArray());
            }
            _SendData = null;

            InterProcessAgglomeration = MPIExtensions.MPIMax(InterProcessAgglomeration);
            var RcvData = SerialisationMessenger.ExchangeData(SendData);

            foreach(var kv in RcvData) {
                int rcvMpiRank = kv.Key;
                var ReceivedAggPairs = kv.Value;

                foreach(var rap in ReceivedAggPairs) {
                    // receive pairs and convert back to local coordinates...
                    long jGlbAggTarget = rap.jGlbTarg;
                    long jGlbAggSource = rap.jGlbSrc;

                    Debug.Assert(CellPart.IsInLocalRange(jGlbAggTarget), "Agglomeration target is expected to be in local cell range.");
                    int jAggTarget = checked((int)(jGlbAggTarget - j0));
                    Debug.Assert(jAggTarget >= 0 && jAggTarget < J);

                    Debug.Assert(!CellPart.IsInLocalRange(jGlbAggSource), "Agglomeration source is expected to be outside of the local cell range.");

                    int jAggSource = g.Parallel.Global2LocalIdx[jGlbAggSource];

                    AggPairs.Add((jAggSource, jAggTarget));
                    TargetCellMpiRank.Add(mpiRank);
                    SourceCellMpiRank.Add(rcvMpiRank);
                }
            }


            return InterProcessAgglomeration;
        }

        void InitCouplingMatrices(int maxDeg) {

            if (maxDeg <= CouplingMtxDegree)
                // nothing to do
                return;

            Basis MaxBasis = new Basis(this.GridDat, maxDeg);
            int mpiRank = this.GridDat.CellPartitioning.MpiRank;

            // compute the coupling matrices
            // =============================

            int Esub = this.AggInfo.AgglomerationPairs.Length;
            var CP = new int[Esub, 2];
            for (int i = 0; i < Esub; i++) {
                var ap = this.AggInfo.AgglomerationPairs[i];
                CP[i, 0] = ap.jCellTarget; // can be an external cell
                CP[i, 1] = ap.jCellSource; // is always on this MPI rank
            }

            int N = MaxBasis.Length;

            //MultidimensionalArray Minv = MultidimensionalArray.Create(Esub, N, N);
            MultidimensionalArray M = MultidimensionalArray.Create(Esub, N, N);
            CouplingMtx = M;
            //CouplingMtx_inv = Minv;

            MaxBasis.GetExtrapolationMatrices(CP, M);
        }


        int CouplingMtxDegree = -1;

        
        

        /// <summary>
        /// coupling matrix.<br/>
        /// - 1st index: agglomeration edge index, i.e. index into the items of <see cref="AgglomerationInfo.AgglomerationEdges"/>; 
        /// - 2nd index: row index. 
        /// - 3rd index: column index. 
        /// </summary>
        MultidimensionalArray CouplingMtx;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="map"></param>
        /// <param name="MaxDegree"></param>
        /// <param name="NoOfVars"></param>
        /// <param name="i0Func">
        /// Start indices for blocks
        /// - 1st argument: local cell index
        /// - 2nd argument: variable index
        /// - return: global index of first DOF
        /// </param>
        /// <param name="NjFunc">
        /// Length of blocks
        /// - 1st argument: local cell index
        /// - 2nd argument: variable index
        /// - return: number of DOFs in cell for respective variable
        /// </param>
        /// <param name="MakeInPlace"></param>
        /// <param name="cm"></param>
        /// <returns></returns>
        public BlockMsrMatrix GetRowManipulationMatrix(UnsetteledCoordinateMapping map,
            int MaxDegree, int NoOfVars, Func<int, int, long> i0Func, Func<int, int, int> NjFunc,
            bool MakeInPlace, CellMask cm) {
            using(new FuncTrace()) {

                if(!object.ReferenceEquals(map.GridDat, this.GridDat)) {
                    throw new ArgumentException("grid object mismatch.");
                }

                int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
                int mpiRank = this.GridDat.CellPartitioning.MpiRank;

                this.InitCouplingMatrices(MaxDegree);

                BlockMsrMatrix CompleteMtx = null;
                for(int AggLevel = 0; AggLevel <= this.AggInfo.MaxLevel; AggLevel++) { // loop over agglomeration level...

                    
                    // alloc matrix for the current agglomeration level
                    // ------------------------------------------------

                    BlockMsrMatrix LevelMtx = new BlockMsrMatrix(map, map);

                    if(!MakeInPlace) {
                        // ++++++++++++
                        // if not an in-place matrix, we have to initialize the identity matrix.
                        // ++++++++++++
                        if(cm == null)
                            cm = CellMask.GetFullMask(this.GridDat);

                        foreach(int j in cm.ItemEnum) { // loop over cells...

                            
                            for(int gamma = 0; gamma < NoOfVars; gamma++) { // loop over DG fields...
                                int N = NjFunc(j, gamma);
                                long i0 = i0Func(j, gamma);

                                for(int n = 0; n < N; n++) { // loop over DOFs for variable 'gamma' in cell 'j'...
                                    LevelMtx[i0 + n, i0 + n] = 1.0;
                                }
                            }
                        }
                    }

                    // accumulate agglomeration entries
                    // --------------------------------

                    int iPair = -1;
                    foreach(var pair in this.AggInfo.AgglomerationPairs) {
                        iPair++;

                        //if(pair.AgglomerationLevel == AggLevel)
                        //    Console.Error.WriteLine($"Rk {mpiRank}, Spc {SpeciesName}, Lv {AggLevel}of{this.AggInfo.MaxLevel}: {this.GridDat.iLogicalCells.GetGlobalID(pair.jCellSource)}->{this.GridDat.iLogicalCells.GetGlobalID(pair.jCellTarget)}");

                        if(pair.OwnerRank4Target == mpiRank) {
                            var Aj = this.CouplingMtx.ExtractSubArrayShallow(iPair, -1, -1);

                            if(pair.AgglomerationLevel != AggLevel)
                                continue;

                            for(int gamma = 0; gamma < NoOfVars; gamma++) { // loop over DG fields...
                                int N = NjFunc(pair.jCellSource, gamma);
                                int N2 = NjFunc(pair.jCellTarget, gamma);
                                if(N != N2) {
                                    throw new NotSupportedException("Different number of DOF in source and target is not supported.");
                                }

                                long i0_Row = i0Func(pair.jCellTarget, gamma);
                                long i0_Col = i0Func(pair.jCellSource, gamma);

                                for(int n = 0; n < N; n++) { // loop over DOFs for variable 'gamma' in cell 'j'...
                                    for(int m = 0; m < N; m++) {
                                        LevelMtx[i0_Row + n, i0_Col + m] = Aj[m, n];
                                    }
                                }
                            }
                        }

                        if(pair.OwnerRank4Source == mpiRank) {
                            // clear agglomeration source

                            if(pair.AgglomerationLevel != AggLevel)
                                continue;

                            for(int gamma = 0; gamma < NoOfVars; gamma++) { // loop over DG fields...
                                int N = NjFunc(pair.jCellSource, gamma);
                                long i0 = i0Func(pair.jCellSource, gamma);

                                for(int n = 0; n < N; n++) { // loop over DOFs for variable 'gamma' in cell 'j'...
                                    if(MakeInPlace)
                                        LevelMtx[i0 + n, i0 + n] = -1.0;
                                    else
                                        LevelMtx[i0 + n, i0 + n] = 0.0;
                                }
                            }
                        }
                    }

                    // combine the matrix for the current agg. level 'LevelMtx' with the complete matrix
                    // ---------------------------------------------------------------------------------

                    if(AggLevel == 0) {
                        CompleteMtx = LevelMtx;
                    } else {
                        if(MakeInPlace) {
                            BlockMsrMatrix P = BlockMsrMatrix.Multiply(LevelMtx, CompleteMtx);
                            CompleteMtx.Acc(1.0, LevelMtx);
                            CompleteMtx.Acc(1.0, P);
                        } else {
                            CompleteMtx = BlockMsrMatrix.Multiply(LevelMtx, CompleteMtx);
                        }
                    }

                }

                // return
                // ======

                return CompleteMtx;
            }
        }

        /// <summary>
        /// In a vector <paramref name="V"/>, this clears all entries which correspond to agglomerated cells.
        /// </summary>
        public void ClearAgglomerated<T>(T V, UnsetteledCoordinateMapping RowMap)
            where T : IList<double> {
            if (V.Count != RowMap.LocalLength)
                throw new ArgumentException();

            var Brow = RowMap.BasisS.ToArray();

            int Nrow = Brow.Length;
            var CellPairs = this.AggInfo.AgglomerationPairs;
            int Esub = CellPairs.Length;

            int MyRank = GridDat.MpiRank;

            for (int ii = 0; ii < Brow.Length; ii++) { // loop over codomain basis
                Basis B = Brow[ii];


                for (int esub = 0; esub < Esub; esub++) {

                    int jCell1 = CellPairs[esub].jCellSource;

                    if (CellPairs[esub].OwnerRank4Source == MyRank) {
                        int N = B.GetLength(jCell1);
                        int i0_j1 = RowMap.LocalUniqueCoordinateIndex(ii, jCell1, 0);
                        for (int n = 0; n < N; n++) {
                            V[i0_j1 + n] = 0.0;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// computes the l2-Norm of the DG coordinates
        /// in the agglomerated cells
        /// </summary>
        /// <remarks>
        /// The l2-Norm of the DG coordinates
        /// is NOT equal to the L2-Norm of the DG-field in the XDG-case,
        /// since the basis is no longer orthonormal and therefore Parcival's equation does not hold.
        /// </remarks>
        public double NormInAgglomerated<T>(T V, UnsetteledCoordinateMapping RowMap)
            where T : IList<double> {
            if (V.Count != RowMap.LocalLength)
                throw new ArgumentException();

            var Brow = RowMap.BasisS.ToArray();

            int Nrow = Brow.Length;
            var CellPairs = this.AggInfo.AgglomerationPairs;
            int Esub = CellPairs.Length;

            double l2Norm = 0;

            for (int ii = 0; ii < Brow.Length; ii++) { // loop over codomain basis
                Basis B = Brow[ii];

                int N = B.Length;

                for (int esub = 0; esub < Esub; esub++) {

                    int jCell1 = CellPairs[esub].jCellSource;
                    {
                        int i0_j1 = RowMap.LocalUniqueCoordinateIndex(ii, jCell1, 0);
                        for (int n = 0; n < N; n++) {
                            l2Norm += V[i0_j1 + n].Pow2();
                        }
                    }
                }
            }

            return l2Norm.MPISum().Sqrt();
        }

        /// <summary>
        /// For a list of DG fields, this method performs a
        /// polynomial extrapolation from agglomeration target cells to agglomeration source cells.
        /// </summary>
        public void Extrapolate(CoordinateMapping DgFields) {
            using (new FuncTrace()) {
                int GAMMA = DgFields.Fields.Count;
                var Brow = DgFields.BasisS.ToArray();
                this.InitCouplingMatrices(Brow.Max(basis => basis.Degree));

                DGField[] DgFlds = DgFields.Fields.ToArray();
                                

                Transceiver trx = null;
                int mpiRank = this.GridDat.CellPartitioning.MpiRank;
                if (this.AggInfo.InterProcessAgglomeration) {
                    trx = new Transceiver(DgFlds);
                }

                
                CellAgglomerator.AgglomerationPair[] AggPairs = this.AggInfo.AgglomerationPairs;

                // loop over agglomeration levels: 
                for (int AggLevel = this.AggInfo.MaxLevel; AggLevel >= 0; AggLevel--) {

                    // MPI exchange
                    if (this.AggInfo.InterProcessAgglomeration) {
                        trx.TransceiveStartImReturn();
                        trx.TransceiveFinish();
                    }


                    for (int iPair = 0; iPair < AggPairs.Length; iPair++) {
                        if (AggLevel != AggPairs[iPair].AgglomerationLevel)
                            continue;

                        // for agglomeration, the source is joined to target, i.e. source-DOFs are removed
                        // for the extrapolation, its target to source
                        int jCellTarget = AggPairs[iPair].jCellTarget; // target
                        int jCellSource = AggPairs[iPair].jCellSource; // source (will be overwritten)

                        if (AggPairs[iPair].OwnerRank4Source != mpiRank)
                            continue;

                        Debug.Assert(jCellSource < this.GridDat.Cells.NoOfLocalUpdatedCells);


                        for (int ii = 0; ii < DgFlds.Length; ii++) { // loop over DG fields 
                            Basis B = Brow[ii];
                            int N = B.Length;

                            IMatrix DgCoord = DgFlds[ii].Coordinates;

                            var M_tmp = CouplingMtx.ExtractSubArrayShallow(new int[] { iPair, 0, 0 }, new int[] { iPair - 1, N - 1, N - 1 });

                            //int i0_j0 = RowMap.LocalUniqueCoordinateIndex(ii, jCellTarget, 0);
                            //int i0_j1 = RowMap.LocalUniqueCoordinateIndex(ii, jCellSource, 0);

                            for (int n = 0; n < N; n++) {
                                double acc0 = 0;
                                for (int m = 0; m < N; m++) {
                                    //acc0 += V[i0_j0 + m] * M_tmp[n, m];
                                    acc0 += DgCoord[jCellTarget, m] * M_tmp[n, m];
                                }

                                DgCoord[jCellSource, n] = acc0;
                            }
                        }
                    }
                }

                // MPI exchange
                if (this.AggInfo.InterProcessAgglomeration) {
                    trx.TransceiveStartImReturn();
                    trx.TransceiveFinish();
                }
            }
        }

        /// <summary>
        /// plots "stays" between agglomerated cells;
        /// </summary>
        public void PlotAgglomerationPairs(StreamWriter writer, MultidimensionalArray CenterOfGravity = null, bool includeDummyPointIfEmpty = false) {
            int D = this.GridDat.SpatialDimension;
            var GridDat = this.GridDat;
            var CenterGlobal = MultidimensionalArray.Create(1, 1, D);


            int I = 20; // number of points to represent the stay

            double[] pt0 = new double[D];
            double[] pt1 = new double[D];

            if (includeDummyPointIfEmpty && this.AggInfo.AgglomerationPairs.Length == 0) {
                writer.WriteLine("-1\t0\t0");
            }

            for (int iPair = 0; iPair < this.AggInfo.AgglomerationPairs.Length; iPair++) {
                int jCell0 = this.AggInfo.AgglomerationPairs[iPair].jCellTarget; // cell which stays
                int jCell1 = this.AggInfo.AgglomerationPairs[iPair].jCellSource; // eliminated cell

                if (CenterOfGravity == null) {
                    GridDat.TransformLocal2Global(GridDat.Cells.GetRefElement(jCell0).Center, jCell0, 1, CenterGlobal, 0);
                    for (int d = 0; d < D; d++) {
                        pt0[d] = CenterGlobal[0, 0, d];
                    }
                    GridDat.TransformLocal2Global(GridDat.Cells.GetRefElement(jCell1).Center, jCell1, 1, CenterGlobal, 0);
                    for (int d = 0; d < D; d++) {
                        pt1[d] = CenterGlobal[0, 0, d];
                    }
                } else {
                    for (int d = 0; d < D; d++) {
                        pt0[d] = CenterOfGravity[jCell0, d];
                        pt1[d] = CenterOfGravity[jCell1, d];
                    }
                }

                for (int i = 0 + 5; i < I; i++) {
                    double alpha = ((double)i) / ((double)(I + 1));

                    writer.Write(iPair);
                    writer.Write("\t");

                    for (int d = 0; d < D; d++) {
                        double x_d = (1.0 - alpha) * pt0[d] + alpha * pt1[d];
                        writer.Write(x_d.ToStringDot());
                        writer.Write("\t");
                    }

                    writer.WriteLine();
                }
            }
        }

        /// <summary>
        /// plots "stays" between agglomerated cells;
        /// </summary>
        public void PlotAgglomerationPairs(string filename, MultidimensionalArray CenterOfGravity = null, bool includeDummyPointIfEmpty = false) {
            using (var writer = new StreamWriter(filename)) {
                PlotAgglomerationPairs(writer, CenterOfGravity, includeDummyPointIfEmpty);
            }
        }



        /// <summary>
        /// Apply the basis agglomeration to a right-hand-side vector.
        /// </summary>
        virtual public void ManipulateRHS<T>(T V, UnsetteledCoordinateMapping RowMap, bool[] RowMapAggSw = null)
            where T : IList<double> //
        {
            MPICollectiveWatchDog.Watch();

            if (RowMapAggSw != null)
                throw new NotImplementedException();

            Basis[] BasisS = RowMap.BasisS.ToArray();

            BlockMsrMatrix AggMtx = this.GetRowManipulationMatrix(RowMap, BasisS.Max(basis => basis.Degree), BasisS.Length,
                 (jCell, iVar) => RowMap.GlobalUniqueCoordinateIndex(iVar, jCell, 0),
                 (jCell, iVar) => BasisS[iVar].GetLength(jCell),
                 true, null);

            double[] tmp = V.ToArray();
            if (object.ReferenceEquals(tmp, V))
                throw new ApplicationException("shallow copy detected");

            AggMtx.SpMV(1.0, tmp, 1.0, V);
        }


        /// <summary>
        /// designed to work together with <see cref="MassMatrixFactory"/>
        /// </summary>
        virtual internal void ManipulateMassMatrixBlocks(MultidimensionalArray MatrixBlocks, Basis B, int[] jSub2Cell, Dictionary<int, int> jCell2jSub) {

            int JSGRD = jSub2Cell.Length;
            int N = B.Length;

            if (MatrixBlocks.Dimension != 3)
                throw new ArgumentException();
            if (MatrixBlocks.GetLength(0) != JSGRD || MatrixBlocks.GetLength(1) != N || MatrixBlocks.GetLength(2) != N)
                throw new ArgumentException();


            this.InitCouplingMatrices(B.Degree);
            var CellPairs = this.AggInfo.AgglomerationPairs;
            int NoOfAgg = CellPairs.Length;

            long i0 = this.GridDat.CellPartitioning.i0;
            int MpiSize = this.GridDat.MpiSize;
            int MpiRank = this.GridDat.MpiRank;
            int Jup = this.GridDat.Cells.NoOfLocalUpdatedCells;
            long[] GidxExt = this.GridDat.Parallel.GlobalIndicesExternalCells;
            Partitioning CellPart = this.GridDat.CellPartitioning;



            // loop over agglomeration levels: 
            for (int AggLevel = 0; AggLevel <= this.AggInfo.MaxLevel; AggLevel++) {

                // MPI exchange of blocks 
                // =======================

                // mass matrix blocks from other processors:
                // key: global cell index
                // value: mass matrix block
                Dictionary<long, MultidimensionalArray> RcvData;
                if (this.AggInfo.InterProcessAgglomeration) {

                    // key: processor rank 'p' <br/>
                    // value: a list of cells which must be sent to rank 'p'
                    Dictionary<int, List<int>> SendLines = new Dictionary<int, List<int>>();

                    foreach (var pair in this.AggInfo.AgglomerationPairs) {
                        if (pair.OwnerRank4Source == MpiRank && pair.OwnerRank4Target != MpiRank && pair.AgglomerationLevel == AggLevel) {
                            // the respective block must be sent to processor 'pair.OwnerRank4Target'

                            List<int> SendLine;
                            if (!SendLines.TryGetValue(pair.OwnerRank4Target, out SendLine)) {
                                SendLine = new List<int>();
                                SendLines.Add(pair.OwnerRank4Target, SendLine);
                            }
                            SendLine.Add(pair.jCellSource);
                        }
                    }


                    // key: destination process rank; value: pairs with (global cell index | mass-matrix block)
                    Dictionary<int, Tuple<long, MultidimensionalArray>[]> SendData = new Dictionary<int, Tuple<long, MultidimensionalArray>[]>();
                    foreach (var kv in SendLines) {
                        int dstRank = kv.Key;
                        List<int> Cells2Send = kv.Value; // local cell indices which must be sent to 'dstRank'

                        var SendData_dstRank = new Tuple<long, MultidimensionalArray>[Cells2Send.Count];
                        for (int k = 0; k < SendData_dstRank.Length; k++) {
                            int jCell = Cells2Send[k]; // local cell index
                            Debug.Assert(jCell < Jup);
                            int jSub = jCell2jSub[jCell];
                            long jCellGlob = jCell + i0; // global cell index (for some locally updated cell)

                            var Block = MatrixBlocks.ExtractSubArrayShallow(jSub, -1, -1);
                            SendData_dstRank[k] = new Tuple<long, MultidimensionalArray>(jCellGlob, Block.CloneAs());
                            Block.Clear();
                        }

                        SendData.Add(dstRank, SendData_dstRank);
                    }

                    IDictionary<int, Tuple<long, MultidimensionalArray>[]> RcvDataTmp = SerialisationMessenger.ExchangeData(SendData);
                    RcvData = new Dictionary<long, MultidimensionalArray>();
                    foreach (var kv in RcvDataTmp) {
                        int iProc = kv.Key;
                        Tuple<long, MultidimensionalArray>[] data = kv.Value;
                        foreach (var t in data) {
                            RcvData.Add(t.Item1, t.Item2);
                        }
                    }

                } else {
                    //Debug.Assert(EsubLoc == Esubtot);
                    RcvData = null;
                }

                // perform agglomeration
                // =====================
                int[] RcvDataCnt = new int[MpiSize];
                for (int iPair = 0; iPair < NoOfAgg; iPair++) {

                    int jCellTarget = CellPairs[iPair].jCellTarget; // agglomeration target
                    int jCellSource = CellPairs[iPair].jCellSource; // agglomeration source 

                    if (jCellTarget < Jup) {
                        Debug.Assert(CellPairs[iPair].OwnerRank4Target == MpiRank);

                        var M_tmp = CouplingMtx.ExtractSubArrayShallow(new int[] { iPair, 0, 0 }, new int[] { iPair - 1, N - 1, N - 1 });

                        MultidimensionalArray Block_j0j0 = MatrixBlocks.ExtractSubArrayShallow(jCell2jSub[jCellTarget], -1, -1);
                        MultidimensionalArray Block_j1j1;
                        if (CellPairs[iPair].OwnerRank4Source == MpiRank) {
                            // agglomeration source is local
                            Block_j1j1 = MatrixBlocks.ExtractSubArrayShallow(jCell2jSub[jCellSource], -1, -1);
                        } else {
                            // agglomeration source was received from another process
                            Debug.Assert(jCellSource >= Jup);
                            long jGlobCell1 = GidxExt[jCellSource - Jup];
                            Block_j1j1 = RcvData[jGlobCell1];
                        }

                        // Block_j0j0 += M^T * Block_j1j1 * M
                        Block_j0j0.Multiply(1.0, M_tmp, Block_j1j1, M_tmp, 1.0, "ij", "li", "lk", "kj");
                    }

                    if (jCellSource < Jup) {
                        // Block_j1j1 = 0
                        var Block_j1j1 = MatrixBlocks.ExtractSubArrayShallow(jCell2jSub[jCellSource], -1, -1);
                        Block_j1j1.Clear();
                    }

                }

            }
        }
    

        class MiniMapping {

            Basis[] Basises;
            UnsetteledCoordinateMapping m_cm;

            public MiniMapping(UnsetteledCoordinateMapping cm) {
                this.Basises = cm.BasisS.ToArray();
                this.MaxDeg = this.Basises.Max(b => b.Degree);
                this.NoOfVars = Basises.Length;
                this.m_cm = cm;
            }

            public int NoOfVars;
            public int MaxDeg;

            public long i0Func(int jCell, int iVar) {
                return m_cm.GlobalUniqueCoordinateIndex(iVar, jCell, 0);
            }

            public int NFunc(int jCell, int iVar) {
                return this.Basises[iVar].GetLength(jCell);
            }
        }


        /// <summary>
        /// applies the agglomeration on a general matrix
        /// </summary>
        /// <param name="Matrix">the matrix that should be manipulated.</param>
        /// <param name="Rhs">the right-hand-side that should be manipulated</param>
        /// <param name="ColMap"></param>
        /// <param name="ColMapAggSw">Turns column agglomeration on/off fore each variable individually; default == null is on. </param>
        /// <param name="RowMap"></param>
        /// <param name="RowMapAggSw">The same shit as for <paramref name="ColMapAggSw"/>, just for rows.</param>
        public void ManipulateMatrixAndRHS<M, T>(M Matrix, T Rhs, UnsetteledCoordinateMapping RowMap, UnsetteledCoordinateMapping ColMap, bool[] RowMapAggSw = null, bool[] ColMapAggSw = null)
            where M : IMutableMatrixEx //
            where T : IList<double> //
        {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                //var mtxS = GetFrameMatrices(Matrix, RowMap, ColMap);

                if (Matrix == null && Rhs == null)
                    // nothing to do
                    return;

                if (TotalNumberOfAgglomerations <= 0)
                    // nothing to do
                    return;

                if (RowMapAggSw != null)
                    throw new NotImplementedException();

                // generate agglomeration sparse matrices
                // ======================================

                int RequireRight;
                if (Matrix == null) {
                    // we don't need multiplication-from-the-right at all
                    RequireRight = 0;
                } else {
                    if (RowMap.EqualsUnsetteled(ColMap) && ArrayTools.ListEquals(ColMapAggSw, RowMapAggSw)) {
                        // we can use the same matrix for right and left multiplication
                        RequireRight = 1;

                    } else {
                        // separate matrix for the multiplication-from-the-right is required
                        RequireRight = 2;
                    }
                }

                BlockMsrMatrix LeftMul, RightMul;
                {
                    MiniMapping rowMini = new MiniMapping(RowMap);
                    LeftMul = this.GetRowManipulationMatrix(RowMap, rowMini.MaxDeg, rowMini.NoOfVars, rowMini.i0Func, rowMini.NFunc, false, null);

                    if (RequireRight == 2) {
                        MiniMapping colMini = new MiniMapping(ColMap);
                        RightMul = this.GetRowManipulationMatrix(ColMap, colMini.MaxDeg, colMini.NoOfVars, colMini.i0Func, colMini.NFunc, false, null);
                    } else if (RequireRight == 1) {
                        RightMul = LeftMul;
                    } else {
                        RightMul = null;
                    }
                }

                // apply the agglomeration to the matrix
                // =====================================

                if (Matrix != null) {
                    BlockMsrMatrix RightMulTr = RightMul.Transpose();

                    BlockMsrMatrix _Matrix;
                    if (Matrix is BlockMsrMatrix) {
                        _Matrix = (BlockMsrMatrix)((object)Matrix);
                    } else {
                        _Matrix = Matrix.ToBlockMsrMatrix(RowMap, ColMap);
                    }

                    var AggMatrix = BlockMsrMatrix.Multiply(LeftMul, BlockMsrMatrix.Multiply(_Matrix, RightMulTr));

                    if (object.ReferenceEquals(_Matrix, Matrix)) {
                        _Matrix.Clear();
                        _Matrix.Acc(1.0, AggMatrix);
                    } else {
                        Matrix.Acc(-1.0, _Matrix); //   das ist so
                        Matrix.Acc(1.0, AggMatrix); //  meagaschlecht !!!!!!
                    }
                }

                // apply the agglomeration to the Rhs
                // ==================================

                if (Rhs != null) {

                    double[] tmp = Rhs.ToArray();
                    if (object.ReferenceEquals(tmp, Rhs))
                        throw new ApplicationException("Flache kopie sollte eigentlich ausgeschlossen sein!?");

                    LeftMul.SpMV(1.0, tmp, 0.0, Rhs);
                }
            }
        }
    }
}
