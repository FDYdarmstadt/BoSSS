using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Sometimes, this provides indeed a correct agglomeration graph ;), 
    /// the result is stored in <see cref="AgglomerationPairs"/>.
    /// If agglomeration fails -- which it does quite regularly -- unleash the <see cref="Katastrophenplot"/> to see the mess!
    /// 
    /// </summary>
    /// <remarks>
    /// Cell agglomeration is used to handle two problems:
    /// first, for the treatment of very small cut cells and, for temporally evolving interfaces, to 
    /// ensure an equal topology of the (agglomerated) XDG cut-cell mesh for all involved temporal levels.
    /// - the issue of small cut cells is described in the paper:
    ///   _Extended discontinuous Galerkin methods for two-phase flows: the spatial discretization; Kummer; IJNMF 109 (2), 2017_. 
    /// - the agglomeration of _newborn_ and _decased_ cells is described in 
    ///   the paper: _Time integration for extended discontinuous Galerkin methods with moving domains; Kummer, Müller, Utz; IJNMF 113 (5), 2018_.
    /// - cell agglomeration is performed per species, this implementation looks only at one species at time;
    /// </remarks>
    public class AgglomerationAlgorithm {

        /// <summary>
        /// Temporary feature; will be removed in future;
        /// Plotting if agglomeration fails.
        /// </summary>
        public static Action<DGField[]> Katastrophenplot;


        /// <summary>
        /// 
        /// </summary>
        public AgglomerationAlgorithm(LevelSetTracker __Tracker, SpeciesId __spId, int CutCellsQuadOrder,
            double AgglomerationThreshold, double[] oldTs__AgglomerationTreshold, double NewbornAndDecasedThreshold,
            bool AgglomerateNewborn, bool AgglomerateDeceased,
            bool ExceptionOnFailedAgglomeration) {


            Tracker = __Tracker;
            spId = __spId;
            this.AgglomerateDeceased = AgglomerateDeceased;
            this.AgglomerateNewborn = AgglomerateNewborn;
            this.AgglomerationThreshold = AgglomerationThreshold;
            this.oldTs__AgglomerationTreshold = oldTs__AgglomerationTreshold;
            this.ExceptionOnFailedAgglomeration = ExceptionOnFailedAgglomeration;
            this.NewbornAndDecasedThreshold = NewbornAndDecasedThreshold;

            var NonAgglomeratedMetrics = Tracker.GetXDGSpaceMetrics(spId, CutCellsQuadOrder, 1).CutCellMetrics;
            CellVolumes = NonAgglomeratedMetrics.CutCellVolumes[spId];
            edgeArea = NonAgglomeratedMetrics.CutEdgeAreas[spId];

            {
                CutCellMetrics[] oldCcm;
                if (AgglomerateNewborn || AgglomerateDeceased) {
                    oldCcm = new CutCellMetrics[oldTs__AgglomerationTreshold.Length];
                    for (int iHistory = 0; iHistory < oldCcm.Length; iHistory++) {
                        oldCcm[iHistory] = Tracker.GetXDGSpaceMetrics(spId, CutCellsQuadOrder, -iHistory).CutCellMetrics;
                    }
                } else {
                    oldCcm = null;
                }

                if ((oldCcm == null) != (oldTs__AgglomerationTreshold == null)) {
                    throw new ArgumentException();
                }

                if (oldCcm != null) {
                    if (oldCcm.Length != oldTs__AgglomerationTreshold.Length)
                        throw new ArgumentException();

                    foreach (double alpha in oldTs__AgglomerationTreshold) {
                        if (alpha < 0.0 || alpha >= 1.0)
                            throw new ArgumentOutOfRangeException();
                    }
                }

                oldCellVolumes = oldCcm != null ? oldCcm.Select(a => a.CutCellVolumes[spId]).ToArray() : null;
            }

            // execute algorithm

            var src = FindAgglomerationSources();
            FindAgglomerationTargets_Mk2(src.AgglomCellsList, src.AgglomCellsBitmask, src.AggCandidates);
        }

        double AgglomerationThreshold;
        double NewbornAndDecasedThreshold;
        double[] oldTs__AgglomerationTreshold;
        bool AgglomerateNewborn;
        bool AgglomerateDeceased;


        MultidimensionalArray CellVolumes;
        
        MultidimensionalArray edgeArea;


        MultidimensionalArray[] oldCellVolumes;

        bool ExceptionOnFailedAgglomeration;




        /// <summary>
        /// 
        /// </summary>
        public LevelSetTracker Tracker {
            get;
            private set;
        }
        
        
        /// <summary>
        /// species, for which the agglomeration graph is determined
        /// </summary>
        public SpeciesId spId {
            get;
            private set;
        }

        GridData grdDat => Tracker.GridDat;

        /// <summary>
        /// 1st pass: of agglomeration algorithm, Identification of agglomeration sources
        /// </summary>
        protected virtual (List<int> AgglomCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates) FindAgglomerationSources(
            ) //
        {

            using (var tracer = new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                //tracer.InfoToConsole = true;
                tracer.Info("agglom newborn? " + AgglomerateNewborn);
                tracer.Info("agglom decased? " + AgglomerateDeceased);
                tracer.Info("AgglomerationThreshold = " + AgglomerationThreshold + ",  NewbornAndDecasedThreshold = " + NewbornAndDecasedThreshold);
                if (oldTs__AgglomerationTreshold != null)
                    tracer.Info("oldTs__AgglomerationTreshold = " + oldTs__AgglomerationTreshold.ToConcatString("", "| ", ";"));

                // init 
                // =======

                
                int[,] Edge2Cell = grdDat.Edges.CellIndices;
                byte[] EdgeTags = grdDat.Edges.EdgeTags;
                var Cell2Edge = grdDat.Cells.Cells2Edges;
                int NoOfEdges = grdDat.Edges.Count;
                int myMpiRank = Tracker.GridDat.MpiRank;
                int Jup = grdDat.Cells.NoOfLocalUpdatedCells;
                int Jtot = grdDat.Cells.Count;
                Partitioning CellPart = Tracker.GridDat.CellPartitioning;
                var GidxExt = Tracker.GridDat.Parallel.GlobalIndicesExternalCells;
                var GidxExt2Lidx = Tracker.GridDat.Parallel.Global2LocalIdx;
                long i0 = CellPart.i0;
                //double[] RefVolumes = grdDat.Grid.RefElements.Select(Kref => Kref.Volume).ToArray();

                if (CellVolumes.GetLength(0) != Jtot)
                    throw new ArgumentException();




                // determine agglomeration source cells
                // ================================
                //var _AccEdgesMask = new BitArray(NoOfEdges);
                var AgglomCellsBitmask = new BitArray(Jup);
                //var _AgglomCellsEdges = new BitArray(NoOfEdges);
                //var _AllowedEdges =  new BitArray(NoOfEdges); _AllowedEdges.SetAll(true); // AllowedEdges.GetBitMask();


                // mask for the cells in which we -- potentially -- want to do agglomeration
                var AggCandidates = Tracker.Regions.GetSpeciesMask(spId).GetBitMaskWithExternal().CloneAs();
                var SpeciesMask = Tracker.Regions.GetSpeciesMask(spId).GetBitMask();

                // pass 1: determine agglomeration sources
                // ---------------------------------------
                List<int> AgglomCellsList = new List<int>();
                {
                    // for the present timestep
                    // - - - - - - - - - - - - - 

                    CellMask suspectsForAgg = Tracker.Regions.GetCutCellMask().Intersect(Tracker.Regions.GetSpeciesMask(spId));
                    foreach (int jCell in suspectsForAgg.ItemEnum) {
                        double totVol = grdDat.Cells.GetCellVolume(jCell);

                        double spcVol = CellVolumes[jCell];
                        double alpha = AgglomerationThreshold;
                        spcVol = Math.Max(spcVol, 0.0);
                        double frac = spcVol / totVol;
                        //if (!(frac >= 0.0 && frac <= 1.00000001))
                        //    throw new Exception("Strange volume fraction: cell " + jCell + ", species" + spId.ToString() + ", volume fraction =" + frac + ".");
                        frac = Math.Min(1.0, Math.Max(0.0, frac));


                        //
                        // NOTE !!!!!!!!!!
                        // Do not exclude empty cells here! Empty cells (volume is zero or negative) must be agglomerated to 
                        // yield a correct matrix structure.
                        //
                        if (frac <= alpha) {
                            // cell 'jCell' should be agglomerated to some other cell
                            AgglomCellsBitmask[jCell] = true;
                            AgglomCellsList.Add(jCell);
                        }
                    }
                }

                int NoTimeLev = oldTs__AgglomerationTreshold != null ? oldTs__AgglomerationTreshold.Length : 0;
                if (NoTimeLev > 0) {
                    // for the previous timestep
                    // - - - - - - - - - - - - - 

                    //ushort[][] PrevRegions = new ushort[NoTimeLev][];
                    //CellMask suspectsForAgg = Tracker.PreviousRegions.GetSpeciesMask(spId);
                    //CellMask[] suspectsForAgg = new CellMask[NoTimeLev];
                    //for(int itl = 0; itl < NoTimeLev; itl++) {
                    //    PrevRegions[itl] = Tracker.RegionsHistory[-itl].LevelSetRegionsCode;
                    //    suspectsForAgg[itl] = Tracker.RegionsHistory[-itl].GetSpeciesMask(spId);
                    //}

                    LevelSetSignCode[] signCodes = Tracker.GetLevelSetSignCodes(spId);
                    int NoOfLevSets = Tracker.LevelSets.Count;

                    for (int iTimeLev = 0; iTimeLev < NoTimeLev; iTimeLev++) {
                        CellMask suspectsForAgg = Tracker.RegionsHistory[-iTimeLev].GetCutCellMask().Intersect(Tracker.RegionsHistory[-iTimeLev].GetSpeciesMask(spId));
                        foreach (int jCell in suspectsForAgg.ItemEnum) {


                            double totVol = grdDat.Cells.GetCellVolume(jCell);
                            double spcVol = oldCellVolumes[iTimeLev][jCell];
                            double alpha = oldTs__AgglomerationTreshold[iTimeLev];
                            spcVol = Math.Max(spcVol, 0.0);
                            double frac = spcVol / totVol;
                            frac = Math.Min(1.0, Math.Max(0.0, frac));

                            if (frac < alpha) {
                                // cell 'jCell' should be agglomerated to some other cell
                                if (!AgglomCellsBitmask[jCell]) {
                                    AgglomCellsBitmask[jCell] = true;
                                    AgglomCellsList.Add(jCell);

                                    //Console.WriteLine("Must agglom cell " + jCell + "#" + Tracker.GetSpeciesName(spId) + " volume frac is " + frac + "on time level " + iTimeLev);

                                }
                            }
                        }
                    }
                }

                if (AgglomerateNewborn) {

                    for (int j = 0; j < Jup; j++) {
                        double vol = grdDat.Cells.GetCellVolume(j);
                        double volNewFrac_j = Math.Max(CellVolumes[j], 0.0) / vol;
                        volNewFrac_j = Math.Min(1.0, Math.Max(0.0, volNewFrac_j));


                        if (volNewFrac_j > NewbornAndDecasedThreshold) {
                            for (int nTs = 0; nTs < oldCellVolumes.Length; nTs++) {

                                double volOldFrac_j = Math.Max(oldCellVolumes[nTs][j], 0.0) / vol;
                                volOldFrac_j = Math.Min(1.0, Math.Max(0.0, volOldFrac_j));
                                if (volOldFrac_j <= NewbornAndDecasedThreshold) {
                                    // cell exists at new time, but not at some old time -> newborn

                                    int jNewbornCell = j;
                                    AggCandidates[jNewbornCell] = false;
                                    if (!AgglomCellsBitmask[jNewbornCell]) {
                                        AgglomCellsList.Add(jNewbornCell);
                                        AgglomCellsBitmask[jNewbornCell] = true;

                                        Console.WriteLine("Must agglom NEWBORN cell " + AgglomCellsBitmask + "#" + Tracker.GetSpeciesName(spId));

                                    }
                                }
                            }
                        }
                    }
                }

                if (AgglomerateDeceased) {

                    for (int j = 0; j < Jup; j++) {
                        double vol = grdDat.Cells.GetCellVolume(j);
                        double volNewFrac_j = Math.Max(CellVolumes[j], 0.0) / vol;
                        volNewFrac_j = Math.Min(1.0, Math.Max(0.0, volNewFrac_j));


                        if (volNewFrac_j <= NewbornAndDecasedThreshold) {
                            for (int nTs = 0; nTs < oldCellVolumes.Length; nTs++) {

                                double volOldFrac_j = Math.Max(oldCellVolumes[nTs][j], 0.0) / vol;
                                volOldFrac_j = Math.Min(1.0, Math.Max(0.0, volOldFrac_j));
                                if (volOldFrac_j > NewbornAndDecasedThreshold) {
                                    // cell does not exist at new time, but at some old time -> decased

                                    int jNewbornCell = j;
                                    AggCandidates[jNewbornCell] = false;
                                    if (!AgglomCellsBitmask[jNewbornCell]) {
                                        AgglomCellsList.Add(jNewbornCell);
                                        AgglomCellsBitmask[jNewbornCell] = true;

                                        Console.WriteLine("Must agglom DEAD cell " + AgglomCellsBitmask + "#" + Tracker.GetSpeciesName(spId));

                                    }
                                }

                            }
                        }

                    }


                }
                //*/


                /*
                if (AgglomerateNewborn || AgglomerateDeceased) {
                    //CellMask oldSpeciesCells = this.Tracker.LevelSetData.PreviousSubGrids[spId].VolumeMask;
                    CellMask newSpeciesCells = Tracker.Regions.GetSpeciesMask(spId);


                    // only accept cells with positive volume (new species cells)
                    BitArray newSpeciesCellsBitmask = newSpeciesCells.GetBitMask().CloneAs();
                    Debug.Assert(newSpeciesCellsBitmask.Count == Jup);
                    for (int j = 0; j < Jup; j++) {
                        double vol = grdDat.Cells.GetCellVolume(j);
                        double volFrac_j = Math.Max(CellVolumes[j], 0.0) / vol;
                        volFrac_j = Math.Min(1.0, Math.Max(0.0, volFrac_j));


                        if (volFrac_j > NewbornAndDecasedThreshold) {
                            Debug.Assert(newSpeciesCellsBitmask[j] == true);
                        } else {
                            newSpeciesCellsBitmask[j] = false;
                        }
                    }
                    newSpeciesCells = new CellMask(grdDat, newSpeciesCellsBitmask);

                    // only accept cells with positive volume (old species cells)
                    BitArray oldSpeciesCellsBitmask = new BitArray(Jup); //  oldSpeciesCells.GetBitMask().CloneAs();
                    //Debug.Assert(oldSpeciesCellsBitmask.Count == Jup);
                    for (int nTs = 0; nTs < oldCellVolumes.Length; nTs++) {
                        MultidimensionalArray _oldCellVolumes = oldCellVolumes[nTs];

                        for (int j = 0; j < Jup; j++) {
                            double vol = grdDat.Cells.GetCellVolume(j);
                            double volFrac_j = Math.Max(_oldCellVolumes[j],0.0) / vol;
                            volFrac_j = Math.Min(1.0, Math.Max(0.0, volFrac_j));

                            if (volFrac_j > NewbornAndDecasedThreshold) {
                                //Debug.Assert(oldSpeciesCellsBitmask[j] == true);
                                oldSpeciesCellsBitmask[j] = true; 
                            } else {
                                //oldSpeciesCellsBitmask[j] = false;
                            }
                        }
                    }
                    var oldSpeciesCells = new CellMask(grdDat, oldSpeciesCellsBitmask);

                    // find newborn and decased
                    CellMask newBorn = newSpeciesCells.Except(oldSpeciesCells);
                    CellMask deceased = oldSpeciesCells.Except(newSpeciesCells);


                    foreach (int jNewbornCell in newBorn.ItemEnum) {
                        AggCandidates[jNewbornCell] = false;
                        if (!AgglomCellsBitmask[jNewbornCell]) {
                            AgglomCellsList.Add(jNewbornCell);
                            AgglomCellsBitmask[jNewbornCell] = true;
                        }
                        //Console.WriteLine("  agglom newborn: " + jNewbornCell);
                    }

                    foreach (int jDeceasedCell in deceased.ItemEnum) {
                        AggCandidates[jDeceasedCell] = false;
                        if (!AgglomCellsBitmask[jDeceasedCell]) {
                            AgglomCellsList.Add(jDeceasedCell);
                            AgglomCellsBitmask[jDeceasedCell] = true;
                        }
                        //Console.WriteLine("  agglom deceased: " + jDeceasedCell);
                    }

                }
                //*/

                return (AgglomCellsList, AgglomCellsBitmask, AggCandidates);

            }
        }


        /// <summary>
        /// 2nd pass of agglomeration algorithm, Identification of agglomeration targets
        /// </summary>
        /// <remarks>
        /// Old version, used until Dec. 2021
        /// </remarks>
        protected virtual void FindAgglomerationTargets(
            List<int> AgglomCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates
            ) {
            using (new FuncTrace()) {

                
                var Cell2Edge = grdDat.Cells.Cells2Edges;
                int[,] Edge2Cell = grdDat.Edges.CellIndices;
                int NoOfEdges = grdDat.Edges.Count;
                byte[] EdgeTags = grdDat.Edges.EdgeTags;
                int myMpiRank = Tracker.GridDat.MpiRank;
                int Jup = grdDat.Cells.NoOfLocalUpdatedCells;

                Partitioning CellPart = Tracker.GridDat.CellPartitioning;
                var GidxExt = Tracker.GridDat.Parallel.GlobalIndicesExternalCells;


                var AgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

                if (edgeArea.GetLength(0) != NoOfEdges)
                    throw new ArgumentException();

                double EmptyEdgeTreshold = 1.0e-10; // edges with a measure below or equal to this threshold are
                //                                     considered to be 'empty', therefore they should not be used for agglomeration;
                //                                     there is, as always, an exception: if all inner edges which belong to a
                //                                     cell that should be agglomerated, the criterion mentioned above must be ignored.



                // pass 2: determine agglomeration targets
                // ---------------------------------------

                var failCells = new List<int>();
                foreach (int jCell in AgglomCellsList) {
                    var Cell2Edge_jCell = Cell2Edge[jCell];
                    //bool[] EdgeIsNonempty = new bool[Cell2Edge_jCell.Length];
                    //int[] jNeigh = new int[Cell2Edge_jCell.Length];
                    //bool[] isAggCandidate = new bool[Cell2Edge_jCell.Length];
                    //bool[] passed1 = new bool[Cell2Edge_jCell.Length];

                    // cell 'jCell' should be agglomerated to some other cell
                    Debug.Assert(AgglomCellsBitmask[jCell] == true);

                    double frac_neigh_max = -1.0;
                    int e_max = -1;
                    int jEdge_max = int.MinValue;
                    int jCellNeigh_max = int.MinValue;

                    int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;

                    //bool print = false;
                    //if (jCell == 29 || jCell == 30) {
                    //    print = true;
                    //    Console.WriteLine("Looking for agglom for cell " + jCell + "#" + Tracker.GetSpeciesName(spId) + " volume is " + CellVolumes[jCell]);
                    //}


                    // determine if there is a non-empty edge which connects cell 'jCell'
                    // to some other cell
                    bool NonEmptyEdgeAvailable = false;
                    for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                        int iEdge = Cell2Edge_jCell[e];
                        int OtherCell;
                        if (iEdge < 0) {
                            // cell 'jCell' is the OUT-cell of edge 'iEdge'
                            OtherCell = 0;
                            iEdge *= -1;
                        } else {
                            OtherCell = 1;
                        }
                        iEdge--;
                        int jCellNeigh = Edge2Cell[iEdge, OtherCell];

                        double EdgeArea_iEdge = edgeArea[iEdge];
                        if (jCellNeigh >= 0 && EdgeArea_iEdge > EmptyEdgeTreshold) {
                            //EdgeIsNonempty[e] = true;
                            NonEmptyEdgeAvailable = true;
                        }
                    }




                    // search for some neighbor cell to agglomerate to:
                    for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                        int iEdge = Cell2Edge_jCell[e];
                        int OtherCell, ThisCell;
                        if (iEdge < 0) {
                            // cell 'jCell' is the OUT-cell of edge 'iEdge'
                            OtherCell = 0;
                            ThisCell = 1;
                            iEdge *= -1;
                        } else {
                            OtherCell = 1;
                            ThisCell = 0;
                        }
                        iEdge--;

                        double EdgeArea_iEdge = edgeArea[iEdge];


                        //_AgglomCellsEdges[iEdge] = true;

                        Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                        int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                        //if (print) {
                        //    Console.WriteLine("  testing with cell " + jCellNeigh);
                        //    Console.WriteLine("    connecting edge area: " + EdgeArea_iEdge);
                        //}
                        //jNeigh[e] = jCellNeigh;
                        if (jCellNeigh < 0 || EdgeTags[iEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG || (EdgeArea_iEdge <= EmptyEdgeTreshold && NonEmptyEdgeAvailable)) {
                            // boundary edge, no neighbour for agglomeration
                            Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell, "sollte aber so sein");
                            continue;
                        }
                        //passed1[e] = true;
                        //isAggCandidate[e] = AggCandidates[jCellNeigh];
                        if (!AggCandidates[jCellNeigh])
                            // not suitable for agglomeration
                            continue;

                        // volume fraction of neighbour cell
                        double spcVol_neigh = CellVolumes[jCellNeigh];
                        double totVol_neigh = grdDat.Cells.GetCellVolume(jCellNeigh);
                        double frac_neigh = spcVol_neigh / totVol_neigh;

                        //if (print)
                        //    Console.WriteLine("    neighbour fraction: " + frac_neigh);

                        // max?
                        if (frac_neigh > frac_neigh_max) {
                            frac_neigh_max = frac_neigh;
                            e_max = e;
                            jCellNeigh_max = jCellNeigh;
                            jEdge_max = iEdge;
                        }
                    }

                    //if (print)
                    //    Console.WriteLine("  selected: " + jCellNeigh_max);

                    if (jCellNeigh_max < 0) {

                        //
                        // no agglomeration target found yet; 
                        //

                        // 2nd try:
                        // Note: at this point, i see no reason, why this second try (in below) should have any different result
                        //       It might be some refactoring artefact - since this algorithm is so critical to many computations,
                        //       i don't dare to change it right now.
                        //       Fk, 14dec21

                        for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                            int iEdge = Cell2Edge_jCell[e];
                            int OtherCell, ThisCell;
                            if (iEdge < 0) {
                                // cell 'jCell' is the OUT-cell of edge 'iEdge'
                                OtherCell = 0;
                                ThisCell = 1;
                                iEdge *= -1;
                            } else {
                                OtherCell = 1;
                                ThisCell = 0;
                            }
                            iEdge--;

                            double EdgeArea_iEdge = edgeArea[iEdge];

                            //_AgglomCellsEdges[iEdge] = true;

                            Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                            int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                            //jNeigh[e] = jCellNeigh;
                            if (jCellNeigh < 0 || EdgeTags[iEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG || (EdgeArea_iEdge <= EmptyEdgeTreshold && NonEmptyEdgeAvailable)) {
                                // boundary edge, no neighbour for agglomeration
                                Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell, "sollte aber so sein");
                                //continue;
                            }
                            //passed1[e] = true;
                            //isAggCandidate[e] = AggCandidates[jCellNeigh];
                            if (jCellNeigh < 0 || !AggCandidates[jCellNeigh])
                                // not suitable for agglomeration
                                continue;

                            // volume fraction of neighbour cell
                            double spcVol_neigh = CellVolumes[jCellNeigh];
                            //double totVol_neigh = RefVolumes[grdDat.Cells.GetRefElementIndex(jCellNeigh)]; 
                            double totVol_neigh = grdDat.Cells.GetCellVolume(jCellNeigh);
                            double frac_neigh = spcVol_neigh / totVol_neigh;

                            // max?
                            if (frac_neigh > frac_neigh_max) {
                                frac_neigh_max = frac_neigh;
                                e_max = e;
                                jCellNeigh_max = jCellNeigh;
                                jEdge_max = iEdge;
                            }
                        }

                        if (jCellNeigh_max < 0) {
                            failCells.Add(jCell);
                        } else {
                            //_AccEdgesMask[jEdge_max] = true;

                            int jCellNeighRank;
                            if (jCellNeigh_max < Jup) {
                                jCellNeighRank = myMpiRank;
                            } else {
                                jCellNeighRank = CellPart.FindProcess(GidxExt[jCellNeigh_max - Jup]);
                            }

                            AgglomerationPairs.Add(new CellAgglomerator.AgglomerationPair() {
                                jCellTarget = jCellNeigh_max,
                                jCellSource = jCell,
                                OwnerRank4Target = jCellNeighRank,
                                OwnerRank4Source = myMpiRank
                            });
                        }
                    } else {
                        //_AccEdgesMask[jEdge_max] = true;

                        int jCellNeighRank;
                        if (jCellNeigh_max < Jup) {
                            // agglomeration target on local processor
                            jCellNeighRank = myMpiRank;
                        } else {
                            // inter-process-agglomeration
                            jCellNeighRank = CellPart.FindProcess(GidxExt[jCellNeigh_max - Jup]);
                        }

                        AgglomerationPairs.Add(new CellAgglomerator.AgglomerationPair() {
                            jCellTarget = jCellNeigh_max,
                            jCellSource = jCell,
                            OwnerRank4Target = jCellNeighRank,
                            OwnerRank4Source = myMpiRank
                        });
                    }
                }

                if (failCells.Count.MPISum() > 0) {
                    PlotFail(CellVolumes, oldCellVolumes, AgglomCellsList, ExceptionOnFailedAgglomeration, failCells);

                }

                // store & return
                // ================
                this.AgglomerationPairs = AgglomerationPairs.Select(pair => (pair.jCellSource, pair.jCellTarget)).ToArray();

            }
        }


        /// <summary>
        /// 2nd pass of agglomeration algorithm, Identification of agglomeration targets
        /// </summary>
        /// <remarks>
        /// Revised algorithm, in use since Dec. 2021
        /// </remarks>
        protected virtual void FindAgglomerationTargets_Mk2(
            List<int> AgglomCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates
            ) {
            using (new FuncTrace()) {


                var Cell2Edge = grdDat.Cells.Cells2Edges;
                int[,] Edge2Cell = grdDat.Edges.CellIndices;
                int NoOfEdges = grdDat.Edges.Count;
                byte[] EdgeTags = grdDat.Edges.EdgeTags;
                int myMpiRank = Tracker.GridDat.MpiRank;
                int Jup = grdDat.Cells.NoOfLocalUpdatedCells;

                Partitioning CellPart = Tracker.GridDat.CellPartitioning;
                var GidxExt = Tracker.GridDat.Parallel.GlobalIndicesExternalCells;


                var AgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

                if (edgeArea.GetLength(0) != NoOfEdges)
                    throw new ArgumentException();

                double EmptyEdgeTreshold = 1.0e-10; // edges with a measure below or equal to this threshold are
                //                                     considered to be 'empty', therefore they should not be used for agglomeration;
                //                                     there is, as always, an exception: if all inner edges which belong to a
                //                                     cell that should be agglomerated, the criterion mentioned above must be ignored.



                // pass 2: determine agglomeration targets
                // ---------------------------------------

                var failCells = new List<int>();
                foreach (int jCell in AgglomCellsList) {
                    var Cell2Edge_jCell = Cell2Edge[jCell];
                    //bool[] EdgeIsNonempty = new bool[Cell2Edge_jCell.Length];
                    //int[] jNeigh = new int[Cell2Edge_jCell.Length];
                    //bool[] isAggCandidate = new bool[Cell2Edge_jCell.Length];
                    //bool[] passed1 = new bool[Cell2Edge_jCell.Length];

                    // cell 'jCell' should be agglomerated to some other cell
                    Debug.Assert(AgglomCellsBitmask[jCell] == true);

                    double frac_neigh_max = -1.0;
                    int e_max = -1;
                    int jEdge_max = int.MinValue;
                    int jCellNeigh_max = int.MinValue;

                    int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;

                  

                    // determine if there is a non-empty edge which connects cell 'jCell' to some other cell
                    bool NonEmptyEdgeAvailable = false;
                    for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                        int iEdge = Cell2Edge_jCell[e];
                        int OtherCell;
                        if (iEdge < 0) {
                            // cell 'jCell' is the OUT-cell of edge 'iEdge'
                            OtherCell = 0;
                            iEdge *= -1;
                        } else {
                            OtherCell = 1;
                        }
                        iEdge--;
                        int jCellNeigh = Edge2Cell[iEdge, OtherCell];

                        double EdgeArea_iEdge = edgeArea[iEdge];
                        if (jCellNeigh >= 0 && EdgeArea_iEdge > EmptyEdgeTreshold) {
                            //EdgeIsNonempty[e] = true;
                            NonEmptyEdgeAvailable = true;
                        }
                    }

                    // for strange reasons, one might encounter (with Saye rules) 
                    // empty cells with non-empty edges...
                    // (Could be problematic for many reasons, but here we just "filter" those cases)
                    if (CellVolumes[jCell] <= 0)
                        NonEmptyEdgeAvailable = false;


                    // search for some neighbor cell to agglomerate to:
                    for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                        int iEdge = Cell2Edge_jCell[e];
                        int OtherCell, ThisCell;
                        if (iEdge < 0) {
                            // cell 'jCell' is the OUT-cell of edge 'iEdge'
                            OtherCell = 0;
                            ThisCell = 1;
                            iEdge *= -1;
                        } else {
                            OtherCell = 1;
                            ThisCell = 0;
                        }
                        iEdge--;

                        double EdgeArea_iEdge = edgeArea[iEdge];


                        //_AgglomCellsEdges[iEdge] = true;

                        Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                        int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                        //if (print) {
                        //    Console.WriteLine("  testing with cell " + jCellNeigh);
                        //    Console.WriteLine("    connecting edge area: " + EdgeArea_iEdge + " NonEmptyEdgeAvailable ? " + NonEmptyEdgeAvailable);
                        //}

                        // check exclusion criteria for neighbour cells: 
                        if (jCellNeigh < 0 // boundary edge
                            || EdgeTags[iEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG // no agglomeration across periodic edge
                            || (EdgeArea_iEdge <= EmptyEdgeTreshold && NonEmptyEdgeAvailable) // edge is empty and there might be another non-empty candidate
                            ) {
                            // +++++++++++++++++++++++++++++++++++++++++++++
                            // Neighbor is ruled out:
                            // boundary edge, periodic edge, or
                            // no neighbor for agglomeration

                            //Console.WriteLine($"    Ignoring: r1? {jCellNeigh < 0}, r2? {EdgeTags[iEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG}, r3? {(EdgeArea_iEdge <= EmptyEdgeTreshold && NonEmptyEdgeAvailable)}");


                            Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell, "sollte aber so sein");
                            continue;
                        }
                        //passed1[e] = true;
                        //isAggCandidate[e] = AggCandidates[jCellNeigh];
                        if (!AggCandidates[jCellNeigh])
                            // not suitable for agglomeration
                            continue;

                        // volume fraction of neighbour cell
                        double spcVol_neigh = CellVolumes[jCellNeigh];
                        double totVol_neigh = grdDat.Cells.GetCellVolume(jCellNeigh);
                        double frac_neigh = spcVol_neigh / totVol_neigh;

                        //if (print)
                        //    Console.WriteLine("    neighbour fraction: " + frac_neigh);

                        // max?
                        if (frac_neigh > frac_neigh_max) {
                            frac_neigh_max = frac_neigh;
                            e_max = e;
                            jCellNeigh_max = jCellNeigh;
                            jEdge_max = iEdge;
                        }
                    }

                                        
                    {
                        if(jCellNeigh_max < 0) {
                            failCells.Add(jCell);
                        }

                        //_AccEdgesMask[jEdge_max] = true;

                        int jCellNeighRank;
                        if (jCellNeigh_max < Jup) {
                            // agglomeration target on local processor
                            jCellNeighRank = myMpiRank;
                        } else {
                            // inter-process-agglomeration
                            jCellNeighRank = CellPart.FindProcess(GidxExt[jCellNeigh_max - Jup]);
                        }

                        AgglomerationPairs.Add(new CellAgglomerator.AgglomerationPair() {
                            jCellTarget = jCellNeigh_max,
                            jCellSource = jCell,
                            OwnerRank4Target = jCellNeighRank,
                            OwnerRank4Source = myMpiRank
                        });
                    }
                }

                if (failCells.Count.MPISum() > 0) {
                    PlotFail(CellVolumes, oldCellVolumes, AgglomCellsList, ExceptionOnFailedAgglomeration, failCells);

                }

                // store & return
                // ================
                this.AgglomerationPairs = AgglomerationPairs.Select(pair => (pair.jCellSource, pair.jCellTarget)).ToArray();

            }
        }

        private void PlotFail(MultidimensionalArray CellVolumes, MultidimensionalArray[] oldCellVolumes, List<int> AgglomCellsList, bool ExceptionOnFailedAgglomeration, List<int> failCells) {
            // ++++++++++++++++++++++++++
            // Error handling / reporting
            // ++++++++++++++++++++++++++

            int Jup = grdDat.Cells.NoOfLocalUpdatedCells;

            Basis b = new Basis(grdDat, 0);
            DGField[] CellVolumesViz = new DGField[1 + (oldCellVolumes != null ? oldCellVolumes.Length : 0)];
            for (int n = -1; n < CellVolumesViz.Length - 1; n++) {
                MultidimensionalArray vol = n < 0 ? CellVolumes : oldCellVolumes[n];
                CellVolumesViz[n + 1] = new SinglePhaseField(b, "VolumeAtTime" + (-n));
                for (int j = 0; j < Jup; j++) {
                    CellVolumesViz[n + 1].SetMeanValue(j, vol[j]);
                }
            }

            DGField AgglomCellsViz = new SinglePhaseField(b, "Cells2Agglom");
            foreach (int j in AgglomCellsList) {
                AgglomCellsViz.SetMeanValue(j, 1);
            }

            DGField FailedViz = new SinglePhaseField(b, "FailedCells");
            BitArray FailCellsMarker = new BitArray(Jup);
            foreach (int j in failCells) {
                FailedViz.SetMeanValue(j, 1);
                FailCellsMarker.Set(j, true);
            }

            if (failCells.Count > 0) {
                CellMask FailCells = new CellMask(grdDat, FailCellsMarker);
                FailCells.SaveToTextFile($"failCells{grdDat.MpiRank}.csv");
            }

            DGField[] LevelSets = Tracker.LevelSets.Select(s => (DGField)s).ToArray();

            if (Katastrophenplot != null)
                Katastrophenplot(CellVolumesViz.Cat(AgglomCellsViz, FailedViz, LevelSets));



            string message = ($"Agglomeration failed - no candidate for agglomeration found for {failCells.Count.MPISum()} cells.");
            if (ExceptionOnFailedAgglomeration)
                throw new Exception(message);
            else
                Console.WriteLine(message);
        }


        /// <summary>
        /// Result of the algorithm
        /// </summary>
        public IEnumerable<(int SourceCell, int TargetCell)> AgglomerationPairs {
            get;
            private set;
        }

        


    }




}
