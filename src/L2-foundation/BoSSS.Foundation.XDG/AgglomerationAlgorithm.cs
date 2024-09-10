using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using BoSSS.Foundation.Comm;

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
        /// defines a group of agglomeration sources to create appropriate agglomeration pairs (source and target)
        /// </summary>
        [Serializable]
        class AgglomerationGroup : IEquatable<AgglomerationGroup> {
            /// <summary>
            /// local cell index of agglomeration target cell (the cell where volume is added) 
            /// </summary>
            public int jCellGroupTarget;

            /// <summary>
            /// MPI rank of the process which owns cell <see cref="OwnerRank4GroupTarget"/>
            /// </summary>
            public int OwnerRank4GroupTarget;

            /// <summary>
            /// list of agglomeration source cells
            /// </summary>
            private List<(int jCell, int rank)> Sources;

            /// <summary>
            /// list of neighbor edges via connected edges
            /// </summary>
            private List<int> NeighborCells;

            public List<int> GetNeighbors { 
                get { return NeighborCells.Distinct().ToList(); }  
            }

            public List<int> GetCells {
                get { return Sources.Select(s => s.jCell).Concat(new[] { jCellGroupTarget }).ToList(); }
            }

            private void AddNeighbors(int jCell) {
                var Cell2Edge = m_grdDat.Cells.Cells2Edges;
                int[,] Edge2Cell = m_grdDat.Edges.CellIndices;
                byte[] EdgeTags = m_grdDat.Edges.EdgeTags;
                var Cell2Edge_jCell = Cell2Edge[jCell];
                int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;

                // create an array for neighbors
                double[,] neighbors = new double[NoOfEdges_4_jCell, 2]; // (jCellNeigh, edgeArea) 

                // Collect neighbors and determine if there is a non-empty edge which connects cell 'jCell' to some other cell
                bool NonEmptyEdgeAvailable = false;
                for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                    var (iEdge, ThisCell, OtherCell) = GetEdgeInfo(Cell2Edge_jCell[e]);

                    int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                    double EdgeArea_iEdge = m_edgeArea[iEdge];
                    Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                    neighbors[e, 0] = (double)jCellNeigh;
                    neighbors[e, 1] = EdgeArea_iEdge;

                    if (jCellNeigh >= 0 && EdgeArea_iEdge > EmptyEdgeTreshold) {
                        NonEmptyEdgeAvailable = true;
                    }
                }

                for (int e = 0; e < NoOfEdges_4_jCell; e++) {
                    double jCellNeigh = neighbors[e, 0];
                    double EdgeArea_iEdge = neighbors[e, 1];

                    if (jCellNeigh >= 0 // not a boundary edge
                     && EdgeTags[e] < GridCommons.FIRST_PERIODIC_BC_TAG // no agglomeration across periodic edge
                     && (EdgeArea_iEdge > EmptyEdgeTreshold || !NonEmptyEdgeAvailable)) // if all edges are empty, add this but otherwise check the criteria
                    {
                        NeighborCells.Add((int)jCellNeigh);
                    }
                }


            }

            private GridData m_grdDat;
            private MultidimensionalArray m_CellVolumes;
            private MultidimensionalArray m_edgeArea;
            private LevelSetTracker m_tracker;

            public AgglomerationGroup(int firstCell, LevelSetTracker Tracker, MultidimensionalArray CellVolumes, MultidimensionalArray edgeArea) {
                m_tracker = Tracker;
                m_grdDat = Tracker.GridDat;
                m_CellVolumes = CellVolumes;
                m_edgeArea = edgeArea;
                jCellGroupTarget = firstCell;
                OwnerRank4GroupTarget = m_grdDat.MpiRank;
                Sources = new List<(int, int)>();
                NeighborCells = new List<int>();
                m_sumFractions = Math.Round(m_CellVolumes[firstCell] / m_grdDat.Cells.GetCellVolume(firstCell), 2);
                AddNeighbors(firstCell);
            }

            public bool IsConnected(int CellNumber) {
                return NeighborCells.Contains(CellNumber);
            }

            public bool IsPart(int CellNumber) {
                bool IsInSources = Sources.Select(p => p.jCell).Contains(CellNumber);
                bool IsTarget = jCellGroupTarget == CellNumber;
                return IsInSources || IsTarget;
            }

            private double m_sumFractions;

            /// <summary>
            /// GetHashCode
            /// </summary>
            /// <returns></returns>
            public override int GetHashCode() {
                int hashCode = 0;
                foreach (var Cell in Sources)
                    hashCode += Cell.GetHashCode();
                return hashCode;
            }

            public void Add(int jCell) {
                if (IsPart(jCell))
                    return;

                Partitioning CellPart = m_tracker.GridDat.CellPartitioning;
                var GidxExt = m_tracker.GridDat.Parallel.GlobalIndicesExternalCells;
                int J = m_grdDat.Cells.NoOfLocalUpdatedCells;


                if (!IsConnected(jCell)) {
                    throw new ArgumentException($"The cell with local index={jCell} has not any connection to the target cell");
                }

                int jCellRank;
                if (jCell >= J) {
                    long jCellGlob = (GidxExt[jCell - J]);
                    jCellRank = CellPart.FindProcess(jCellGlob);
                } else {
                    jCellRank = m_tracker.GridDat.MpiRank;
                }

                double newVolume =  Math.Round(m_CellVolumes[jCell] / m_grdDat.Cells.GetCellVolume(jCell), 2);
                m_sumFractions = m_sumFractions + newVolume;
                Sources.Add((jCell, jCellRank));
                AddNeighbors(jCell);
                //Console.WriteLine($"Added jcell={jCell} with fraction {Math.Round(m_CellVolumes[jCell] / m_grdDat.Cells.GetCellVolume(jCell), 2)}");
            }

            public List<CellAgglomerator.AgglomerationPair> GetAggPairs {
                get {
                    var AggPairs = new List<CellAgglomerator.AgglomerationPair>();
                    AggPairs.Add(new CellAgglomerator.AgglomerationPair() {
                        jCellTarget = jCellGroupTarget,
                        jCellSource = jCellGroupTarget,
                        OwnerRank4Target = OwnerRank4GroupTarget,
                        OwnerRank4Source = OwnerRank4GroupTarget,
                        posTarget = m_grdDat.Cells.GetCenter(jCellGroupTarget),
                        fracTarget = m_CellVolumes[jCellGroupTarget],
                        AgglomerationLevel = 0 //OwnerRank4GroupTarget == Source.rank ? 0 : 1
                    });


                    foreach (var Source in Sources) {
                        AggPairs.Add(new CellAgglomerator.AgglomerationPair() {
                            jCellTarget = jCellGroupTarget,
                            jCellSource = Source.jCell,
                            OwnerRank4Target = OwnerRank4GroupTarget,
                            OwnerRank4Source = Source.rank,
                            posTarget = m_grdDat.Cells.GetCenter(jCellGroupTarget),
                            fracTarget = m_CellVolumes[jCellGroupTarget],
                            AgglomerationLevel = 0 //OwnerRank4GroupTarget == Source.rank ? 0 : 1
                        });
                    }
                    return AggPairs;
                }
            }

            public double SumFractions { get => m_sumFractions; }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="other"></param>
            /// <returns></returns>
            public bool Equals(AgglomerationGroup other) {
                if (object.ReferenceEquals(this, other)) {
                    return true;
                } else if (this.GetHashCode() != other.GetHashCode()) {
                    return false;
                }

                if (this.Sources.Count() != other.Sources.Count())
                    return false;

                bool result = true;

                for (int k = 0; k < this.Sources.Count(); k++)
                    result = result && this.Sources.Take(k).Equals(other.Sources.Take(k));

                return result;
            }

            /// <summary>
            /// Returns a string listing the chain target cell and the agglomeration pairs in the chain 
            /// </summary>
            public override string ToString() {
                string str = $"AggGroup to {jCellGroupTarget} [rnk {OwnerRank4GroupTarget}]) with total frac={SumFractions} \n";

                foreach (var Cell in GetAggPairs)
                    str += $"- ({Cell.jCellSource} [rnk {Cell.OwnerRank4Source}] -> {Cell.jCellTarget} [rnk {Cell.OwnerRank4Target}], Level {Cell.AgglomerationLevel}) \n";

                str += "NeighborCells: ";
                foreach (var Cell in NeighborCells) {
                    str += Cell.ToString() + ", ";
                }

                return str;
            }

        }

        /// <summary>
        /// Temporary feature; will be removed in future;
        /// Plotting if agglomeration fails.
        /// </summary>
        public static Action<DGField[], string> Katastrophenplot;

        /// <summary>
        /// If agglomeration plots are demanded. A flag for debugging purposes in the agglomeration algorithm.
        /// </summary>
        public static bool PlotAgglomeration = false;

        /// <summary>
        /// 
        /// </summary>
        public AgglomerationAlgorithm(LevelSetTracker __Tracker, SpeciesId __spId, int CutCellsQuadOrder,
            double AgglomerationThreshold, double[] oldTs__AgglomerationTreshold, double NewbornAndDecasedThreshold,
            bool AgglomerateNewborn, bool AgglomerateDeceased,
            bool ExceptionOnFailedAgglomeration, string _tag) {


            Tracker = __Tracker;
            spId = __spId;
            Tag = string.IsNullOrEmpty(_tag) ? "" : _tag + "_";  
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
            //FindAgglomerationTargets_Mk3(src.AgglomCellsList, src.AgglomCellsBitmask, src.aggCandidates);

            FindAgglomerationTargets_Mk3(src.AgglomCellsList, src.AgglomCellsBitmask, src.AggCandidates, src.NewbornCellsList, src.VanishingCellsList);
        }

        double AgglomerationThreshold;
        double NewbornAndDecasedThreshold;
        double[] oldTs__AgglomerationTreshold;
        bool AgglomerateNewborn;
        bool AgglomerateDeceased;
        string Tag;

        private List<CellAgglomerator.AgglomerationPair> m_AggPairsOnNeighborPairs;
        private List<CellAgglomerator.AgglomerationPair> m_AggPairs;
        private List<int> m_failCells;
        private List<int> m_CellsNeedChainAgglomeration;
        private List<AgglomerationGroup> m_aggGroups;

        MultidimensionalArray CellVolumes;
        MultidimensionalArray edgeArea;
        MultidimensionalArray[] oldCellVolumes;

        bool ExceptionOnFailedAgglomeration;
        bool m_globalAgglomerationMapping = true; //to force agglomeration mapping such that every iteration is done by a processor (see e.g., group agglomeration), practically not needed for most of the cases
        string marker = ""; //marker for debugging purposes

        /// <summary>
        /// edges with a measure below or equal to this threshold are
        ///  considered to be 'empty', therefore they should not be used for agglomeration;
        ///  there is, as always, an exception: if all inner edges which belong to a
        ///  cell that should be agglomerated, the criterion mentioned above must be ignored.
        /// </summary>
        static double EmptyEdgeTreshold = 1.0e-10;

        /// <summary>
        /// reference level set tracker
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

        /// <summary>
        /// reference grid data
        /// </summary>
        GridData grdDat => Tracker.GridDat;

        /// <summary>
        /// returns the unique pairs on this this processor and neighbors pairs (not all pairs of the neighbor procs)
        /// </summary>
        private List<CellAgglomerator.AgglomerationPair> m_AggPairsWithExtNeighborPairs
        {
            get => m_AggPairsOnNeighborPairs.ToList().Union(m_AggPairs).Distinct().ToList();
        }

        /// <summary>
        /// 1st pass: of agglomeration algorithm, Identification of agglomeration sources
        /// </summary>
        /// <returns>
        /// - `AgglomCellsList`: all agglomeration sources
        /// - `AgglomCellsBitmask`: the same as `AgglomCellsList`, just in bit-mask form
        /// - `AggCandidates`: all cells which are allowed as **potential** agglomeration **targets*
        /// </returns>
        protected virtual (List<int> AgglomCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates, List<int> NewbornCellsList, List<int> VanishingCellsList) FindAgglomerationSources(
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
                List<int> VanishingCellsList = new List<int>();
                List<int> NewbornCellsList = new List<int>();

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

                // Topological changes (only valid for dynamic meshes)
                // Determine newborn cells
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
                                    NewbornCellsList.Add(jNewbornCell);
                                    AggCandidates[jNewbornCell] = false;
                                    if (!AgglomCellsBitmask[jNewbornCell]) {
                                        AgglomCellsList.Add(jNewbornCell);
                                        AgglomCellsBitmask[jNewbornCell] = true;

                                        Console.WriteLine("Must agglom NEWBORN cell " + jNewbornCell + "#" + Tracker.GetSpeciesName(spId) + " on RANK-" + myMpiRank);

                                    }
                                }
                            }
                        }
                    }
                }

                // Determine vanishing cells
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

                                    int jVanishingCell = j;
                                    VanishingCellsList.Add(jVanishingCell);
                                    AggCandidates[jVanishingCell] = false;
                                    if (!AgglomCellsBitmask[jVanishingCell]) {
                                        AgglomCellsList.Add(jVanishingCell);
                                        VanishingCellsList.Add(jVanishingCell);
                                        AgglomCellsBitmask[jVanishingCell] = true;

                                        Console.WriteLine("Must agglom DEAD cell " + jVanishingCell + "#" + Tracker.GetSpeciesName(spId) + " on RANK-" + myMpiRank);

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

                return (AgglomCellsList, AgglomCellsBitmask, AggCandidates, NewbornCellsList, VanishingCellsList);

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


                Console.WriteLine("Count of AgglomCellsList:" + AgglomCellsList.Count + " AgglomCellsBitmask: " + AgglomCellsBitmask.Count + " AggCandidates: " + AggCandidates.Count);

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

                //double EmptyEdgeTreshold = 1.0e-10; // edges with a measure below or equal to this threshold are
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
                int NoFailedCells = failCells.Count.MPISum();
                if (NoFailedCells > 0 || PlotAgglomeration) {
                    double[] volFrac = new double[Jup];
                    int Dim = grdDat.Cells.GetCenter(0).Dim;
                    int[] pairIdentification = new int[Jup];
                    int[] pairColor = new int[Jup];
                    Vector[] aggDirection = new Vector[Jup];

                    for (int j = 0; j < Jup; j++) {
                        double totVol = grdDat.Cells.GetCellVolume(j);
                        double spcVol = CellVolumes[j];
                        volFrac[j] = spcVol / totVol;
                        aggDirection[j] = new Vector(Dim);
                    }

                    int k = 1;
                    foreach (var pair in AgglomerationPairs) {
                        if (pair.jCellSource < Jup) {
                            pairIdentification[pair.jCellSource] = pair.jCellTarget;
                            pairColor[pair.jCellSource] = k;
                            Vector direction = new Vector(Dim);
                            direction = grdDat.Cells.GetCenter(pair.jCellTarget) - grdDat.Cells.GetCenter(pair.jCellSource);
                            aggDirection[pair.jCellSource] = direction;
                        }
                        if ((pair.jCellTarget < Jup)) {
                            pairColor[pair.jCellTarget] = k; //can override when it is chain
                        }
                    }

                    int[] AggTargets = AgglomerationPairs.Select(p => p.jCellTarget).Distinct().ToArray();

                    if (NoFailedCells > 0)
                        PlotFail(false, CellVolumes, oldCellVolumes, AgglomCellsList, m_failCells, AggCandidates, AggTargets, pairIdentification, pairColor, aggDirection, volFrac, $"{Tag}AgglomerationKatastrophe{spId.ToString()}");

                        if (PlotAgglomeration)
                        PlotFail(false, CellVolumes, oldCellVolumes, AgglomCellsList, m_failCells, AggCandidates, AggTargets, pairIdentification, pairColor, aggDirection, volFrac, $"{Tag}AgglomerationGraphOf{spId.ToString()}");

                }

                // store & return
                // ================
                this.AgglomerationPairs = AgglomerationPairs.Select(pair => (pair.jCellSource, pair.jCellTarget)).ToArray();

            }
        }

        // depreciated
        //private int DoAggSourceCellsMPIexchangeToTheirNeighbors(List<int> AggSourceCells, ref List<int> AggSourceCellsOnExtNeighborPairs) {
        //    int InterProcessAgglomeration = 0;

        //    int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
        //    Partitioning CellPart = grdDat.CellPartitioning;
        //    long j0 = CellPart.i0;
        //    long[] GidxExt = grdDat.Parallel.GlobalIndicesExternalCells;

        //    Dictionary<int, List<int>> _SendData = new Dictionary<int, List<int>>(Math.Max(1, (int)Math.Round(((double)(AggSourceCells.Count)) * 1.03)));

        //    foreach (int jAggSource in AggSourceCells) {

        //        // neighbors of source cell (to avoid duplication we loop over only the neighbors of source cells in a pair)
        //        var Cell2Edge = grdDat.Cells.Cells2Edges;
        //        var Cell2Edge_jCell = Cell2Edge[jAggSource];
        //        int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;
        //        int[,] Edge2Cell = grdDat.Edges.CellIndices;

        //        // loop over faces/neighbor cells...
        //        for (int e = 0; e < NoOfEdges_4_jCell; e++) {
        //            int iEdge = Cell2Edge_jCell[e];
        //            int OtherCell, ThisCell;
        //            if (iEdge < 0) {
        //                // cell 'jCell' is the OUT-cell of edge 'iEdge'
        //                OtherCell = 0;
        //                ThisCell = 1;
        //                iEdge *= -1;
        //            } else {
        //                OtherCell = 1;
        //                ThisCell = 0;
        //            }
        //            iEdge--;
        //            int jCellNeigh = Edge2Cell[iEdge, OtherCell];
        //            Debug.Assert(Edge2Cell[iEdge, ThisCell] == jAggSource);

        //            // is neighbor in another process? (  if not, no need to send it)
        //            if (jCellNeigh >= J) {
        //                InterProcessAgglomeration++;
        //                //find rank of jCellNeigh
        //                long jGlbNeigh = (GidxExt[jCellNeigh - J]);
        //                int jGlbNeighRank = CellPart.FindProcess(jGlbNeigh);

        //                List<int> SendDataList;
        //                if (!_SendData.TryGetValue(jGlbNeighRank, out SendDataList)) {
        //                    SendDataList = new List<int>();
        //                    _SendData.Add(jGlbNeighRank, SendDataList);
        //                }

        //                // source cells are always in local processor
        //                int jAggSourceGlob = (int)j0 + jAggSource;

        //                if (!SendDataList.Contains(jAggSourceGlob))
        //                    SendDataList.Add(jAggSourceGlob);

        //            }
        //        }
        //    }

        //    Dictionary<int, int[]> SendData = new Dictionary<int, int[]>();
        //    foreach (var kv in _SendData) {
        //        SendData.Add(kv.Key, kv.Value.ToArray());
        //    }
        //    _SendData = null;

        //    InterProcessAgglomeration = MPIExtensions.MPIMax(InterProcessAgglomeration);

        //    if (InterProcessAgglomeration > 0) {
        //        var RcvData = SerialisationMessenger.ExchangeData(SendData);

        //        foreach (var kv in RcvData) {
        //            int rcvMpiRank = kv.Key;
        //            var ReceivedAggSources = kv.Value;

        //            // receive pairs and convert back to local coordinates...
        //            foreach (var jGlbAggSource in ReceivedAggSources) {

        //                Debug.Assert(CellPart.FindProcess((long)jGlbAggSource) == rcvMpiRank);
        //                Debug.Assert(!CellPart.IsInLocalRange(jGlbAggSource), $"The Agglomeration source received is expected to be outside of the local cell range. proc-{grdDat.MpiRank}");

        //                // jGlbAggSource to jAggSource (local)
        //                int jAggSource = grdDat.Parallel.Global2LocalIdx[jGlbAggSource];
        //                AggSourceCellsOnExtNeighborPairs.Add(jAggSource);
        //            }
        //        }
        //    }
        //    AggSourceCellsOnExtNeighborPairs = AggSourceCellsOnExtNeighborPairs.Distinct().ToList();
        //    return InterProcessAgglomeration;
        //}

        /// <summary>
        /// Exchange the local agglomeration pairs with their neighbors (only the ranks that is adjacent to the cells in a pair)
        /// </summary>
        /// <param name="AggPairs">Local pairs <see cref="CellAgglomerator.AgglomerationPair" >See agglomeration pairs</see></param>
        /// <param name="AggPairsOnExtNeighborPairs">Incoming pairs</param>
        /// <returns></returns>
        public int DoAggPairsMPIexchangeForGhostCells(List<CellAgglomerator.AgglomerationPair> AggPairs, ref List<CellAgglomerator.AgglomerationPair> AggPairsOnExtNeighborPairs) {
            int InterProcessAgglomeration = 0;
            int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
            Partitioning CellPart = grdDat.CellPartitioning;
            long j0 = CellPart.i0;
            long[] GidxExt = grdDat.Parallel.GlobalIndicesExternalCells;

            Dictionary<int, List<CellAgglomerator.AgglomerationPair>> _SendData = new Dictionary<int, List<CellAgglomerator.AgglomerationPair>>(Math.Max(1, (int)Math.Round(((double)(AggPairs.Count)) * 1.03)));

            foreach (var AggPair in AggPairs) {
                int jAggSource = AggPair.jCellSource;
                int jAggTarget = AggPair.jCellTarget;

                // neighbors of source cell (to avoid duplication we loop over only the neighbors of source cells in a pair)
                var Cell2Edge = grdDat.Cells.Cells2Edges;
                var Cell2Edge_jCell = Cell2Edge[jAggSource];
                int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;
                int[,] Edge2Cell = grdDat.Edges.CellIndices;

                // loop over faces/neighbor cells...
                for (int e = 0; e < NoOfEdges_4_jCell; e++) {
                    var (iEdge, ThisCell, OtherCell) = GetEdgeInfo(Cell2Edge_jCell[e]);
                    int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                    Debug.Assert(Edge2Cell[iEdge, ThisCell] == jAggSource);

                    // is neighbor in another process? (  if not no need to change it)
                    if (jCellNeigh >= J) {
                        InterProcessAgglomeration++;
                        //find rank of jCellNeigh
                        long jGlbNeigh = (GidxExt[jCellNeigh - J]);
                        int jGlbNeighRank = CellPart.FindProcess(jGlbNeigh);

                        List<CellAgglomerator.AgglomerationPair> SendDataList;
                        if (!_SendData.TryGetValue(jGlbNeighRank, out SendDataList)) {
                            SendDataList = new List<CellAgglomerator.AgglomerationPair>();
                            _SendData.Add(jGlbNeighRank, SendDataList);
                        }

                        // source cells are always in local processor
                        int jAggSourceGlob = (int)j0 + jAggSource;
                        int jGlbTarg;

                        //target cell can be on another processor 
                        if (jAggTarget >= J) {
                            jGlbTarg = (int)GidxExt[jAggTarget - J];
                        } else {
                            jGlbTarg = (int)j0 + jAggTarget;
                        }

                        SendDataList.Add(new CellAgglomerator.AgglomerationPair() {
                            jCellSource = jAggSourceGlob,
                            jCellTarget = jGlbTarg,
                            OwnerRank4Source = AggPair.OwnerRank4Source,
                            OwnerRank4Target = AggPair.OwnerRank4Target,
                            posTarget = AggPair.posTarget,
                            fracTarget = AggPair.fracTarget
                        });

                    }
                }
            }

            Dictionary<int, CellAgglomerator.AgglomerationPair[]> SendData = new Dictionary<int, CellAgglomerator.AgglomerationPair[]>(_SendData.Count);
            foreach (var kv in _SendData) {
                SendData.Add(kv.Key, kv.Value.ToArray());
            }
            _SendData = null;

            InterProcessAgglomeration = MPIExtensions.MPIMax(InterProcessAgglomeration);

            if (InterProcessAgglomeration > 0) {
                var RcvData = SerialisationMessenger.ExchangeData(SendData);
                foreach (var kv in RcvData) {
                    int rcvMpiRank = kv.Key;
                    var ReceivedAggPairs = kv.Value;

                    foreach (var rap in ReceivedAggPairs) {
                        // receive pairs and convert back to local coordinates...
                        long jGlbAggTarget = rap.jCellTarget;
                        long jGlbAggSource = rap.jCellSource;

                        //Debug.Assert(CellPart.IsInLocalRange(jGlbAggTarget), "Agglomeration target is expected to be in local cell range.");
                        //int jGlbAggSourceInt = checked((int)(jGlbAggSource - j0));
                        //Debug.Assert(jGlbAggSourceInt >= 0 && jGlbAggSourceInt < J);

                        Debug.Assert(!CellPart.IsInLocalRange(jGlbAggSource), $"Agglomeration source is expected to be outside of the local cell range. proc-{grdDat.MpiRank}");

                        int jAggSource = grdDat.Parallel.Global2LocalIdx[jGlbAggSource]; //checked((int)(jGlbAggSource - j0));

                        int jAggTarget = -1; //by default assign this to invalid number (-1) to indicate that it is not known to receiving proc
                        if (CellPart.IsInLocalRange(jGlbAggTarget)) {                           // is the target cell in this proc
                            jAggTarget = (int)(jGlbAggTarget - j0);
                        } else if (grdDat.Parallel.Global2LocalIdx.ContainsKey(jGlbAggTarget)) { // is the target cell in external/ghost cell list
                            jAggTarget = grdDat.Parallel.Global2LocalIdx[jGlbAggTarget];
                        }

                        AggPairsOnExtNeighborPairs.Add(new CellAgglomerator.AgglomerationPair() {
                            jCellSource = jAggSource,
                            jCellTarget = jAggTarget,
                            OwnerRank4Source = rap.OwnerRank4Source,
                            OwnerRank4Target = rap.OwnerRank4Target,
                            posTarget = rap.posTarget,
                            fracTarget = rap.fracTarget,
                        });

                    }
                }
            }
            AggPairsOnExtNeighborPairs = AggPairsOnExtNeighborPairs.Distinct().ToList();
            return InterProcessAgglomeration;
        }

        /// <summary>
        /// 2nd pass of agglomeration algorithm, Identification of agglomeration targets
        /// </summary>
        /// <remarks>
        /// Revised algorithm, in use since Dec. 2021
        /// </remarks>
        protected virtual void FindAgglomerationTargets_Mk2(
            List<int> AgglomSourceCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates
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
                m_AggPairs = AgglomerationPairs;

                //exchange the source cell information to be able to know about external/ghost cells (needed only once)
                BitArray AggSourcesWithExternalCell = new BitArray(grdDat.iLogicalCells.Count, false);

                foreach (var jCell in AgglomSourceCellsList)
                    AggSourcesWithExternalCell[jCell] = true;

                AggSourcesWithExternalCell.MPIExchange(grdDat);

                //make sure source cells cannot be target cells
                for (int j = 0; j < Jup; j++) {
                    if (AggSourcesWithExternalCell[j] == true) {
                        AggCandidates[j] = false;
                    }
                }

                //exchange the candidate information to be able to know about external/ghost cells (needed only once)
                AggCandidates.MPIExchange(grdDat);

                if (edgeArea.GetLength(0) != NoOfEdges)
                    throw new ArgumentException();

                double EmptyEdgeTreshold = 1.0e-10; // edges with a measure below or equal to this threshold are
                //                                     considered to be 'empty', therefore they should not be used for agglomeration;
                //                                     there is, as always, an exception: if all inner edges which belong to a
                //                                     cell that should be agglomerated, the criterion mentioned above must be ignored.



                // pass 2: determine agglomeration targets
                // ---------------------------------------
                var failCells = new List<int>();
                var CellsNeedChainAgglomeration = new List<int>();

                // First check if a target can be found among neighbors
                #region DirectAgg
                // first check if source cells can find a "proper" target (excluding vanishing and newborn cells)
                foreach (int jCell in AgglomSourceCellsList) {
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

                    // create an array for neighbors
                    double[,] neighbors = new double[NoOfEdges_4_jCell, 3]; // (jCellNeigh, iEdge, edgeArea) 

                    // Collect neighbors and determine if there is a non-empty edge which connects cell 'jCell' to some other cell
                    bool NonEmptyEdgeAvailable = false;
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
                        int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                        double EdgeArea_iEdge = edgeArea[iEdge];
                        Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                        neighbors[e, 0] = (double)jCellNeigh;
                        neighbors[e, 1] = (double)iEdge;
                        neighbors[e, 2] = EdgeArea_iEdge;

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
                        int jCellNeigh = (int)neighbors[e, 0];
                        int iEdge = (int)neighbors[e, 1];
                        double EdgeArea_iEdge = neighbors[e, 2];

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


                            //Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell, "sollte aber so sein");
                            continue;
                        }

                        if (!AggCandidates[jCellNeigh]) {
                            continue; // not suitable for agglomeration
                        }

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

                        if (jCellNeigh_max < 0) { //is(not) it a proper target?
                            bool isAggChainPossible = false;

                            // check neighbors if they can be both agg. source
                            for (int e = 0; e < NoOfEdges_4_jCell; e++) {
                                if (AggSourcesWithExternalCell[(int)neighbors[e, 0]]) {
                                    isAggChainPossible = true; // there is a possibility to form a chain, which means that the neighbor cell could be carrying element to the final target)
                                    break;
                                }
                            }

                            if (isAggChainPossible) {
                                CellsNeedChainAgglomeration.Add(jCell);
                            } else {
                                failCells.Add(jCell); // jcell has no possible targets to be agglomerated
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
                                OwnerRank4Source = myMpiRank,
                                AgglomerationLevel = 0
                            });
                        }
                    }
                }
                #endregion

                // If there are remaining cells, try to create agg. chains (the level of chains would be handled in CellAgglomerator.cs)
                #region ChainAgg
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                // in certain cases, agglomeration sources can require to create a chain/cluster (a set of cut cells that should be agglomerated to a target) 
                var ChainAgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

                int iteration_threeshold = (int)Math.Pow(Math.Max(CellsNeedChainAgglomeration.Count, 2), 3); // 
                iteration_threeshold = iteration_threeshold.MPIMax();
                int ii = 0;

                // Exchange already paired cells (needs to be updated in every change in the loop below)
                var AggPairsOnExtNeighborPairs = new List<CellAgglomerator.AgglomerationPair>();
                var InterProcessChainAgglomeration = DoAggPairsMPIexchangeForGhostCells(AgglomerationPairs, ref AggPairsOnExtNeighborPairs);
                m_AggPairsOnNeighborPairs = AggPairsOnExtNeighborPairs;

                // Create agglomeration chains
                // m_AgglomerationChains = new List<AgglomerationChain>();

                int ChainCountMax = CellsNeedChainAgglomeration.Count;
                ChainCountMax = ChainCountMax.MPIMax();

                bool anyUpdate = true;

                // Global-MPI level
                while (ChainCountMax > 0 && anyUpdate && ii < iteration_threeshold) { // as long as a cell needs to be agglomerated or not exceeding the threshold
                    anyUpdate = false;
                    var LoopChainAgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

                    List<(double SpeciesFrac, double Dist, int AggLevel, int jCell, int jCellNeigh, int targetCell, int targetRank)> weightedEdges = new List<(double, double, int, int, int, int, int)>();
                    List<int> ImmediateConnectionCells = new List<int>(CellsNeedChainAgglomeration.Count);

                    // Any cell to be agglomerated in a chain on the local proc (a parallel variation of Kruskal's algorithm for Minimum Spanning Forest adapted to our context
                    if (CellsNeedChainAgglomeration.Count > 0) {

                        //create edge information for possible connections, remind that we could only create a chain if there is a neighbor that is already paired with some cell.  
                        foreach (int jCell in CellsNeedChainAgglomeration) {
                            var Cell2Edge_jCell = Cell2Edge[jCell];
                            int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;

                            List<(double SpeciesFrac, double Dist, int AggLevel, int jCell, int jCellNeigh, int targetCell, int targetRank)> weightedEdgesjCell = new List<(double, double, int, int, int, int, int)>(NoOfEdges_4_jCell);

                            // Collect neighbors and determine if there is a possible source cell through which 'jCell' can be connected to a proper target
                            int NeighborAggSourceCells = 0;
                            int PairedNeighborAggSourceCells = 0;
                            for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                                bool IsPossibleTarget = false;
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
                                int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                                double EdgeArea_iEdge = edgeArea[iEdge];
                                Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                                if ((jCellNeigh >= 0 && EdgeTags[iEdge] < GridCommons.FIRST_PERIODIC_BC_TAG) || CellVolumes[jCell] <= 0) {
                                    IsPossibleTarget = true;

                                    if (AggSourcesWithExternalCell[jCellNeigh])
                                        NeighborAggSourceCells++;

                                } else {
                                    continue;
                                }

                                // assume maximum distance as default value for cases that there is no additional info
                                double Distance = double.MaxValue;

                                // Searching for a pair is needed to ensure the neighbor cell has a "real" target, which will be connected in CellAgglomerator.cs later
                                var TargetPairs = m_AggPairsWithExtNeighborPairs.Where(p => p.jCellSource == jCellNeigh);

                                if (TargetPairs.Count() > 1) {
                                    throw new ArgumentException("Agglomeration fail: a cell is mapped to two different target cells");
                                }

                                // only already paired cells can form a target for the chain
                                if (TargetPairs.Any() && IsPossibleTarget) { //check if a neighbor is agglomerated
                                    var TargetPair = TargetPairs.First();
                                    int TargetCellOfTargetPair = TargetPair.jCellTarget; //looking for the final target
                                    int possibleTarget; // due to the mpi boundaries we may select the neighbor as target instead of the final target
                                    int targetRank = -1; //by default assigned an invalid value to ensure code works

                                    // To turn off the direct agg to the final target, assign the below to -1
                                     TargetCellOfTargetPair = -1; // normally the code is able to choose a target that forms higher level agg. chains. However, it is decided to be handled in CellAgglomerator.cs

                                    // if TargetCellOfTargetPair is one of ext/ghost cells or already on local proc
                                    bool IsTargetCellKnown = TargetCellOfTargetPair >= 0 && TargetPairs.First().OwnerRank4Target >= 0;

                                    // We must also ensure that source cell is known to the owner of the target as well (otherwise it will lead to local indexing errors)
                                    bool IsSourceKnownToOwnerOfTargetCell = grdDat.Cells.CellNeighbours[jCell].Contains(TargetCellOfTargetPair) || TargetPairs.First().OwnerRank4Target == myMpiRank;
                                    int AggLevel;
                                    if (IsTargetCellKnown && IsSourceKnownToOwnerOfTargetCell) {
                                        possibleTarget = TargetCellOfTargetPair; //if so, directly pair jCell with the final target
                                        targetRank = TargetPair.OwnerRank4Target;
                                        AggLevel = TargetPairs.First().AgglomerationLevel;
                                    } else {
                                        possibleTarget = jCellNeigh;  //if not, pair jCell with jNeigh (which is also TargetPair.jCellSource). This cell will behave like a "carrying" cell which ultimately would lead to the target cell
                                        targetRank = TargetPair.OwnerRank4Source;
                                        AggLevel = TargetPairs.First().AgglomerationLevel + 1;
                                        Debug.Assert(jCellNeigh == TargetPair.jCellSource);
                                    }

                                    // Create the AggLevel info


                                    double SpeciesFrac = Math.Round(CellVolumes[possibleTarget] / grdDat.Cells.GetCellVolume(possibleTarget), 2);

                                    Vector posTarget = grdDat.Cells.GetCenter(possibleTarget);
                                    Vector posSource = grdDat.Cells.GetCenter(jCell);
                                    Distance = Vector.Dist(posTarget, posSource);
                                    weightedEdgesjCell.Add((SpeciesFrac, Distance, AggLevel, jCell, jCellNeigh, possibleTarget, targetRank));
                                    PairedNeighborAggSourceCells++;
                                }
                            }

                            bool IsAllPossibleNeighborCellsPaired = (PairedNeighborAggSourceCells == NeighborAggSourceCells) && (PairedNeighborAggSourceCells > 0);

                            // If all possible neighbors are already paired, we can choose one of them directly
                            if (IsAllPossibleNeighborCellsPaired) {
                                var ConnectionEdge = weightedEdgesjCell.Where(p => p.jCell == jCell).OrderByDescending(p => p.SpeciesFrac).ThenBy(p => p.Dist).ToList().First();

                                if (ConnectionEdge.targetRank >= 0) {
                                    var newPair = new CellAgglomerator.AgglomerationPair() {
                                        jCellTarget = ConnectionEdge.targetCell,
                                        //jNeighborSource = ConnectionEdge.jCellNeigh,
                                        jCellSource = ConnectionEdge.jCell,
                                        OwnerRank4Target = ConnectionEdge.targetRank,
                                        OwnerRank4Source = myMpiRank,
                                        AgglomerationLevel = ConnectionEdge.AggLevel
                                    };

                                    //var targetChains = this.m_AgglomerationChains.Where(p => p.jCellChainTarget == ConnectionEdge.targetCell);

                                    //// Create agglomeration chain
                                    //if (targetChains.Any()) {
                                    //    if (targetChains.Count() > 1)
                                    //        throw new Exception("A cell is associated with two different chains");

                                    //    var targetChain = targetChains.First();
                                    //    targetChain.Add(newPair);
                                    //} else {
                                    //    var firstPair = m_AggPairsWithExtNeighborPairs.Where(p => p.jCellTarget == ConnectionEdge.targetCell && p.jCellSource == ConnectionEdge.jCellNeigh).First();
                                    //    var newChain = new AgglomerationChain(firstPair, grdDat, CellVolumes);
                                    //    newChain.Add(newPair);
                                    //    this.m_AgglomerationChains.Add(newChain);
                                    //}

                                    LoopChainAgglomerationPairs.Add(newPair);
                                    ImmediateConnectionCells.Add(ConnectionEdge.jCell);

                                    continue;
                                }

                            } else {
                                weightedEdges.AddRange(weightedEdgesjCell);
                            }

                            Debug.Assert(NeighborAggSourceCells < 1, $"No possible chain target for {jCell} on proc-{ilPSP.Environment.MPIEnv.MPI_Rank}");
                        }

                        // discard already connected cells
                        foreach (int DirectConnected in ImmediateConnectionCells) {
                            CellsNeedChainAgglomeration.Remove(DirectConnected);
                            anyUpdate = true;
                        }

                        // sort remaining edges
                        weightedEdges = weightedEdges.OrderByDescending(p => p.SpeciesFrac).ThenBy(p => p.Dist).ToList(); // (descending, ascending)

                        //weightedEdges.SaveToTextFileDebugUnsteady("weightedEdges", ".txt");

                        // if any target available in case a target needed?
                        Debug.Assert(weightedEdges.Any() || ImmediateConnectionCells.Any() || !CellsNeedChainAgglomeration.Any(), "Cell agglomeration failed." +
                                " There are cells that cannot be connected any target cells. (Cycle between cells to be agglomerated");

                        // choose the first edge and add the corresponding agg. pair
                        if (weightedEdges.Any() && CellsNeedChainAgglomeration.Any()) {
                            var aggConnectionEdge = weightedEdges.First();


                            // AddToList
                            if (aggConnectionEdge.targetRank > -1) {
                                var SourceCell2Edge_jCell = Cell2Edge[aggConnectionEdge.jCell];

                                var newPair = new CellAgglomerator.AgglomerationPair() {
                                    jCellTarget = aggConnectionEdge.targetCell,
                                    jCellSource = aggConnectionEdge.jCell,
                                    OwnerRank4Target = aggConnectionEdge.targetRank,
                                    OwnerRank4Source = myMpiRank,
                                    AgglomerationLevel = aggConnectionEdge.AggLevel

                                };

                                LoopChainAgglomerationPairs.Add(newPair);

                                //var targetChains = this.m_AgglomerationChains.Where(p => p.jCellChainTarget == aggConnectionEdge.targetCell);

                                //if (targetChains.Any()) {
                                //    if (targetChains.Count() > 1)
                                //        throw new Exception("A cell is associated with two different chains");

                                //    var targetChain = targetChains.First();
                                //    targetChain.Add(newPair);


                                //} else {
                                //    var firstPair = m_AggPairsWithExtNeighborPairs.Where(p => p.jCellTarget == aggConnectionEdge.targetCell && p.jCellSource == aggConnectionEdge.jCellNeigh).First();
                                //    var newChain = new AgglomerationChain(firstPair, grdDat, CellVolumes);


                                //    newChain.Add(newPair);
                                //    this.m_AgglomerationChains.Add(newChain);

                                //}

                                CellsNeedChainAgglomeration.Remove(aggConnectionEdge.jCell);
                                anyUpdate = true;
                            }
                        }
                    }

                    #region update lists and variables for the while loop
                    // Exchange the new pairs with neighbor cells and their processors
                    if (InterProcessChainAgglomeration > 0) {
                        // add new neighbor pairs to external pairs lists 
                        InterProcessChainAgglomeration = DoAggPairsMPIexchangeForGhostCells(LoopChainAgglomerationPairs, ref AggPairsOnExtNeighborPairs);
                        anyUpdate = anyUpdate.MPIOr();
                        //AggPairsOnExtNeighborPairs.SaveToTextFileDebugUnsteady("e_AggPairsOnExtNeighborPairs", ".txt");
                    }

                    // Add loop chain aggs to the lists
                    ChainAgglomerationPairs.AddRange(LoopChainAgglomerationPairs);
                    m_AggPairs.AddRange(LoopChainAgglomerationPairs);

                    ii++;

                    ChainCountMax = CellsNeedChainAgglomeration.Count;
                    ChainCountMax = ChainCountMax.MPIMax();
                    #endregion
                }
                #endregion

                if (ChainAgglomerationPairs.Count() > 0)
                    //m_AgglomerationChains.SaveToTextFileDebugUnsteady("aggChains_" + Tag, ".txt");

                    // If there is still cells waiting for agglomeration, this means that agg. failed
                    #region AgglomerationKatastrophe
                    if (CellsNeedChainAgglomeration.Count > 0) {
                        Console.WriteLine($"## Chain Agglomeration is failed on proc-{ilPSP.Environment.MPIEnv.MPI_Rank} ##");
                    }

                // Save the data for debugging purposes
                if (ChainCountMax > 0) {
                    ChainAgglomerationPairs.SaveToTextFileDebugUnsteady("agg_ChainAgglomerationPairs_" + Tag, ".txt");
                    m_AggPairs.SaveToTextFileDebugUnsteady("agg_AggPairs_" + Tag, ".txt");
                    AggPairsOnExtNeighborPairs.SaveToTextFileDebugUnsteady("agg_AggPairsOnExtNeighborPairs_ " + Tag, ".txt");
                }

                failCells.AddRange(CellsNeedChainAgglomeration);
                int NoFailedCells = failCells.Count.MPISum();
                if (NoFailedCells > 0 || PlotAgglomeration) {
                    int[] pairIdentification = new int[Jup];
                    int[] pairColor = new int[Jup];
                    Vector[] aggDirection = new Vector[Jup];
                    double[] volFrac = new double[Jup];
                    int Dim = grdDat.Cells.GetCenter(0).Dim;


                    for (int j = 0; j < Jup; j++) {
                        double totVol = grdDat.Cells.GetCellVolume(j);
                        double spcVol = CellVolumes[j];
                        volFrac[j] = spcVol / totVol;
                        aggDirection[j] = new Vector(Dim);
                    }

                    int k = (int)CellPart.i0 + 1;
                    foreach (var pair in AgglomerationPairs) {
                        pairIdentification[pair.jCellSource] = pair.jCellTarget;
                        pairColor[pair.jCellSource] = k;
                        if (pair.jCellTarget < Jup) {
                            pairColor[pair.jCellTarget] = k; //can override when it is chain
                            Vector direction = new Vector(Dim);
                            direction = grdDat.Cells.GetCenter(pair.jCellTarget) - grdDat.Cells.GetCenter(pair.jCellSource);
                            aggDirection[pair.jCellSource] = direction;
                        } else {
                            int jAggTarget = (int)grdDat.Parallel.GetGlobalCellIndex(pair.jCellTarget);
                            pairIdentification[pair.jCellSource] = -jAggTarget;
                            Vector direction = new Vector(Dim);
                            direction = grdDat.Cells.GetCenter(pair.jCellTarget) - grdDat.Cells.GetCenter(pair.jCellSource);
                            aggDirection[pair.jCellSource] = direction;
                        }

                        k++;
                    }

                    int[] AggTargets = AgglomerationPairs.Select(p => p.jCellTarget).Distinct().ToArray();

                    if (NoFailedCells > 0)
                        PlotFail(false, CellVolumes, oldCellVolumes, AgglomSourceCellsList, m_failCells, AggCandidates, AggTargets, pairIdentification, pairColor, aggDirection, volFrac, $"{Tag}AgglomerationKatastrophe{spId.ToString()}");


                    if (PlotAgglomeration)
                        PlotFail(false, CellVolumes, oldCellVolumes, AgglomSourceCellsList, m_failCells, AggCandidates, AggTargets, pairIdentification, pairColor, aggDirection, volFrac, $"{Tag}AgglomerationGraphOf{spId.ToString()}");

                }
                #endregion



                //int[] intAggCandidates = new int[Jup];
                //for (int j=0; j < Jup; j++) {
                //    if (AggCandidates[j])
                //        intAggCandidates[j] = 1;
                //}
                //intAggCandidates.SaveToTextFileDebugUnsteadyNumbered("intAggCandidates",".txt");

                // store & return;
                // ================
                this.AgglomerationPairs = AgglomerationPairs.Select(pair => (pair.jCellSource, pair.jCellTarget)).ToArray();
            }
        }

        /// <summary>
        /// 3rd pass of agglomeration algorithm, Identification of agglomeration targets
        /// </summary>
        /// <remarks>
        /// Revised algorithm for chains, in use since Sept. 2023
        /// </remarks>
        protected virtual void FindAgglomerationTargets_Mk3(
            List<int> AgglomSourceCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates, List<int>NewbornCellsList, List<int> VanishingCellsList
            ) {
            using (new FuncTrace()) {
                int myMpiRank = Tracker.GridDat.MpiRank;
                Partitioning CellPart = Tracker.GridDat.CellPartitioning;
                int Jup = grdDat.Cells.NoOfLocalUpdatedCells;

                // The incoming AggCandidates represent newborn or vanishing cells
                List<int> CellsRequireTopologyChanges = NewbornCellsList.Union(VanishingCellsList).ToList();

                // Initiate new agglomeration pair list
                var AgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();
                m_AggPairs = AgglomerationPairs;

                //exchange the source cell information to be able to know about external/ghost cells (needed only once)
                BitArray AggSourcesWithExternalCell = new BitArray(grdDat.iLogicalCells.Count, false);

                foreach (var jCell in AgglomSourceCellsList)
                    AggSourcesWithExternalCell[jCell] = true;

                AggSourcesWithExternalCell.MPIExchange(grdDat);

                //make sure source cells cannot be target cells
                for (int j = 0; j < Jup; j++) {
                    if (AggSourcesWithExternalCell[j] == true) {
                        AggCandidates[j] = false;
                    }
                }

                //exchange the candidate information to be able to know about external/ghost cells (needed only once)
                AggCandidates.MPIExchange(grdDat);

                // pass 2: determine agglomeration targets
                // ---------------------------------------
                m_failCells = new List<int>();
                m_CellsNeedChainAgglomeration = new List<int>();

                // First check if a target can be found among immediate neighbors
                SearchDirectAgglomeration_Mk3(AgglomSourceCellsList, AggCandidates, AggSourcesWithExternalCell);

                // If there are remaining cells, try to create agg. chains
                m_AggPairsOnNeighborPairs = new List<CellAgglomerator.AgglomerationPair>();
                SearchChainAgglomeration_Mk3(m_CellsNeedChainAgglomeration, AggCandidates, AggSourcesWithExternalCell);

                // If there are still cells waiting for agglomeration, this means that we can perhaps create agg groups with source cells
                m_aggGroups = new List<AgglomerationGroup>();
                SearchGroupAgglomeration_Mk3(m_CellsNeedChainAgglomeration, AggCandidates, AggSourcesWithExternalCell, CellsRequireTopologyChanges);

                // In case of still failed cases
                #region AgglomerationKatastrophe               
                var NoFailedCells = m_failCells.Count.MPISum();
                if (NoFailedCells > 0 || PlotAgglomeration) {
                    int[] pairIdentification = new int[Jup];
                    int[] pairColor = new int[Jup];
                    Vector[] aggDirection = new Vector[Jup];
                    double[] volFrac = new double[Jup];
                    int Dim = grdDat.Cells.GetCenter(0).Dim;


                    for (int j = 0; j < Jup; j++) {
                        double totVol = grdDat.Cells.GetCellVolume(j);
                        double spcVol = CellVolumes[j];
                        volFrac[j] = spcVol / totVol;
                        aggDirection[j] = new Vector(Dim);
                    }

                    int k = (int)CellPart.i0 + 1;
                    foreach (var pair in AgglomerationPairs) {
                        pairIdentification[pair.jCellSource] = pair.jCellTarget;
                        pairColor[pair.jCellSource] = k;
                        if (pair.jCellTarget < Jup) {
                            pairColor[pair.jCellTarget] = k; //can override when it is chain
                            Vector direction = new Vector(Dim);
                            direction = grdDat.Cells.GetCenter(pair.jCellTarget) - grdDat.Cells.GetCenter(pair.jCellSource);
                            aggDirection[pair.jCellSource] = direction;
                        } else {
                            int jAggTarget = (int)grdDat.Parallel.GetGlobalCellIndex(pair.jCellTarget);
                            pairIdentification[pair.jCellSource] = -jAggTarget;
                            Vector direction = new Vector(Dim);
                            direction = grdDat.Cells.GetCenter(pair.jCellTarget) - grdDat.Cells.GetCenter(pair.jCellSource);
                            aggDirection[pair.jCellSource] = direction;
                        }

                        k++;
                    }

                    int[] AggTargets = AgglomerationPairs.Select(p => p.jCellTarget).Distinct().ToArray();

                    if (NoFailedCells > 0)
                        PlotFail(false, CellVolumes, oldCellVolumes, AgglomSourceCellsList, m_failCells, AggCandidates, AggTargets, pairIdentification, pairColor, aggDirection, volFrac, $"{Tag}AgglomerationKatastrophe{spId.ToString()}");


                    if (PlotAgglomeration)
                        PlotFail(false, CellVolumes, oldCellVolumes, AgglomSourceCellsList, m_failCells, AggCandidates, AggTargets, pairIdentification, pairColor, aggDirection, volFrac, $"{Tag}AgglomerationGraphOf{spId.ToString()}");

                }
                #endregion


                // store & return;
                // ================
                //this.AgglomerationPairs = AgglomerationPairs.Select(pair => (pair.jCellSource, pair.jCellTarget)).ToArray();
                this.AgglomerationPairsWithRanks = AgglomerationPairs.ToArray();
                AgglomerationPairsWithRanks.SaveToTextFileDebugUnsteadyNumbered("aggPairs",".txt");
                //DoAggPairsMPIexchangeForGhostCells(AgglomerationPairs, ref m_AggPairsOnNeighborPairs);
                // this.AgglomerationPairsWithRanks.m_AggPairsWithExtNeighborPairs.Where(p => p.OwnerRank4Source == myMpiRank || p.OwnerRank4Target == myMpiRank).ToArray(););
            }
        }

        /// <summary>
        /// searches chaing agglomeration targets
        /// </summary>
        /// <param name="CellsNeedChainAgglomeration"></param>
        /// <param name="AggCandidates"></param>
        /// <param name="AggSourcesWithExternalCell"></param>
        /// <exception cref="ArgumentException"></exception>
        protected virtual void SearchChainAgglomeration_Mk3(List<int> CellsNeedChainAgglomeration, BitArray AggCandidates,
            BitArray AggSourcesWithExternalCell
            ) {
            var Cell2Edge = grdDat.Cells.Cells2Edges;
            int[,] Edge2Cell = grdDat.Edges.CellIndices;
            int NoOfEdges = grdDat.Edges.Count;
            byte[] EdgeTags = grdDat.Edges.EdgeTags;
            int myMpiRank = Tracker.GridDat.MpiRank;
            int Jup = grdDat.Cells.NoOfLocalUpdatedCells;

            Partitioning CellPart = Tracker.GridDat.CellPartitioning;
            var GidxExt = Tracker.GridDat.Parallel.GlobalIndicesExternalCells;

            // MPI Watch
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            // Get already assigned pairs
            var AgglomerationPairs = m_AggPairs;

            // in certain cases, agglomeration sources can require to create a chain/cluster (a set of cut cells that should be agglomerated to a target) 
            var ChainAgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

            // Exchange already paired cells (needs to be updated in every change in the loop below)
            var AggPairsOnExtNeighborPairs = m_AggPairsOnNeighborPairs;
            var InterProcessChainAgglomeration = DoAggPairsMPIexchangeForGhostCells(AgglomerationPairs, ref AggPairsOnExtNeighborPairs);

            int ChainCountMax = CellsNeedChainAgglomeration.Count;
            ChainCountMax = ChainCountMax.MPIMax();

            bool anyUpdate = true;

            // Global-MPI level
            while (ChainCountMax > 0 && anyUpdate) { // as long as a cell needs to be agglomerated and a possibility
                anyUpdate = false;
                var LoopChainAgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

                List<(double SpeciesFrac, double Dist, int AggLevel, int jCell, int jCellNeigh, int targetCell, int targetRank)> weightedEdges = new List<(double, double, int, int, int, int, int)>();
                List<int> ImmediateConnectionCells = new List<int>(CellsNeedChainAgglomeration.Count);
                bool isThisProcAuthorized = true;

                // Any cell to be agglomerated in a chain on the local proc (a parallel variation of Kruskal's algorithm for Minimum Spanning Forest adapted to our context
                if (CellsNeedChainAgglomeration.Count > 0) {

                    //create edge information for possible connections, remind that we could only create a chain if there is a neighbor that is already paired with some cell.  
                    foreach (int jCell in CellsNeedChainAgglomeration) {
                        var Cell2Edge_jCell = Cell2Edge[jCell];
                        int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;

                        List<(double SpeciesFrac, double Dist, int AggLevel, int jCell, int jCellNeigh, int targetCell, int targetRank)> weightedEdgesjCell = new List<(double, double, int, int, int, int, int)>(NoOfEdges_4_jCell);

                        // Collect neighbors and determine if there is a possible source cell through which 'jCell' can be connected to a proper target
                        int NeighborAggSourceCells = 0;
                        int PairedNeighborAggSourceCells = 0;
                        for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                            bool IsPossibleTarget = false;
                            var (iEdge, ThisCell, OtherCell) = GetEdgeInfo(Cell2Edge_jCell[e]);
                            int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                            double EdgeArea_iEdge = edgeArea[iEdge];
                            Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                            if (jCellNeigh >= 0 && (EdgeTags[iEdge] < GridCommons.FIRST_PERIODIC_BC_TAG || CellVolumes[jCell] <= 0)) {
                                IsPossibleTarget = true;

                                if (AggSourcesWithExternalCell[jCellNeigh])
                                    NeighborAggSourceCells++;

                            } else {
                                continue;
                            }

                            // assume maximum distance as default value for cases that there is no additional info
                            double Distance = double.MaxValue;

                            // Searching for a pair is needed to ensure the neighbor cell has a "real" target, which will be connected in CellAgglomerator.cs later
                            var TargetPairs = m_AggPairsWithExtNeighborPairs.Where(p => p.jCellSource == jCellNeigh);

                            if (TargetPairs.Count() > 1) {
                                throw new ArgumentException("Agglomeration fail: a cell is mapped to two different target cells");
                            }

                            // only already paired cells can form a target for the chain
                            if (TargetPairs.Any() && IsPossibleTarget) { //check if a neighbor is agglomerated
                                var TargetPair = TargetPairs.First();
                                int TargetCellOfTargetPair = TargetPair.jCellTarget; //looking for the final target
                                int possibleTarget; // due to the mpi boundaries we may select the neighbor as target instead of the final target
                                int targetRank = -1; //by default assigned an invalid value to ensure code works

                                // To turn off the direct agg to the final target, assign the below to -1
                                // TargetCellOfTargetPair = -1; // normally the code is able to choose a target that forms higher level agg. chains. However, it is decided to be handled in CellAgglomerator.cs

                                // if TargetCellOfTargetPair is one of ext/ghost cells or already on local proc
                                bool IsTargetCellKnown = TargetCellOfTargetPair >= 0 && TargetPair.OwnerRank4Target >= 0;

                                // We must also ensure that source cell is known to the owner of the target as well (otherwise it will lead to local indexing errors)
                                bool IsSourceKnownToOwnerOfTargetCell = grdDat.Cells.CellNeighbours[jCell].Contains(TargetCellOfTargetPair) || TargetPair.OwnerRank4Target == myMpiRank;
                                int AggLevel;
                                if (IsTargetCellKnown && IsSourceKnownToOwnerOfTargetCell) {
                                    possibleTarget = TargetCellOfTargetPair; //if so, directly pair jCell with the final target
                                    targetRank = TargetPair.OwnerRank4Target;
                                    AggLevel = TargetPair.AgglomerationLevel;
                                } else {
                                    possibleTarget = jCellNeigh;  //if not, pair jCell with jNeigh (which is also TargetPair.jCellSource). This cell will behave like a "carrying" cell which ultimately would lead to the target cell
                                    targetRank = TargetPair.OwnerRank4Source;
                                    AggLevel = TargetPair.AgglomerationLevel;
                                    Debug.Assert(jCellNeigh == TargetPair.jCellSource);
                                }

                                // Create the weighting info (using the real target with lvl=0)
                                double SpeciesFrac = TargetPair.fracTarget;
                                Vector posTarget = TargetPair.posTarget;

                                Vector posSource = grdDat.Cells.GetCenter(jCell);
                                Distance = Vector.Dist(posTarget, posSource);
                                weightedEdgesjCell.Add((SpeciesFrac, Distance, AggLevel, jCell, jCellNeigh, possibleTarget, targetRank));
                                PairedNeighborAggSourceCells++;
                            }
                        }

                        bool IsAllPossibleNeighborCellsPaired = (PairedNeighborAggSourceCells == NeighborAggSourceCells) && (PairedNeighborAggSourceCells > 0);

                        // If all possible neighbors are already paired, we can choose one of them directly
                        if (IsAllPossibleNeighborCellsPaired) {
                            var ConnectionEdge = weightedEdgesjCell
                                .Where(p => p.jCell == jCell)
                                .OrderBy(p => p.Dist)
                                .ThenByDescending(p => p.SpeciesFrac)
                                .ToList()
                                .First();

                            if (ConnectionEdge.targetRank >= 0) {
                                //Get the target info from the target pair
                                var TargetPair = m_AggPairsWithExtNeighborPairs.Where(p => p.jCellSource == ConnectionEdge.jCellNeigh).First();
                                double SpeciesFrac = TargetPair.fracTarget;
                                Vector posTarget = TargetPair.posTarget;

                                //Add new pair
                                var newPair = new CellAgglomerator.AgglomerationPair() {
                                    jCellTarget = ConnectionEdge.targetCell,
                                    posTarget = posTarget,
                                    fracTarget = SpeciesFrac,
                                    //jNeighborSource = ConnectionEdge.jCellNeigh,
                                    jCellSource = ConnectionEdge.jCell,
                                    OwnerRank4Target = ConnectionEdge.targetRank,
                                    OwnerRank4Source = myMpiRank,
                                    AgglomerationLevel = ConnectionEdge.AggLevel
                                };

                                LoopChainAgglomerationPairs.Add(newPair);
                                ImmediateConnectionCells.Add(ConnectionEdge.jCell);

                                continue;
                            }
                        } else {
                            weightedEdges.AddRange(weightedEdgesjCell);
                        }

                    }

                    // discard already connected cells
                    foreach (int DirectConnected in ImmediateConnectionCells) {
                        CellsNeedChainAgglomeration.Remove(DirectConnected);
                        anyUpdate = true;
                    }

                    // sort remaining edges
                    weightedEdges = weightedEdges.OrderBy(p => p.Dist).ThenByDescending(p => p.SpeciesFrac).ToList(); // (ascending, descending)

                }

                // check if this proc has the globally best edge (only if required)
                if (m_globalAgglomerationMapping && InterProcessChainAgglomeration > 0)
                {
                    var bestEdge = weightedEdges.FirstOrDefault();

                    //convert the local number to the global and then create a new tuple with it
                    var globalCellNumberForBestEdge = grdDat.CellPartitioning.i0 + (long)bestEdge.jCell;
                    var bestEdgeWithGlobalIndex = new { Dist = bestEdge.Dist, SpeciesFrac = bestEdge.SpeciesFrac, GlobalCellNumber = globalCellNumberForBestEdge };
                    var allBestEdges = bestEdgeWithGlobalIndex.MPIAllGatherO();
                    var bestGlobalEdge = allBestEdges.OrderBy(p => p.Dist).ThenByDescending(p => p.SpeciesFrac).ThenBy(p => p.GlobalCellNumber).First();
                    //isThisProcAuthorized = bestGlobalEdge.GlobalCellNumber == globalCellNumberForBestEdge ? true : false;
                }
                //weightedEdges.SaveToTextFileDebugUnsteady("weightedEdges", ".txt");

                // if any target available in case a target needed?
                Debug.Assert((weightedEdges.Any() || ImmediateConnectionCells.Any()) && !CellsNeedChainAgglomeration.Any(), "Cell agglomeration failed." +
                        " There are cells that cannot be connected any target cells. (Cycle between cells to be agglomerated");

                // choose the first edge and add the corresponding agg. pair
                if (weightedEdges.Any() && CellsNeedChainAgglomeration.Any() && isThisProcAuthorized)
                {
                    var aggConnectionEdge = weightedEdges.First();


                    // AddToList
                    if (aggConnectionEdge.targetRank > -1)
                    {
                        //Get the target info from the target pair
                        var TargetPair = m_AggPairsWithExtNeighborPairs.Where(p => p.jCellSource == aggConnectionEdge.jCellNeigh).First();
                        double SpeciesFrac = TargetPair.fracTarget;
                        Vector posTarget = TargetPair.posTarget;

                        //Add new pair
                        //var SourceCell2Edge_jCell = Cell2Edge[aggConnectionEdge.jCell];
                        var newPair = new CellAgglomerator.AgglomerationPair()
                        {
                            jCellTarget = aggConnectionEdge.targetCell,
                            posTarget = posTarget,
                            fracTarget = SpeciesFrac,
                            jCellSource = aggConnectionEdge.jCell,
                            OwnerRank4Target = aggConnectionEdge.targetRank,
                            OwnerRank4Source = myMpiRank,
                            AgglomerationLevel = aggConnectionEdge.AggLevel

                        };

                        LoopChainAgglomerationPairs.Add(newPair);
                        CellsNeedChainAgglomeration.Remove(aggConnectionEdge.jCell);
                        anyUpdate = true;
                    }
                }

                #region update lists and variables for the while loop
                ChainCountMax = CellsNeedChainAgglomeration.Count;

                // Exchange the new pairs with neighbor cells and their processors
                if (InterProcessChainAgglomeration > 0) {
                    // add new neighbor pairs to external pairs lists 
                    InterProcessChainAgglomeration = DoAggPairsMPIexchangeForGhostCells(LoopChainAgglomerationPairs, ref AggPairsOnExtNeighborPairs);
                    anyUpdate = anyUpdate.MPIOr();
                    ChainCountMax = ChainCountMax.MPIMax();
                    //AggPairsOnExtNeighborPairs.SaveToTextFileDebugUnsteady("e_AggPairsOnExtNeighborPairs", ".txt");
                }

                // Add loop chain aggs to the lists
                ChainAgglomerationPairs.AddRange(LoopChainAgglomerationPairs);
                m_AggPairs.AddRange(LoopChainAgglomerationPairs);
                #endregion
            }

            if (ChainAgglomerationPairs.Any() && PlotAgglomeration) //for debugging purposes
                ChainAgglomerationPairs.SaveToTextFileDebugUnsteady(marker  + "ChainAgglomerationPairs" + Tag + spId.ToString(), ".txt");

        }

        /// <summary>
        /// searches direct agglomeration targets
        /// </summary>
        /// <param name="AgglomSourceCellsList"></param>
        /// <param name="AggCandidates"></param>
        /// <param name="AggSourcesWithExternalCell"></param>
        /// <exception cref="ArgumentException"></exception>
        protected virtual void SearchDirectAgglomeration_Mk3(
            List<int> AgglomSourceCellsList, BitArray AggCandidates, BitArray AggSourcesWithExternalCell
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

                var AgglomerationPairs = m_AggPairs;

                if (edgeArea.GetLength(0) != NoOfEdges)
                    throw new ArgumentException();

                // pass 2: determine agglomeration targets
                // ---------------------------------------
                var failCells = m_failCells;
                var CellsNeedChainAgglomeration = m_CellsNeedChainAgglomeration;

                // First check if a target can be found among neighbors
                #region DirectAgg
                // first check if source cells can find a "proper" target (excluding vanishing and newborn cells)
                foreach (int jCell in AgglomSourceCellsList) {
                    var Cell2Edge_jCell = Cell2Edge[jCell];

                    double frac_neigh_max = -1.0;
                    int jCellNeigh_max = int.MinValue;

                    int NoOfEdges_4_jCell = Cell2Edge_jCell.Length;

                    // create an array for neighbors
                    double[,] neighbors = new double[NoOfEdges_4_jCell, 3]; // (jCellNeigh, iEdge, edgeArea) 

                    // Collect neighbors and determine if there is a non-empty edge which connects cell 'jCell' to some other cell
                    bool NonEmptyEdgeAvailable = false;
                    for (int e = 0; e < NoOfEdges_4_jCell; e++) { // loop over faces/neighbour cells...
                        var (iEdge, ThisCell, OtherCell) = GetEdgeInfo(Cell2Edge_jCell[e]);
                        int jCellNeigh = Edge2Cell[iEdge, OtherCell];
                        double EdgeArea_iEdge = edgeArea[iEdge];
                        Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell);

                        neighbors[e, 0] = (double)jCellNeigh;
                        neighbors[e, 1] = (double)iEdge;
                        neighbors[e, 2] = EdgeArea_iEdge;

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
                        int jCellNeigh = (int)neighbors[e, 0];
                        int iEdge = (int)neighbors[e, 1];
                        double EdgeArea_iEdge = neighbors[e, 2];

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

                            //Debug.Assert(Edge2Cell[iEdge, ThisCell] == jCell, "sollte aber so sein");
                            continue;
                        }

                        if (!AggCandidates[jCellNeigh]) {
                            continue; // not suitable for agglomeration
                        }

                        // volume fraction of neighbour cell
                        double spcVol_neigh = CellVolumes[jCellNeigh];
                        double totVol_neigh = grdDat.Cells.GetCellVolume(jCellNeigh);
                        double frac_neigh = spcVol_neigh / totVol_neigh;

                        // max?
                        if (frac_neigh > frac_neigh_max) {
                            frac_neigh_max = frac_neigh;
                            jCellNeigh_max = jCellNeigh;
                        }
                    }


                    {
                        if (jCellNeigh_max < 0) { //is(not) it a proper target?
                            bool isAggChainPossible = false;

                            // check neighbors if they can be both agg. source
                            for (int e = 0; e < NoOfEdges_4_jCell; e++) {
                                int jCellNeigh = (int)neighbors[e, 0];

                                //Exclude boundary cells
                                if (jCellNeigh < 0)
                                    continue;

                                if (AggSourcesWithExternalCell[jCellNeigh]) {
                                    isAggChainPossible = true; // there is a possibility to form a chain, which means that the neighbor cell could be carrying element to the final target)
                                    break;
                                }
                            }

                            if (isAggChainPossible) {
                                CellsNeedChainAgglomeration.Add(jCell);
                            } else {
                                failCells.Add(jCell); // jcell has no possible targets to be agglomerated (it is a single isolated cell)
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
                                posTarget = grdDat.Cells.GetCenter(jCellNeigh_max),
                                fracTarget = frac_neigh_max,
                                OwnerRank4Source = myMpiRank,
                                AgglomerationLevel = 0
                            });
                        }
                    }
                }
                #endregion

                if (AgglomerationPairs.Any() && PlotAgglomeration) //for debugging purposes
                    AgglomerationPairs.SaveToTextFileDebugUnsteady("DirectAgglomerationPairs" + Tag +  spId.ToString(), ".txt");

            }
        }

        /// <summary>
        /// creates agglomeration groups from cut cells isolated from uncut cells, with the biggests being targets
        /// </summary>
        /// <param name="CellsNeedChainAgglomeration"></param>
        /// <param name="AggCandidates"></param>
        /// <param name="AggSourcesWithExternalCell"></param>
        /// <param name="CellsRequireTopologyChanges"></param>
        protected virtual void SearchGroupAgglomeration_Mk3(List<int> CellsNeedChainAgglomeration, BitArray AggCandidates, BitArray AggSourcesWithExternalCell, List<int> CellsRequireTopologyChanges) {
            var possibleGroupTargets = m_CellsNeedChainAgglomeration.Except(CellsRequireTopologyChanges).ToList(); //cells requiring topology change cannot be a target at any point
            var AgglomerationPairs = m_AggPairs; // Get already assigned pairs

            marker = "grp";

            if (m_globalAgglomerationMapping)
            {
                //create global groups (costly as the remaining cells are usually close to each other but results in unique agg. mapping)
                while (possibleGroupTargets.Count.MPIMax() > 0) { //these cells have at least one adjacent source cell, in contrast to failedCells, which do not have any neighbor to form agglomeration. So, let's group them.

                    //choose the local biggest cell
                    var cellWithMaxVolume = possibleGroupTargets
                            .Select(p => (CellNumber: p, CellVolume: CellVolumes[p]))
                            .OrderByDescending(p => p.CellVolume)
                            .FirstOrDefault(); //due to MPIAllGatherO trivial or non-trivial every processor needs to send data

                    int localCellNumberWithMaxVolume = cellWithMaxVolume.CellNumber;

                    //convert the local number to the global and then create a new tuple with it
                    var globalCellNumberWithMaxVolume = grdDat.CellPartitioning.i0 + (long)localCellNumberWithMaxVolume;
                    var cellWithMaxVolumeGlobal = new { GlobalCellNumber = globalCellNumberWithMaxVolume, CellVolume = cellWithMaxVolume.CellVolume };

                    //Gather all the biggest cells from each processor and choose the global biggest (cells can have the same volume so use global cell index as a second criterion)
                    var cellWithMaxVolumeArray = cellWithMaxVolumeGlobal.MPIAllGatherO();
                    var cellWithGlobalMaxVolume = cellWithMaxVolumeArray.OrderByDescending(c => c.CellVolume).ThenBy(c => c.GlobalCellNumber).First();

                    //check if this processor has the biggest cell
                    if (cellWithGlobalMaxVolume.GlobalCellNumber == globalCellNumberWithMaxVolume && possibleGroupTargets.Count != 0) {

                        //mark the cell as candidate
                        var aggGroupLoop = new AgglomerationGroup(localCellNumberWithMaxVolume, Tracker, CellVolumes, edgeArea);
                        m_aggGroups.Add(aggGroupLoop);

                        m_CellsNeedChainAgglomeration.Remove(localCellNumberWithMaxVolume);
                        bool anyUpdates = true;

                        while (anyUpdates)
                        {
                            var intersection = aggGroupLoop.GetNeighbors.Intersect(m_CellsNeedChainAgglomeration);
                            anyUpdates = intersection.Any();
                            foreach (int cell in intersection)
                            {
                                aggGroupLoop.Add(cell);
                                m_CellsNeedChainAgglomeration.Remove(cell);
                            }
                        }

                        // add groups into the pairs list
                        AgglomerationPairs.AddRange(aggGroupLoop.GetAggPairs);
                        //if (aggGroupLoop.SumFractions < AgglomerationThreshold) {  //mark unsuccessful groups as failure and add them the fail list
                        //    m_failCells.AddRange(aggGroupLoop.GetCells);
                        //}
                    } 

                    //re-run chain agglomerations to attach source cells from other processors
                    SearchChainAgglomeration_Mk3(m_CellsNeedChainAgglomeration, AggCandidates, AggSourcesWithExternalCell);

                    possibleGroupTargets = m_CellsNeedChainAgglomeration.Except(CellsRequireTopologyChanges).ToList();
                }


            } else {
                // create only local groups
                while (possibleGroupTargets.Count > 0)
                { //these cells have at least one adjacent source cell, in contrast to failedCells, which do not have any neighbor to form agglomeration. So, let's group them.
                    int cellNumberWithMaxVolume = possibleGroupTargets
                            .Select(p => new { CellNumber = p, CellVolume = CellVolumes[p] })
                            .OrderByDescending(p => p.CellVolume)
                            .First()
                            .CellNumber;

                    var aggGroupLoop = new AgglomerationGroup(cellNumberWithMaxVolume, Tracker, CellVolumes, edgeArea);
                    m_aggGroups.Add(aggGroupLoop);

                    m_CellsNeedChainAgglomeration.Remove(cellNumberWithMaxVolume);
                    bool anyUpdates = true;

                    while (anyUpdates)
                    {
                        var intersection = aggGroupLoop.GetNeighbors.Intersect(m_CellsNeedChainAgglomeration);
                        anyUpdates = intersection.Any();
                        foreach (int cell in intersection)
                        {
                            aggGroupLoop.Add(cell);
                            m_CellsNeedChainAgglomeration.Remove(cell);
                        }
                    }
                    possibleGroupTargets = m_CellsNeedChainAgglomeration.Except(CellsRequireTopologyChanges).ToList();
                }

                // add groups into the pairs list
                foreach (var group in m_aggGroups)
                {
                    AgglomerationPairs.AddRange(group.GetAggPairs);
                    //if (group.SumFractions < AgglomerationThreshold)
                    //{  //mark unsuccessful groups as failure and add them the fail list
                    //    m_failCells.AddRange(group.GetCells);
                    //}
                }

            }

            // So far if a cell that requires topology change cannot be mapped, this means a topologicalFailure
            m_failCells.AddRange(m_CellsNeedChainAgglomeration);
            var topologicalFailures = m_failCells.Intersect(CellsRequireTopologyChanges);
            if (topologicalFailures.Any())
                topologicalFailures.SaveToTextFileDebugUnsteady("failedTopologicalCells", ".txt");

            // Debugging info
            if (m_aggGroups.Any() && PlotAgglomeration)
            {
                m_aggGroups.SaveToTextFileDebugUnsteady("m_aggGroups" + Tag + spId.ToString(), ".txt");

            }
            AgglomerationPairs.RemoveAll(p => p.jCellSource == p.jCellTarget);   //remove self mapped elements due to the agglomeration groups
        }

        /// <summary>
        /// plot failure as tecplot
        /// </summary>
        /// <param name="ExceptionOnFailedAgglomeration"></param>
        /// <param name="CellVolumes"></param>
        /// <param name="oldCellVolumes"></param>
        /// <param name="AgglomCellsList"></param>
        /// <param name="failCells"></param>
        /// <param name="aggCandidates"></param>
        /// <param name="aggTargets"></param>
        /// <param name="pairIdentification"></param>
        /// <param name="pairColor"></param>
        /// <param name="aggDirection"></param>
        /// <param name="volFrac"></param>
        /// <param name="Tag"></param>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        /// <exception cref="Exception"></exception>
        private void PlotFail(bool ExceptionOnFailedAgglomeration, MultidimensionalArray CellVolumes, MultidimensionalArray[] oldCellVolumes, List<int> AgglomCellsList, List<int> failCells, BitArray aggCandidates, int[] aggTargets, int[] pairIdentification, int[] pairColor, Vector[] aggDirection, double[] volFrac, string Tag = "") {
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
                FailCells.SaveToTextFile($"aggFailCells{grdDat.MpiRank}.csv");
            }


            DGField CellNumbers = new SinglePhaseField(b, "CellNumbers");
            CellNumbers.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {
                        result[j, k] = j0 + j;
                    }
                }
            }, new CellQuadratureScheme());


            DGField CellNumbersGlob = new SinglePhaseField(b, "CellNumbersGlob");
            CellNumbersGlob.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                for (int j = 0; j < Len; j++) {
                    long gidx = grdDat.CellPartitioning.i0 + j0 + j;
                    for (int k = 0; k < K; k++) {
                        result[j, k] = gidx;
                    }
                }
            }, new CellQuadratureScheme());


            DGField GlobalID = new SinglePhaseField(b, "GlobalID");
            GlobalID.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                for (int j = 0; j < Len; j++) {
                    long gid = this.grdDat.iLogicalCells.GetGlobalID(j + j0); ;
                    for (int k = 0; k < K; k++) {
                        result[j, k] = gid;
                    }
                }
            }, new CellQuadratureScheme());
            int my_rank = ilPSP.Environment.MPIEnv.MPI_Rank;

            DGField MPIRank = new SinglePhaseField(b, "MPIRank");
            MPIRank.ProjectField(1.0, delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {
                        result[j, k] = my_rank;
                    }
                }
            }, new CellQuadratureScheme());

            DGField[] LevelSets = Tracker.LevelSets.Select(s => (DGField)s).ToArray();

            DGField LvSetDist = new SinglePhaseField(b, "LvSetDist");
            LvSetDist.AccLevelSetDist(1, Tracker, 1);

            DGField Identity = new SinglePhaseField(b, "Identity");
            Identity.Clear();
            Identity.AccConstant(1);

            DGField AggCandidates = new SinglePhaseField(b, "aggCandidate");
            AggCandidates.Clear();
            for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                if (aggCandidates[j]) AggCandidates.SetMeanValue(j, 1);
            }

            DGField AggTargets = new SinglePhaseField(b, "aggTarget");
            AggTargets.Clear();
            foreach (int jCell in aggTargets) {
                if(jCell < grdDat.Cells.NoOfLocalUpdatedCells) AggTargets.SetMeanValue(jCell, 1);
            }

            DGField PairIdentification = new SinglePhaseField(b, "pairTargetNo");
            PairIdentification.Clear();
            for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                PairIdentification.SetMeanValue(j, pairIdentification[j]);
            }

            DGField PairColor = new SinglePhaseField(b, "pairColor");
            PairColor.Clear();
            for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                PairColor.SetMeanValue(j, pairColor[j]);
            }

            if (aggDirection[0].Dim > 0) {
                DGField AggDirX = new SinglePhaseField(b, "aggDirX");
                AggDirX.Clear();
                for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                    AggDirX.SetMeanValue(j, aggDirection[j][0]);
                }
                CellVolumesViz = CellVolumesViz.Cat(AggDirX);
            } 
            if (aggDirection[0].Dim > 1) {
                DGField AggDirY = new SinglePhaseField(b, "aggDirY");
                AggDirY.Clear();
                for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                    AggDirY.SetMeanValue(j, aggDirection[j][1]);
                }
                CellVolumesViz = CellVolumesViz.Cat(AggDirY);
            } 
            if (aggDirection[0].Dim > 2) {
                DGField AggDirZ = new SinglePhaseField(b, "aggDirZ");
                AggDirZ.Clear();
                for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                    AggDirZ.SetMeanValue(j, aggDirection[j][2]);
                }
                CellVolumesViz = CellVolumesViz.Cat(AggDirZ);
            }
            if (aggDirection[0].Dim > 3 || aggDirection[0].Dim < 1) {
                    throw new ArgumentOutOfRangeException("Agglomeration fail plot can only support 1D,2D and 3D geometries");
            }
          
            DGField VolFrac = new SinglePhaseField(b, "volFrac");
            VolFrac.Clear();
            for (int j = 0; j < grdDat.Cells.NoOfLocalUpdatedCells; j++) {
                VolFrac.SetMeanValue(j, volFrac[j]);
            }

            if (Katastrophenplot != null)
                Katastrophenplot(CellVolumesViz.Cat(MPIRank,AgglomCellsViz, AggCandidates, AggTargets, FailedViz, LevelSets, LvSetDist, CellNumbers, CellNumbersGlob, GlobalID, PairIdentification, PairColor, VolFrac),Tag);


            if (failCells.Count.MPISum() > 0) { 
               string message = ($"!!! Warning: Agglomeration failed - no candidate for agglomeration found for {failCells.Count.MPISum()} cells. ");
                if (ExceptionOnFailedAgglomeration)
                    throw new Exception(message);
                else
                    Console.WriteLine(message);
            }
        }

        /// <summary>
        /// returns the cell indices for an edge
        /// </summary>
        /// <param name="iEdge">local edge index</param>
        /// <returns></returns>
        static (int iEdge, int ThisCell, int OtherCell) GetEdgeInfo(int iEdge)
        {
            int OtherCell, ThisCell;
            if (iEdge < 0)
            {
                // Cell is the OUT-cell of edge
                OtherCell = 0;
                ThisCell = 1;
                iEdge *= -1;
            } else
            {
                OtherCell = 1;
                ThisCell = 0;
            }
            iEdge--;  // Adjust edge index for 0-based indexing
            return (iEdge, ThisCell, OtherCell);
        }


        /// <summary>
        /// Result of the algorithm
        /// </summary>
        public IEnumerable<(int SourceCell, int TargetCell)> AgglomerationPairs {
            get;
            private set;
        }

        /// <summary>
        /// Result of the algorithm
        /// </summary>
        public IEnumerable<CellAgglomerator.AgglomerationPair> AgglomerationPairsWithRanks {
            get;
            private set;
        }


    }




}
