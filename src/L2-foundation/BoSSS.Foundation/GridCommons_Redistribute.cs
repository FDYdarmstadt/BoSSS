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
using BoSSS.Foundation.IO;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Kraypis;
using ilPSP.HilbertCurve;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Foundation.Grid.Classic {

    partial class GridCommons {

        /// <summary>
        /// Driver method for grid redistribution; this includes 
        /// - computing a new partition (e.g. by <see cref="ComputePartitionMETIS"/>, <see cref="ComputePartitionHilbert"/>) 
        ///   *or* loading a predefined partition <see cref="PredefinedGridPartitioning"/> from the database
        /// - application of this partition to this grid, by a call to <see cref="RedistributeGrid(int[])"/>
        /// </summary>
        public void Redistribute(IDatabaseInfo db, GridPartType method, string PartOptions) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                int size;
                int rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (size <= 0)
                    return;

                //// invalid from now on
                //m_GlobalId2CellIndexMap = null;

                tr.Info($"Grid repartitioning method: {method}");
                tr.Info($"Grid repartitioning options: {PartOptions}");

                int[] part;
                switch (method) {
                    case GridPartType.METIS: {
                            int.TryParse(PartOptions, out int noOfPartitioningsToChooseFrom);
                            noOfPartitioningsToChooseFrom = Math.Max(1, noOfPartitioningsToChooseFrom);
                            part = ComputePartitionMETIS(noOfPartitioningsToChooseFrom: noOfPartitioningsToChooseFrom);
#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);
                            break;
                        }


                    case GridPartType.ParMETIS: {
                            int.TryParse(PartOptions, out int noOfRefinements);

                            part = ComputePartitionParMETIS();
#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);

                            for (int i = 0; i < noOfRefinements; i++) {
                                part = ComputePartitionParMETIS(refineCurrentPartitioning: true);
#if DEBUG
                                CheckPartitioning(part);
#endif
                                RedistributeGrid(part);
                            }
                            break;
                        }

                    case GridPartType.clusterHilbert: {
                            part = ComputePartitionHilbert(Functype: 0);
#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);
                            break;
                        }

                    case GridPartType.Hilbert: {
                            part = ComputePartitionHilbert(Functype: 1);
#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);
                            break;
                        }

                    case GridPartType.none:
                        break;

                    case GridPartType.Predefined: {
                            if (PartOptions == null || PartOptions.Length <= 0)
                                //throw new ArgumentException("'" + GridPartType.Predefined.ToString() + "' requires, as an option, the name of the Partition.", "PartOptions");
                                PartOptions = size.ToString();

                            if (!m_PredefinedGridPartitioning.ContainsKey(PartOptions)) {
                                StringWriter stw = new StringWriter();
                                for (int i = 0; i < m_PredefinedGridPartitioning.Count; i++) {
                                    stw.Write("'" + m_PredefinedGridPartitioning.Keys[i] + "'");
                                    if (i < (m_PredefinedGridPartitioning.Count - 1))
                                        stw.Write(", ");
                                }

                                throw new ArgumentException("Grid Partitioning with name '" + PartOptions + "' is unknown; known are: " + stw.ToString() + ";");
                            }

                            Console.WriteLine("redistribution according to " + PartOptions);

                            var partHelp = m_PredefinedGridPartitioning[PartOptions];
                            part = partHelp.CellToRankMap;
                            if (part == null) {
                                var cp = this.CellPartitioning;
                                var iom = db?.Controller?.DBDriver;
                                part = iom.LoadVector<int>(partHelp.Guid, ref cp).ToArray();

                                
                            }

#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);
                        }
                        break;

                    case GridPartType.OtherSession: {
                            Guid SessionID;
                            try {
                                SessionID = Guid.Parse(PartOptions);
                            } catch (Exception ex) {
                                throw new ArgumentException($"Exception parsing session ID from grid partitioning options {PartOptions}: {ex}");
                            }
                            tr.Info("Grid partition from Session: " + SessionID);
                            var sess = db.Sessions.Single(sx => sx.ID == SessionID);
                            var ts = sess.Timesteps.Last();
                            tr.Info("From timestep: " + ts);
                            if (ts.GridID != this.ID)
                                throw new ArgumentException($"Grid mismatch in restoring partitioning of session {SessionID}: last grid is {ts.GridID}, this sessions grid is {this.ID}");

                            var iom = db?.Controller?.DBDriver;
                            var MPIrankField = iom.LoadFields(ts, this.GridData, new string[] { "MPIrank" }).Single();
                            MPIrankField.GetExtremalValues(out double _minRank, out double _maxRank);
                            int minRank = Convert.ToInt32(_minRank), maxRank = Convert.ToInt32(_maxRank);
                            if (minRank < 0 || maxRank != this.Size - 1)
                                throw new ArgumentException($"`MPIrank` from previous session is in the range of {minRank} to {maxRank}, but expecting 0 to {this.Size - 1}");

                            int J = this.GridData.CellPartitioning.LocalLength;
                            part = new int[J];
                            double diff = 0;
                            int myRank = this.MyRank;
                            for (int j = 0; j < J; j++) {
                                
                                part[j] = Convert.ToInt32(MPIrankField.GetMeanValue(j));
                                int diff_j = part[j] - myRank;
                                diff += diff_j*diff_j;
                            }
                            diff = diff.MPISum().Sqrt();
                            tr.Info("l2 distance to current partition: " + diff);
#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);
                        }
                        break;

                    default:
                        throw new NotImplementedException();
                }
            }
        }

        private static void CheckPartitioning(int[] part) {
            int MpiRank, MpiSize;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MpiRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out MpiSize);


            int LocMin = 0;
            int LocMax = 0;
            for (int j = part.Length - 1; j >= 0; j--) {
                LocMin = Math.Min(LocMin, part[j]);
                LocMax = Math.Max(LocMax, part[j]);
            }

            int[] MinMax = (new int[] { -LocMin, LocMax }).MPIMax();
            int Min = -MinMax[0];
            int Max = MinMax[1];

            if (Min < 0) {
                throw new ApplicationException("Illegal MPI grid partition: minimum processor ranks is " + Min + ";");
            }
            if (Max >= MpiSize) {
                throw new ApplicationException("MPI grid partition not usable: specifies " + (Max + 1) + " processors, but currently running on " + MpiSize + " processors.");
            }

            int[] CellsPerRank = new int[MpiSize];
            for (int i = 0; i < part.Length; i++) {
                CellsPerRank[part[i]]++;
            }
            CellsPerRank = CellsPerRank.MPISum();

            List<int> Problems = new List<int>();
            for (int rnk = 0; rnk < MpiSize; rnk++) {
                if (CellsPerRank[rnk] == 0) {
                    Problems.Add(rnk);
                }
            }

            if (Problems.Count > 0) {
                using (var str = new StringWriter()) {
                    foreach (int r in Problems) {
                        str.Write(r);
                        str.Write(" ");
                    }
                    throw new ApplicationException("Illegal MPI partition: zero cells on MPI ranks: " + str.ToString());
                }
            }
        }

        /// <summary>
        /// Computes a grid partitioning (which cell should be on which processor)
        /// using the serial METIS library -- work is only done on MPi rank 0.
        /// </summary>
        /// <param name="ListofLocalCellWeights">
        /// If not null, defines the weight associated with each cell on the current process
        /// </param>
        /// <param name="noOfPartitioningsToChooseFrom">
        /// Tells METIS to compute
        /// </param>
        /// <returns>
        /// - Index: local cell index, 
        /// - content: MPI Processor rank
        /// This is the suggestion
        /// of ParMETIS for the grid partitioning:
        /// For each local cell index, the returned array contains the MPI
        /// process rank where the cell should be placed.
        /// </returns>
        public int[] ComputePartitionMETIS(IList<int[]> ListofLocalCellWeights = null, int noOfPartitioningsToChooseFrom = 1) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                int size = this.Size;
                int rank = this.MyRank;

                tr.Info("Number of cell Weights: " + (ListofLocalCellWeights?.Count ?? 0));


                if (size == 1) {
                    return new int[NoOfUpdateCells];
                }
                // int.MaxValue is 2^31 and fairly a big number with order of the magnitude 10^9. This exception might be wrong texted (!)
                if (this.NumberOfCells_l > int.MaxValue) {
                    throw new Exception(String.Format(
                        "Grid contains more than {0} cells and can thus not be partitioned using METIS. Use ParMETIS instead.",
                        int.MaxValue));
                }
                int J = (rank == 0) ? this.NumberOfCells : 0;

                // Setup communication; all send to rank 0
                SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
                if (rank != 0) {
                    sms.SetCommPath(0);
                }
                sms.CommitCommPaths();

                // Assemble adjacency lists on rank 0
                IEnumerable<Neighbour>[] neighboursGlobal = new IEnumerable<Neighbour>[J];
                {
                    IEnumerable<Neighbour>[] neighboursLocal = GetCellNeighbourship(IncludeBcCells: false, FilterPeriodicDuplicities:true).Take(NoOfUpdateCells).ToArray();
                    if (rank == 0) {
                        long localOffset = m_CellPartitioning.GetI0Offest(rank);
                        int localLength = m_CellPartitioning.GetLocalLength(rank);

                        for (int i = 0; i < localLength; i++) {
                            neighboursGlobal[localOffset + i] = neighboursLocal[i];
                        }
                    } else {
                        sms.Transmit(0, neighboursLocal);
                    }

                    while (sms.GetNext(out int senderRank, out IEnumerable<Neighbour>[] neighbours)) {
                        long localOffset = m_CellPartitioning.GetI0Offest(senderRank);
                        int localLength = m_CellPartitioning.GetLocalLength(senderRank);

                        if (neighbours.Length != localLength) {
                            throw new Exception();
                        }

                        for (int i = 0; i < localLength; i++) {
                            neighboursGlobal[localOffset + i] = neighbours[i];
                        }
                    }
                }

                if (ListofLocalCellWeights == null) {
                    ListofLocalCellWeights = new List<int[]>();
                    int L = this.CellPartitioning.LocalLength;
                    var dummy = new int[L];
                    dummy.SetAll(1);

                    ListofLocalCellWeights.Add(dummy);
                }

                var GlobalCellWeightsCollection = new List<int[]>();

                foreach (var cellWeightsLocal in ListofLocalCellWeights) {
                    // Gather global weights on rank 0
                    int[] cellWeightsGlobal = null;
                    if (cellWeightsLocal != null) {
                        cellWeightsGlobal = new int[J];
                        if (rank == 0) {
                            long localOffset = m_CellPartitioning.GetI0Offest(rank);
                            int localLength = m_CellPartitioning.GetLocalLength(rank);

                            for (int i = 0; i < localLength; i++) {
                                cellWeightsGlobal[localOffset + i] = cellWeightsLocal[i];
                            }
                        } else {
                            sms.Transmit(0, cellWeightsLocal);
                        }

                        while (sms.GetNext(out int senderRank, out int[] cellWeights)) {
                            long localOffset = m_CellPartitioning.GetI0Offest(senderRank);
                            int localLength = m_CellPartitioning.GetLocalLength(senderRank);

                            if (cellWeights.Length != localLength) {
                                throw new Exception("Length of cell weights exceeds Length of cell partitioning");
                            }

                            for (int i = 0; i < localLength; i++) {
                                cellWeightsGlobal[localOffset + i] = cellWeights[i];
                            }
                        }
                    }
                    GlobalCellWeightsCollection.Add(cellWeightsGlobal);
                }

                int[] globalResult = new int[J];
                if (rank == 0) {
                    int[] xadj = new int[J + 1];
                    List<int> adjncy = new List<int>(J * m_RefElements[0].NoOfFaces);
                    for (int j = 0; j < J; j++) {
                        var cNj = neighboursGlobal[j];
                        int E = cNj.Count();

                        for (int e = 0; e < E; e++) {
                            var NN = cNj.ElementAt(e);

                            if (NN.Neighbour_GlobalIndex >= 0
                                && !NN.IsPeriodicNeighbour) {
                                adjncy.Add((int)NN.Neighbour_GlobalIndex);
                            }
                        }
                        xadj[j + 1] = adjncy.Count;
                    }

                    Debug.Assert(
                        GlobalCellWeightsCollection.All(w => w.Length == GlobalCellWeightsCollection.First().Length),
                        "All cell weights arrays must have the same length!");

                    int[] cellWeightsFlattened = new int[GlobalCellWeightsCollection.Sum(c => c.Length)];
                    int index = 0;
                    for (int iCell = 0; iCell < GlobalCellWeightsCollection[0].Length; iCell++) {
                        for (int iConstraint = 0; iConstraint < GlobalCellWeightsCollection.Count; iConstraint++) {
                            cellWeightsFlattened[index] = GlobalCellWeightsCollection[iConstraint][iCell];
                            index++;
                        }
                    }

                    // Call Metis
                    int nparts = size;
                    int ncon = GlobalCellWeightsCollection.Count; // each vertex has ncon number of constraints
                    MPI_Comm wrld = csMPI.Raw._COMM.WORLD;

                    // Equal distribution of balance constraints (default)
                    float[] tpwgts = new float[ncon * nparts];
                    for (int i = 0; i < tpwgts.Length; i++) {
                        tpwgts[i] = 1.0f / nparts;
                    }

                    // Default imbalance tolerance (5%)
                    float[] ubvec = new float[ncon];
                    for (int i = 0; i < ubvec.Length; i++) {
                        ubvec[i] = 1.05F;
                    }

                    Debug.Assert((cellWeightsFlattened == null) == (ListofLocalCellWeights == null));
                    int objval = -1; // Value of the objective function at return time

                    tr.Info($"noOfPartitioningsToChooseFrom = {noOfPartitioningsToChooseFrom}");

                    int[] Options = new int[METIS.METIS_NOPTIONS];
                    Options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = noOfPartitioningsToChooseFrom; // 5 cuts
                    Options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
                    Options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
                    Options[(int)METIS.OptionCodes.METIS_OPTION_SEED] = 0;
                    Options[(int)METIS.OptionCodes.METIS_OPTION_OBJTYPE] = 0; // corresponds to Edge-cut minimization

                    METIS.ReturnCodes status = (METIS.ReturnCodes)METIS.PARTGRAPHKWAY(
                        nvtxs: ref J,
                        ncon: ref ncon,
                        xadj: xadj,
                        adjncy: adjncy.ToArray(),
                        vwgt: cellWeightsFlattened, // if null, METIS assumes all have weight 1
                        vsize: null, // No information about communication size
                        adjwgt: null, // No edge weights
                        nparts: ref nparts,
                        tpwgts: tpwgts, // equal partitioning of constrain
                        ubvec: ubvec,
                        options: Options,
                        objval: ref objval,
                        part: globalResult);

                    if (status != METIS.ReturnCodes.METIS_OK) {
                        throw new Exception(String.Format(
                            "Error partitioning the mesh. METIS reported {0}",
                            status));
                    }

                    int[] CountCheck = new int[size];
                    int J2 = this.NumberOfCells;
                    for (int i = 0; i < J2; i++) {
                        CountCheck[globalResult[i]]++;
                    }
                    for (int rnk = 0; rnk < size; rnk++) {
                        if (CountCheck[rnk] <= 0) {
                            throw new ApplicationException("METIS produced illegal partitioning - 0 cells on process " + rnk + ".");
                        }
                    }
                    globalResult = SortPartitioning(globalResult);
                }


                int[] localLengths = new int[size];
                for (int p = 0; p < localLengths.Length; p++) {
                    localLengths[p] = this.CellPartitioning.GetLocalLength(p);
                }
                int[] localResult = globalResult.MPIScatterv(localLengths);

                int[] CellsPerRank = new int[size];
                foreach (var rnk in localResult)
                    CellsPerRank[rnk]++;
                CellsPerRank = CellsPerRank.MPISum();
                tr.Info($"Cells per rank: " + CellsPerRank.ToConcatString("[", ", ", "]"));


                return localResult;
            }
        }

        /// <summary>
        /// Sort the new partitioning w.r.t. current one to minimize the cost of re-distribution.
        /// Because partitioning algorithms do not take the current partitioning into account
        /// and assigns an arbitrary processor to partitions.
        /// </summary>
        /// <returns>
        /// Global list for partitioning
        /// For each local cell index, the returned array contains the MPI
        /// process rank where the cell should be placed.
        /// - Index: cell index, 
        /// - content: MPI Processor rank
        /// </returns>
        /// <param name="part"> Global list for partitioning
        /// - Index: cell index, 
        /// - content: MPI Processor rank
        ///</param>  
        internal int[] SortPartitioning(int[] part) {
            if (m_Size < 2)
                    return part;

            int[] sortedPart = new int[part.Length];
            List<int> mapping = new List<int>(m_Size);

            // Find a suitable partition for each processor(p)
            for (int p = 0; p < m_Size; p++){
                // The info about the old partitioning can be accessed by using offsets
                // Get the offset and length of the processor to mask cells with respect to (w.r.t) old partitioning
                long localOffset = m_CellPartitioning.GetI0Offest(p);
                int localLength = m_CellPartitioning.GetLocalLength(p);

                // Sort the input cell-packages w.r.t their occurrences in the old partitioning
                var occurrences = part.GetSubVector((int)localOffset, localLength).GroupBy(v => v)
                  .Select(g => new { Value = g.Key, Count = g.Count() })
                  .OrderByDescending(x => x.Count)
                  .ToList();

                // Make sure that you supply each processor as a possible target for the partitions
                for (int pp = 0; pp < m_Size; pp++) {
                    // As long as they are not assigned to anyone ...
                    if (!mapping.Contains(pp) && !occurrences.Any(x => x.Value == pp)) {
                        // ... we can add this processor(pp) even though they do not have any cells in the current partitioning. 
                        occurrences.Add(new { Value = pp, Count = 0 });
                    }
                }

                // Find the most repeated rank w.r.t. old partitioning
                foreach (var processor in occurrences) {
                    if (!mapping.Contains((int)processor.Value)) { // make sure it is not already assigned to another partition
                        mapping.Add((int)processor.Value); // content: old rank, index: new rank
                        break;
                    }
                }
            }

            // Check if every rank is assigned to mapping
            for (int p = 0; p < m_Size; p++) {
                if (!mapping.Contains(p))
                    throw new Exception("Cannot sort the distribution of ranks");
            }

            // Assign the new partitioning w.r.t. mapping
            for (int j = 0; j < part.Length; j++) {
                sortedPart[j] = mapping.IndexOf(part[j]); // Index of the list represents the new rank
            }

            return sortedPart;
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <returns>
        ///  - Index: local cell index
        ///  - content: MPI Processor rank;
        /// This is the suggestion
        /// of ParMETIS for the grid partitioning:
        /// For each local cell index, the returned array contains the MPI processor rank
        /// where the cell should be placed.<br/>
        /// may be null, if no repartitioning is required at all.
        /// </returns>
        internal int[] SMPRedistribution() {
            using (new FuncTrace()) {
                return null;
                //throw new NotImplementedException("todo");
                /*
                Comm.Master cm = ctx.CommMaster;

                if (cm.Size <= 1 || cm.SMPSize <= 1) {
                    return null; // nothing to do
                }

                int J = NoOfUpdateCells;
                //int E = GridSimplex.NoOfEdges;

                int[] _MyOwnNeighbourProcesses; // a list of all processes which exchange info with this process
                {
                    // this section is very expensive - however, it is uses only at startup,
                    // i don't care; (Die Arbeitszeit von Programmierern ist teuer; spare sie auf Kosten der Rechenzeit. Eric S. Raymond, The Art of UNIX Programming)

                    long[] CellNeighbours_oneCol = ArrayTools.Resize(CellNeighbours, false); // temporary resize the array
                    for (int i = 0; i < CellNeighbours_oneCol.Length; i++)
                        if (CellNeighbours_oneCol[i] < 0)
                            // Neighbour not defined - we reset it to a cell which is definetly
                            // on this process, so this won't effect the neighbourhood calc.
                            // (An undefined, i.e. negative entry will lead to an exeption during
                            // globalIndices.EvaluatePermutation(...), a few lines below).
                            CellNeighbours_oneCol[i] = GlobalID[0];

                    // transform GlobalID's into global cell indices
                    Permutation _globalID = new Permutation(GlobalID, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    Permutation globalIndices = _globalID.Invert();

                    long[] CellNeighbours_ByIndex = new long[CellNeighbours_oneCol.Length];
                    globalIndices.EvaluatePermutation(CellNeighbours_oneCol, CellNeighbours_ByIndex);
                    CellNeighbours_oneCol = null;


                    List<int> NeighProcesses = new List<int>();

                    Partition GridPerm = new Partition(J, MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                    for (int i = 0; i < CellNeighbours_ByIndex.Length; i++) {
                        int neig_cell_procrank = GridPerm.FindProcess(CellNeighbours_ByIndex[i]);

                        if (!NeighProcesses.Contains(neig_cell_procrank) && neig_cell_procrank != cm.MyRank)
                            NeighProcesses.Add(neig_cell_procrank);
                    }

                    _MyOwnNeighbourProcesses = NeighProcesses.ToArray();
                }


                // define SMP rank;
                // for each MPI process, the SMP node index
                // index: MPI rank; content: SMP rank;
                int[] SMPRank = new int[cm.Size];
                int NoOfSMPs = -1;
                {
                    NoOfSMPs = cm.SMPSize;
                    for (int r = 0; r < cm.Size; r++)
                        SMPRank[r] = cm.SMPRank(r);
                }

                // 1st index: processor rank p
                // 2nd index: only a list
                // content: all mpi processor ranks that communicate with p
                int[][] NeighbourProcesses = null;
                {
                    SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
                    if (cm.MyRank > 0)
                        sms.SetCommPath(0);
                    sms.CommitCommPaths();

                    if (cm.MyRank > 0) {
                        //Console.WriteLine("rank: " + cm.MyRank + ": sending to 0");
                        sms.Transmitt(0, _MyOwnNeighbourProcesses);
                    }

                    int[] nighs;
                    int rcvRnk;
                    sms.GetNext<int[]>(out rcvRnk, out nighs);

                    if (cm.MyRank == 0) {

                        NeighbourProcesses = new int[cm.Size][];
                        NeighbourProcesses[0] = _MyOwnNeighbourProcesses;

                        while (nighs != null) {
                            NeighbourProcesses[rcvRnk] = nighs;
                            //Console.WriteLine("rank: 0: received from " + rcvRnk);

                            sms.GetNext(out rcvRnk, out nighs);
                        }

                        for (int rnk = 0; rnk < cm.Size; rnk++) {

                            string msg = "";
                            for (int i = 0; i < NeighbourProcesses[rnk].Length; i++)
                                msg = msg + NeighbourProcesses[rnk][i] + " ";

                            Logger.Info("Before SMP redist.: Proc " + rnk + " communicates with " + msg);
                        }

                    } else {
                        if (nighs != null)
                            throw new ApplicationException("should not happen.");
                    }

                    sms.Dispose();
                }

                // compute number of inter-smp - communication
                {
                    //throw new Exception("todo");
                }


                // find smp distribution using metis on processor 0
                int[] nodesPart = null;
                bool Sucess = true;
                {

                    if (cm.MyRank == 0) {

                        // define METIS input
                        int n = cm.Size;
                        int[] xadjncy = new int[n + 1];
                        List<int> _adjncy = new List<int>();

                        int cnt = 0;
                        for (int i = 0; i < cm.Size; i++) {
                            xadjncy[i] = cnt;
                            foreach (int neigh_rnk in NeighbourProcesses[i]) {
                                _adjncy.Add(neigh_rnk);
                                cnt++;
                            }
                        }
                        xadjncy[cm.Size] = cnt;
                        int[] adjncy = _adjncy.ToArray();

                        //Debugger.Break();

                        float[] tpwgts = new float[cm.SMPSize];
                        for (int smpRnk = 0; smpRnk < cm.SMPSize; smpRnk++) {
                            tpwgts[smpRnk] = (float)((double)cm.MPIProcessesPerSMP(smpRnk) / (double)cm.Size);
                        }


                        // call METIS
                        int numflag = 0;
                        int wgtflag = 0;
                        int nparts = NoOfSMPs;
                        int[] options = new int[5];
                        nodesPart = new int[n];
                        int edgecut = -1;


                        Logger.Info("Trying 'METIS_WPartGraphKway'...");
                        METIS.WPartGraphKway(ref n, xadjncy, adjncy, null, null, ref wgtflag, ref numflag,
                            ref nparts, tpwgts, options, ref edgecut, nodesPart);

                        if (!CheckPartition(cm, nodesPart)) {

                            Logger.Info("failed. Trying 'METIS_WPartGraphRecursive'...");
                            METIS.WPartGraphRecursive(ref n, xadjncy, adjncy, null, null, ref wgtflag, ref numflag,
                                ref nparts, tpwgts, options, ref edgecut, nodesPart);

                            if (!CheckPartition(cm, nodesPart)) {
                                Logger.Info("failed. Trying 'METIS_WPartGraphVKway'...");
                                METIS.WPartGraphVKway(ref n, xadjncy, adjncy, null, null, ref wgtflag, ref numflag,
                                    ref nparts, tpwgts, options, ref edgecut, nodesPart);


                                if (!CheckPartition(cm, nodesPart)) {

                                    Sucess = false;
                                }
                            }
                        }

                        //// check METIS output
                        //{
                        //    int[] res = new int[NoOfSMPs];
                        //    for (int mpiRnk = 0; mpiRnk < cm.Size; mpiRnk++)
                        //        res[nodesPart[mpiRnk]]++;
                        //    for (int smpRnk = 0; smpRnk < NoOfSMPs; smpRnk++) {
                        //        ctx.IOMaster.tracer.Message(ht, "METIS: Processes for SMP node #" + smpRnk + ": " + res[smpRnk]);
                        //        ctx.IOMaster.tracer.Message(ht, "Available cores for SMP node #" + smpRnk + ": " + cm.MPIProcessesPerSMP(smpRnk));
                        //    }
                        //    for (int smpRnk = 0; smpRnk < NoOfSMPs; smpRnk++)
                        //        if (res[smpRnk] != cm.MPIProcessesPerSMP(smpRnk))
                        //            throw new ApplicationException("METIS output is strange.");
                        //}
                    }

                    Sucess = Comm.Master._Bcast(Sucess, 0, csMPI.Raw._COMM.WORLD);
                    if (!Sucess) {
                        Logger.Info("FAILED to find a balanced SMP distribution with METIS. - NO SMP redistribution.");
                        return null;
                    }
                }

                // for each MPI rank, the MPI rank to which the grid data should be send;
                // index: MPI rank
                int[] mpiNodes = null;
                bool noChange = true;
                {
                    if (cm.MyRank == 0) {
                        mpiNodes = new int[cm.Size];
                        bool[] RankOccupied = new bool[cm.Size];
                        ArrayTools.SetAll(mpiNodes, -1);

                        for (int rnk = 0; rnk < cm.Size; rnk++) {
                            if (SMPRank[rnk] == nodesPart[rnk]) {
                                // this rank does not need to be moved - it
                                // stays on the same computer
                                mpiNodes[rnk] = rnk;
                                RankOccupied[rnk] = true;
                            }
                        }

                        for (int rnk = 0; rnk < cm.Size; rnk++) {
                            if (mpiNodes[rnk] < 0) {
                                // this rank needs to be moved - search a new rank

                                noChange = false;

                                bool found = false;
                                for (int targ_rnk = 0; targ_rnk < cm.Size; targ_rnk++) {
                                    if (RankOccupied[targ_rnk])
                                        continue; // target Rank (targ_rank) allready occupied
                                    if (SMPRank[targ_rnk] != nodesPart[rnk])
                                        continue; // target Rank (targ_rank) on wrong computer

                                    mpiNodes[rnk] = targ_rnk;
                                    RankOccupied[targ_rnk] = true;
                                    found = true;
                                    break;
                                }

                                if (!found)
                                    // unable to find some target rank
                                    throw new ApplicationException("fatal error in algorithm");
                            }
                        }
                    }

                    noChange = Comm.Master._Bcast(noChange, 0, csMPI.Raw._COMM.WORLD);
                    mpiNodes = Comm.Master._Bcast(mpiNodes, 0, csMPI.Raw._COMM.WORLD);
                }

                int[] part = null;
                if (!noChange) {
                    part = new int[J];

                    int targRnk = mpiNodes[cm.MyRank];
                    for (int j = 0; j < J; j++)
                        part[j] = targRnk;
                }


                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                return part;
              */
            }
        }


        private delegate METIS.ReturnCodes ParMETISAction(
            int[] vtxdist,
            int[] xadj,
            int[] adjncy,
            int[] vwgt,
            int[] adjwgt,
            ref int wgtflag,
            ref int numflag,
            ref int ncon,
            ref int nparts,
            float[] tpwgts,
            float[] ubvec,
            int[] options,
            ref int edgecut,
            int[] partitioningResult,
            MPI_Comm MPIComm);

        /// <summary>
        /// Computes a grid partitioning (which cell should be on which processor)
        /// by calling ParMETIS.
        /// </summary>
        /// <param name="cellWeights">
        /// If not null, defines the (list of) weights associated with each
        /// cell on this process. If multiple weights are present per cell,
        /// this implies that there are multiple balance constraints to be
        /// obeyed by ParMETIS
        /// </param>
        /// <param name="refineCurrentPartitioning">
        /// Refines the current partitioning instead of starting from scratch
        /// </param>
        /// <returns>
        /// Index: local cell index, content: MPI Processor rank;<br/>
        /// This is the suggestion
        /// of ParMETIS for the grid partitioning:
        /// For each local cell index, the returned array contains the MPI
        /// process rank where the cell should be placed.
        /// </returns>
        public int[] ComputePartitionParMETIS(IList<int[]> cellWeights = null, bool refineCurrentPartitioning = false) {
            using (new FuncTrace()) {
                int size = this.Size;
                int rank = this.MyRank;
                if (size > 1) {
                    // vtxdist array
                    var CPart = this.CellPartitioning;
                    int[] currentCellPartitioning = new int[size + 1];
                    for (int p = 0; p < currentCellPartitioning.Length; p++) {
                        currentCellPartitioning[p] = (int)CPart.GetI0Offest(p);
                    }

                    // Assemble adjacency lists
                    var CellNeighs = GetCellNeighbourship(IncludeBcCells: false, FilterPeriodicDuplicities:true);
                    int J = NoOfUpdateCells;
                    int[] xadj = new int[J + 1];
                    List<int> adjncyL = new List<int>(J * m_RefElements[0].NoOfFaces);
                    for (int j = 0; j < J; j++) {
                        var Cl = Cells[j];
                        var cNj = CellNeighs[j];
                        int E = cNj.Count();

                        for (int e = 0; e < E; e++) {
                            var NN = cNj.ElementAt(e);

                            if (NN.Neighbour_GlobalIndex >= 0
                                && !NN.IsPeriodicNeighbour) {
                                adjncyL.Add((int)NN.Neighbour_GlobalIndex);
                            }
                        }
                        xadj[j + 1] = adjncyL.Count;
                    }

                    // Prepare ParMETIS parameters
                    ParMETISAction parmetisAction;
                    if (refineCurrentPartitioning) {
                        parmetisAction = ParMETIS.V3_RefineKway;
                    } else {
                        parmetisAction = ParMETIS.V3_PartKway;
                    }


                    if (cellWeights == null) {
                        // Cell weights null causes problems with ParMETIS
                        cellWeights = new List<int[]>() { new int[NoOfUpdateCells] };
                        cellWeights.Single().SetAll(1);
                    }

                    Debug.Assert(
                        cellWeights.All(w => w.Length == cellWeights[0].Length),
                        "All cell weights arrays must have the same length!");

                    int[] cellWeightsFlattened = new int[cellWeights.Sum(c => c.Length)];
                    int index = 0;
                    for (int iCell = 0; iCell < cellWeights[0].Length; iCell++) {
                        for (int iConstraint = 0; iConstraint < cellWeights.Count; iConstraint++) {
                            cellWeightsFlattened[index] = cellWeights[iConstraint][iCell];
                            index++;
                        }
                    }

                    int nparts = size;
                    int ncon = cellWeights.Count; // Number of balance constraints
                    MPI_Comm wrld = csMPI.Raw._COMM.WORLD;
                    int wgtflag = 2; // Cell weights only
                    int numflag = 0; // 0 -> use C-style ordering

                    // Equal distribution of balance constraints (default)
                    float[] tpwgts = new float[ncon * nparts];
                    for (int i = 0; i < tpwgts.Length; i++) {
                        tpwgts[i] = 1.0f / nparts;
                    }

                    // Default imbalance tolerance (5%)
                    float[] ubvec = new float[ncon];
                    for (int i = 0; i < ubvec.Length; i++) {
                        ubvec[i] = 1.05F;
                    }

                    // Call ParMETIS
                    int edgecut = 0;
                    int[] result = new int[J + 1];
                    METIS.ReturnCodes status = parmetisAction(
                        vtxdist: currentCellPartitioning,
                        xadj: xadj,
                        adjncy: adjncyL.ToArray(),
                        vwgt: cellWeightsFlattened,
                        adjwgt: null, // No edge weights
                        wgtflag: ref wgtflag,
                        numflag: ref numflag,
                        ncon: ref ncon,
                        nparts: ref nparts,
                        tpwgts: tpwgts,
                        ubvec: ubvec,
                        options: new int[] { 0, 0, 0 }, // use default
                        edgecut: ref edgecut,
                        partitioningResult: result,
                        MPIComm: wrld);

                    if (status != METIS.ReturnCodes.METIS_OK) {
                        throw new Exception(String.Format(
                            "Error partitioning the mesh. ParMETIS reported {0}",
                            status));
                    }

                    Array.Resize(ref result, J);
                    return result;
                } else {
                    int J = NoOfUpdateCells;
                    return new int[J];
                }
            }
        }

        /*
        public int[] ComputePartitionHilbert(int[] cellCosts = null) {
            // Step 1: Compute some global indexing pattern for all cells
            // according to space-filling curve
            long globalNumberOfCells = long.MaxValue;
            int[] localToGlobalIndexMap = null; // Maps local cell index to global index in clusterHilbert curve

            // Step 2: Compute numbers of cells per process (TODO: obey weighting!)
            int size = this.Size;
            int minNoOfCellsPerProcess = (int)(globalNumberOfCells / size);
            int noOfSurplusCells = (int)(globalNumberOfCells % size);

            int[] numbersOfCellsPerProcess = new int[size];
            numbersOfCellsPerProcess.SetAll(minNoOfCellsPerProcess);
            for (int i = 0; i < noOfSurplusCells; i++) {
                numbersOfCellsPerProcess[i]++;
            }

            // Step 3: Determine target rank for each local cell
            int[] partitioning = new int[NoOfUpdateCells];
            for (int i = 0; i < NoOfUpdateCells; i++) {
                long globalIndex = localToGlobalIndexMap[i];

                int targetRank = 0;
                long firstIndexForTargetRank = 0;
                while (globalIndex >= firstIndexForTargetRank + numbersOfCellsPerProcess[targetRank]) {
                    firstIndexForTargetRank += numbersOfCellsPerProcess[targetRank];
                    targetRank++;
                }
                partitioning[i] = targetRank;
            }

            return partitioning;
        }*/

        /*
        private double[] GetShortestDistance() {
            int D = this.SpatialDimension;
            IEnumerable<Neighbour>[] bla = this.GetCellNeighbourship(true);
            var GlobalBB = this.GetGridBoundingBox();
            int cellInd = 0;
            double[] ShortestGlobalDistance = GlobalBB.Max;

            foreach (IEnumerable<Neighbour> enumNb in bla) {
                //Loop over all Cells
                Cell C_O = this.Cells[cellInd];
                double[] CenterC_O = new double[D];

                //Barycentre of Origin
                for (int d = 0; d < D; d++) {
                    double NoOfNodes = C_O.TransformationParams.NoOfRows;
                    double center = 0;
                    for (int k = 0; k < NoOfNodes; k++) {
                        center += C_O.TransformationParams[k, d];
                    }
                    CenterC_O[d] = center / (NoOfNodes);
                }

                foreach (Neighbour nb in enumNb) {
                    //Loop over Neighbours
                    long NIndex = nb.Neighbour_GlobalIndex;
                    Cell C_N = this.Cells[NIndex - this.CellPartitioning.i0];
                    double[] CenterC_N = new double[D];

                    //Barycentre of Neighbour
                    for (int d = 0; d < D; d++) {
                        double NoOfNodes = C_N.TransformationParams.NoOfRows;
                        double center = 0;
                        for (int k = 0; k < NoOfNodes; k++) {
                            center += C_N.TransformationParams[k, d];
                        }
                        CenterC_N[d] = center / (NoOfNodes);

                        //Save shortest distance
                        double distance = Math.Abs(CenterC_O[d] - CenterC_N[d]);
                        ShortestGlobalDistance[d] = Math.Min(ShortestGlobalDistance[d], distance);
                    }
                }
                cellInd++;
            }
            return ShortestGlobalDistance.MPIMin();
        }
        */

        private BoundingBox GetGridBoundingBox() {
            int D = this.SpatialDimension;
            var BB = new BoundingBox(D);

            long J = this.CellPartitioning.LocalLength;
            for (int j = 0; j < J; j++) {
                Cell Cj = this.Cells[j];
                BB.AddPoints(Cj.TransformationParams);
            }

            BB.Min = BB.Min.MPIMin();
            BB.Max = BB.Max.MPIMax();

            return BB;
        }

        /*
        #region for testing and debugging
        private List<long> VerdachtsCellen = new List<long>();
        public List<long> GetZellsOfChangingProc {
            get { return VerdachtsCellen.MPIBroadcast(0); }
        }
        private ulong[] m_HilbertIdx;
        public ulong[] GetLocHilbertIdcs {
            get { return m_HilbertIdx; }
        }
        #endregion
        */

        /// <summary>
        /// Computes a grid partitioning (which cell should be on which processor) based on a Hilbertcurve of maximum order (64 bit>nBit*nDim).
        /// Note: The identification of barycenter only works for rectangular-shaped cells.
        /// </summary>
        public int[] ComputePartitionHilbert(IList<int[]> localcellCosts = null, int Functype = 0, bool adjustRefinement=false) {

            //CONSTRUCTIONSITE: GetShortestDistance will be called from here, the return argument will be used to calculate a refinement order
            //This refinement order has to be stored once and is available in here during runtime
            
            //Functype: Constraint mapping (0) or direct Costmapping (1)
            //Notice: Functype=1 will lead to bad behavior, when using Clusters

            int[] local_Rank_RedistributionList;

            if (this.Size > 1) {
                int D = this.SpatialDimension;
                var GlobalBB = this.GetGridBoundingBox();

                int J = this.CellPartitioning.LocalLength;
                
                ulong[] discreteCenter = new ulong[D];
                ulong[] local_HilbertIndex = new ulong[J];
                long[] local_CellIndex = new long[J];
                int MaxRef = 64;

                for (int j = 0; j < J; j++) {
                    Cell Cj = this.Cells[j];
                    int NoOfNodes = Cj.TransformationParams.NoOfRows;
                    //Compute Barycenter of rectangular cells
                    for (int d = 0; d < D; d++) {
                        double center = 0;
                        for (int k = 0; k < NoOfNodes; k++) {
                            center += Cj.TransformationParams[k, d];
                        }

                        center = center / ((double)NoOfNodes); // ''center of gravity'' for coordinate direction 'd'
                        double centerTrf = (center - GlobalBB.Min[d]) * (1.0 / (GlobalBB.Max[d] - GlobalBB.Min[d])) * Math.Pow(2, MaxRef / D);
                        centerTrf = Math.Round(centerTrf);
                        if (centerTrf < 0)
                            centerTrf = 0;
                        if (centerTrf > ulong.MaxValue)
                            centerTrf = ulong.MaxValue;
                        discreteCenter[d] = (ulong)centerTrf;
                        //Debugger.Break();
                    }
                    //Derive Hilbertindex from Barycenter
                    ulong iH = HilbertCurve.hilbert_c2i(MaxRef / D, discreteCenter);
                    local_HilbertIndex[j] = iH;
                    local_CellIndex[j] = j + this.CellPartitioning.i0;
                }

                //Gather all stuff for computation on rank==0
                int numberofcells = this.NumberOfCells;
                int[] CellsPerRank = new int[this.Size];
                for (int r = 0; r < CellsPerRank.Length; r++) {
                    CellsPerRank[r] = this.CellPartitioning.GetLocalLength(r);
                }
                long[] CellIndex = local_CellIndex.MPIGatherv(CellsPerRank);
                ulong[] HilbertIndex = local_HilbertIndex.MPIGatherv(CellsPerRank);
                //m_HilbertIdx = local_HilbertIndex;
                List<int[]> cellCosts = new List<int[]>();
                if (localcellCosts != null) {
                    cellCosts = new List<int[]>() { new int[numberofcells] };
                    foreach (int[] cellCostmap in localcellCosts) {
                        cellCosts.Add(cellCostmap.MPIGatherv(CellsPerRank));
                    }
                    cellCosts.RemoveAt(0);
                }
                int[] RankIndex = new int[numberofcells];

                //Sequential Part
                if (MyRank == 0) {
                    //catch case localcellsCosts is empty
                    if (localcellCosts == null) {
                        cellCosts = new List<int[]>() { new int[numberofcells] };
                        switch (Functype) {
                            case 0:
                                cellCosts.Single().SetAll(10);
                                break;
                            case 1:
                                cellCosts.Single().SetAll(1);
                                break;
                        }
                    }

                    ulong[] HilbertIndex_tmp = HilbertIndex.CloneAs<ulong[]>();
                    Array.Sort(HilbertIndex, CellIndex);

                    
                    int numproc = this.Size;
                    int[] CostClustermap = new int[numberofcells];

                    switch (Functype) {
                    case 0:
                        //Distribution of MPI-Rank along Hilbertcurve considering multiple Constraints with equal weight
                        //Assuming: Constraintmap only contains two types = (1, a number>1); Constraints do not overlap, but in sum cover whole grid
                        Debug.Assert(cellCosts.Count >= 1);
                        int[,] buckets = new int[numproc, cellCosts.Count];
                        int[] CellsPerCluster = new int[cellCosts.Count];
                        int CheckCount = 0;
                        int[] CheckInt = new int[numberofcells];

                        //unite all constraints in CostClustermap
                        for (int k = 0; k < cellCosts.Count; k++) {
                            int CountTargetCells = 0;
                            for (int i = 0; i < cellCosts.ElementAt(k).Length; i++) {
                                Debug.Assert(cellCosts.ElementAt(k).Length == numberofcells);
                                if (cellCosts.ElementAt(k)[i] > 1) {
                                    CountTargetCells++;
                                    CostClustermap[i] = k;
                                    CheckInt[i] += 1;
                                }
                                if (CheckInt[i] > 1) {
                                    throw new ArgumentException("Clusters are overlapping!");
                                }
                            }
                            CellsPerCluster[k] = CountTargetCells;
                            CheckCount += CountTargetCells;
                        }

                        if (CheckCount != numberofcells)
                            throw new ArgumentException("There are unmasked Cells in cellCosts!");
                        Array.Sort(HilbertIndex_tmp, CostClustermap); //Would be better to unify with Array.Sort(HilbertIndex, CellIndex)

                        //fill buckets, buckets contain number of cells to distribute per rank
                        for (int iCluster = 0; iCluster < CellsPerCluster.Length; iCluster++) {
                            int N = CellsPerCluster[iCluster];
                            int TestSum = 0;
                            for (int iProc = 0; iProc < numproc; iProc++) {
                                buckets[iProc, iCluster] = N * (iProc + 1) / numproc - N * iProc / numproc;
                                TestSum += buckets[iProc, iCluster];
                            }
                            Debug.Assert(CellsPerCluster[iCluster] == TestSum);

                        }

                        //Go through HilbertCurve and map rank->cell
                        int rank = 0;
                        for (int cell = 0; cell < CostClustermap.Length;) {
                            if (buckets[rank, CostClustermap[cell]] != 0) {
                                buckets[rank, CostClustermap[cell]]--;
                                RankIndex[cell] = rank;
                                cell++;
                            } else {
                                rank = (rank + 1) % numproc;
                                //Console.WriteLine("switch at: {0}, from rank{1}->{2}", cell, oldrank, rank);
                                //Console.WriteLine("proc{0},cluster{2}: {1}", oldrank, buckets[oldrank, CostClustermap[cell]], CostClustermap[cell]);
                                //Console.WriteLine("proc{0},cluster{2}: {1}", rank, buckets[rank, CostClustermap[cell]], CostClustermap[cell]);
                                //VerdachtsCellen.Add(CellIndex[cell]);
                            }
                        }
                        break;
                    case 1:
                        //direct cost mapping
                        CostClustermap = cellCosts[0];
                        if (cellCosts.Count > 1) {
                            Console.WriteLine("More than one cell cost cluster chosen. You can use cluster Hilbert instead. Proceeding with equal cell costs ...");
                            var tmp = new int[cellCosts[0].Length];
                            tmp.SetAll(1);
                            CostClustermap = tmp;
                        }

                        Array.Sort(HilbertIndex_tmp, CostClustermap);
                        int CellCostSum = 0;
                        for (int cell = 0; cell < CostClustermap.Length; cell++)
                            CellCostSum += CostClustermap[cell];
                       
                        int MPIrank = -1;
                        int CostofProc = 0;
                        int DistributedCost = 0;
                        int TestSumC = 0;
                        for (int iCell = 0; iCell < numberofcells; iCell++) {
                            if (DistributedCost >= CostofProc) {
                                MPIrank++;
                                CostofProc = CellCostSum * (MPIrank + 1) / numproc - CellCostSum * MPIrank / numproc;
                                TestSumC += CostofProc;
                                DistributedCost = 0;
                            }
                            RankIndex[iCell] = MPIrank;
                            DistributedCost += CostClustermap[iCell];
                        }

                        Debug.Assert(TestSumC == CellCostSum);
                        Debug.Assert(MPIrank == this.Size-1);
                        break;
                    default:
                        throw new NotImplementedException();
                    }
                    Array.Sort(CellIndex, RankIndex);
                }
                //RankIndex = SortPartitioning(RankIndex);
                //Scatter Rank-Array for local Process
                local_Rank_RedistributionList = RankIndex.MPIScatterv(CellsPerRank);
            } else {
                local_Rank_RedistributionList = new int[NoOfUpdateCells];
            }
            return local_Rank_RedistributionList;
        }

        private bool CheckPartitioning(Master cm, int[] nodesPart) {

            int NoOfSMPs = cm.SMPSize;

            int[] res = new int[NoOfSMPs];
            for (int mpiRnk = 0; mpiRnk < cm.Size; mpiRnk++)
                res[nodesPart[mpiRnk]]++;

            for (int smpRnk = 0; smpRnk < NoOfSMPs; smpRnk++) {
                Logger.Info("METIS: Processes for SMP node #" + smpRnk + ": " + res[smpRnk]);
                Logger.Info("Available cores for SMP node #" + smpRnk + ": " + cm.MPIProcessesPerSMP(smpRnk));
            }

            for (int smpRnk = 0; smpRnk < NoOfSMPs; smpRnk++)
                if (res[smpRnk] != cm.MPIProcessesPerSMP(smpRnk))
                    return false;

            return true;

        }

        /// <summary>
        /// redistributes this grid, i.e. sends cells to different processors
        /// </summary>
        /// <param name="part">
        /// MPI processor rank for each cell; index: local cell index;
        /// </param>
        public void RedistributeGrid(int[] part) {
            int Size = this.Size;
            int MyRank = this.MyRank;

            CheckPartitioning(part);
            
            // partition is no longer valid anymore!
            InvalidateGridData();
            m_CellPartitioning = null;

            //
            int J = NoOfUpdateCells;

            if (part.Length < J)
                throw new ArgumentException();

            // send cells to other processors
            // ==============================


            // count number of items for each processor
            int[] ItmCnt = new int[Size];
            for (int j = 0; j < J; j++)
                ItmCnt[part[j]]++;

            // create messenger
            SerialisationMessenger sm = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
            for (int p = 0; p < Size; p++)
                if (ItmCnt[p] > 0 && p != MyRank)
                    sm.SetCommPath(p);
            sm.CommitCommPaths();

            // sort cells according to processors (part notes which proc. will get a cell)
            List<Cell>[] exch = new List<Cell>[Size];
            for (int p = 0; p < Size; p++)
                if (ItmCnt[p] > 0 || p == MyRank) {
                    exch[p] = new List<Cell>(ItmCnt[p]);
                }
            for (int j = 0; j < J; j++) {
                exch[part[j]].Add(Cells[j]);
            }

            // start transmission
            for (int p = 0; p < Size; p++) {
                if (exch[p] != null && p != MyRank)
                    sm.Transmit(p, exch[p]);
            }

            // receive cells from other processors
            var myown = exch[MyRank];
            List<Cell> o;
            int prc;
            while (sm.GetNext<List<Cell>>(out prc, out o)) {
                myown.AddRange(o);
            }

            // write back data to arrays
            this.Cells = myown.ToArray();
        }

    }
}
