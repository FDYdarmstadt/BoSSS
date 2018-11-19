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
        public void Redistribute(IDatabaseDriver iom, GridPartType method, string PartOptions) {
            using (new FuncTrace()) {
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                int size;
                int rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                if (size <= 0)
                    return;

                //// invalid from now on
                //m_GlobalId2CellIndexMap = null;

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

                    case GridPartType.Hilbert: {
                            part = ComputePartitionHilbert(Functype:0);
#if DEBUG
                            CheckPartitioning(part);
#endif
                            RedistributeGrid(part);
                            break;
                        }

                    case GridPartType.directHilbert: {
                            part = ComputePartitionHilbert(Functype:1);
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
                                part = iom.LoadVector<int>(partHelp.Guid, ref cp).ToArray();
                            }

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
        /// <param name="cellWeightsLocal">
        /// If not null, defines the weight associted with each cell on the current process
        /// </param>
        /// <param name="noOfPartitioningsToChooseFrom">
        /// Tells METIS to compute
        /// </param>
        /// <returns>
        /// Index: local cell index, content: MPI Processor rank;<br/>
        /// This is the suggestion
        /// of ParMETIS for the grid partitioning:
        /// For each local cell index, the returned array contains the MPI
        /// process rank where the cell should be placed.
        /// </returns>
        public int[] ComputePartitionMETIS(int[] cellWeightsLocal = null, int noOfPartitioningsToChooseFrom = 1) {
            using (new FuncTrace()) {
                int size = this.Size;
                int rank = this.MyRank;


                if (size == 1) {
                    return new int[NoOfUpdateCells];
                }

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
                    IEnumerable<Neighbour>[] neighboursLocal = GetCellNeighbourship(IncludeBcCells: false).Take(NoOfUpdateCells).ToArray();
                    if (rank == 0) {
                        int localOffset = m_CellPartitioning.GetI0Offest(rank);
                        int localLength = m_CellPartitioning.GetLocalLength(rank);

                        for (int i = 0; i < localLength; i++) {
                            neighboursGlobal[localOffset + i] = neighboursLocal[i];
                        }
                    } else {
                        sms.Transmitt(0, neighboursLocal);
                    }

                    while (sms.GetNext(out int senderRank, out IEnumerable<Neighbour>[] neighbours)) {
                        int localOffset = m_CellPartitioning.GetI0Offest(senderRank);
                        int localLength = m_CellPartitioning.GetLocalLength(senderRank);

                        if (neighbours.Length != localLength) {
                            throw new Exception();
                        }

                        for (int i = 0; i < localLength; i++) {
                            neighboursGlobal[localOffset + i] = neighbours[i];
                        }
                    }
                }

                // Gather global weights on rank 0
                int[] cellWeightsGlobal = null;
                if (cellWeightsLocal != null) {
                    cellWeightsGlobal = new int[J];
                    if (rank == 0) {
                        int localOffset = m_CellPartitioning.GetI0Offest(rank);
                        int localLength = m_CellPartitioning.GetLocalLength(rank);

                        for (int i = 0; i < localLength; i++) {
                            cellWeightsGlobal[localOffset + i] = cellWeightsLocal[i];
                        }
                    } else {
                        sms.Transmitt(0, cellWeightsLocal);
                    }

                    while (sms.GetNext(out int senderRank, out int[] cellWeights)) {
                        int localOffset = m_CellPartitioning.GetI0Offest(senderRank);
                        int localLength = m_CellPartitioning.GetLocalLength(senderRank);

                        if (cellWeights.Length != localLength) {
                            throw new Exception();
                        }

                        for (int i = 0; i < localLength; i++) {
                            cellWeightsGlobal[localOffset + i] = cellWeights[i];
                        }
                    }
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

                    // Call METIS
                    int nparts = size;
                    Debug.Assert((cellWeightsGlobal == null) == (cellWeightsLocal == null));
                    int ncon = 1;  // One weight per vertex/cell
                    int objval = -1; // Value of the objective function at return time

                    int[] Options = new int[METIS.METIS_NOPTIONS];
                    Options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = noOfPartitioningsToChooseFrom; // 5 cuts
                    Options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
                    Options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)

                    METIS.ReturnCodes status = (METIS.ReturnCodes)METIS.PartGraphKway(
                        nvtxs: ref J,
                        ncon: ref ncon,
                        xadj: xadj,
                        adjncy: adjncy.ToArray(),
                        vwgt: cellWeightsGlobal, // if null, METIS assumes all have weight 1
                        vsize: null, // No information about communication size
                        adjwgt: null, // No edge weights
                        nparts: ref nparts,
                        tpwgts: null, // No weights for partition constraints
                        ubvec: null, // No imbalance tolerance for constraints
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
                }

                int[] localLengths = new int[size];
                for (int p = 0; p < localLengths.Length; p++) {
                    localLengths[p] = this.CellPartitioning.GetLocalLength(p);
                }
                int[] localResult = globalResult.MPIScatterv(localLengths);

                return localResult;
            }
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
                    var CellNeighs = GetCellNeighbourship(IncludeBcCells: false);
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
            int[] localToGlobalIndexMap = null; // Maps local cell index to global index in Hilbert curve

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

        private BoundingBox GetGridBoundingBox() {
            int D = this.SpatialDimension;
            var BB = new BoundingBox(D);

            int J0 = this.CellPartitioning.i0;
            int JE = this.CellPartitioning.iE;
            for (int j = J0; j < JE; j++) {
                Cell Cj = this.Cells[j - J0];
                BB.AddPoints(Cj.TransformationParams);
            }

            BB.Min = BB.Min.MPIMin();
            BB.Max = BB.Max.MPIMax();

            return BB;
        }

        /// <summary>
        /// Computes a grid partitioning (which cell should be on which processor) based on a Hilbertcurve of maximum order (64 bit>nBit*nDim).
        /// </summary>
        public int[] ComputePartitionHilbert(IList<int[]> localcellCosts = null, int Functype = 0) {
#if DEBUG
            System.Threading.Thread.Sleep(5000);
#endif
            //Functype: Constraint mapping (0) or direct Costmapping (1)
            //Notice: Functype=1 will lead to bad behavior, when using Clusters

            if (this.Size > 1) {
                int D = this.SpatialDimension;
                var GlobalBB = this.GetGridBoundingBox();

                int J0 = this.CellPartitioning.i0;
                int JE = this.CellPartitioning.iE;

                ulong[] discreteCenter = new ulong[D];
                ulong[] local_HilbertIndex = new ulong[JE - J0];
                int[] local_CellIndex = new int[JE - J0];

                for (int j = J0; j < JE; j++) {
                    Cell Cj = this.Cells[j - J0];
                    int NoOfNodes = Cj.TransformationParams.NoOfRows;
                    //Compute Barycenter
                    for (int d = 0; d < D; d++) {
                        double center = 0;
                        for (int k = 0; k < NoOfNodes; k++) {
                            center += Cj.TransformationParams[k, d];
                        }

                        center = center / ((double)NoOfNodes); // ''center of gravity'' for coordinate direction 'd'
                        double centerTrf = (center - GlobalBB.Min[d]) * (1.0 / (GlobalBB.Max[d] - GlobalBB.Min[d])) * Math.Pow(2, 64 / D);
                        //double centerTrf = (center - GlobalBB.Min[d]) * (1.0 / (GlobalBB.Max[d] - GlobalBB.Min[d])) * ((double)long.MaxValue);
                        centerTrf = Math.Round(centerTrf);
                        if (centerTrf < 0)
                            centerTrf = 0;
                        if (centerTrf > ulong.MaxValue)
                            centerTrf = ulong.MaxValue;
                        discreteCenter[d] = (ulong)centerTrf;
                        //Debugger.Break();
                    }
                    //Derive Hilbertindex from Barycenter
                    ulong iH = HilbertCurve.hilbert_c2i(64 / D, discreteCenter);
                    local_HilbertIndex[j - J0] = iH;
                    local_CellIndex[j - J0] = j;
                }

                //Gather all stuff for computation on rank==0
                int numberofcells = this.NumberOfCells;
                int[] CellsPerRank = new int[this.Size];
                for (int r = 0; r < CellsPerRank.Length; r++) {
                    CellsPerRank[r] = this.CellPartitioning.GetLocalLength(r);
                }
                
                int[] CellIndex = local_CellIndex.MPIGatherv(CellsPerRank);
                ulong[] HilbertIndex = local_HilbertIndex.MPIGatherv(CellsPerRank);
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
                            for (int field = 0; field < CellsPerCluster.Length; field++) {
                                int counter = CellsPerCluster[field] % numproc;
                                int base_size = CellsPerCluster[field] / numproc; //sauberer: Math.Methode hierfür finden
                                for (int tRank = 0; tRank < buckets.GetLength(0); tRank++)
                                    buckets[tRank, field] = base_size;
                                // Distribution of rest, if distribution is uneven
                                while (counter != 0) {
                                    counter--;
                                    buckets[counter, field]++;
                                }
                            }

                            //Go through HilbertCurve and map cell->rank
                            int rank = 0;
                            for (int cell = 0; cell < CostClustermap.Length;) {
                                if (buckets[rank, CostClustermap[cell]] != 0) {
                                    buckets[rank, CostClustermap[cell]]--;
                                    RankIndex[cell] = rank;
                                    cell++;
                                } else {
                                    rank = (rank + 1) % (numproc);
                                }
                            }
                            break;
                        case 1:
                            //direct cost mapping
                            if (cellCosts.Count > 1)
                                throw new ArgumentOutOfRangeException("Only one CellCost map allowed for directHilbert! Select Hilbert if you want to use Clusters!");
                            CostClustermap = cellCosts[0];
                            Array.Sort(HilbertIndex_tmp, CostClustermap);
                            int CellCostSum = 0;
                            for (int cell = 0; cell < CostClustermap.Length; cell++)
                                CellCostSum += CostClustermap[cell];
                            int CostPerRank = CellCostSum / numproc;
                            Debug.Assert(CostPerRank > 0);
                            int MPIrank = 0;
                            int CostCount = 0;
                            for (int cell = 0; cell < numberofcells; cell++) {
                                if (MPIrank == numproc - 1) {
                                    RankIndex[cell] = MPIrank;
                                } else {
                                    RankIndex[cell] = MPIrank;
                                    CostCount += CostClustermap[cell];
                                    if (CostCount >= CostPerRank) {
                                        MPIrank++;
                                        CostCount = 0;
                                    }
                                }
                            }
                            Debug.Assert(MPIrank <= numproc - 1);
                            break;
                        default:
                            throw new NotImplementedException();
                    }
                    Array.Sort(CellIndex, RankIndex);
                }
                //Scatter Rank-Array for local Process
                int[] local_Rank_RedistributionList = RankIndex.MPIScatterv(CellsPerRank);
                return local_Rank_RedistributionList;
            } else {
                int[] local_Rank_RedistributionList = new int[NoOfUpdateCells];
                return local_Rank_RedistributionList;
            }
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
                    sm.Transmitt(p, exch[p]);
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
