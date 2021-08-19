using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Foundation.ConstrainedDGprojection {

    /// <summary>
    /// Projection of a DG field onto a continuous subspace of the DG space:
    /// Here, the projection operation is split up into local patches.
    /// This scales much better with larger grids, but is only an approximation to the L2-projection 
    /// and likely produce higher oscillations.
    /// </summary>
    public class ConstrainedDgField_Patchwise : ConstrainedDGField {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="b"><see cref="ConstrainedDGField.Basis"/></param>
        /// <param name="__domainLimit">
        /// <see cref="ConstrainedDGField.domainLimit"/>
        /// </param>
        /// <param name="NoOfPatchesPerProcess">
        /// Number of patches which should be used on the current MPI process.
        /// if 0 or negative, determined automatically.
        /// </param>
        public ConstrainedDgField_Patchwise(Basis b, CellMask __domainLimit, int NoOfPatchesPerProcess = -1) : base(b, __domainLimit) {
            m_grd = (GridData)b.GridDat;
            ProjectDGField_patchwise_setup(NoOfPatchesPerProcess);
        }

        GridData m_grd;

        /// <summary>
        /// release of internal solvers
        /// </summary>
        public override void Dispose() {
            foreach(var s in localSeaming)
                s.Dispose();
            localSeaming.Clear();
            foreach(var s in localPatches)
                s.Dispose();
            localPatches.Clear();
            if(globalSeaming != null)
                globalSeaming.Dispose();
            globalSeaming = null;


        }


        /// <summary>
        /// Threshold for the patch-wise projection; patches are determined
        /// that the size of the linear system stays below this number.
        /// </summary>
        int maxNoOfCoordinates = 10000;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="NoOfPatchesPerProcess">
        /// number of patches on current MPI process
        /// </param>
        void ProjectDGField_patchwise_setup(int NoOfPatchesPerProcess) {
            using(new FuncTrace()) {
                MPICollectiveWatchDog.Watch();

                if(diagnosticOutput)
                    Console.WriteLine("starting patch-wise procedure: No of (local) cells {0}", domainLimit.NoOfItemsLocally);

                SubGrid maskSG = new SubGrid(domainLimit);
                EdgeMask innerEM = maskSG.InnerEdgesMask;

                int J = m_grd.CellPartitioning.LocalLength;


                // continuity projection between processes
                // =======================================
                //EdgeMask fixedInterProcBoundary = null;
                if(m_grd.MpiSize > 1) {

                    int mpiRank = this.m_grd.MpiRank;

                    int[][] edgeSendLists = this.m_grd.Edges.EdgeSendLists;
                    //List<int> ownedInterProcEdges = new List<int>();
                    int[][] edgeInsertLists = this.m_grd.Edges.EdgeInsertLists;
                    List<int> interProcEdges = new List<int>();
                    for(int proc = 0; proc < this.m_grd.MpiSize; proc++) {
                        if(edgeSendLists[proc] != null) {
                            //ownedInterProcEdges.AddRange(edgeSendLists[proc]);
                            interProcEdges.AddRange(edgeSendLists[proc]);
                        }
                        if(edgeInsertLists[proc] != null) {
                            interProcEdges.AddRange(edgeInsertLists[proc]);
                        }
                    }

                    int[,] cellInd = this.m_grd.Edges.CellIndices;
                    BitArray localInterProcBA = new BitArray(J);
                    BitArray fixedInterProcBoundaryBA = new BitArray(m_grd.Edges.Count);
                    foreach(int edge in interProcEdges) {
                        int cell1 = cellInd[edge, 0];
                        int cell2 = cellInd[edge, 1];
                        Debug.Assert(!(cell1 < J && cell2 < J), "both cells stored locally: no interproc edge!");
                        if(cell1 < J) {
                            if(domainLimit.Contains(cell1)) {
                                localInterProcBA[cell1] = true;
                                fixedInterProcBoundaryBA[edge] = true;
                            }
                        } else {
                            if(domainLimit.Contains(cell2)) {
                                localInterProcBA[cell2] = true;
                                fixedInterProcBoundaryBA[edge] = true;
                            }
                        }
                    }
                    CellMask localInterProcPatch = new CellMask(this.m_grd, localInterProcBA);

                    // add neighbours
                    int NoLocalNeigh = 2; // Gets the neighbor and the neighbor's neighbor
                    CellMask localInterProcNeigh = localInterProcPatch;
                    BitArray localInterProcNeighBA = new BitArray(J);
                    for(int n = 0; n < NoLocalNeigh; n++) {
                        foreach(int j in localInterProcNeigh.ItemEnum) {
                            localInterProcNeighBA[j] = true;
                            int[] neigh = m_grd.Cells.CellNeighbours[j];
                            foreach(int jNeigh in neigh) {
                                if(jNeigh >= 0 && jNeigh < J && domainLimit.Contains(jNeigh))
                                    localInterProcNeighBA[jNeigh] = true;
                            }

                        }
                        localInterProcNeigh = new CellMask(domainLimit.GridData, localInterProcNeighBA);
                    }


                    //this.ProjectDGField_global(domainLimit, localInterProcNeigh);
                    globalSeaming = new ConstrainedProjectionInternal(this, domainLimit, localInterProcNeigh, false);
                }

                // determine patches per process
                // =============================

                int NoPatches = 1;
                if(NoOfPatchesPerProcess > 0) {
                    NoPatches = NoOfPatchesPerProcess;
                } else {
                    int NoOfCoordOnProc = domainLimit.NoOfItemsLocally * Basis.Length;
                    if(NoOfCoordOnProc > maxNoOfCoordinates) {
                        NoPatches = (NoOfCoordOnProc / maxNoOfCoordinates) + 1;
                    }
                }
                Console.WriteLine("No of local patches: {0}", NoPatches);

                // divide projection domain into non-overlapping patches
                var m_BlockingStrategy = new METISBlockingStrategy() {
                    NoOfPartsOnCurrentProcess = NoPatches
                };
                var blocking = m_BlockingStrategy.GetBlocking(m_grd, domainLimit);

                // define patches
                List<CellMask> patches = new List<CellMask>();
                foreach(List<int> block in blocking) {

                    BitArray patchBA = new BitArray(J);
                    foreach(int j in block) {
                        patchBA[j] = true;
                    }
                    CellMask patch = new CellMask(m_grd, patchBA);
                    patches.Add(patch);
                }


                // determine interpatch domain 
                // ===========================
                EdgeMask interPatchEM = innerEM;
                foreach(CellMask patch in patches) {
                    SubGrid maskPatch = new SubGrid(patch);
                    EdgeMask innerPatch = maskPatch.InnerEdgesMask;
                    interPatchEM = interPatchEM.Except(innerPatch);
                }
                Console.WriteLine("inter patch EM No of edges: {0}", interPatchEM.NoOfItemsLocally);

                BitArray interPatchBA = new BitArray(J);
                foreach(int j in interPatchEM.ItemEnum) {
                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    if(cell1 >= 0 && cell1 < J)
                        interPatchBA[cell1] = true;
                    int cell2 = m_grd.Edges.CellIndices[j, 1];
                    if(cell2 >= 0 && cell2 < J)
                        interPatchBA[cell2] = true;
                }
                CellMask interPatchCM = new CellMask(domainLimit.GridData, interPatchBA);
                // add neighbours
                int NoNeigh = 2;
                CellMask interPatchNeigh = interPatchCM;
                BitArray interPatchNeighBA = new BitArray(J);
                for(int n = 0; n < NoNeigh; n++) {
                    foreach(int j in interPatchNeigh.ItemEnum) {
                        interPatchNeighBA[j] = true;
                        int[] neigh = m_grd.Cells.CellNeighbours[j];
                        foreach(int jNeigh in neigh) {
                            if(jNeigh >= 0 && jNeigh < J && domainLimit.Contains(jNeigh))
                                interPatchNeighBA[jNeigh] = true;
                        }

                    }
                    interPatchNeigh = new CellMask(domainLimit.GridData, interPatchNeighBA);
                }

                // continuity projection between patches within one process
                // ========================================================
                if(NoPatches > 1) {

                    if(diagnosticOutput) {
                        Console.WriteLine("======================");
                        Console.WriteLine("project local merging patch: No of cells {0}", interPatchNeigh.NoOfItemsLocally);
                    }

                    //this.ProjectDGFieldOnPatch(domainLimit, interPatchNeigh);
                    localSeaming.Add(new ConstrainedProjectionInternal(this, domainLimit, interPatchNeigh, true));
                }


                // constrained projection on all patches
                // =====================================
                int pC = 0;
                foreach(CellMask patch in patches) {
                    //this.ProjectDGFieldOnPatch(domainLimit, patch);
                    localPatches.Add(new ConstrainedProjectionInternal(this, domainLimit, patch, true));
                    pC++;
                }
            }
        }
        
        List<ConstrainedProjectionInternal> localPatches = new List<ConstrainedProjectionInternal>();
        List<ConstrainedProjectionInternal> localSeaming = new List<ConstrainedProjectionInternal>();
        ConstrainedProjectionInternal globalSeaming;

        /// <summary>
        /// Projects some DG field <paramref name="orgDGField"/> onto the internal, continuous representation.
        /// </summary>
        /// <param name="orgDGField">
        /// input; unchanged on exit
        /// </param>
        public override void ProjectDGField(ConventionalDGField orgDGField) {
            SetDGCoordinatesOnce(orgDGField, domainLimit);
            UpdateInternalProjection(csMPI.Raw._COMM.WORLD);

            ProjectDGField_SinglePass();

            UpdateInternalProjection(csMPI.Raw._COMM.WORLD);
        }


        void ProjectDGField_SinglePass() {
            if(globalSeaming != null) {
                globalSeaming.PerformProjection();
            }

            foreach(var s in localSeaming)
                s.PerformProjection();
            foreach(var s in localPatches)
                s.PerformProjection();
        }






        /*
        /// <summary>
        /// 
        /// </summary>
        /// <param name="NoOfPatchesPerProcess">
        /// number of patches on current MPI process
        /// </param>
        void ProjectDGField_patchwise_setup(int NoOfPatchesPerProcess = 0) {
            using(new FuncTrace()) {
                MPICollectiveWatchDog.Watch();

                if(diagnosticOutput)
                    Console.WriteLine("starting patch-wise procedure: No of (local) cells {0}", domainLimit.NoOfItemsLocally);

                SubGrid maskSG = new SubGrid(domainLimit);
                EdgeMask innerEM = maskSG.InnerEdgesMask;

                int J = m_grd.CellPartitioning.LocalLength;

                //List<DGField> returnFields = new List<DGField>();

                UpdateInternalProjection(csMPI.Raw._COMM.WORLD);

                // continuity projection between processes
                // =======================================
                //EdgeMask fixedInterProcBoundary = null;
                if(m_grd.MpiSize > 1) {

                    int mpiRank = this.m_grd.MpiRank;

                    int[][] edgeSendLists = this.m_grd.Edges.EdgeSendLists;
                    //List<int> ownedInterProcEdges = new List<int>();
                    int[][] edgeInsertLists = this.m_grd.Edges.EdgeInsertLists;
                    List<int> interProcEdges = new List<int>();
                    for(int proc = 0; proc < this.m_grd.MpiSize; proc++) {
                        if(edgeSendLists[proc] != null) {
                            //ownedInterProcEdges.AddRange(edgeSendLists[proc]);
                            interProcEdges.AddRange(edgeSendLists[proc]);
                        }
                        if(edgeInsertLists[proc] != null) {
                            interProcEdges.AddRange(edgeInsertLists[proc]);
                        }
                    }

                    int[,] cellInd = this.m_grd.Edges.CellIndices;
                    BitArray localInterProcBA = new BitArray(J);
                    BitArray fixedInterProcBoundaryBA = new BitArray(m_grd.Edges.Count);
                    foreach(int edge in interProcEdges) {
                        int cell1 = cellInd[edge, 0];
                        int cell2 = cellInd[edge, 1];
                        Debug.Assert(!(cell1 < J && cell2 < J), "both cells stored locally: no interproc edge!");
                        if(cell1 < J) {
                            if(domainLimit.Contains(cell1)) {
                                localInterProcBA[cell1] = true;
                                fixedInterProcBoundaryBA[edge] = true;
                            }
                        } else {
                            if(domainLimit.Contains(cell2)) {
                                localInterProcBA[cell2] = true;
                                fixedInterProcBoundaryBA[edge] = true;
                            }
                        }
                    }
                    CellMask localInterProcPatch = new CellMask(this.m_grd, localInterProcBA);
                    //fixedInterProcBoundary = new EdgeMask(this.m_grd, fixedInterProcBoundaryBA);
                    //Console.WriteLine("No of edges on inter process boundary: {0}", fixedInterProcBoundary.NoOfItemsLocally);
                    //innerEM = innerEM.Except(fixedInterProcBoundary);

                    // add neighbours
                    int NoLocalNeigh = 2; // Gets the neighbor and the neighbor's neighbor
                    CellMask localInterProcNeigh = localInterProcPatch;
                    BitArray localInterProcNeighBA = new BitArray(J);
                    for(int n = 0; n < NoLocalNeigh; n++) {
                        foreach(int j in localInterProcNeigh.ItemEnum) {
                            localInterProcNeighBA[j] = true;
                            int[] neigh = m_grd.Cells.CellNeighbours[j];
                            foreach(int jNeigh in neigh) {
                                if(jNeigh >= 0 && jNeigh < J && domainLimit.Contains(jNeigh))
                                    localInterProcNeighBA[jNeigh] = true;
                            }

                        }
                        localInterProcNeigh = new CellMask(domainLimit.GridData, localInterProcNeighBA);
                    }

                    // plot patches for debugging
                    //SinglePhaseField interProcField = new SinglePhaseField(m_Basis, "interProc");
                    //interProcField.AccConstant(1.0, localInterProcNeigh);
                    //returnFields.Add(interProcField);

                    //SubGrid localInterProcSbgrd = new SubGrid(localInterProcPatch);
                    //EdgeMask lipPatchInnerEM = localInterProcSbgrd.InnerEdgesMask;
                    //BitArray lipPatchInnerBA = new BitArray(lipPatchInnerEM.GetBitMask().Length);
                    //foreach (int edge in lipPatchInnerEM.ItemEnum) {
                    //    if (ownedInterProcEdges.Contains(edge))
                    //        lipPatchInnerBA[edge] = true;
                    //    if (!interProcEdges.Contains(edge))
                    //        lipPatchInnerBA[edge] = true;
                    //}
                    //lipPatchInnerEM = new EdgeMask(this.m_grd, lipPatchInnerBA);

                    //EdgeMask lipPatchBoundaryEM = localInterProcSbgrd.BoundaryEdgesMask;
                    //EdgeMask fixedBoundaryMask = lipPatchBoundaryEM.Intersect(innerEM);

                    //if(diagnosticOutput) {
                    //    Console.WriteLine("======================");
                    //    Console.WriteLine("project on interProc patch: No of local/wExternal cells {0}/{1}", localInterProcNeigh.NoOfItemsLocally, localInterProcNeigh.NoOfItemsLocally_WithExternal);
                    //}
                    this.ProjectDGField_global(domainLimit, localInterProcNeigh);
                    //if(diagnosticOutput) {
                    //    double jumpNorm = CheckLocalProjection(localInterProcNeigh, true);
                    //    Console.WriteLine("L2 jump norm = {0}", jumpNorm);
                    //}
                }

                // determine patches per process
                // =============================

                int NoPatches = 1;
                if(NoOfPatchesPerProcess > 0) {
                    NoPatches = NoOfPatchesPerProcess;
                } else {
                    int NoOfCoordOnProc = domainLimit.NoOfItemsLocally * m_Basis.Length;
                    if(NoOfCoordOnProc > maxNoOfCoordinates) {
                        NoPatches = (NoOfCoordOnProc / maxNoOfCoordinates) + 1;
                    }
                }
                Console.WriteLine("No of local patches: {0}", NoPatches);

                // divide projection domain into non-overlapping patches
                var m_BlockingStrategy = new METISBlockingStrategy() {
                    NoOfPartsOnCurrentProcess = NoPatches
                };
                var blocking = m_BlockingStrategy.GetBlocking(m_grd, domainLimit);

                // define patches
                List<CellMask> patches = new List<CellMask>();
                foreach(List<int> block in blocking) {

                    BitArray patchBA = new BitArray(J);
                    foreach(int j in block) {
                        patchBA[j] = true;
                    }
                    CellMask patch = new CellMask(m_grd, patchBA);
                    patches.Add(patch);
                }


                // determine interpatch domain 
                // ===========================
                EdgeMask interPatchEM = innerEM;
                foreach(CellMask patch in patches) {
                    SubGrid maskPatch = new SubGrid(patch);
                    EdgeMask innerPatch = maskPatch.InnerEdgesMask;
                    interPatchEM = interPatchEM.Except(innerPatch);
                }
                Console.WriteLine("inter patch EM No of edges: {0}", interPatchEM.NoOfItemsLocally);

                BitArray interPatchBA = new BitArray(J);
                foreach(int j in interPatchEM.ItemEnum) {
                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    if(cell1 >= 0 && cell1 < J)
                        interPatchBA[cell1] = true;
                    int cell2 = m_grd.Edges.CellIndices[j, 1];
                    if(cell2 >= 0 && cell2 < J)
                        interPatchBA[cell2] = true;
                }
                CellMask interPatchCM = new CellMask(domainLimit.GridData, interPatchBA);
                // add neighbours
                int NoNeigh = 2;
                CellMask interPatchNeigh = interPatchCM;
                BitArray interPatchNeighBA = new BitArray(J);
                for(int n = 0; n < NoNeigh; n++) {
                    foreach(int j in interPatchNeigh.ItemEnum) {
                        interPatchNeighBA[j] = true;
                        int[] neigh = m_grd.Cells.CellNeighbours[j];
                        foreach(int jNeigh in neigh) {
                            if(jNeigh >= 0 && jNeigh < J && domainLimit.Contains(jNeigh))
                                interPatchNeighBA[jNeigh] = true;
                        }

                    }
                    interPatchNeigh = new CellMask(domainLimit.GridData, interPatchNeighBA);
                }

                // continuity projection between patches within one process
                // ========================================================
                if(NoPatches > 1) {

                    if(diagnosticOutput) {
                        Console.WriteLine("======================");
                        Console.WriteLine("project local merging patch: No of cells {0}", interPatchNeigh.NoOfItemsLocally);
                    }

                    //if(fixedInterProcBoundary != null) {
                    //    BitArray fixedBoundaryBA = new BitArray(m_grd.Edges.Count);
                    //    foreach(int edg in fixedInterProcBoundary.ItemEnum) {
                    //        int cell1 = m_grd.Edges.CellIndices[edg, 0];
                    //        if(cell1 >= 0 && cell1 < J && interPatchNeigh.Contains(cell1))
                    //            fixedBoundaryBA[edg] = true;
                    //        int cell2 = m_grd.Edges.CellIndices[edg, 1];
                    //        if(cell2 >= 0 && cell2 < J && interPatchNeigh.Contains(cell2))
                    //            fixedBoundaryBA[edg] = true;
                    //    }
                    //    EdgeMask fixedBoundary = new EdgeMask(m_grd, fixedBoundaryBA);
                    //    Console.WriteLine("inter patch neighbor fixed boundary EM No of cells: {0}", fixedBoundary.NoOfItemsLocally);
                    //    this.ProjectDGFieldOnPatch(domainLimit, interPatchNeigh);
                    //} else {
                    //    this.ProjectDGFieldOnPatch(domainLimit, interPatchNeigh);
                    //}
                    this.ProjectDGFieldOnPatch(domainLimit, interPatchNeigh);

                    //if(diagnosticOutput) {
                    //    double jumpNorm = CheckLocalProjection(interPatchNeigh);
                    //    Console.WriteLine("L2 jump norm = {0}", jumpNorm);
                    //}

                }


                // constrained projection on all patches
                // =====================================
                int pC = 0;
                foreach(CellMask patch in patches) {

                    //EdgeMask fixedBoundary = (fixedInterProcBoundary != null) ? interPatchEM.Union(fixedInterProcBoundary) : interPatchEM;

                    //BitArray fixedBoundaryBA = new BitArray(m_grd.Edges.Count);
                    //foreach(int edg in fixedBoundary.ItemEnum) {
                    //    int cell1 = m_grd.Edges.CellIndices[edg, 0];
                    //    if(cell1 >= 0 && cell1 < J && patch.Contains(cell1))
                    //        fixedBoundaryBA[edg] = true;
                    //    int cell2 = m_grd.Edges.CellIndices[edg, 1];
                    //    if(cell2 >= 0 && cell2 < J && patch.Contains(cell2))
                    //        fixedBoundaryBA[edg] = true;
                    //}
                    //fixedBoundary = new EdgeMask(m_grd, fixedBoundaryBA);

                    //if(diagnosticOutput) {
                    //    Console.WriteLine("======================");
                    //    Console.WriteLine("project patch {0}: No of cells {1}", pC, patch.NoOfItemsLocally);
                    //    Console.WriteLine("local patch fixed boundary EM No of cells: {0}", fixedBoundary.NoOfItemsLocally);
                    //}
                    //this.ProjectDGFieldOnPatch(patch, fixedBoundary);
                    //if(diagnosticOutput) {
                    //    double jumpNorm = CheckLocalProjection(patch);
                    //    Console.WriteLine("L2 jump norm = {0}", jumpNorm);
                    //}
                    this.ProjectDGFieldOnPatch(domainLimit, patch);
                    pC++;
                }

                //// plot patches for debugging
                //SinglePhaseField patchField = new SinglePhaseField(m_Basis, "Patches");
                //int pColor = 1;
                //foreach(CellMask patch in patches) {
                //    patchField.AccConstant(pColor, patch);
                //    pColor++;
                //}

                //returnFields.Add(patchField);
                //return returnFields;
            }

            */
    }
    
}
