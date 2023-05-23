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
using System.Text;

namespace BoSSS.Foundation.ConstrainedDGprojection {

    /// <summary>
    /// Projection of a DG field onto a continuous subspace of the DG space:
    /// Here, the projection operation is split up into local patches.
    /// This scales much better with larger grids, but is only an approximation to the L2-projection 
    /// and likely produce higher oscillations.
    /// </summary>
    public class ConstrainedDgField_Patchwise : ConstrainedDGFieldMk3 {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="b"><see cref="ConstrainedDGFieldMk3.Basis"/></param>
        /// <param name="__domainLimit">
        /// <see cref="ConstrainedDGFieldMk3.domainLimit"/>
        /// </param>
        /// <param name="NoOfPatchesPerProcess">
        /// Number of patches which should be used on the current MPI process.
        /// if 0 or negative, determined automatically.
        /// </param>
        public ConstrainedDgField_Patchwise(Basis b, CellMask __domainLimit, int NoOfPatchesPerProcess = -1) : base(b, __domainLimit) {
            m_grd = (GridData)b.GridDat;
            ProjectDGField_patchwise_setup(NoOfPatchesPerProcess);
        }

        readonly GridData m_grd;

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
        const int maxNoOfCoordinates = 20000;

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
                    //SubGrid maskPatch = new SubGrid(patch);
                    //EdgeMask innerPatch = maskPatch.InnerEdgesMask;
                    EdgeMask innerPatch = patch.GetAllInnerEdgesMask();
                    interPatchEM = interPatchEM.Except(innerPatch);
                }
                Console.WriteLine("inter patch EM No of edges: {0}", interPatchEM.NoOfItemsLocally);

                BitArray interPatchBA = new BitArray(J);
                foreach(int edg in interPatchEM.ItemEnum) {
                    int cell1 = m_grd.Edges.CellIndices[edg, 0];
                    if(cell1 >= 0 && cell1 < J)
                        interPatchBA[cell1] = true;
                    int cell2 = m_grd.Edges.CellIndices[edg, 1];
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
            using(var tr = new FuncTrace()) {
                SetDGCoordinatesOnce(orgDGField, domainLimit);

                double updateNorm = 0;
                double relUpdateNorm = 0;
                bool CheckTerm(int iIter) {
                    if(TotNoOfSolvers <= 1 && iIter >= 1)
                        return false;
                    if(iIter < MinIterations)
                        return true;
                    if(iIter >= MaxIterations)
                        return false;

                    if(updateNorm < AbsTerminationThreshold)
                        return false;
                    if(relUpdateNorm < RelTerminationThreshold)
                        return false;



                    return true;
                }


                for(int i = 0; CheckTerm(i); i++) {

                    UpdateInternalProjection(csMPI.Raw._COMM.WORLD);

                    updateNorm = ProjectDGField_SinglePass();

                    UpdateInternalProjection(csMPI.Raw._COMM.WORLD);

                    double ResNorm = base.internalProjection.L2Error(orgDGField, this.domainLimit);
                    relUpdateNorm = updateNorm / Math.Max(ResNorm, 1e-30); // avoid division by zero
                    //tr.Info
                    Console.WriteLine("Iteration " + (i + 1) + ", delta to DG = " + ResNorm + "  abs change = " + updateNorm + " rel change = " + relUpdateNorm);
                }
            }
        }

        /// <summary>
        /// Minimum required iterations for <see cref="ProjectDGField"/>
        /// </summary>
        public int MinIterations = 3;

        /// <summary>
        /// Maximum allowed iterations for <see cref="ProjectDGField"/>
        /// </summary>
        public int MaxIterations = 40;

        /// <summary>
        /// Relative threshold for the termination of the iterations for <see cref="ProjectDGField"/>
        /// </summary>
        public double RelTerminationThreshold = 1e-10;


        /// <summary>
        /// Absolute threshold for the termination of the iterations for <see cref="ProjectDGField"/>
        /// </summary>
        public double AbsTerminationThreshold = 1e-15;


        List<double[]> OnbSystem = new List<double[]>();

        void AddCoord() {
            var vN = base.m_Coordinates.CloneAs();

            foreach(var v in OnbSystem) {
                double prj = v.InnerProd(vN);
                vN.AccV(-prj, v);
            }

            double alpha = vN.L2Norm();
            if(alpha > 1.0e-8) {
                vN.ScaleV(1 / alpha);

                for(int i = 0; i < OnbSystem.Count; i++) {
                    double test = OnbSystem[i].InnerProd(vN);
                    if(test > 1e-6) 
                        Console.WriteLine("orthonormality test failed: " + test);
                }

                OnbSystem.Add(vN);
            }
        }

        double Minimize() {
            double prj = 0;
            m_Coordinates.ClearEntries();
            for(int i = 0; i < OnbSystem.Count; i++) {
                prj = m_Coordinates0.InnerProd(OnbSystem[i]);
                m_Coordinates.AccV(prj, OnbSystem[i]);
            }
            return prj;
        }



        double ProjectDGField_SinglePass() {
            //double kryCh = 0;
            double GlbChange = 0.0, LocChange = 0.0;
            if(globalSeaming != null) {
                GlbChange += globalSeaming.PerformProjection().Pow2();

                //AddCoord();
                //kryCh += Minimize().Pow2();
            }

            foreach(var s in localSeaming)
                LocChange += s.PerformProjection().Pow2();
            foreach(var s in localPatches)
                LocChange += s.PerformProjection().Pow2(); // Note: since the domains are disjoint, we have almost |x + y|^2 = |x]^2 + |y]^2
            LocChange = LocChange.MPISum();

            //AddCoord();
            //kryCh += Minimize().Pow2();

            double Totchange = (GlbChange + LocChange).Sqrt();
            return Totchange;
        }

        int TotNoOfSolvers {
            get {
                int Sum = 0;
                if(globalSeaming != null)
                    Sum++;
                Sum += localPatches.Count;
                Sum += localSeaming.Count;
                return Sum;
            }
        }


    }
    
}
