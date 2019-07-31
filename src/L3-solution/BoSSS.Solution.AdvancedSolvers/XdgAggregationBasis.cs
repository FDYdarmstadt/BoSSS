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
using System.Linq;
using System.Text;
using ilPSP;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using ilPSP.Utils;
using BoSSS.Foundation;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using System.Collections;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Solution.AdvancedSolvers {
    public class XdgAggregationBasis : AggregationGridBasis {

        /// <summary>
        /// XDG basis on original grid
        /// </summary>
        public XDGBasis XDGBasis {
            get;
            private set;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="xb">
        /// XDG basis on original grid
        /// </param>
        /// <param name="parentBasis">
        /// basis on parent grid
        /// </param>
        /// <param name="ag">
        /// aggregation grid level.
        /// </param>
        /// <param name="inj">
        /// injection operators.
        /// </param>
        internal XdgAggregationBasis(XDGBasis xb, XdgAggregationBasis parentBasis, AggregationGridData ag, MultidimensionalArray[] inj)
            : base(xb.NonX_Basis, parentBasis, ag, inj) //
        {
            using(new FuncTrace()) {
                this.XDGBasis = xb;
                this.XCompositeBasis = new MultidimensionalArray[base.AggGrid.iLogicalCells.NoOfLocalUpdatedCells][];
            }
        }

        public void Update(MultiphaseCellAgglomerator Agglomerator) {
            using(new FuncTrace()) {
                //if(!this.XDGBasis.IsSubBasis(mmf.MaxBasis))
                //    throw new ArgumentException();

                var LsTrk = this.XDGBasis.Tracker;
                var CCBit = LsTrk.Regions.GetCutCellMask().GetBitMask();

                int N = this.DGBasis.Length;
                int JAGG = base.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                UpdateSpeciesMapping(Agglomerator);



                // mass matrix in the aggregate cell, before we orthonormalize:
                var AggCellMMb4Ortho = MultidimensionalArray.Create(N, N);

                for(int jagg = 0; jagg < JAGG; jagg++) {
                    int NoOfSpc_jagg = this.NoOfSpecies[jagg];
                    int[] AggCell = base.AggGrid.iLogicalCells.AggregateCellToParts[jagg];
                    int K = AggCell.Length;

                    int[,] sim = this.SpeciesIndexMapping[jagg];
                    Debug.Assert((sim == null) || (sim.GetLength(0) == NoOfSpc_jagg));
                    Debug.Assert((sim == null) || (sim.GetLength(1) == K));


                    if(sim == null) {
                        // nothing to do.
                        // +++++++++++++++++++++++++++++++++++++++++++++++

                        this.XCompositeBasis[jagg] = null;
                    } else {
                        // the cut-cell may require some special treatment
                        // +++++++++++++++++++++++++++++++++++++++++++++++


                        for(int iSpc_agg = 0; iSpc_agg < NoOfSpc_jagg; iSpc_agg++) { // loop over all species in aggregate cell

                            bool emptyCell = false; // true: at least one of the base cells which form 
                            //                         aggregate cell 'jagg' is empty with respect to the 'iSpc_agg'--th species.
                            //                         => this will alter the projection operator.
                            for(int k = 0; k < K; k++) {
                                if(sim[iSpc_agg, k] < 0) {
                                    emptyCell = true;
                                    break;
                                }
                            }
#if DEBUG
                            {
                                bool NonEmpty = false;
                                for(int k = 0; k < K; k++) {
                                    if(sim[iSpc_agg, k] >= 0) {
                                        NonEmpty = true;
                                        break;
                                    }
                                }
                                Debug.Assert(NonEmpty); // there should be at least one non-empty base cell (for species 'iSpc_agg')
                                //                         in aggregate cell 'jagg'
                            }
#endif

                            if(this.XCompositeBasis[jagg] == null || (this.XCompositeBasis[jagg].Length != NoOfSpc_jagg))
                                this.XCompositeBasis[jagg] = new MultidimensionalArray[NoOfSpc_jagg];

                            if(emptyCell) {

                                // compute mass matrix in aggregate cell 'jagg' for species index 'iSpc_agg'
                                // -------------------------------------------------------------------------

                                AggCellMMb4Ortho.Clear();
                                for(int k = 0; k < K; k++) { // loop over the base cells in the aggregate cell
                                    if(sim[iSpc_agg, k] >= 0) {
                                        var ExPolMtx = base.CompositeBasis[jagg].ExtractSubArrayShallow(k, -1, -1);
                                        AggCellMMb4Ortho.Multiply(1.0, ExPolMtx, ExPolMtx, 1.0, "lm", "im", "il");
                                    }
                                }


                                // change to orthonormal basis
                                // ---------------------------
                                MultidimensionalArray B = MultidimensionalArray.Create(N, N);
                                AggCellMMb4Ortho.SymmetricLDLInversion(B, default(double[]));


                                if(this.XCompositeBasis[jagg][iSpc_agg] == null)
                                    this.XCompositeBasis[jagg][iSpc_agg] = MultidimensionalArray.Create(K, N, N);

                                var X_ExPolMtx = this.XCompositeBasis[jagg][iSpc_agg];
                                var NonX_ExPolMtx = CompositeBasis[jagg];
                                X_ExPolMtx.Allocate(NonX_ExPolMtx.Lengths); // should not reallocate if lengths stay the same;

                                X_ExPolMtx.Multiply(1.0, NonX_ExPolMtx, B, 0.0, "imn", "imk", "kn");
                            } else {
                                // no empty cell with respect to species 'iSpc_agg' in aggregate cell
                                // --> non-XDG projector can be used.
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                this.XCompositeBasis[jagg][iSpc_agg] = null;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// All used species.
        /// Index: enumeration over species.
        /// </summary>
        public SpeciesId[] UsedSpecies {
            get;
            private set;
        }

        /// <summary>
        /// defines how species in the agglomerated cells
        /// map to species in the base grid cells.
        /// </summary>
        private void UpdateSpeciesMapping(MultiphaseCellAgglomerator Agglomerator) {
            Debug.Assert(object.ReferenceEquals(Agglomerator.Tracker, this.XDGBasis.Tracker));
            var LsTrk = this.XDGBasis.Tracker;
            int GlobalNoOfSpc = LsTrk.SpeciesIdS.Count;
            var agg = base.AggGrid;
            int[][] compCells = agg.iLogicalCells.AggregateCellToParts;
            int JAGG = compCells.Length;

            this.UsedSpecies = Agglomerator.SpeciesList.ToArray();

            if(NoOfSpecies == null) {
                NoOfSpecies = new int[JAGG];
                SpeciesIndexMapping = new int[JAGG][,];
            } else {
                Debug.Assert(NoOfSpecies.Length == JAGG);
                Array.Clear(NoOfSpecies, 0, JAGG);
            }

            if(AggCellsSpecies == null)
                AggCellsSpecies = new SpeciesId[JAGG][];
            //if (CoarseToFineSpeciesIndex == null && agg.MgLevel > 0)
            //    CoarseToFineSpeciesIndex = new int[JAGG][,];


            Dictionary<SpeciesId,BitArray> agglomeratedCells = new Dictionary<SpeciesId, BitArray>();
            foreach(var SpId in Agglomerator.SpeciesList) {
                agglomeratedCells.Add(SpId, Agglomerator.GetAgglomerator(SpId).AggInfo.SourceCells.GetBitMaskWithExternal());
            }
            
            // temp buffer
            int[] _NoOfSpecies = new int[agg.jCellCoarse2jCellFine[0].Length];
            ReducedRegionCode[] _RRcs = new ReducedRegionCode[agg.iLogicalCells.AggregateCellToParts[0].Length];
            SpeciesId[] allPresentSpecies = new SpeciesId[LsTrk.TotalNoOfSpecies];

            // loop over all aggregate (multigrid) cells...
            for(int jagg = 0; jagg < JAGG; jagg++) {
                int[] compCell = compCells[jagg];
                int K = compCell.Length;

                // buffer sizes
                if(K > _NoOfSpecies.Length) {
                    int newLen = (int)Math.Ceiling(K * 1.2);
                    Array.Resize(ref _NoOfSpecies, newLen);
                    Array.Resize(ref _RRcs, newLen);
                }

                // check no of species:
                int MaxNoOfSpecies = 0; // number of species in the composite cell
                for(int k = 0; k < K; k++) { // loop over original cells in composite cell
                    int jCell = compCell[k];
                    _NoOfSpecies[k] = LsTrk.Regions.GetNoOfSpecies(jCell, out _RRcs[k]);
                    MaxNoOfSpecies = Math.Max(_NoOfSpecies[k], MaxNoOfSpecies);

                }
                
                if(MaxNoOfSpecies == 1) {
                    // only one species -- use single-phase implementation in underlying class

                    this.SpeciesIndexMapping[jagg] = null;
                    this.NoOfSpecies[jagg] = 1;
                    //this.XCompositeBasis[jagg] = null;
                    this.AggCellsSpecies[jagg] = null;
                    //this.CoarseToFineSpeciesIndex[jagg] = null;
                } else {
                    //LsTrk.ContainesSpecies

                    int w = 0; // counter for found species
                    
                    if(this.SpeciesIndexMapping[jagg] == null || this.SpeciesIndexMapping[jagg].GetLength(0) != w) {
                        this.SpeciesIndexMapping[jagg] = new int[0, compCell.Length];
                    }
                    ArrayTools.SetAll(this.SpeciesIndexMapping[jagg], -897645);

                    
                    // collect all species which are present in the composite cell 'jagg':
                    // -------------------------------------------------------------------
                    for(int k = 0; k < K; k++) { // loop over all base cells in the aggregate cell 'jagg'
                        var rrc = _RRcs[k];
                        int jCell = compCell[k]; // base cell index

                        for(int iSpc = _NoOfSpecies[k] - 1; iSpc >= 0; iSpc--) {
                            SpeciesId spId = LsTrk.GetSpeciesIdFromIndex(rrc, iSpc);
                            bool isPresent = LsTrk.Regions.IsSpeciesPresentInCell(spId, jCell);
                            bool isAgglomerated = agglomeratedCells.ContainsKey(spId) ? agglomeratedCells[spId][jCell] : false;
                            bool isUsed = Array.IndexOf(this.UsedSpecies, spId) >= 0;

                            if(isUsed && isPresent && !isAgglomerated) {
                                bool bfound = false;

                                int z;
                                for(z = 0; z < w; z++) { // loop over species found so far...
                                    if(allPresentSpecies[z] == spId) {
                                        bfound = true;
                                        break;
                                    }
                                }
                                if(!bfound) {
                                    // species not found -> add to list
                                    allPresentSpecies[w] = spId;
                                    z = w;
                                    w++;
                                }

                                // z    is the species index in the composite cell
                                // iSpc is the species index in 'jCell'
                                Resize2DArray(ref this.SpeciesIndexMapping[jagg], w, K, -897645);
                                this.SpeciesIndexMapping[jagg][z, k] = iSpc;
                            }
                        }
                    }

                    this.AggCellsSpecies[jagg] = allPresentSpecies.GetSubVector(0, w);
                    this.NoOfSpecies[jagg] = w;
                    Debug.Assert(this.NoOfSpecies[jagg] == this.SpeciesIndexMapping[jagg].GetLength(0));
                }
                Debug.Assert(this.NoOfSpecies[jagg] >= 0);
                Debug.Assert(this.NoOfSpecies[jagg] <= this.UsedSpecies.Length);
            }
        }


        void Resize2DArray<T>(ref T[,] A, int l1, int l2, T defVal) {
            if(A.GetLength(0) != l1 || A.GetLength(1) != l2) {
                T[,] O = A;
                A = new T[l1, l2];
                A.SetAll(defVal);
                for(int i = Math.Min(O.GetLength(0), A.GetLength(0)) - 1; i >= 0; i--)
                    for(int j = Math.Min(O.GetLength(1), A.GetLength(1)) - 1; j >= 0; j--)
                        A[i, j] = O[i, j];

                //for(int i = O.GetLength(0); i < A.GetLength(0); i++)
                //    for(int j = O.GetLength(1); j < A.GetLength(1); j++)
                //        A[i, j] = defVal;
            }

        }


        /// <summary>
        /// Number of species per composite cell; 
        ///  - index: composite cell index;
        /// </summary>
        int[] NoOfSpecies;

       

        /// <summary>
        /// Mapping from species indices in composite/aggregate cells to species index in base grid.
        ///  - 1st index: composite cell index;
        ///  - 2nd index: species index in composite/aggregate cell.
        ///  - 3rd index: enumeration of base grid cells in the composite cell.
        /// </summary>
        /// <remarks>
        /// If (<see cref="SpeciesIndexMapping"/>[j] == null), this indicates that 
        /// the single-phase implementation in underlying class should be used.
        /// </remarks>
        internal int[][,] SpeciesIndexMapping;


        /// <summary>
        /// species in the composite/aggregate cells;<br/>
        ///  - 1st index: aggregate cell index.
        ///  - 2nd index: species index in the composite/aggregate cell.
        /// </summary>
        internal SpeciesId[][] AggCellsSpecies;


        //int[][,] CoarseToFineSpeciesIndex;


        /// <summary>
        /// Returns the species index of species <paramref name="spid"/> in
        /// composite/aggregate cell <paramref name="jAgg"/> .
        /// </summary>
        public int GetSpeciesIndex(int jAgg, SpeciesId spid) {
            if(AggCellsSpecies[jAgg] != null) {
                return Array.IndexOf(AggCellsSpecies[jAgg], spid);
            } else {
                int[] BaseCells = this.AggGrid.iLogicalCells.AggregateCellToParts[jAgg];
                var LsTrk = this.XDGBasis.Tracker;

                int iSpc = LsTrk.Regions.IsSpeciesPresentInCell(spid, BaseCells[0]) ? 0 : -1;
#if DEBUG
                if(iSpc >= 0) {
                    foreach(int j in BaseCells) {
                        Debug.Assert(LsTrk.Regions.GetSpeciesIndex(spid, j) == iSpc);
                    }
                } else {
                    foreach(int j in BaseCells) {
                        Debug.Assert(LsTrk.Regions.IsSpeciesPresentInCell(spid, BaseCells[0]) == false);
                    }
                }
#endif

                return iSpc;
            }
        }

        /// <summary>
        /// Number of species in composite/aggregate cell <paramref name="jAgg"/>.
        /// </summary>
        override public int GetNoOfSpecies(int jAgg) {
            return NoOfSpecies[jAgg];
        }


        /// <summary>
        /// local vector-space dimension.
        /// </summary>
        override public int LocalDim {
            get {
                int L = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells * this.XDGBasis.MaximalLength;
                Debug.Assert(base.DGBasis.Length * this.XDGBasis.Tracker.TotalNoOfSpecies == this.XDGBasis.MaximalLength);
                return L;
            }
        }


        public override void ProlongateToFullGrid<T, V>(T FullGridVector, V AggGridVector) {

            var fullMapping = new UnsetteledCoordinateMapping(this.XDGBasis);
            if(FullGridVector.Count != fullMapping.LocalLength)
                throw new ArgumentException("mismatch in vector length", "FullGridVector");
            int L = this.LocalDim;
            if(AggGridVector.Count != L)
                throw new ArgumentException("mismatch in vector length", "AggGridVector");


            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int Nmax = this.XDGBasis.MaximalLength;
            int N = this.DGBasis.Length;
            Debug.Assert(Nmax % N == 0);


            var FulCoords = new double[N];
            var AggCoords = new double[N];


            for(int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over all composite cells...
                int[] agCl = agCls[jAgg];
                int K = agCl.Length;

                
                if(this.SpeciesIndexMapping[jAgg] == null) {
                    Debug.Assert(this.XCompositeBasis[jAgg] == null);

                    int i0 = jAgg * Nmax; // index offset into 'AggGridVector'
                    for(int n = 0; n < N; n++)
                        AggCoords[n] = AggGridVector[n + i0];
                    
                    for(int k = 0; k < K; k++) { // loop over the cells wich form the aggregated cell...
                        int jCell = agCl[k];
                        int j0 = fullMapping.LocalUniqueCoordinateIndex(0, jCell, 0);

                        MultidimensionalArray Trf = base.CompositeBasis[jAgg].ExtractSubArrayShallow(k, -1, -1);

                        //Trf.Solve(FulCoords, AggCoords);
                        Trf.gemv(1.0, AggCoords, 0.0, FulCoords);

                        for(int n = 0; n < N; n++) {
                            FullGridVector[j0 + n] = FulCoords[n];
                        }
                    }

                } else {
                    int NoSpc = this.NoOfSpecies[jAgg];
                    int[,] sim = SpeciesIndexMapping[jAgg];
                    Debug.Assert(sim.GetLength(0) == NoSpc);
                    Debug.Assert(sim.GetLength(1) == K);

                    for(int iSpcAgg = 0; iSpcAgg < NoSpc; iSpcAgg++) {

                        int i0 = jAgg * Nmax + iSpcAgg*N; // index offset into 'AggGridVector'
                        for(int n = 0; n < N; n++)
                            AggCoords[n] = AggGridVector[n + i0];


                        for(int k = 0; k < K; k++) { // loop over the cells wich form the aggregated cell...
                            int jCell = agCl[k];
                            int iSpcBase = sim[iSpcAgg, k];
                            if(iSpcBase < 0)
                                continue;

                            int j0 = fullMapping.LocalUniqueCoordinateIndex(0, jCell, iSpcBase * N);

                            MultidimensionalArray Trf;
                            if(this.XCompositeBasis[jAgg] == null || this.XCompositeBasis[jAgg][iSpcAgg] == null)
                                Trf = base.CompositeBasis[jAgg].ExtractSubArrayShallow(k, -1, -1);
                            else
                                Trf = this.XCompositeBasis[jAgg][iSpcAgg].ExtractSubArrayShallow(k, -1, -1);

                            Trf.gemv(1.0, AggCoords, 0.0, FulCoords);

                            for(int n = 0; n < N; n++) {
                                FullGridVector[j0 + n] = FulCoords[n];
                            }
                        }
                    }
                }
            }
        }


        MultidimensionalArray[][] XCompositeBasis;


        public override void RestictFromFullGrid<T, V>(T FullGridVector, V AggGridVector) {

            var fullMapping = new UnsetteledCoordinateMapping(this.XDGBasis);
            if(FullGridVector.Count != fullMapping.LocalLength)
                throw new ArgumentException("mismatch in vector length", "FullGridVector");
            int L = this.LocalDim;
            if(AggGridVector.Count != L)
                throw new ArgumentException("mismatch in vector length", "AggGridVector");


            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int Nmax = this.XDGBasis.MaximalLength;
            int N = this.DGBasis.Length;
            Debug.Assert(Nmax % N == 0);


            var Buffer = MultidimensionalArray.Create(agCls.Max(cc => cc.Length), N);
            var AggCoords = MultidimensionalArray.Create(N);



            for(int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over all composite cells...
                int[] agCl = agCls[jAgg];
                int K = agCl.Length;

                MultidimensionalArray FulCoords = (K * N == Buffer.GetLength(0)) ? Buffer : Buffer.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, N - 1 });

                int i0Agg = jAgg * Nmax, i0Full;


                if(this.SpeciesIndexMapping[jAgg] == null) {
                    // default branch:
                    // use restiction/prolongation from un-cut basis
                    // +++++++++++++++++++++++++++++++++++++++++++++

                    for(int k = 0; k < K; k++) { // loop over the cells wich form the aggregated cell...
                        int jCell = agCl[k];
                        i0Full = jCell * Nmax;
                        for(int n = 0; n < N; n++) {
                            FulCoords[k, n] = FullGridVector[n + i0Full];
                        }
                    }

                    MultidimensionalArray Trf = base.CompositeBasis[jAgg];
                    AggCoords.Clear();
                    AggCoords.Multiply(1.0, Trf, FulCoords, 0.0, "n", "kmn", "km");
                    
                    int _i0 = jAgg * N;
                    for(int n = 0; n < N; n++)
                        AggGridVector[n + i0Agg] = AggCoords[n];

                } else {
                    int NoSpc = this.NoOfSpecies[jAgg];
                    int[,] sim = SpeciesIndexMapping[jAgg];
                    Debug.Assert(sim.GetLength(0) == NoSpc);
                    Debug.Assert(sim.GetLength(1) == K);

                    for(int iSpcAgg = 0; iSpcAgg < NoSpc; iSpcAgg++) {
                        for(int k = 0; k < K; k++) { // loop over the cells wich form the aggregated cell...
                            int jCell = agCl[k];
                            int iSpcBase = sim[iSpcAgg, k];
                             
                            if(iSpcBase < 0) {
                                for(int n = 0; n < N; n++)
                                    FulCoords[k, n] = 0;
                            } else {
                                i0Full = fullMapping.LocalUniqueCoordinateIndex(0, jCell, iSpcBase * N);
                           
                                for(int n = 0; n < N; n++)
                                    FulCoords[k, n] = FullGridVector[i0Full + n];
                            }
                        }

                        MultidimensionalArray Trf;
                        if(this.XCompositeBasis[jAgg] == null || this.XCompositeBasis[jAgg][iSpcAgg] == null)
                            Trf = base.CompositeBasis[jAgg];
                        else
                            Trf = this.XCompositeBasis[jAgg][iSpcAgg];

                        AggCoords.Clear();
                        AggCoords.Multiply(1.0, Trf, FulCoords, 0.0, "n", "kmn", "km");

                        for(int n = 0; n < N; n++) {
                            AggGridVector[i0Agg + iSpcAgg * N + n] = AggCoords[n];
                        }
                    }
                }
            }
        }

        public override int MaximalLength {
            get {
                return this.XDGBasis.MaximalLength;
            }
        }

        public override int MinimalLength {
            get {
                return this.XDGBasis.MinimalLength;
            }
        }

        public override int GetLength(int jCell, int p) {
            return this.NoOfSpecies[jCell] * base.GetLength(jCell, p);
        }


        public override int GetMaximalLength(int p) {
            return this.XDGBasis.Tracker.TotalNoOfSpecies * base.GetMaximalLength(p);
        }

        public override int GetMinimalLength(int p) {
            return base.GetMinimalLength(p);
        }


        public override void GetRestrictionMatrix(BlockMsrMatrix RST, MultigridMapping mgMap, int iF) {
            if(!object.ReferenceEquals(mgMap.AggBasis[iF], this))
                throw new ArgumentException();
            
            //MsrMatrix RST = new MsrMatrix(mgMap.Partitioning, mgMap.ProblemMapping);
            int JAGG = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;

            int[] degrees = mgMap.DgDegree;
            int NoFld = degrees.Length;

            int N_rest;// = new int[NoFld];
            int N_full;// = new int[NoFld];
            {
                BoSSS.Foundation.Basis b = mgMap.ProblemMapping.BasisS[iF];
                BoSSS.Foundation.Basis NxB = (b is BoSSS.Foundation.XDG.XDGBasis) ? ((BoSSS.Foundation.XDG.XDGBasis)b).NonX_Basis : b;
                N_rest = NxB.Polynomials[0].Where(poly => poly.AbsoluteDegree <= degrees[iF]).Count();
                N_full = NxB.Length;
            }

          
            int mgMap_Offset = mgMap.Partitioning.i0;

            for(int jAgg = 0; jAgg < JAGG; jAgg++) {
                int[] AgCell = this.AggGrid.iLogicalCells.AggregateCellToParts[jAgg];
                int K = AgCell.Length;


                int NoSpc_Agg = this.NoOfSpecies[jAgg]; // number of species in aggregate cell
                for(int iSpc_Agg = 0; iSpc_Agg < NoSpc_Agg; iSpc_Agg++) { // loop over all species in aggregate cell
                    { // loop over DG fields in the mapping
                        int NROW = N_rest;
                        int NCOL = N_full;
                        Debug.Assert(mgMap.AggBasis[iF].GetLength(jAgg, degrees[iF]) == N_rest * NoSpc_Agg);
                        int i0Agg_Loc = mgMap.LocalUniqueIndex(iF, jAgg, N_rest * iSpc_Agg);
                        int i0Agg = i0Agg_Loc + mgMap_Offset;
                        Debug.Assert(i0Agg >= mgMap.Partitioning.i0);
                        Debug.Assert(i0Agg < mgMap.Partitioning.iE);
                        
                        if(this.SpeciesIndexMapping[jAgg] == null) {
                            Debug.Assert(NoSpc_Agg == 1);
                            Debug.Assert(NoSpc_Agg == 1);


                            // default branch:
                            // use restriction/prolongation from un-cut basis
                            // +++++++++++++++++++++++++++++++++++++++++++++

                            MultidimensionalArray Trf = base.CompositeBasis[jAgg];
                            
                            for(int k = 0; k < K; k++) { // loop over the cells which form the aggregated cell...
                                int jCell = AgCell[k];
                                int i0Full = mgMap.ProblemMapping.GlobalUniqueCoordinateIndex(iF, jCell, 0);
                                var Block = Trf.ExtractSubArrayShallow(k, -1, -1);


                                for(int nRow = 0; nRow < NROW; nRow++) {
                                    for(int nCol = 0; nCol < NCOL; nCol++) {
                                        RST[i0Agg + nRow, i0Full + nCol] = Block[nCol, nRow];
                                    }
                                }
                            }
                        } else {
                            int NoSpc = this.NoOfSpecies[jAgg];
                            int[,] sim = SpeciesIndexMapping[jAgg];
                            Debug.Assert(sim.GetLength(0) == NoSpc);
                            Debug.Assert(sim.GetLength(1) == K);

                            MultidimensionalArray Trf;
                            if(this.XCompositeBasis[jAgg] == null || this.XCompositeBasis[jAgg][iSpc_Agg] == null)
                                Trf = base.CompositeBasis[jAgg];
                            else
                                Trf = this.XCompositeBasis[jAgg][iSpc_Agg];
                            
                            for(int k = 0; k < K; k++) { // loop over the cells which form the aggregated cell...
                                int jCell = AgCell[k];
                                int iSpcBase = sim[iSpc_Agg, k];

                                if(iSpcBase < 0) {
                                    //for(int n = 0; n < N; n++)
                                    //    FulCoords[k, n] = 0;
                                } else {
                                    int i0Full = mgMap.ProblemMapping.GlobalUniqueCoordinateIndex(iF, jCell, iSpcBase * N_full);
                                    var Block = Trf.ExtractSubArrayShallow(k, -1, -1);

                                    for(int nRow = 0; nRow < NROW; nRow++) {
                                        for(int nCol = 0; nCol < NCOL; nCol++) {
                                            RST[i0Agg + nRow, i0Full + nCol] = Block[nCol, nRow];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }


            }
            
            //return RST;
        }

        /// <summary>
        /// for XDG, the cell mode index <paramref name="n"/> may not be equal
        /// in the full and the aggregated grid. This method performs the transformation.
        /// </summary>
        override internal int N_Murks(int j, int n, int N) {
            Debug.Assert(j >= 0 && j < this.DGBasis.GridDat.iLogicalCells.Count);
            Debug.Assert(n >= 0 && n < this.XDGBasis.GetLength(j));

            int[,] sim = this.SpeciesIndexMapping[j];
            if(sim == null) {
                Debug.Assert(this.NoOfSpecies[j] == 1);
                return n;
            } else {
                int N0 = N / this.NoOfSpecies[j];
                Debug.Assert(N0*this.NoOfSpecies[j] == N);
                int iSpc = n / N0;
                int nx = n - iSpc * N0;
                Debug.Assert(sim.GetLength(1) == 1);
                int iSpcQ = sim[iSpc, 0];
                return iSpcQ * N0 + nx;
            }
        }

        override internal bool ReqModeIndexTrafo {
            get {
                return true;
            }
        }

        public override int[] ModeIndexForDegree(int j, int p, int Pmax) {
            if(this.NoOfSpecies[j] == 1) {
                return base.ModeIndexForDegree(j, p, Pmax);
            } else {
                int N = this.GetLength(j, Pmax);
                int NSpc = this.NoOfSpecies[j];
                int N0 = N / NSpc;

                int[] Nmode = base.ModeIndexForDegree(j, p, Pmax);

                int[] R = Nmode;
                for(int iSpc = 1; iSpc < NSpc; iSpc++) {
                    int[] NmodeX = Nmode.CloneAs();
                    for(int k = 0; k < NmodeX.Length; k++) {
                        NmodeX[k] += iSpc * N0;
                    }
                    R = ArrayTools.Cat(R, NmodeX);
                }

                return R;
            }
        }
    }


}
